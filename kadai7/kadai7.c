#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define NX 50 
#define NY 40

#define LX 100.0
#define LY 80.0

#define MAX_ITER 100000
#define EPS      1.0e-3

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    double dx = LX / NX;
    double dy = LY / NY;
    double beta  = dx / dy;
    double beta2 = beta * beta;

    static double psi[NX+1][NY+1];
    static double psi_old[NX+1][NY+1];
    static int    is_solid[NX+1][NY+1]; 

    int i, j, iter;
    double maxdiff;

    for (i = 0; i <= NX; i++) {
        for (j = 0; j <= NY; j++) {
            psi[i][j] = 0.0;
            is_solid[i][j] = 0;
        }
    }

    int iB = (int)round(30.0 / dx);
    int iE = (int)round(60.0 / dx);
    int jC = (int)round(40.0 / dy);

    for (i = iB; i <= iE; i++) {
        for (j = 0; j <= jC; j++) {
            is_solid[i][j] = 1;
            psi[i][j] = 0.0;      
        }
    }

    for (i = 0; i <= NX; i++) {
        psi[i][0] = 0.0;
    }


    for (i = 0; i <= NX; i++) {
        psi[i][NY] = (double)NY;
    }

    for (j = 0; j <= NY; j++) {
        psi[0][j]  = (double)j;
        psi[NX][j] = (double)j;
    }

    printf("差分法ラプラス方程式: NX=%d, NY=%d, dx=%g, dy=%g\n",
           NX, NY, dx, dy);
    printf("障害物: x=[30,60], y=[0,40]\n");


    for (iter = 1; iter <= MAX_ITER; iter++) {

        for (i = 0; i <= NX; i++) {
            for (j = 0; j <= NY; j++) {
                psi_old[i][j] = psi[i][j];
            }
        }

        maxdiff = 0.0;

        for (i = 1; i < NX; i++) {
            for (j = 1; j < NY; j++) {

                if (is_solid[i][j]) continue;
                if (j == 0 || j == NY) continue;
                if (i == 0 || i == NX) continue;

                double psi_new =
                    ( psi[i+1][j] + psi[i-1][j]
                    + beta2 * (psi[i][j+1] + psi[i][j-1]) )
                    / (2.0 * (1.0 + beta2));

                double diff = fabs(psi_new - psi[i][j]);
                if (diff > maxdiff) maxdiff = diff;

                psi[i][j] = psi_new;
            }
        }

        if (iter % 1000 == 0) {
            printf("iter = %6d, maxdiff = %e\n", iter, maxdiff);
        }

        if (maxdiff < EPS) {
            printf("収束しました: iter = %d, maxdiff = %e\n", iter, maxdiff);
            break;
        }
    }

    if (iter > MAX_ITER) {
        printf("MAX_ITER まで計算しても収束しませんでした。\n");
    }

    FILE *fp = fopen("psi.dat", "w");
    if (!fp) {
        printf("psi.dat を開けませんでした。\n");
        return 1;
    }

    for (j = 0; j <= NY; j++) {
        double y = j * dy;
        for (i = 0; i <= NX; i++) {
            double x = i * dx;
            fprintf(fp, "%g %g %g\n", x, y, psi[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("結果を psi.dat に出力しました。\n");

    const char *script_name = "kadai7.gp";
    FILE *gp = fopen(script_name, "w");
    if (!gp) {
        printf("gnuplot スクリプト %s を作成できませんでした。\n", script_name);
        return 1;
    }

    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "set view map\n");
    fprintf(gp, "set contour base\n");
    fprintf(gp, "unset surface\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set object 1 rect from 30,0 to 60,40 fc rgb 'gray' fs solid 0.4 border lc rgb 'black'\n");
    fprintf(gp, "set cntrparam levels discrete 0,1,2,3,4,5,6,7,8,9,10,20,30\n");
    fprintf(gp, "splot 'psi.dat' using 1:2:3 with lines\n");
    fclose(gp);

    char cmd[256];
#ifdef _WIN32
    snprintf(cmd, sizeof(cmd), "cmd /c gnuplot -persist \"%s\"", script_name);
#else
    snprintf(cmd, sizeof(cmd), "gnuplot -persist \"%s\"", script_name);
#endif
    int ret = system(cmd);
    if (ret != 0) {
        printf("kadai7 の gnuplot 実行に失敗しました (return=%d)。\n", ret);
        printf("コマンド: %s\n", cmd);
    } else {
        printf("kadai7 の結果を gnuplot で表示しました。\n");
    }

    return 0;
}
