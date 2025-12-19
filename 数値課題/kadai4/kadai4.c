#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define N     100      
#define R     0.4   
#define NSTEP 5000

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    int i, n;
    double dx = 1.0 / N;
    double dt = R * dx * dx; 
    double u_old[N + 1];
    double u_new[N + 1];
    double x;

    for (i = 0; i <= N; i++) {
        x = i * dx;
        if (x <= 0.5) {
            u_old[i] = 2.0 * x;   
        } else {
            u_old[i] = 2.0 * (1.0 - x);  
        }
    }

    u_old[0] = 0.0;
    u_old[N] = 0.0;

    FILE *fp0 = fopen("heat_t0.dat", "w");
    if (fp0 == NULL) {
        printf("ファイル heat_t0.dat を開けませんでした。\n");
        return 1;
    }
    for (i = 0; i <= N; i++) {
        x = i * dx;
        fprintf(fp0, "%lf %lf\n", x, u_old[i]);
    }
    fclose(fp0);

    int output_interval = 1000;
    char fname[64];

    for (n = 0; n < NSTEP; n++) {

        for (i = 1; i < N; i++) {
            u_new[i] = u_old[i]
                      + R * (u_old[i + 1] - 2.0 * u_old[i] + u_old[i - 1]);
        }

        u_new[0] = 0.0;
        u_new[N] = 0.0;

        for (i = 0; i <= N; i++) {
            u_old[i] = u_new[i];
        }

        if ((n + 1) % output_interval == 0) {
            snprintf(fname, sizeof(fname), "heat_t%04d.dat", n + 1);
            FILE *fp_mid = fopen(fname, "w");
            if (fp_mid != NULL) {
                for (i = 0; i <= N; i++) {
                    x = i * dx;
                    fprintf(fp_mid, "%lf %lf\n", x, u_old[i]);
                }
                fclose(fp_mid);
                printf("t ステップ %d の分布を %s に出力しました。\n", n + 1, fname);
            } else {
                printf("中間ファイル %s を開けませんでした。\n", fname);
            }
        }
    }

    FILE *fp = fopen("heat_tend.dat", "w");
    if (fp == NULL) {
        printf("ファイル heat_tend.dat を開けませんでした。\n");
        return 1;
    }

    for (i = 0; i <= N; i++) {
        x = i * dx;
        fprintf(fp, "%lf %lf\n", x, u_old[i]);
    }

    fclose(fp);

    printf("t=0 の分布:      heat_t0.dat\n");
    printf("t=t_end の分布:  heat_tend.dat\n");
    printf("r = dt/dx^2 = %g (0 <= r <= 0.5 で安定)\n", R);
    printf("N = %d, dx = %g, dt = %g, NSTEP = %d\n", N, dx, dt, NSTEP);

    {
        const char *script_name = "kadai4.gp";
        FILE *gp = fopen(script_name, "w");
        if (!gp) {
            printf("gnuplot スクリプト %s を作成できませんでした。\n", script_name);
            return 1;
        }

        fprintf(gp, "set xlabel 'dimensionless distance x'\n");
        fprintf(gp, "set ylabel 'dimensionless temperature u(x,t)'\n");
        fprintf(gp, "set grid\n");

        fprintf(gp, "plot \\\n");
        fprintf(gp, "  'heat_t0.dat' using 1:2 with lines lw 2 lc rgb 'black' title 't=0'");

        int output_interval = 1000;
        int first_mid = output_interval;
        for (int step = first_mid; step <= NSTEP; step += output_interval) {
            double t = step * dt;
            fprintf(gp,
                    ", \\\n  'heat_t%04d.dat' using 1:2 with lines lw 1 lc rgb 'blue' title 't=%g'",
                    step, t);
        }

        fprintf(gp,
                ", \\\n  'heat_tend.dat' using 1:2 with lines lw 2 lc rgb 'red' title 't=end'\n");
        fclose(gp);

        char cmd[256];
    #ifdef _WIN32
        snprintf(cmd, sizeof(cmd), "cmd /c gnuplot -persist \"%s\"", script_name);
    #else
        snprintf(cmd, sizeof(cmd), "gnuplot -persist \"%s\"", script_name);
    #endif

        int ret = system(cmd);
        if (ret != 0) {
            printf("gnuplot 実行に失敗しました (return=%d)。\n", ret);
            printf("コマンド:\n%s\n", cmd);
        } else {
            printf("グラフを gnuplot で表示しました。\n");
        }
    }

    return 0;
}
