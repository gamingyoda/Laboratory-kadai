#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define NTHETA 720
#define PI 3.14159265358979323846

typedef struct {
    double a;
    double xi0;
    double eta0;
    const char *label;
} CaseParam;

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    CaseParam cases[9] = {
        /* (1) a = 1, xi0 = 0,   eta0 = 0.1 */
        {1.0,  0.0,  0.1, "1"},
        /* (2) a = 1, xi0 = 0,   eta0 = 0.2 */
        {1.0,  0.0,  0.2, "2"},
        /* (3) a = 1, xi0 = 0,   eta0 = 0.3 */
        {1.0,  0.0,  0.3, "3"},
        /* (4) a = 1, xi0 = -0.1, eta0 = 0   */
        {1.0, -0.1,  0.0, "4"},
        /* (5) a = 1, xi0 = -0.2, eta0 = 0   */
        {1.0, -0.2,  0.0, "5"},
        /* (6) a = 1, xi0 = -0.3, eta0 = 0   */
        {1.0, -0.3,  0.0, "6"},
        /* (7) a = 1, xi0 = -0.1, eta0 = 0.1 */
        {1.0, -0.1,  0.1, "7"},
        /* (8) a = 1, xi0 = -0.1, eta0 = 0.2 */
        {1.0, -0.1,  0.2, "8"},
        /* (9) a = 1, xi0 = -0.1, eta0 = 0.3 */
        {1.0, -0.1,  0.3, "9"}
    };

    for (int c = 0; c < 9; c++) {
        double a    = cases[c].a;
        double xi0  = cases[c].xi0;
        double eta0 = cases[c].eta0;

        double r = sqrt((a - xi0) * (a - xi0) + eta0 * eta0);

        char fname[64];
        snprintf(fname, sizeof(fname), "joukowski_%s.dat", cases[c].label);

        FILE *fp = fopen(fname, "w");
        if (fp == NULL) {
            printf("ファイル %s を開けませんでした。\n", fname);
            continue;
        }

        char fname_circle[64];
        snprintf(fname_circle, sizeof(fname_circle), "circle_%s.dat", cases[c].label);
        FILE *fp_circle = fopen(fname_circle, "w");
        if (fp_circle == NULL) {
            printf("ファイル %s を開けませんでした。\n", fname_circle);
            fclose(fp);
            continue;
        }

        fprintf(fp, "# Joukowski airfoil case %s\n", cases[c].label);
        fprintf(fp, "# a = %g, xi0 = %g, eta0 = %g, r = %g\n",
                a, xi0, eta0, r);
        fprintf(fp, "# x   y   x/a   y/a\n");

        fprintf(fp_circle, "# circle in zeta-plane case %s\n", cases[c].label);
        fprintf(fp_circle, "# xi0 = %g, eta0 = %g, r = %g\n", xi0, eta0, r);
        fprintf(fp_circle, "# xi   eta   xi/a   eta/a\n");

        for (int k = 0; k <= NTHETA; k++) {
            double theta = 2.0 * PI * k / NTHETA;

            double xi  = xi0  + r * cos(theta);
            double eta = eta0 + r * sin(theta);

            double xi_nd  = xi  / a;
            double eta_nd = eta / a;
            fprintf(fp_circle, "%lf %lf %lf %lf\n", xi, eta, xi_nd, eta_nd);

            double denom = xi * xi + eta * eta;

            if (denom < 1.0e-12) {
                continue;
            }

            double x = xi  * (1.0 + (a * a) / denom);
            double y = eta * (1.0 - (a * a) / denom);

            double x_nd = x / a;
            double y_nd = y / a;

            fprintf(fp, "%lf %lf %lf %lf\n", x, y, x_nd, y_nd);
        }

        fclose(fp);
        fclose(fp_circle);

        double r_nd = r / a;
        printf("(%s): a=%g, xi0=%g, eta0=%g, r=%g, r/a=%g\n",
               cases[c].label, a, xi0, eta0, r, r_nd);
        printf("         データを %s, %s に出力しました。\n",
               fname, fname_circle);
    }

    printf("\n出力された .dat ファイルをプロットして円と翼型を確認してください。\n");
    printf("例 (gnuplot):\n");
    printf("  plot 'circle_1.dat' using 3:4 with lines, \\\n");
    printf("       'joukowski_1.dat' using 3:4 with lines\n");

    for (int c = 0; c < 9; c++) {
        char script_name[64];
        snprintf(script_name, sizeof(script_name), "kadai5_%s.gp", cases[c].label);

        FILE *gp = fopen(script_name, "w");
        if (!gp) {
            printf("gnuplot スクリプト %s を作成できませんでした。\n", script_name);
            continue;
        }

        double a    = cases[c].a;
        double xi0  = cases[c].xi0;
        double eta0 = cases[c].eta0;
        double r    = sqrt((a - xi0) * (a - xi0) + eta0 * eta0);
        double xi0_nd  = xi0 / a;
        double eta0_nd = eta0 / a;
        double r_nd    = r    / a;

        fprintf(gp, "set xlabel 'x, xi'\n");
        fprintf(gp, "set ylabel 'iy, ieta'\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange[-2.5:2.5]\n");
        fprintf(gp, "set yrange[-2.5:2.5]\n");
        fprintf(gp, "set title 'a=%g, xi0=%g, eta0=%g'\n", a, xi0, eta0);

        fprintf(gp,
                "set label 1 sprintf('x_c/a = %.3f', %g) at 0.6, 2.1\n",
                xi0_nd, xi0_nd);
        fprintf(gp,
                "set label 2 sprintf('y_c/a = %.3f', %g) at 0.6, 1.8\n",
                eta0_nd, eta0_nd);
        fprintf(gp,
                "set label 3 sprintf('r/a   = %.3f', %g) at 0.6, 1.5\n",
                r_nd, r_nd);

        fprintf(gp, "plot \\\n");
        fprintf(gp,
                "  'circle_%s.dat' using 3:4 with lines lw 2 lc rgb 'black' notitle, \\\n",
                cases[c].label);
        fprintf(gp,
                "  'joukowski_%s.dat' using 3:4 with lines lw 2 lc rgb 'blue'  notitle\n",
                cases[c].label);

        fclose(gp);

        char cmd[256];
    #ifdef _WIN32
        snprintf(cmd, sizeof(cmd), "cmd /c gnuplot -persist \"%s\"", script_name);
    #else
        snprintf(cmd, sizeof(cmd), "gnuplot -persist \"%s\"", script_name);
    #endif

        int ret = system(cmd);
        if (ret != 0) {
            printf("gnuplot 実行に失敗しました (script=%s, return=%d)。\n",
                   script_name, ret);
            printf("コマンド:\n%s\n", cmd);
        } else {
            printf("(%s) の円と翼型のグラフを gnuplot で表示しました。\n",
                   cases[c].label);
        }
    }

    return 0;
}
