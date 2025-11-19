#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define DT    0.01
#define TMAX  5.0
#define STEPS ((int)(TMAX/DT))

#define N_R      5
#define N_THETA  8
#define N_INIT   (N_R * N_THETA)

#define R_MAX  6.0
#define PI 3.14159265358979323846

void rhs(int sys, double x, double y, double *dxdt, double *dydt)
{
    double r = sqrt(x * x + y * y);
    double r2 = r * r;

    if (sys == 1) {
        /* 
        (1)
        dx/dt = y + x/√(x^2+y^2) * (1 - (x^2+y^2))
        dy/dt = -x + y/√(x^2+y^2) * (1 - (x^2+y^2)) 
        */
        double coef = (1.0 - r2) / r;
        *dxdt = y + x * coef;
        *dydt = -x + y * coef;

    } else if (sys == 2) {
        /*
        (2)
        dx/dt = -y + x(x^2+ y^2 -1)
        dy/dt =  x + y(x^2+ y^2 -1) 
        */
        double coef = (r2 - 1.0);
        *dxdt = -y + x * coef;
        *dydt =  x + y * coef;

    } else if (sys == 3) {
        /* 
        (3)
        dx/dt =  y + x√(x^2+y^2) (x^2 + y^2 -1)^2
        dy/dt = -x + y√(x^2+y^2) (x^2 + y^2 -1)^2 
        */
        double coef = r * (r2 - 1.0) * (r2 - 1.0);
        *dxdt =  y + x * coef;
        *dydt = -x + y * coef;
    } else {
        *dxdt = 0.0;
        *dydt = 0.0;
    }
}

void rk4_step(int sys, double *x, double *y, double dt)
{
    double k1x, k1y;
    double k2x, k2y;
    double k3x, k3y;
    double k4x, k4y;

    double xx, yy;

    rhs(sys, *x, *y, &k1x, &k1y);

    xx = *x + 0.5 * dt * k1x;
    yy = *y + 0.5 * dt * k1y;
    rhs(sys, xx, yy, &k2x, &k2y);

    xx = *x + 0.5 * dt * k2x;
    yy = *y + 0.5 * dt * k2y;
    rhs(sys, xx, yy, &k3x, &k3y);

    xx = *x + dt * k3x;
    yy = *y + dt * k3y;
    rhs(sys, xx, yy, &k4x, &k4y);

    *x += dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    *y += dt / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
}

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    double r_init[N_R] = {
        0.5, 0.9, 1.0, 1.1, 1.5
    };

    double theta_init[N_THETA] = {
        0.0,
        2.0 * PI / N_THETA * 1,
        2.0 * PI / N_THETA * 2,
        2.0 * PI / N_THETA * 3,
        2.0 * PI / N_THETA * 4,
        2.0 * PI / N_THETA * 5,
        2.0 * PI / N_THETA * 6,
        2.0 * PI / N_THETA * 7
    };

    for (int sys = 1; sys <= 3; sys++) {

        char fname[64];
        snprintf(fname, sizeof(fname), "kadai3-2_sys%d.dat", sys);
        FILE *fp = fopen(fname, "w");
        if (fp == NULL) {
            printf("ファイル %s を開けませんでした。\n", fname);
            continue;
        }

        double x[N_INIT];
        double y[N_INIT];
        double r_current[N_INIT];
        double theta_current[N_INIT];

        int n = 0;
        for (int i_r = 0; i_r < N_R; i_r++) {
            for (int i_th = 0; i_th < N_THETA; i_th++) {
                double r0 = r_init[i_r];
                double th0 = theta_init[i_th];
                x[n] = r0 * cos(th0);
                y[n] = r0 * sin(th0);
                r_current[n] = r0;
                theta_current[n] = th0;
                n++;
            }
        }

        fprintf(fp, "# t");
        for (int i = 0; i < N_INIT; i++) {
            fprintf(fp, "  x%d  y%d", i + 1, i + 1);
        }
        fprintf(fp, "\n");

        double t = 0.0;

        for (int step = 0; step < STEPS; step++) {

            fprintf(fp, "%lf", t);

            for (int i = 0; i < N_INIT; i++) {

                if (r_current[i] > R_MAX) {
                    fprintf(fp, " %lf %lf", x[i], y[i]);
                    continue;
                }

                fprintf(fp, " %lf %lf", x[i], y[i]);

                rk4_step(sys, &x[i], &y[i], DT);

                r_current[i] = sqrt(x[i] * x[i] + y[i] * y[i]);
                theta_current[i] = atan2(y[i], x[i]);
            }

            fprintf(fp, "\n");
            t += DT;
        }

        fclose(fp);
        printf("(%d) の位相平面データを %s に出力しました。\n", sys, fname);
    }

    printf("\n出力された kadai3-2_sys*.dat を使って位相図を描写\n");

    for (int sys = 1; sys <= 3; sys++) {

        char script_name[64];
        snprintf(script_name, sizeof(script_name), "kadai3-2_sys%d.gp", sys);

        FILE *gp = fopen(script_name, "w");
        if (!gp) {
            printf("((%d) の gnuplot スクリプト %s を作成できませんでした。\n", sys, script_name);
            continue;
        }

        if (sys == 1) {
            fprintf(gp,
                "set title '(1) dx/dt = y + x/sqrt(x**2+y**2)*(1-(x**2+y**2)), "
                "dy/dt = -x + y/sqrt(x**2+y**2)*(1-(x**2+y**2))'\n");
        } else if (sys == 2) {
            fprintf(gp,
                "set title '(2) dx/dt = -y + x*(x**2+y**2-1), "
                "dy/dt = x + y*(x**2+y**2-1)'\n");
        } else if (sys == 3) {
            fprintf(gp,
                "set title '(3) dx/dt = y + x*sqrt(x**2+y**2)*(x**2+y**2-1)**2, "
                "dy/dt = -x + y*sqrt(x**2+y**2)*(x**2+y**2-1)**2'\n");
        }

        fprintf(gp, "set xrange[-2:2]\n");
        fprintf(gp, "set yrange[-2:2]\n");
        fprintf(gp, "plot \\\n");

        int first = 1;
        for (int i_r = 0; i_r < N_R; i_r++) {
            int color_index = i_r + 1;

            for (int i_th = 0; i_th < N_THETA; i_th++) {
                int idx   = i_r * N_THETA + i_th;
                int col_x = 2 + 2 * idx;
                int col_y = 3 + 2 * idx;

                if (!first) {
                    fprintf(gp, ", \\\n");
                }
                first = 0;

                if (i_th == 0) {
                    fprintf(
                        gp,
                        "'kadai3-2_sys%d.dat' using %d:%d with lines linecolor %d "
                        "title 'r=%.2f'",
                        sys, col_x, col_y,
                        color_index,
                        r_init[i_r]
                    );
                } else {
                    fprintf(
                        gp,
                        "'kadai3-2_sys%d.dat' using %d:%d with lines linecolor %d notitle",
                        sys, col_x, col_y,
                        color_index
                    );
                }
            }
        }

        fprintf(gp, "\n");
        fclose(gp);

        char cmd[256];
#ifdef _WIN32
        snprintf(cmd, sizeof(cmd), "cmd /c gnuplot -persist \"%s\"", script_name);
#else
        snprintf(cmd, sizeof(cmd), "gnuplot -persist \"%s\"", script_name);
#endif

        int ret = system(cmd);
        if (ret != 0) {
            printf("(%d) の gnuplot 実行に失敗しました (return=%d)。\n", sys, ret);
            printf("実行しようとしたコマンド:\n%s\n", cmd);
        } else {
            printf("(%d) の位相図を gnuplot で表示しました。\n", sys);
        }
    }

    return 0;
}
