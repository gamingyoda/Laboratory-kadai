#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define DT    0.001
#define STEPS 200000

typedef struct {
    double a1, a2, a3;
    double c1, c2, c3;
} Params;

int set_params(int type, Params *p)
{
    switch (type) {
    case 1: /* スパイラルアトラクタ */
        p->a1 = 0.19;
        p->a2 = -1.9644;
        p->a3 = 6.31222;
        p->c1 = -2.77;
        p->c2 = -2.61;
        p->c3 = -13.91;
        return 1;
    case 2: /* ダブルスクロールアトラクタ */
        p->a1 = 0.19;
        p->a2 = -1.964;
        p->a3 = 6.312;
        p->c1 = -2.69;
        p->c2 = -2.383;
        p->c3 = -14.049;
        return 1;
    case 3: /* ダブルスクリューアトラクタ */
        p->a1 = 0.604;
        p->a2 = -6.3047;
        p->a3 = 3.1264;
        p->c1 = -1.284;
        p->c2 = 5.2757;
        p->c3 = 7.1022;
        return 1;
    case 4: /* スパローアトラクタ */
        p->a1 = -1.164;
        p->a2 = -0.827824;
        p->a3 = -1.306011;
        p->c1 = 1.2482;
        p->c2 = -8.532633;
        p->c3 = 16.98886;
        return 1;
    default:
        return 0;
    }
}

void rhs(double x, double y, double z, const Params *p,
         double *dxdt, double *dydt, double *dzdt)
{
    double h = 0.5 * (fabs(x - 1.0) - fabs(x + 1.0));

    *dxdt = p->c1 * x + y + p->c1 * h;
    *dydt = p->c2 * x + z + p->c2 * h;
    *dzdt = (p->c3 + p->a3) * x + p->a2 * y + p->a1 * z + p->c3 * h;
}

void rk4_step(double *x, double *y, double *z, double dt, const Params *p)
{
    double k1x, k1y, k1z;
    double k2x, k2y, k2z;
    double k3x, k3y, k3z;
    double k4x, k4y, k4z;
    double xx, yy, zz;

    rhs(*x, *y, *z, p, &k1x, &k1y, &k1z);

    xx = *x + 0.5 * dt * k1x;
    yy = *y + 0.5 * dt * k1y;
    zz = *z + 0.5 * dt * k1z;
    rhs(xx, yy, zz, p, &k2x, &k2y, &k2z);

    xx = *x + 0.5 * dt * k2x;
    yy = *y + 0.5 * dt * k2y;
    zz = *z + 0.5 * dt * k2z;
    rhs(xx, yy, zz, p, &k3x, &k3y, &k3z);

    xx = *x + dt * k3x;
    yy = *y + dt * k3y;
    zz = *z + dt * k3z;
    rhs(xx, yy, zz, p, &k4x, &k4y, &k4z);

    *x += dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    *y += dt / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
    *z += dt / 6.0 * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);
}

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    for (int type = 1; type <= 4; type++) {

        Params p;
        if (!set_params(type, &p)) {
            continue;
        }

        const char *fname = NULL;
        const char *title = NULL;
        double xrange_min, xrange_max, yrange_min, yrange_max;

        switch (type) {
        case 1:
            fname = "spiral.dat";
            title = "(1) Spiral attractor";
            xrange_min = -6.0; xrange_max = 6.0;
            yrange_min = -3.0; yrange_max = 3.0;
            break;
        case 2:
            fname = "dscroll.dat";
            title = "(2) Double scroll attractor";
            xrange_min = -6.0; xrange_max = 6.0;
            yrange_min = -3.0; yrange_max = 3.0;
            break;
        case 3:
            fname = "dscrew.dat";
            title = "(3) Double screw attractor";
            xrange_min = -70.0; xrange_max = 70.0;
            yrange_min = -70.0; yrange_max = 70.0;
            break;
        case 4:
            fname = "sparrow.dat";
            title = "(4) Sparrow attractor";
            xrange_min = -2.0; xrange_max = 2.0;
            yrange_min = -2.0; yrange_max = 2.0;
            break;
        }

        FILE *fp = fopen(fname, "w");
        if (!fp) {
            printf("ファイル %s を開けませんでした。\n", fname);
            continue;
        }

        double x = 0.1;
        double y = 0.0;
        double z = 0.0;
        int skip = 2000;

        fprintf(fp, "# type = %d\n", type);
        fprintf(fp, "# x0 = %g, y0 = %g, z0 = %g\n", x, y, z);
        fprintf(fp, "# projection: ");
        if (type == 3) {
            fprintf(fp, "horizontal = y, vertical = z\n");
            fprintf(fp, "# y   z\n");
        } else if (type == 4) {
            fprintf(fp, "horizontal = y, vertical = x\n");
            fprintf(fp, "# y   x\n");
        } else {
            fprintf(fp, "horizontal = z, vertical = x\n");
            fprintf(fp, "# z   x\n");
        }

        for (int i = 0; i < STEPS; i++) {
            rk4_step(&x, &y, &z, DT, &p);
            if (i < skip) continue;

            if (type == 3) {
                fprintf(fp, "%lf %lf\n", y, z);
            } else if (type == 4) {
                fprintf(fp, "%lf %lf\n", y, x);
            } else {
                fprintf(fp, "%lf %lf\n", z, x);
            }
        }
        fclose(fp);
        printf("(%d) のデータを %s に出力しました。\n", type, fname);

        char script_name[64];
        snprintf(script_name, sizeof(script_name), "kadai3-3_type%d.gp", type);

        FILE *gp = fopen(script_name, "w");
        if (!gp) {
            printf("gnuplot スクリプト %s を作成できませんでした。\n", script_name);
            continue;
        }

        fprintf(gp, "set title '%s'\n", title);
        if (type == 3) {
            fprintf(gp, "set xlabel 'y'\n");
            fprintf(gp, "set ylabel 'z'\n");
        } else if (type == 4) {
            fprintf(gp, "set xlabel 'y'\n");
            fprintf(gp, "set ylabel 'x'\n");
        } else {
            fprintf(gp, "set xlabel 'z'\n");
            fprintf(gp, "set ylabel 'x'\n");
        }
        fprintf(gp, "set xrange[%g:%g]\n", xrange_min, xrange_max);
        fprintf(gp, "set yrange[%g:%g]\n", yrange_min, yrange_max);
        fprintf(gp, "plot '%s' using 1:2 with lines notitle\n", fname);
        fclose(gp);

        char cmd[256];
#ifdef _WIN32
        snprintf(cmd, sizeof(cmd), "cmd /c gnuplot -persist \"%s\"", script_name);
#else
        snprintf(cmd, sizeof(cmd), "gnuplot -persist \"%s\"", script_name);
#endif
        int ret = system(cmd);
        if (ret != 0) {
            printf("(%d) の gnuplot 実行に失敗しました (return=%d)。\n", type, ret);
            printf("コマンド:\n%s\n", cmd);
        } else {
            printf("(%d) のグラフを gnuplot で表示しました。\n", type);
        }
    }

    return 0;
}
