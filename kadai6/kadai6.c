#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define NPSI   9
#define NTHETA 720
#define NSTEPS 8000
#define HSTEP  0.01
#define PI 3.14159265358979323846
#define VIEW_MIN -5.0
#define VIEW_MAX  5.0

void velocity(double xi, double eta,
              double xi0, double eta0, double r,
              double U, double alpha, double Gamma,
              double cosA, double sinA,
              double *ux, double *uy)
{
    double xr = xi  - xi0;
    double yr = eta - eta0;

    double rho2 = xr * xr + yr * yr;

    if (rho2 < 1.0e-8) {
        *ux = 0.0;
        *uy = 0.0;
        return;
    }

    double inv_re = xr / rho2;
    double inv_im = -yr / rho2;

    double inv2_re = inv_re * inv_re - inv_im * inv_im;
    double inv2_im = 2.0 * inv_re * inv_im;

    double term1_re = U * cosA;
    double term1_im = -U * sinA;

    double tmp_re =  cosA * inv2_re - sinA * inv2_im;
    double tmp_im =  cosA * inv2_im + sinA * inv2_re;
    double term2_re = -U * r * r * tmp_re;
    double term2_im = -U * r * r * tmp_im;

    double coefG = Gamma / (2.0 * PI);
    double term3_re =  coefG * inv_im;
    double term3_im =  coefG * (-inv_re);

    double v_re = term1_re + term2_re + term3_re;
    double v_im = term1_im + term2_im + term3_im;

    *ux = v_re;
    *uy = -v_im;
}

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    /* ====== ここでパラメータを指定 ====== */
    double a         = 1.0;     
    double xi0       = -0.10;   
    double eta0      =  0.25; 
    double U         = 1.0;    
    double alpha_deg = 15.0;   

    double alpha = alpha_deg * PI / 180.0;
    double cosA  = cos(alpha);
    double sinA  = sin(alpha);

    double r = sqrt((a - xi0)*(a - xi0) + eta0*eta0);
    printf("円の半径 r = %g\n", r);

    double beta = atan2(-eta0, a - xi0);
    double Gamma = -4.0 * PI * r * U * sin(alpha + beta);

    printf("beta (rad) = %g\n", beta);
    printf("Gamma      = %g\n", Gamma);
    printf("Gamma/(2π a U) = %g\n", Gamma / (2.0 * PI * a * U));

    FILE *fp_airfoil = fopen("airfoil.dat", "w");
    if (!fp_airfoil) {
        printf("airfoil.dat を開けませんでした\n");
        return 1;
    }

    fprintf(fp_airfoil, "# Airfoil (Joukowski)\n");
    fprintf(fp_airfoil, "# a=%g, xi0=%g, eta0=%g, r=%g\n", a, xi0, eta0, r);
    fprintf(fp_airfoil, "# Gamma=%g, alpha=%g[deg]\n", Gamma, alpha_deg);

    for (int n = 0; n <= NTHETA; n++) {
        double theta = 2.0 * PI * n / NTHETA;

        double xi_circle  = xi0 + r * cos(theta);
        double eta_circle = eta0 + r * sin(theta);

        double denom = xi_circle*xi_circle + eta_circle*eta_circle;
        if (denom < 1e-12) continue;

        double x = xi_circle  * (1.0 + (a*a)/denom);
        double y = eta_circle * (1.0 - (a*a)/denom);

        fprintf(fp_airfoil, "%lf %lf\n", x, y);
    }
    fclose(fp_airfoil);

    double L = 8.0 * r;
    double S = 3.0 * r;

    char fname[64];

    for (int k = 0; k < NPSI; k++) {
        double s = ((double)k/(NPSI-1) - 0.5) * 2.0 * S;

        double xi  = xi0 - L*cosA + s*(-sinA);
        double eta = eta0 - L*sinA + s*( cosA);

        snprintf(fname, sizeof(fname), "stream_%02d.dat", k);
        FILE *fp = fopen(fname, "w");
        if (!fp) {
            printf("%s を開けませんでした\n", fname);
            return 1;
        }

        bool entered_view = false;

        for (int n = 0; n < NSTEPS; n++) {

            double xr   = xi - xi0;
            double yr   = eta - eta0;
            double rho2 = xr*xr + yr*yr;
            if (rho2 < 1.0e-8) break;
            double denom = xi*xi + eta*eta;
            if (denom < 1e-12) break;

            double x = xi  * (1.0 + (a*a)/denom);
            double y = eta * (1.0 - (a*a)/denom);

            double nx = x / a;
            double ny = y / a;
            int inside_view = (nx >= VIEW_MIN && nx <= VIEW_MAX &&
                               ny >= VIEW_MIN && ny <= VIEW_MAX);

            if (inside_view) {
                entered_view = true;
                fprintf(fp, "%lf %lf\n", x, y);
            } else if (entered_view) {
                break;
            }

            double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y;
            double tx,ty;

            velocity(xi,eta,xi0,eta0,r,U,alpha,Gamma,cosA,sinA,&k1x,&k1y);

            tx = xi + 0.5*HSTEP*k1x;
            ty = eta+ 0.5*HSTEP*k1y;
            velocity(tx,ty,xi0,eta0,r,U,alpha,Gamma,cosA,sinA,&k2x,&k2y);

            tx = xi + 0.5*HSTEP*k2x;
            ty = eta+ 0.5*HSTEP*k2y;
            velocity(tx,ty,xi0,eta0,r,U,alpha,Gamma,cosA,sinA,&k3x,&k3y);

            tx = xi + HSTEP*k3x;
            ty = eta+ HSTEP*k3y;
            velocity(tx,ty,xi0,eta0,r,U,alpha,Gamma,cosA,sinA,&k4x,&k4y);

            xi  += (HSTEP/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
            eta += (HSTEP/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
        }

        fclose(fp);
        printf("%s を出力\n", fname);
    }

    {
        char script_name[64];
        snprintf(script_name, sizeof(script_name), "kadai6.gp");

        double xi_over_a   = xi0 / a;
        double eta_over_a  = eta0 / a;
        double r_over_a    = r   / a;
        double G_over_2pUa = Gamma / (2.0 * PI * a * U);

        FILE *gp = fopen(script_name, "w");
        if (!gp) {
            printf("gnuplot スクリプト %s を作成できませんでした。\n", script_name);
        } else {
            fprintf(gp, "set xlabel 'x/a'\n");
            fprintf(gp, "set ylabel 'iy/a'\n");
            fprintf(gp, "set grid\n");
            fprintf(gp, "set xrange [-5:5]\n");
            fprintf(gp, "set yrange [-5:5]\n");
            fprintf(gp, "set size ratio -1\n");

            fprintf(gp,
                    "set label 1 sprintf('{/Symbol x}/a = %.3f, {/Symbol h}/a = %.3f, r/a = %.3f', %.15g, %.15g, %.15g) at graph 0.98,0.10 right\n",
                    xi_over_a, eta_over_a, r_over_a,
                    xi_over_a, eta_over_a, r_over_a);
            fprintf(gp,
                    "set label 2 sprintf('{/Symbol G}/(2{/Symbol p} a U_0) = %.3f, {/Symbol a} = %.1f deg', %.15g, %.15g) at graph 0.98,0.05 right\n",
                    G_over_2pUa, alpha_deg,
                    G_over_2pUa, alpha_deg);

            fprintf(gp, "plot \\\n");
            fprintf(gp,
                    "  'airfoil.dat' using ($1/%g):($2/%g) with lines lt -1 notitle", a, a);
            for (int k = 0; k < NPSI; k++) {
                fprintf(gp,
                        ", \\\n  'stream_%02d.dat' using ($1/%g):($2/%g) with lines lw 1 notitle",
                        k, a, a);
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
                printf("kadai6 の gnuplot 実行に失敗しました (return=%d)。\n", ret);
                printf("実行しようとしたコマンド:\n%s\n", cmd);
            } else {
                printf("kadai6 の翼型と流線のグラフを gnuplot で表示しました。\n");
            }
        }
    }

    return 0;
}
