#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define LX 20.0
#define LY 4.0

#define NX 160
#define NY 32

#define RE   4.0

#define DT   5.0e-4 
#define TMAX 10.0 

#define MAX_SOR_IT 3000
#define SOR_EPS    1.0e-4
#define OMEGA_SOR  1.5

#define STEADY_EPS 1.0e-5

static double psi     [NX+1][NY+1];
static double omega   [NX+1][NY+1]; 
static double psi_old [NX+1][NY+1]; 
static int    is_solid[NX+1][NY+1];

double dx, dy;        
double beta, beta2;   

void setup_solid_region(void);
void set_initial_conditions(void);
void apply_bc(double psi[][NY+1], double omega[][NY+1]);
void update_vorticity(double psi[][NY+1], double omega[][NY+1], double dt);
void solve_poisson_sor(double psi[][NY+1], double omega[][NY+1]);
double max_abs_diff(double a[][NY+1], double b[][NY+1]);
void print_psi_range(const char *msg);

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    dx    = LX / NX;
    dy    = LY / NY;
    beta  = dx / dy;
    beta2 = beta * beta;

    printf("2D 非定常・非圧縮・粘性流 (バックステップ)\n");
    printf("NX=%d, NY=%d, dx=%g, dy=%g, Re=%g\n", NX, NY, dx, dy, RE);

    setup_solid_region();

    set_initial_conditions();
    print_psi_range("初期状態");

    double t = 0.0;
    int    step = 0;
    int    reached_steady = 0;

    while (t < TMAX) {
        step++;
        t += DT;

        for (int i = 0; i <= NX; i++) {
            for (int j = 0; j <= NY; j++) {
                psi_old[i][j] = psi[i][j];
            }
        }

        update_vorticity(psi, omega, DT);

        solve_poisson_sor(psi, omega);

        apply_bc(psi, omega);

        if (step % 200 == 0) {
            print_psi_range("time step");
        }

        double diff = max_abs_diff(psi, psi_old);
        if (step % 500 == 0) {
            printf("step=%6d, t=%8.4f, max|Δψ|=%e\n", step, t, diff);
        }

        if (!isnan(diff) && diff < STEADY_EPS) {
            printf("定常に到達したとみなします: step=%d, t=%g\n", step, t);
            reached_steady = 1;
            break;
        }

        if (fabs(psi[NX/2][NY/2]) > 1.0e6) {
            printf("ψ が非常に大きくなったため、発散と判断して終了します。\n");
            break;
        }
    }

    if (!reached_steady) {
        printf("TMAX=%g まで計算しましたが、完全な定常には達していない可能性があります。\n", TMAX);
    }

    FILE *fp = fopen("field.dat", "w");
    if (!fp) {
        printf("field.dat を開けませんでした。\n");
        return 1;
    }
    for (int j = 0; j <= NY; j++) {
        double y = j * dy;
        for (int i = 0; i <= NX; i++) {
            double x = i * dx;
            int solid = is_solid[i][j];
            fprintf(fp, "%g %g %g %g %d\n", x, y, psi[i][j], omega[i][j], solid);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    const char *script_name = "kadai8.gp";
    FILE *gp = fopen(script_name, "w");
    if (!gp) {
        printf("%s を作成できませんでした。\n", script_name);
        return 1;
    }
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "set view map\n");
    fprintf(gp, "set contour base\n");
    fprintf(gp, "unset surface\n");
    fprintf(gp, "set cntrparam levels discrete 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0\n");
    fprintf(gp, "set object 1 rect from 0,0 to 4,1 fc rgb 'gray' fs solid 0.4 border lc rgb 'black'\n");
    fprintf(gp, "splot 'field.dat' using 1:2:3 with lines\n");
    fclose(gp);
    char cmd[256];
#ifdef _WIN32
    snprintf(cmd, sizeof(cmd), "cmd /c gnuplot -persist \"%s\"", script_name);
#else
    snprintf(cmd, sizeof(cmd), "gnuplot -persist \"%s\"", script_name);
#endif
    int ret = system(cmd);
    if (ret != 0) {
        printf("kadai8 の gnuplot 実行に失敗しました (return=%d)。\n", ret);
    } else {
        printf("kadai8 の結果を gnuplot で表示しました。\n");
    }

    return 0;
}

void setup_solid_region(void)
{
    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            is_solid[i][j] = 0;
        }
    }


    double H     = 1.0;      
    double Lstep = 4.0; 

    int i_max = (int)round(Lstep / dx);
    int j_max = (int)round(H     / dy);

    for (int i = 0; i <= i_max; i++) {
        for (int j = 0; j <= j_max; j++) {
            is_solid[i][j] = 1; 
        }
    }
}

void set_initial_conditions(void)
{
    double U0 = 1.0; 

    for (int i = 0; i <= NX; i++) {
        double x = i * dx;
        (void)x; 

        for (int j = 0; j <= NY; j++) {
            double y = j * dy;

            if (is_solid[i][j]) {
                psi[i][j]   = 0.0;
                omega[i][j] = 0.0;
            } else {
                psi[i][j]   = U0 * y;
                omega[i][j] = 0.0;
            }
        }
    }

    apply_bc(psi, omega);
}

void apply_bc(double psi[][NY+1], double omega[][NY+1])
{
    double U0 = 1.0; 

    for (int i = 0; i <= NX; i++) {
        if (is_solid[i][0]) {
            psi[i][0]   = 0.0;
            omega[i][0] = 0.0;
        } else {
            psi[i][0] = 0.0;  
            if (!is_solid[i][1]) {
                omega[i][0] = -2.0*(psi[i][1] - psi[i][0])/(dy*dy);
            } else {
                omega[i][0] = 0.0;
            }
        }
    }

    for (int i = 0; i <= NX; i++) {
        if (is_solid[i][NY]) {
            psi[i][NY]   = 0.0;
            omega[i][NY] = 0.0;
        } else {
            psi[i][NY] = U0 * LY;

            if (!is_solid[i][NY-1]) {
                omega[i][NY] =
                    -2.0*(psi[i][NY-1] - psi[i][NY])/(dy*dy);
            } else {
                omega[i][NY] = 0.0;
            }
        }
    }

    for (int j = 0; j <= NY; j++) {

        if (is_solid[0][j]) {
            psi[0][j]   = 0.0;
            omega[0][j] = 0.0;
        } else {
            psi[0][j]   = psi[1][j];
            omega[0][j] = omega[1][j];
        }
    }

    for (int j = 0; j <= NY; j++) {
        if (is_solid[NX][j]) {
            psi[NX][j]   = 0.0;
            omega[NX][j] = 0.0;
        } else {
            psi[NX][j]   = psi[NX-1][j];
            omega[NX][j] = omega[NX-1][j];
        }
    }

    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            if (is_solid[i][j]) {
                psi[i][j]   = 0.0;
                omega[i][j] = 0.0;
            }
        }
    }
}


void update_vorticity(double psi[][NY+1], double omega[][NY+1], double dt)
{
    static double omega_new[NX+1][NY+1];

    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            omega_new[i][j] = omega[i][j];
        }
    }

    for (int i = 1; i < NX; i++) {
        for (int j = 1; j < NY; j++) {

            if (is_solid[i][j]) continue; 

            double dpsidx = (psi[i+1][j] - psi[i-1][j]) / (2.0*dx);
            double dpsidy = (psi[i][j+1] - psi[i][j-1]) / (2.0*dy);
            double dwdx   = (omega[i+1][j] - omega[i-1][j]) / (2.0*dx);
            double dwdy   = (omega[i][j+1] - omega[i][j-1]) / (2.0*dy);

            double J = dpsidx * dwdy - dpsidy * dwdx;

            double lap =
                (omega[i+1][j] - 2.0*omega[i][j] + omega[i-1][j]) / (dx*dx)
              + (omega[i][j+1] - 2.0*omega[i][j] + omega[i][j-1]) / (dy*dy);

            double rhs = J + (1.0/RE) * lap;

            omega_new[i][j] = omega[i][j] + dt * rhs;
        }
    }

    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            omega[i][j] = omega_new[i][j];
        }
    }
}


void solve_poisson_sor(double psi[][NY+1], double omega[][NY+1])
{
    for (int it = 0; it < MAX_SOR_IT; it++) {

        double maxdiff = 0.0;

        for (int i = 1; i < NX; i++) {
            for (int j = 1; j < NY; j++) {

                if (is_solid[i][j]) continue;  

                double old = psi[i][j];


                double term =
                    psi[i+1][j] + psi[i-1][j]
                  + psi[i][j+1] + psi[i][j-1];

                double rhs = -omega[i][j] * dx * dx; 

                double psi_new = (term - rhs) / 4.0;

                psi_new = old + OMEGA_SOR * (psi_new - old);

                double diff = fabs(psi_new - old);
                if (!isnan(diff) && diff > maxdiff) maxdiff = diff;

                psi[i][j] = psi_new;
            }
        }

        if (maxdiff < SOR_EPS) {
            break;
        }
    }
}

double max_abs_diff(double a[][NY+1], double b[][NY+1])
{
    double m = 0.0;
    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            double d = fabs(a[i][j] - b[i][j]);
            if (isnan(d)) continue;
            if (d > m) m = d;
        }
    }
    return m;
}

void print_psi_range(const char *msg)
{
    double psi_min = DBL_MAX;
    double psi_max = -DBL_MAX;
    int    nan_count = 0;

    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            if (isnan(psi[i][j])) {
                nan_count++;
                continue;
            }
            if (psi[i][j] < psi_min) psi_min = psi[i][j];
            if (psi[i][j] > psi_max) psi_max = psi[i][j];
        }
    }

    printf("psi range (%s): [%g, %g], NaN count = %d\n",
           msg, psi_min, psi_max, nan_count);
}
