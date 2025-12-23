/*  
1次元の衝撃波管問題(Sod問題)
1次元Euler方程式 + FDS（Roe-FDS）+ Runge–Kutta

γ=1.4
左：
ρ=1.0,u=0.0,p=1.0
右：
ρ=0.125,u=0.0,p=0.1
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define GAMMA 1.4

#define NX 1000
#define X0 0.0
#define X1 1.0
#define XDIFF 0.5
#define CFL 0.5
#define T_END 0.20

#define RHO_L  1.0
#define U_L    0.0
#define P_L    1.0

#define RHO_R  0.125
#define U_R    0.0
#define P_R    0.1

#define WL_INIT  ((Primitive){RHO_L, U_L, P_L})
#define WR_INIT  ((Primitive){RHO_R, U_R, P_R})


typedef struct {double rho, rho_u, E;} Conservative;
typedef struct {double rho, u, p;} Primitive;

static Primitive cons_to_prim(Conservative U)
{
    Primitive W;
    W.rho = U.rho;
    W.u   = U.rho_u / U.rho;
    W.p   = (GAMMA - 1.0) * (U.E - 0.5 * U.rho * W.u * W.u);
    return W;
}

static Conservative prim_to_cons(Primitive W)
{
    Conservative U;
    U.rho   = W.rho;
    U.rho_u = W.rho * W.u;
    U.E     = W.p / (GAMMA - 1.0) + 0.5 * W.rho * W.u * W.u;
    return U;
}

static Conservative Flux(Primitive W)
{
    Conservative F;
    F.rho   = W.rho * W.u;
    F.rho_u = W.rho * W.u * W.u + W.p;
    F.E     = W.u * ((W.p / (GAMMA - 1.0) + 0.5 * W.rho * W.u * W.u) + W.p);
    return F;
}

static double ROE_Flux(Primitive WL, Primitive WR, Conservative *Fout)
{
    Conservative FL = Flux(WL);
    Conservative FR = Flux(WR);
    Conservative UL = prim_to_cons(WL);
    Conservative UR = prim_to_cons(WR);

    double u_t = (sqrt(WL.rho)*WL.u + sqrt(WR.rho)*WR.u) / (sqrt(WL.rho)+sqrt(WR.rho));
    double H_t = (sqrt(WL.rho)*(UL.E + WL.p)/WL.rho + sqrt(WR.rho)*(UR.E + WR.p)/WR.rho)/(sqrt(WL.rho)+sqrt(WR.rho));
    double a_t = sqrt((GAMMA-1)*(H_t - 0.5 * u_t * u_t));
    
    double dr = WR.rho - WL.rho;
    double du = WR.u - WL.u;  
    double dp = WR.p - WL.p;

    double alpha2 = dr - dp/(a_t*a_t);
    double alpha1 = (dp - sqrt(WL.rho*WR.rho)*a_t*du)/(2*a_t*a_t);
    double alpha3 = (dp + sqrt(WL.rho*WR.rho)*a_t*du)/(2*a_t*a_t);

    double lamda1 = u_t - a_t;
    double lamda2 = u_t;
    double lamda3 = u_t + a_t;

    double delta = 0.1 * a_t;
    double a1 = (fabs(lamda1) >= delta) ? fabs(lamda1) : (lamda1*lamda1 + delta*delta)/(2.0*delta); 
    double a2 = (fabs(lamda2) >= delta) ? fabs(lamda2) : (lamda2*lamda2 + delta*delta)/(2.0*delta); 
    double a3 = (fabs(lamda3) >= delta) ? fabs(lamda3) : (lamda3*lamda3 + delta*delta)/(2.0*delta);

    double d0 = a1*alpha1*1.0 + a2*alpha2*1.0 + a3*alpha3*1.0;
    double d1 = a1*alpha1*(u_t - a_t) + a2*alpha2*u_t + a3*alpha3*(u_t + a_t);
    double d2 = a1*alpha1*(H_t - u_t*a_t) + a2*alpha2*0.5*u_t*u_t + a3*alpha3*(H_t + u_t*a_t);
    
    Fout->rho   = 0.5 * (FL.rho   + FR.rho  ) - 0.5 * d0;
    Fout->rho_u = 0.5 * (FL.rho_u + FR.rho_u) - 0.5 * d1;
    Fout->E     = 0.5 * (FL.E     + FR.E    ) - 0.5 * d2;

    return 0.0;
}

static void zero_slope(Conservative *U, int n_end)
{
    for (int i=0; i<2; i++) {
        U[i] = U[2];
        U[n_end-1-i] = U[n_end-3];
    }
}

static void calculete_rhs(Conservative *U, Conservative *rhs, double dx, int n_end)
{
    Conservative *F = malloc(sizeof(Conservative)*(n_end-1));
    for (int i=1; i<=n_end-3; i++) {
        Primitive WL = cons_to_prim(U[i]);
        Primitive WR = cons_to_prim(U[i+1]);
        ROE_Flux(WL, WR, &F[i]);
    }
    for (int i=2; i<=n_end-3; i++) {
        rhs[i].rho   = -(F[i].rho   - F[i-1].rho  ) / dx;
        rhs[i].rho_u = -(F[i].rho_u - F[i-1].rho_u) / dx;
        rhs[i].E     = -(F[i].E     - F[i-1].E    ) / dx;
    }
    free(F);
}

static double calculate_dt(const Conservative *U, double dx, int n_end)
{
    const double aL = sqrt(GAMMA * (P_L / RHO_L));
    const double aR = sqrt(GAMMA * (P_R / RHO_R));
    double smax = fmax(fabs(U_L) + aL, fabs(U_R) + aR);

    for (int i = 2; i <= n_end - 3; i++) { 
        Primitive W = cons_to_prim(U[i]);
        if (W.rho <= 0.0 || W.p <= 0.0) continue;

        const double a = sqrt(GAMMA * W.p / W.rho);
        const double s = fabs(W.u) + a;
        if (s > smax) smax = s;
    }

    return CFL * dx / smax;
}

static void write_and_plot(const Conservative *U, int n_end, int ng, double dx)
{
    FILE *fp = fopen("sod1d.dat", "w");
    if (!fp) { perror("sod1d.dat"); return; }

    fprintf(fp, "# x rho u p\n");
    for (int i=ng; i<n_end-ng; i++) {
        const double x = X0 + ((i-ng) + 0.5)*dx;
        Primitive W = cons_to_prim(U[i]);
        fprintf(fp, "%.10f %.10f %.10f %.10f\n", x, W.rho, W.u, W.p);
    }
    fclose(fp);

    const char *gpname = "sod1d-1.gp";
    FILE *gp = fopen(gpname, "w");
    if (!gp) { perror("sod1d-1.gp"); return; }

    fprintf(gp, "set grid\n");
    fprintf(gp, "set multiplot layout 3,1 title 'Sod 1D: Roe-FDS + TVD RK2 (t=%.3f, NX=%d)'\n", T_END, NX);

    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'rho'\n");
    fprintf(gp, "plot 'sod1d.dat' using 1:2 with lines title 'rho'\n");

    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'u'\n");
    fprintf(gp, "plot 'sod1d.dat' using 1:3 with lines title 'u'\n");

    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'p'\n");
    fprintf(gp, "plot 'sod1d.dat' using 1:4 with lines title 'p'\n");

    fprintf(gp, "unset multiplot\n");
    fclose(gp);

    char cmd[256];
#ifdef _WIN32
    snprintf(cmd, sizeof(cmd), "cmd /c gnuplot -persist \"%s\"", gpname);
#else
    snprintf(cmd, sizeof(cmd), "gnuplot -persist \"%s\"", gpname);
#endif
    int ret = system(cmd);
    if (ret != 0) {
        fprintf(stderr, "gnuplot 実行に失敗しました (return=%d)。gnuplot が入っているか確認してください。\n", ret);
    }
}

int main(void)
{
    const int ng = 2;
    const int n_end = NX + 2*ng;
    double dx = (X1 - X0)/NX;

    Conservative *U  = (Conservative*)malloc(sizeof(Conservative)*n_end);
    Conservative *U1 = (Conservative*)malloc(sizeof(Conservative)*n_end);
    Conservative *rhs= (Conservative*)malloc(sizeof(Conservative)*n_end);

    Primitive WL = WL_INIT;
    Primitive WR = WR_INIT;

    for (int i=0; i<n_end; i++) {
        double x = X0 + ( (i-ng) + 0.5 )*dx;
        Primitive W = (x < XDIFF) ? WL : WR;
        U[i] = prim_to_cons(W);
    }
    zero_slope(U, n_end);

    double t = 0.0;
    int step = 0;
    
    while (t < T_END) {
        double dt = calculate_dt(U, dx, n_end);
        if (t + dt > T_END) dt = T_END - t;

        //stage 1
        zero_slope(U, n_end);
        calculete_rhs(U, rhs, dx, n_end);
        for (int i=0; i<n_end; i++) {
            U1[i].rho   = U[i].rho   + dt * rhs[i].rho;
            U1[i].rho_u = U[i].rho_u + dt * rhs[i].rho_u;
            U1[i].E     = U[i].E     + dt * rhs[i].E;
        }

        //stage 2
        zero_slope(U1, n_end);
        calculete_rhs(U1, rhs, dx, n_end);
        for (int i=0; i<n_end; i++) {
            U[i].rho   = 0.5 * (U[i].rho   + U1[i].rho   + dt * rhs[i].rho);
            U[i].rho_u = 0.5 * (U[i].rho_u + U1[i].rho_u + dt * rhs[i].rho_u);
            U[i].E     = 0.5 * (U[i].E     + U1[i].E     + dt * rhs[i].E);
        }

        t += dt;
        step++;
        printf("Step: %d, Time: %.6f, dt: %.6f\n", step, t, dt);
    }

    write_and_plot(U, n_end, ng, dx);

    free(U);
    free(U1);
    free(rhs);

    return 0;
}
