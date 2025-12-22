/*  
2次元の衝撃波管問題(Sod問題)
2次元Euler方程式 + FDS（Roe-FDS）+ Runge–Kutta

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
#define NY 1000
#define X0 0.0
#define X1 1.0
#define Y0 0.0
#define Y1 1.0

#define XDIFF 0.5
#define CFL 0.5
#define T_END 0.20

#define RHO_L  1.0
#define U_L    0.0
#define P_L    1.0

#define RHO_R  0.125
#define U_R    0.0
#define P_R    0.1

#define WL_INIT  ((Primitive){RHO_L, U_L, 0.0, P_L})
#define WR_INIT  ((Primitive){RHO_R, U_R, 0.0, P_R})


typedef struct { double rho, rho_u, rho_v, E; } Conservative;
typedef struct { double rho, u, v, p; } Primitive;

static Primitive cons_to_prim(Conservative U)
{
    Primitive W;
    W.rho = U.rho;
    W.u   = U.rho_u / U.rho;
    W.v   = U.rho_v / U.rho;
    W.p   = (GAMMA - 1.0) * (U.E - 0.5 * U.rho * (W.u*W.u + W.v*W.v));
    return W;
}

static Conservative prim_to_cons(Primitive W)
{
    Conservative U;
    U.rho  = W.rho;
    U.rho_u = W.rho * W.u;
    U.rho_v = W.rho * W.v;
    U.E    = W.p/(GAMMA-1.0) + 0.5 * W.rho * (W.u*W.u + W.v*W.v);
    return U;
}

static Conservative FluxX(Primitive W)
{
    const double E = W.p/(GAMMA-1.0) + 0.5*W.rho*(W.u*W.u + W.v*W.v);
    Conservative F;
    F.rho  = W.rho * W.u;
    F.rho_u = W.rho * W.u * W.u + W.p;
    F.rho_v = W.rho * W.u * W.v;
    F.E    = W.u * (E + W.p);
    return F;
}

static Conservative FluxY(Primitive W)
{
    const double E = W.p/(GAMMA-1.0) + 0.5*W.rho*(W.u*W.u + W.v*W.v);
    Conservative G;
    G.rho  = W.rho * W.v;
    G.rho_u = W.rho * W.u * W.v;
    G.rho_v = W.rho * W.v * W.v + W.p;
    G.E    = W.v * (E + W.p);
    return G;
}

static void ROE_Flux_X(Primitive WL, Primitive WR, Conservative *Fout)
{
    Conservative FL = FluxX(WL);
    Conservative FR = FluxX(WR);
    Conservative UL = prim_to_cons(WL);
    Conservative UR = prim_to_cons(WR);

    double u_t = (sqrt(WL.rho)*WL.u + sqrt(WR.rho)*WR.u) / (sqrt(WL.rho)+sqrt(WR.rho));
    double v_t = (sqrt(WL.rho)*WL.v + sqrt(WR.rho)*WR.v) / (sqrt(WL.rho)+sqrt(WR.rho));
    double H_t = (sqrt(WL.rho)*(UL.E + WL.p)/WL.rho + sqrt(WR.rho)*(UR.E + WR.p)/WR.rho)/(sqrt(WL.rho)+sqrt(WR.rho));
    double a_t = sqrt((GAMMA-1)*(H_t - 0.5 * (u_t * u_t + v_t * v_t)));
    
    double dr = WR.rho - WL.rho;
    double du = WR.u - WL.u;  
    double dv = WR.v - WL.v;  
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
    double d3 = a1*alpha1*v_t + a2*alpha2*v_t + a3*alpha3*v_t;

    Fout->rho   = 0.5 * (FL.rho   + FR.rho  ) - 0.5 * d0;
    Fout->rho_u = 0.5 * (FL.rho_u + FR.rho_u) - 0.5 * d1;
    Fout->rho_v  = 0.5 * (FL.rho_v  + FR.rho_v ) - 0.5 * d3;
    Fout->E     = 0.5 * (FL.E     + FR.E    ) - 0.5 * d2;
}

static void ROE_Flux_Y(Primitive WL, Primitive WR, Conservative *Fout)
{
    Conservative GL = FluxY(WL);
    Conservative GR = FluxY(WR);
    Conservative UL = prim_to_cons(WL);
    Conservative UR = prim_to_cons(WR);

    double u_t = (sqrt(WL.rho)*WL.u + sqrt(WR.rho)*WR.u) / (sqrt(WL.rho)+sqrt(WR.rho));
    double v_t = (sqrt(WL.rho)*WL.v + sqrt(WR.rho)*WR.v) / (sqrt(WL.rho)+sqrt(WR.rho));
    double H_t = (sqrt(WL.rho)*(UL.E + WL.p)/WL.rho + sqrt(WR.rho)*(UR.E + WR.p)/WR.rho)/(sqrt(WL.rho)+sqrt(WR.rho));
    double a_t = sqrt((GAMMA-1)*(H_t - 0.5 * (u_t * u_t + v_t * v_t)));
    
    double dr = WR.rho - WL.rho;
    double du = WR.u - WL.u;  
    double dv = WR.v - WL.v;  
    double dp = WR.p - WL.p;

    double alpha2 = dr - dp/(a_t*a_t);
    double alpha1 = (dp - sqrt(WL.rho*WR.rho)*a_t*dv)/(2*a_t*a_t);
    double alpha3 = (dp + sqrt(WL.rho*WR.rho)*a_t*dv)/(2*a_t*a_t);

    double lamda1 = v_t - a_t;
    double lamda2 = v_t;
    double lamda3 = v_t + a_t;

    double delta = 0.1 * a_t;
    double a1 = (fabs(lamda1) >= delta) ? fabs(lamda1) : (lamda1*lamda1 + delta*delta)/(2.0*delta); 
    double a2 = (fabs(lamda2) >= delta) ? fabs(lamda2) : (lamda2*lamda2 + delta*delta)/(2.0*delta); 
    double a3 = (fabs(lamda3) >= delta) ? fabs(lamda3) : (lamda3*lamda3 + delta*delta)/(2.0*delta);

    double d0 = a1*alpha1*1.0 + a2*alpha2*1.0 + a3*alpha3*1.0;
    double d1 = a1*alpha1*(H_t - v_t*a_t) + a2*alpha2*0.5*v_t*v_t + a3*alpha3*(H_t + v_t*a_t);
    double d2 = a1*alpha1*v_t + a2*alpha2*v_t + a3*alpha3*v_t;
    double d3 = a1*alpha1*(v_t - a_t) + a2*alpha2*v_t + a3*alpha3*(v_t + a_t);
    
    Fout->rho   = 0.5 * (GL.rho   + GR.rho  ) - 0.5 * d0;
    Fout->rho_u = 0.5 * (GL.rho_u + GR.rho_u) - 0.5 * d1;
    Fout->rho_v  = 0.5 * (GL.rho_v  + GR.rho_v ) - 0.5 * d2;
    Fout->E     = 0.5 * (GL.E     + GR.E    ) - 0.5 * d3;
}

static inline double minmod(double a, double b)
{
    if (a*b <= 0.0) return 0.0;
    return (fabs(a) < fabs(b)) ? a : b;
}

static void zero_slope(Conservative *U, int nx_end, int ny_end, int ng)
{
    for (int j = 0; j < ny_end; j++) {
        for (int i = 0; i < ng; i++) {
            U[j * nx_end + i] = U[j * nx_end + ng];
            U[j * nx_end + (nx_end - 1 - i)] = U[j * nx_end + (nx_end - 1 - ng)];
        }
    }

    for (int i = 0; i < nx_end; i++) {
        for (int j = 0; j < ng; j++) {
            U[j * nx_end + i] = U[ng * nx_end + i];
            U[(ny_end - 1 - j) * nx_end + i] = U[(ny_end - 1 - ng) * nx_end + i];
        }
    }

}

static void calculate_rhs(Conservative *U, Conservative *rhs, double dx, double dy, int nx_end, int ny_end)
{
    for (int j = 0; j < ny_end; j++) {
        for (int i = 0; i < nx_end; i++) {
            rhs[j * nx_end + i] = (Conservative){0, 0, 0, 0};
        }
    }

    Conservative *Fx = (Conservative*)malloc(sizeof(Conservative) * (nx_end - 1) * ny_end);
    Conservative *Gy = (Conservative*)malloc(sizeof(Conservative) * nx_end * (ny_end - 1));

    // mimod

    for (int j = 0; j < ny_end; j++) {
        for (int i = 0; i < nx_end - 1; i++) {
            Primitive WL = cons_to_prim(U[j * nx_end + i]);
            Primitive WR = cons_to_prim(U[j * nx_end + i + 1]);
            ROE_Flux_X(WL, WR, &Fx[j * (nx_end - 1) + i]);
        }
    }

    for (int j = 0; j < ny_end - 1; j++) {
        for (int i = 0; i < nx_end; i++) {
            Primitive WL = cons_to_prim(U[j * nx_end + i]);
            Primitive WR = cons_to_prim(U[(j + 1) * nx_end + i]);
            ROE_Flux_Y(WL, WR, &Gy[j * nx_end + i]);
        }
    }

    for (int j = 1; j < ny_end - 1; j++) {
        for (int i = 1; i < nx_end - 1; i++) {
            rhs[j * nx_end + i].rho   = -(Fx[j * (nx_end - 1) + i].rho   - Fx[j * (nx_end - 1) + i - 1].rho) / dx -(Gy[j * nx_end + i].rho   - Gy[(j - 1) * nx_end + i].rho) / dy;
            rhs[j * nx_end + i].rho_u = -(Fx[j * (nx_end - 1) + i].rho_u - Fx[j * (nx_end - 1) + i - 1].rho_u) / dx -(Gy[j * nx_end + i].rho_u - Gy[(j - 1) * nx_end + i].rho_u) / dy;
            rhs[j * nx_end + i].rho_v  = -(Fx[j * (nx_end - 1) + i].rho_v  - Fx[j * (nx_end - 1) + i - 1].rho_v) / dx -(Gy[j * nx_end + i].rho_v - Gy[(j - 1) * nx_end + i].rho_v) / dy;
            rhs[j * nx_end + i].E     = -(Fx[j * (nx_end - 1) + i].E     - Fx[j * (nx_end - 1) + i - 1].E) / dx -(Gy[j * nx_end + i].E     - Gy[(j - 1) * nx_end + i].E) / dy;
        }
    }
    free(Fx);
    free(Gy);
}

static double calculate_dt(const Conservative *U, double dx, double dy, int nx_end, int ny_end)
{
    double smax = 0.0;

    for (int j = 0; j < ny_end; j++) {
        for (int i = 0; i < nx_end; i++) {
            Primitive W = cons_to_prim(U[j * nx_end + i]);
            if (W.rho <= 0.0 || W.p <= 0.0) continue;

            const double a = sqrt(GAMMA * W.p / W.rho);
            const double s = fabs(W.u) + fabs(W.v) + a;
            if (s > smax) smax = s;
        }
    }

    return CFL / (smax * (1.0/dx + 1.0/dy));
}

static void write_and_plot(const Conservative *U, int nx_end, int ny_end, int ng, double dx, double dy)
{
    FILE *fp = fopen("sod2d.dat", "w");
    if (!fp) { perror("sod2d.dat"); return; }

    fprintf(fp, "# x y rho u v p\n");
    for (int j=ng; j<ny_end-ng; j++) {
        for (int i=ng; i<nx_end-ng; i++) {
            const double x = X0 + ((i-ng) + 0.5)*dx;
            const double y = Y0 + ((j-ng) + 0.5)*dy;
            Primitive W = cons_to_prim(U[j * nx_end + i]);
            fprintf(fp, "%.10f %.10f %.10f %.10f %.10f %.10f\n", x, y, W.rho, W.u, W.v, W.p);
        }
    }
    fclose(fp);

    const char *gpname = "sod2d.gp";
    FILE *gp = fopen(gpname, "w");
    if (!gp) { perror("sod2d.gp"); return; }

    fprintf(gp, "set grid\n");
    fprintf(gp, "set pm3d map\n");
    fprintf(gp, "set title 'Sod 2D: Roe-FDS + TVD RK2 (t=%.3f, NX=%d, NY=%d)'\n", T_END, NX, NY);
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set cblabel 'Density rho'\n");
    fprintf(gp, "splot 'sod2d.dat' using 1:2:3 with pm3d notitle\n");
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
    const int nx_end = NX + 2*ng;
    const int ny_end = NY + 2*ng;
    double dx = (X1 - X0)/NX;
    double dy = (Y1 - Y0)/NY;

    Conservative *U  = (Conservative*)malloc(sizeof(Conservative)*nx_end*ny_end);
    Conservative *U1 = (Conservative*)malloc(sizeof(Conservative)*nx_end*ny_end);
    Conservative *rhs= (Conservative*)malloc(sizeof(Conservative)*nx_end*ny_end);

    Primitive WL = WL_INIT;
    Primitive WR = WR_INIT;

    for (int j=0; j<ny_end; j++) {
        for (int i=0; i<nx_end; i++) {
            double x = X0 + ( (i-ng) + 0.5 )*dx;
            double y = Y0 + ( (j-ng) + 0.5 )*dy;
            Primitive W = (x < XDIFF) ? WL : WR;
            U[j * nx_end + i] = prim_to_cons(W);
        }
    }
    zero_slope(U, nx_end, ny_end, ng);

    double t = 0.0;
    int step = 0;
    
    while (t < T_END) {
        double dt = calculate_dt(U, dx, dy, nx_end, ny_end);
        if (t + dt > T_END) dt = T_END - t;

        //stage 1
        zero_slope(U, nx_end, ny_end, ng);
        calculate_rhs(U, rhs, dx, dy, nx_end, ny_end);
        for (int j=0; j<ny_end; j++) {
            for (int i=0; i<nx_end; i++) {
                int idx = j * nx_end + i;
                U1[idx].rho   = U[idx].rho   + dt * rhs[idx].rho;
                U1[idx].rho_u = U[idx].rho_u + dt * rhs[idx].rho_u;
                U1[idx].rho_v  = U[idx].rho_v  + dt * rhs[idx].rho_v;
                U1[idx].E     = U[idx].E     + dt * rhs[idx].E;
            }
        }

        //stage 2
        zero_slope(U1, nx_end, ny_end, ng);
        calculate_rhs(U1, rhs, dx, dy, nx_end, ny_end);
        for (int j=0; j<ny_end; j++) {
            for (int i=0; i<nx_end; i++) {
                int idx = j * nx_end + i;
                U[idx].rho   = 0.5 * (U[idx].rho   + U1[idx].rho   + dt * rhs[idx].rho);
                U[idx].rho_u = 0.5 * (U[idx].rho_u + U1[idx].rho_u + dt * rhs[idx].rho_u);
                U[idx].rho_v  = 0.5 * (U[idx].rho_v  + U1[idx].rho_v  + dt * rhs[idx].rho_v);
                U[idx].E     = 0.5 * (U[idx].E     + U1[idx].E     + dt * rhs[idx].E);
            }
        }
        t += dt;
        step++;
        printf("Step: %d, Time: %.6f, dt: %.6f\n", step, t, dt);
    }

    write_and_plot(U, nx_end, ny_end, ng, dx, dy);

    free(U);
    free(U1);
    free(rhs);

    return 0;
}