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

typedef struct {double rho, rho_u, E;} Conservative;
typedef struct {double rho, u, p;} Primitive;

static Primitive P_to_C(Conservative U)
{
    Primitive W;
    W.rho = U.rho;
    W.u   = U.rho_u / U.rho;
    W.p   = (GAMMA - 1.0) * (U.E - 0.5 * U.rho * W.u * W.u);
    return W;
}

static Conservative C_to_P(Primitive W)
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
    Conservative UL = C_to_P(WL);
    Conservative UR = C_to_P(WR);

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


int main(void)
{
    return 0;
}
