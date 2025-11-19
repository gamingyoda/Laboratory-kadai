#include <stdio.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif

#define G      9.81
#define L      1.0 
#define DT     0.01  // 時間刻み
#define STEPS  10000

// ここを弄る
const double theta0_list[] = { 
    0.5,  1.0,  1.5,  2.0,  2.5, 3.0,
};

enum {
    N_PARAM_THETA = sizeof(theta0_list) / sizeof(theta0_list[0]),
};

#define N_PARAM N_PARAM_THETA

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    FILE *fp;
    int i, j;
    double t;

    double theta[N_PARAM];
    double omega[N_PARAM];

    double k1_theta, k1_omega;
    double k2_theta, k2_omega;
    double k3_theta, k3_omega;
    double k4_theta, k4_omega;
    double theta2, omega2;
    double theta3, omega3;
    double theta4, omega4;

    t = 0.0;
    for (j = 0; j < N_PARAM; j++) {
        theta[j] = theta0_list[j];
        omega[j] = 0.0;
    }

    fp = fopen("kadai3-1.dat", "w");
    if (fp == NULL) {
        printf("kadai3-1.datを開けません\n");
        return 1;
    }

    fprintf(fp, "# t");
    for (j = 0; j < N_PARAM; j++) {
        fprintf(fp, "  theta%d(rad)  omega%d(rad/s)", j+1, j+1);
    }
    fprintf(fp, "\n");

    for (i = 0; i < STEPS; i++) {

        fprintf(fp, "%lf", t);

        for (j = 0; j < N_PARAM; j++) {

            k1_theta = omega[j];
            k1_omega = - (G / L) * sin(theta[j]);

            theta2 = theta[j] + 0.5 * DT * k1_theta;
            omega2 = omega[j] + 0.5 * DT * k1_omega;
            k2_theta = omega2;
            k2_omega = - (G / L) * sin(theta2);

            theta3 = theta[j] + 0.5 * DT * k2_theta;
            omega3 = omega[j] + 0.5 * DT * k2_omega;
            k3_theta = omega3;
            k3_omega = - (G / L) * sin(theta3);

            theta4 = theta[j] + DT * k3_theta;
            omega4 = omega[j] + DT * k3_omega;
            k4_theta = omega4;
            k4_omega = - (G / L) * sin(theta4);

            theta[j] = theta[j] + (DT / 6.0) *
                       (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
            omega[j] = omega[j] + (DT / 6.0) *
                       (k1_omega + 2.0 * k2_omega + 2.0 * k3_omega + k4_omega);

            fprintf(fp, " %lf %lf", theta[j], omega[j]);
        }

        fprintf(fp, "\n");
        t = t + DT;
    }

    fclose(fp);

    {
        char cmd[512];
        int offset = 0;

        offset += snprintf(cmd + offset, sizeof(cmd) - offset,
                           "gnuplot -persist -e \"plot ");

        for (int j = 0; j < N_PARAM; j++) {
            int col_theta = 2 + 2 * j;
            int col_omega = 3 + 2 * j;
            char legend[64];
            snprintf(legend, sizeof(legend), "theta0=%.2f", theta0_list[j]);

            offset += snprintf(cmd + offset, sizeof(cmd) - offset,
                               "'kadai3-1.dat' using %d:%d with lines title '%s'%s",
                               col_theta, col_omega, legend,
                               (j != N_PARAM - 1) ? ", " : "");
        }
        snprintf(cmd + offset, sizeof(cmd) - offset, "\"");

        system(cmd);
    }

    return 0;
}