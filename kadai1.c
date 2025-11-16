#include <stdio.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif

#define N 5

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    double x[N] = {3.0, 3.1, 3.2, 3.3, 3.4};
    double y[N] = {1.09861, 1.13140, 1.16315, 1.19392, 1.22378};

    double h = x[1] - x[0];

    double diff[N][N];
    int i, j, k;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            diff[i][j] = 0.0;
        }
    }

    for (i = 0; i < N; i++) diff[i][0] = y[i];

    for (k = 1; k < N; k++) {
        for (i = 0; i < N - k; i++) {
            diff[i][k] = diff[i + 1][k - 1] - diff[i][k - 1];
        }
    }

    double X = 3.08;
    double u = (X - x[0]) / h;

    double sum = diff[0][0];
    for (k = 1; k < N; k++) {
        double numer = 1.0;
        for (j = 0; j < k; j++) {
            numer *= (u - j); /* (u)(u-1)(u-2)... */
        }
        double fact = 1.0;
        for (j = 1; j <= k; j++) fact *= j; 
        sum += (numer / fact) * diff[0][k];
    }

    printf("前進差分表:\n");
    for (k = 0; k < N; k++) {
        printf("Δ^%d: ", k);
        for (i = 0; i < N - k; i++) {
            printf("%f ", diff[i][k]);
        }
        printf("\n");
    }

    printf("\nX = %.2f, h = %.1f, u = %.2f\n", X, h, u);
    printf("Newton 前進公式による補間 ln(%.2f) = %f\n", X, sum);

    return 0;
}
