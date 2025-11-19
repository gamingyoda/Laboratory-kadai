#include <stdio.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif

#define PI 3.14159265358979323846

double f(double x)
{
    return sin(x);
}

double trap(double a, double b, int n)
{
    double h = (b - a) / n;
    double s = 0.5 * (f(a) + f(b));
    int i;
    for (i = 1; i < n; i++) {
        s += f(a + i * h);
    }
    return s * h;
}

int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    double a = 0.0;
    double b = PI / 2.0;
    double true_val = 1.0;

    int nlist[] = {2, 6, 12};
    int cnt = sizeof(nlist) / sizeof(nlist[0]);

    printf("台形公式による積分: I = ∫_0^{π/2} sin x dx\n");
    printf("真値 = %.15f\n\n", true_val);
    printf("n       近似値                 絶対誤差\n");
    printf("------  --------------------  ------------\n");

    int idx;
    for (idx = 0; idx < cnt; idx++) {
        int n = nlist[idx];
        double approx = trap(a, b, n);
        double err = fabs(approx - true_val);

        printf("%6d  %20.10f  %12.5e\n", n, approx, err);
    }

    return 0;
}

