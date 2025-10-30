#include <stdio.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif

#define PI 3.14159265358979323846

// 被積分関数
static inline double f1(double x) {
    return 1.0 / (x + 1.0);
}
static inline double f2(double x) {
    return sqrt( (1.0 + cos(x)) / 2.0 );
}

// 複合 Simpson 1/3 則
double simpson13(double (*f)(double), double a, double b, int n) {
    if (n % 2 != 0) {
        fprintf(stderr, "[simpson13] n must be even. got %d\n", n);
        return NAN;
    }
    const double h = (b - a) / n;
    double s_odd = 0.0, s_even = 0.0;
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 2) s_odd += f(x);
        else       s_even += f(x);
    }
    return (h/3.0) * (f(a) + f(b) + 4.0*s_odd + 2.0*s_even);
}

// 複合 Simpson 3/8 則
double simpson38(double (*f)(double), double a, double b, int n) {
    if (n % 3 != 0) {
        fprintf(stderr, "[simpson38] n must be a multiple of 3. got %d\n", n);
        return NAN;
    }
    const double h = (b - a) / n;
    double sum = f(a) + f(b);
    double s3 = 0.0, s2 = 0.0; 
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 3 == 0) s2 += f(x);
        else            s3 += f(x);
    }
    return (3.0*h/8.0) * (sum + 3.0*s3 + 2.0*s2);
}

int main(void) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    const double I1_true = log(2.0);                 
    const double I2_true = 2.0;

    const double a1 = 0.0, b1 = 1.0;
    const double a2 = 0.0, b2 = PI; 

    int n13_list[] = {4, 8, 12};   // 1/3 則は偶数
    int n38_list[] = {6, 9, 12};   // 3/8 則は 3 の倍数
    int m13 = sizeof(n13_list)/sizeof(n13_list[0]);
    int m38 = sizeof(n38_list)/sizeof(n38_list[0]);

    printf("=== Simpson 1/3 則（複合） ===\n");
    printf("%6s | %-12s | %-20s | %-20s\n", "n", "積分値", "絶対誤差 (I1=ln2)", "絶対誤差 (I2=2)");
    for (int i = 0; i < m13; ++i) {
        int n = n13_list[i];
        double I1 = simpson13(f1, a1, b1, n);
        double I2 = simpson13(f2, a2, b2, n);
        printf("%6d | %-12.10f | %-20.10e | %-20.10e\n",
               n, I1, fabs(I1 - I1_true), fabs(I2 - I2_true));
    }

    printf("\n=== Simpson 3/8 則（複合） ===\n");
    printf("%6s | %-12s | %-20s | %-20s\n", "n", "積分値", "絶対誤差 (I1=ln2)", "絶対誤差 (I2=2)");
    for (int i = 0; i < m38; ++i) {
        int n = n38_list[i];
        double I1 = simpson38(f1, a1, b1, n);
        double I2 = simpson38(f2, a2, b2, n);
        printf("%6d | %-12.10f | %-20.10e | %-20.10e\n",
               n, I1, fabs(I1 - I1_true), fabs(I2 - I2_true));
    }

    return 0;
}
