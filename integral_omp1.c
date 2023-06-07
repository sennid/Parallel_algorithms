#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>

double f(double x)
{
        return sqrt(4 - x * x);
}

int main(int argc, char *argv[])
{
        int N, size, sectionSize;
        double a, b;
        N = 1000000000;
        a = 0, b = 2;
        size = omp_get_max_threads();
        sectionSize = (int)(N / size);
        double h = (b - a) / N;
        int i = 0;
        double I_k = 0;
#pragma omp parallel reduction (+: I_k)
{
        double x;
        int i;
#pragma omp for schedule (static, sectionSize)
        for (i = 1; i <= N - 1; i++)
        {
                x = a + i * h;
                I_k += h * f(x);
        }
}
        I_k += h * (f(a) + f(b)) / 2; //добавка из метода трапеций
        printf("Итоговое значение интеграла: %7f\n", I_k);
        return 0;

}
