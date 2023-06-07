#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
        int N = 100000;
        double T = 0.000001;
        double a = 0;
        double b = 1;
        double l = 1;
        int c = 1;
        double h = l / N;
        int size, i, j;
        size = omp_get_max_threads();
        double tau = 0.3 * h * h;
        double *section;
        double *newsection;
        int sectionSize = N / size;
        newsection = (double*)malloc((N + 1) * sizeof(double));
        section = (double*)malloc((N + 1) * sizeof(double));
        for (i = 0; i < N + 1; i++)
                section[i] = 0;
        section[0] = a;
        section[N] = b;
        for (i = 0; i < (int)(T / tau); i++)
        {
#pragma omp parallel private(j)
{
#pragma omp for schedule (dynamic, sectionSize)
        for (j = 1; j < N; j++)
        {
                int n;
                n = omp_get_thread_num();
                newsection[j] = section[j] + tau * c * c / (h * h) * (section[j - 1] - 2 * section[j] + section[j + 1]);
        }
#pragma omp for schedule (dynamic, sectionSize)
        for (j = 1; j < N; j++)
                section[j] = newsection[j];
}
        }
        double x = -h;
        FILE * fout = fopen("out.txt", "w");
        for (i = 0; i <= N; i++)
        {
                x += h;
                fprintf(fout, "%lf %lf\n", x, section[i]);
        }
        free(section);
        free(newsection);
        return 0;
}
