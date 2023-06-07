#include <stdio.h>
#include <mpi.h>
#include <math.h>

double f(double x)
{
        return sqrt(4 - x*x);
}

int main(int argc, char *argv[]){
        int N, rank, size;
        double a, b;

        N = 2000000000;
        a = 0, b = 2;

        double h = (b - a) / N; // длина шага

        MPI_Init(&argc, &argv);

        double I_k = 0;

        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        double x;
        int sectionSize = (int)(N / size);
        int i = 0;
        for (i = rank; i <= N; i += size)
        {                                       //распределяем рассчеты между процессами
                x = a + i * h;
                I_k += h * f(x);
        }
        if (N % size != 0 && rank + 1 + size * sectionSize <= N)
                I_k += h * f(size * sectionSize + rank + 1);

        if(rank == 0)
        {                                        // 0-й процесс будет считать общий интеграл. Он принимает значения I_k, посчитанные другими процессами
                double I = I_k;
                int j;
                for(j = 1; j < size; j++)
                {
                        double I_j = 0;
                        MPI_Recv(&I_j, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        I += I_j;
                }
                I -= h * (f(a)+f(b))/2;

                printf("Итоговое значение интеграла: %7f\n",I);      // печатаем итоговое значение
        }
        else
        {
                MPI_Send(&I_k, 1, MPI_DOUBLE, 0,  1, MPI_COMM_WORLD);           //ненулевые процессы отправляет посчитанные I_k 0-му процессу
        }

        MPI_Finalize();
        return 0;
}
