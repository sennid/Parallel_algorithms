#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char* argv[]){
        int N = 100000;
        double T = 0.000001;
        double a = 0;
        double b = 1;
        double l = 1;
        int c = 1;
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        double h = l / N;
        double tau = 0.3 * h * h;
        double *section;
        double *newsection;
        int sectionSize = N / size;
        if (rank == size - 1 && N % size != 0)
                sectionSize += N % size;
        newsection = (double*)malloc((sectionSize + 2) * sizeof(double));
        section = (double*)malloc((sectionSize + 2) * sizeof(double));
        int i;
        int m;
        FILE* fout;
        for (i = 0; i < sectionSize + 2; i++)
                section[i] = 0;
        if (rank == 0)
                section[1] = a;
        if (rank == size - 1)
        {
                section[sectionSize + 1] = b;
                printf("b = %lf", section[sectionSize]);
        }
        printf("rg: %d, h: %lf, T: %lf, tau: %10.3e T/tau:  %d\n", rank, h, T, tau, (int)(T / tau));
        for (i = 0; i < (int)(T / tau); i++)
        {
                if (rank % 2 == 0)
                {
                        if (rank > 0)
                        {
                                MPI_Send(section + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
                        //      printf("rank: %d send", rank);
                                MPI_Recv(section, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //      printf("rank: %d recv section[0]: %lf", rank, section[0]);
                        }
                        if (rank < size - 1)
                        {
                                MPI_Send(section + sectionSize, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                        //      printf("rank: %d send\n", rank);
                                MPI_Recv(section + sectionSize + 1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //      printf("rank: %d recv section[-1]: %lf\n", rank, section[sectionSize + 1]);
                        }
                }
                else
                {
                        if (rank > 0)
                        {
                                MPI_Recv(section, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //      printf("rank: %d recv section[0]: %lf", rank, section[0]);
                                MPI_Send(section + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
                        //      printf("rank: %d send\n", rank);
                        }
                        if (rank < size - 1)
                        {
                                MPI_Recv(section + sectionSize + 1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //      printf("rank: %d recv section[-1]: %lf", rank, section[sectionSize + 1]);
                                MPI_Send(section + sectionSize, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                        //      printf("rank: %d send\n", rank);
                        }

                }
                int left = 0;
                if (rank == 0)
                        left = 1;
                for (m = 1 + left; m <= sectionSize; m++)
                {
                        newsection[m] = section[m] + tau * c * c / (h * h) * (section[m - 1] - 2 * section[m] + section[m + 1]);
                //      printf("rg: %d, newsection[%d]: %lf, iteration: %d\n", rank, m, newsection[m], i);
                }
                for (m = 1 + left; m <= sectionSize; m++)
                        section[m] = newsection[m];

        }

        double x = -h;
        int right = 0;
        if (rank == size - 1)
                right = 1;
        if (rank == 0)
        {
                double *ux;
                int add = 0;
                ux = (double*)malloc((N + 1) * sizeof(double));;
                fout = fopen("out.txt", "w");
                for (i = 0; i < sectionSize; i++)
                        ux[i] = section[i + 1];
                for (i = 1; i <= size - 1; i++)
                {
                        if (i == size - 1)
                                add += N % size + 1;
                        MPI_Recv(ux + sectionSize * i, sectionSize + add, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                for (i = 0; i <= N; i++)
                {
                        x += h;
                        fprintf(fout, "%lf %lf\n", x, ux[i]);
                }
                fclose(fout);
                free(ux);
        }
        else
                MPI_Send(section + 1, sectionSize + right, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        free(section);
        free(newsection);
        MPI_Finalize();
        return 0;
}
