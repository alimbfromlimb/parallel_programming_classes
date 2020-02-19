#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

double function(double arg) {
    return 1/(1+arg^2);
}

double trap_area(double x0, double x1) {
    f0 = function(x0);
    f1 = function(x1);
    return ((f0+f1)/2)*(x1-x0);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Please provide the program with the number of segments and the number of processes\n");
        exit(0);
    } else {
        int N = argv[1];
        int size = argv[2];
    }

    int i;
    int myrank;

    int batch_size = (int) floor(N/size);
    float residual = (float) N/size;
    if (residual != 0) {size += 1;}

    double *buf = (double *) calloc(N + 1, sizeof(double));
    if (buf == NULL) { 
        printf("Memory for buf not allocated.\n"); 
        exit(0); 
    } 
    
    double *integral = (double *) calloc(size, sizeof(double));
    if (integral == NULL) { 
        printf("Memory for integral not allocated.\n"); 
        exit(0); 
    }

    double integral_0 = 0;
    double integral_ = 0;

    double begin_seq, end_seq;
    double begin_par, end_par;
     
    MPI_Status Status;
    // N processes start living after Init command
    MPI_Init(&argc, &argv);
    // size присваивается число, равное количеству процессов
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // каждому процессу присваивается уникальный номер
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // тут они все работают одновременно
    if (myrank == 0) {       
        for (i = 0; i < N + 1; i++) {
            buf[i] = i/N;
        }

        // считает всё самостоятельно
        begin_seq = MPI_Wtime();
        for (i = 0; i < N; i++) {
            integral_0 += trap_area(buf[i], buf[i+1]);
        }        
        end_seq = MPI_Wtime();

        // считает параллельно
        begin_par = MPI_Wtime();
        for (i = 1; i < size; i++) {
            // отсылка сообщения процессам с номерами i
            MPI_Send(buf[i*batch_size], batch_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
            MPI_Send(integral[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        }
        // пусть сам хоть что-то посчитает
        for (i = 0; i < batch_size - 1; i++) {
            integral[0] += trap_area(buf[i], buf[i+1]);
            integral_ += integral[0];
        }
        for (i = 1; i < size; i++) {
            // принятие сообщений процессов с номерами i
            MPI_Recv(integral[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Status);
            integral_ += integral[i];
        }
        end_par = Wtime();
    }
    if (myrank != 0) {
        // каждый из процессов получает сообщение от процесса с номером 0
        MPI_Recv(buf[myrank*batch_size], batch_size, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, &Status);
        MPI_Recv(integral[myrank], 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, &Status);
        for (i = 0; i < batch_size - 1; i++) {
            integral[myrank] += trap_area(buf[i], buf[i+1]);
        }
        MPI_Send(integral[myrank], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    printf("J_seq = %e, J_par = %e\n", integral_0, integral_);
    printf("Time_seq = %e, Time_par = \n", end_seq - begin_seq, end_par - begin_par);
    MPI_Finalize();
    free(buf);    
    free(integral);
    return 0;
}