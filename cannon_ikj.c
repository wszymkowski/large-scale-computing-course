#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


//const int n = 10, PP = 2, P = 4;
#define n 4002
#define PP 3
#define P 9

float temp[n / PP];
float A[n / PP][n / PP], B[n / PP][n / PP], C[n / PP][n / PP];


int main(int argc, char *argv[]) {
        FILE *file;
        FILE *file_out;

        int my_rank, ncpus;
        int row, col, mod = 0;
        int data_received = -1;
        int tag;
        int end;
        int from_rank, to_rank;

        MPI_Status statRecv[2];
        MPI_Request reqSend[2], reqRecv[2];

        double startwtime1, startwtime2, endtime;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Status status;

        if(my_rank == 0){
                printf("obliczenia metoda Cannona dla tabilicy %d x %d elementow \n", n, n);
        }
        if(my_rank == 0){
                startwtime1 = MPI_Wtime();
        }

        // wczytanie danych przez proces rank=0
        if(my_rank == 0){
                file = fopen("liczby4002.txt","r");
                if(file == NULL){
                        printf("Error while opening file.\n");
                        end = 1;
                        MPI_Bcast(&end, 1, MPI_INT, 0, MPI_COMM_WORLD);
                        MPI_Finalize();
                        exit(0);
                }else{
                        end = 0;
                        MPI_Bcast(&end, 1, MPI_INT, 0, MPI_COMM_WORLD);
                }
        }else{
                MPI_Bcast(&end, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if(end){
                        MPI_Finalize();
                        exit(0);
                }
        }

        if(ncpus != P){ 
                if(my_rank == 0){
                        printf("wywolano obliczenia iloczynu macierzy metoda cannona na %d procesach - uruchom mpiexec -n %d matrixmult\n",ncpus, P);
                } 
                MPI_Finalize();
                exit (0);
       }

        
        

        //wczytanie fragmentu danych przez proces 0 do tablic dwuwymiarowych A i B oraz rozesłanie pozostałych fragmentów do innych procesów
        if(my_rank == 0){
                for(int i = 0; i < n * PP; i++){
                        for(int j = 0; j < (n / PP); j++){
                                fscanf(file, "%f", &temp[j]);
                        }
                        to_rank = (i / n) * PP + (i % PP);
                        if(to_rank == 0){
                                for(int k = 0; k < (n / PP); k++){
                                        A[i / PP][k] = temp[k];
                                        B[i / PP][k] = temp[k];
                                }
                        }else{
                                MPI_Send(temp, (n / PP), MPI_FLOAT, to_rank, (i / PP) % (n / PP), MPI_COMM_WORLD);
                        }
                }
                fclose(file);
        }else{                          // odbiór fragmentów danych przez pozostale procesy
                for(int i = 0; i < (n / PP); i++){
                        MPI_Recv(temp, (n / PP), MPI_FLOAT, 0, i, MPI_COMM_WORLD, &status);
                        for(int j = 0; j < (n / PP); j++){
                                A[i][j] = temp[j];
                                B[i][j] = temp[j];
                        }
                }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        
        //przygotowanie tablicy wynikowej
        row = my_rank / PP;
        col = my_rank % PP;

        for(int i = 0; i < n / PP; i++){
                for(int j = 0; j < n / PP; j++){
                        C[i][j] = 0;
                }
        }

        if(my_rank == 0){
                startwtime2 = MPI_Wtime();
        }

        //obliczenia iloczynu macierzy zgodnie z algorytmem Cannona 
        int numprocs_sqrt = sqrt((double)ncpus);

        int row_id = row;
        int col_id = col;
        int first_in_row = row_id * numprocs_sqrt;

        int right = ((col_id + 1) % numprocs_sqrt) + first_in_row; 
        int left = ((col_id - 1 + numprocs_sqrt) % numprocs_sqrt) + first_in_row;

        int down = col_id + ((row_id + 1) % numprocs_sqrt) * numprocs_sqrt;
        int up = col_id + ((row_id - 1 + numprocs_sqrt) % numprocs_sqrt) * numprocs_sqrt;
        

        for(int kk = 0; kk < PP; kk++){
                for(int i = 0; i < n/PP; i++){
                        for(int k = 0; k < n/PP; k++){
                                for(int j = 0; j < n/PP; j++){
                                        C[i][j] += A[i][k] * B[k][j];
                                }
                        }
                }
                
                MPI_Isend(A, n*n / PP / PP, MPI_FLOAT, left, 0, MPI_COMM_WORLD, reqSend);
                MPI_Irecv(A, n*n / PP / PP, MPI_FLOAT, right, 0, MPI_COMM_WORLD, reqRecv);
                MPI_Isend(B, n*n / PP / PP, MPI_FLOAT, up, 0, MPI_COMM_WORLD, &reqSend[1]);
                MPI_Irecv(B, n*n / PP / PP, MPI_FLOAT, down, 0, MPI_COMM_WORLD, &reqRecv[1]);
                MPI_Wait(reqRecv, statRecv);
                MPI_Wait(&reqRecv[1], &statRecv[1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if(my_rank == 0){
                endtime = MPI_Wtime();
                printf("Total time: %f\n", endtime - startwtime1);
                printf("Calculation time: %f\n", endtime - startwtime2);
        }

        if(my_rank == 0){ //jeśli rank == 0 to zapisujemy wynik
                file_out = fopen("cannon_result.txt", "w");
                for (int i = 0; i < n * PP; i++) {
                        from_rank = (i / n) * PP + (i % PP);
                        if(from_rank == 0){
                                for(int j = 0; j < (n / PP); j++){
                                        fprintf(file_out, "%10.1f", C[(i / PP) % (n / PP)][j]);
                                }
                        }else{
                                MPI_Recv(temp, (n / PP), MPI_FLOAT, from_rank, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                for(int j = 0; j < (n / PP); j++){
                                        fprintf(file_out, "%10.1f", temp[j]);
                                }
                        }
                        if(i % PP == PP - 1){
                                fprintf(file_out, "\n");
                        }
                }
                fclose(file_out);
        }else{                          //jesli rank != 0 to wysylamy swoj fragment 
                for(int i = 0; i < n/PP; i++){
                        for(int j = 0; j < n/PP; j++){
                                temp[j] = C[i][j];
                        }
                        tag = (row) * n + i * PP + col;
                        MPI_Isend(temp, (n / PP), MPI_FLOAT, 0, tag, MPI_COMM_WORLD, reqSend);
                }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
}

