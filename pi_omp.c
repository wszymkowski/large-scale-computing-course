#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <omp.h>


int main(int argc,char* argv[]){
    int steps, rank, numprocs, i;
    double mypi, pi, h, sum, x;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double start_time, end_time;
	
    while (1){
		if (rank == 0) {
			printf("Podaj liczbe krokow: (0 konczy prace) ");
			steps=1000000000;
			start_time = MPI_Wtime();
		}

		//start_time = MPI_Wtime();

		MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (steps == 0){
			printf("Koniec pracy\n");
			break;
		} 
	
		h   = 1.0 / (double) steps;
		sum = 0.0;

        omp_set_num_threads(8);
	    #pragma omp parallel for shared(numprocs, rank, steps, h), private(i, x), reduction(+:sum)
        
		for (i = rank + 1; i <= steps; i += numprocs) {
			x = h * ((double)i - 0.5);
			sum += 4.0 / (1.0 + x*x);
		}
		mypi = h * sum;

		//end_time = MPI_Wtime();
		
		MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (rank == 0){
			end_time = MPI_Wtime();
			printf("    pi:  %f\n", pi);
			printf("    czas obliczen: %f\n", end_time - start_time);
		}
    }
    MPI_Finalize();
    return 0;
}
