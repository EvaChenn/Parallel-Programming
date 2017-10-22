#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>


int swap = 0;

void cmp_swap(float *a, float *b)
{
	float temp;
	if(*a > *b)
	{
		temp = *b;
		*b = *a;
		*a = temp;
		swap = 1;
	}
}

int isOdd(int a)
{

	if(a%2 != 0)
		return 1;
	else
		return 0;
}

int main(int argc, char **argv){

	assert(argc == 4);

	float temp, receiveTemp, *num;
	int i,total;
	int rank,  world_size;
	int base, basePlus, firstIndex;
	int swap_sum=0, swap_zero=0;

	long compu_time=0,commu_time=0, io_time=0;
	double round=0;
	int even_i, odd_i;

	MPI_File fh;
	MPI_Status status;

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	total = atoi(argv[1]);

	// base: 	number of array that each process manage
	// basePlus: 	when (size of array) % (number of process) !=0, 
	// 		then add these remainders to 0~basePlus ranked process
	// firstIndex: 	real index in the input file


	/****	Count base, firstIndex; Read input file	****/
	// number of process < size of array
	if(rank < total)
	{
		// number of process == 1
		if(world_size == 1){
			firstIndex = 0;
			base = total;
		}
		else{

			base = total / world_size;
			basePlus = total % world_size;

			if(rank < basePlus){
				base++;
				firstIndex = rank*base;
			}
			else{
				firstIndex = basePlus*(base+1) + base*(rank-basePlus);
			}
		}

		MPI_File_seek(fh, firstIndex*4, MPI_SEEK_SET);
		num = (float*) malloc(base* sizeof(float));

		MPI_File_read(fh, num, base, MPI_FLOAT, &status);

	}
	// when (number of process > size of array) && (rank > size of array)
	else
		base = 0;


	if(isOdd(firstIndex)){
		odd_i = 0;
		even_i = 1;
	}
	else{
		odd_i = 1;
		even_i = 0;
	}


	/**** 	Sort ****/
	while(swap_zero!=2){	
		round++;
		swap = swap_sum = 0;


		// Even phase
		for(i=even_i; i<base; i+=2){

			// last index in num[] 
			if(i == (base-1)  )
			{
				// do nothing with last real index 
				if( (rank+1)<world_size && (rank+1)<total)
				{
					temp = num[i];

					// Send num[local rightmost] to rank+1 with Tag = real index of local rightmost
					// Recv num[rank+1's local leftmost] form rank+1 with Tag = real index of rank+1's local leftmost

					MPI_Sendrecv(&temp, 1, MPI_FLOAT, rank+1, firstIndex+i,
							&receiveTemp, 1, MPI_FLOAT, rank+1, firstIndex+i+1,
							MPI_COMM_WORLD,&status);

					if(receiveTemp < temp)
					{
						num[i] = receiveTemp;
						swap = 1;
					}
				}
			}
			else 
				cmp_swap(&num[i], &num[i+1]);
		}

		
		if(rank!=0 && even_i && base){

			temp = num[0];
			MPI_Sendrecv(&temp, 1, MPI_FLOAT, rank-1, firstIndex,
					&receiveTemp, 1, MPI_FLOAT, rank-1, firstIndex-1,
					MPI_COMM_WORLD,&status);

			if(receiveTemp > temp){
				num[0] = receiveTemp;
				swap = 1;
			}
		}

		MPI_Allreduce(&swap, &swap_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
		if(swap_sum == 0)
			swap_zero++;
		else
			swap_zero = 0;

		swap = swap_sum = 0;

		if(swap_zero == 2)
			break;


		// Odd phase
		for(i=odd_i; i<base; i+=2){
			if(i == (base-1) ){
				if( (rank+1)<world_size && (rank+1)<total){
					temp = num[i];

					// Send num[local rightmost] to rank+1 with Tag = real index of local rightmost
					// Recv num[rank+1's local leftmost] form rank+1 with Tag = real index of rank+1's local leftmost	
					MPI_Sendrecv(&temp, 1, MPI_FLOAT, rank+1, firstIndex+i,
							&receiveTemp, 1, MPI_FLOAT, rank+1, firstIndex+i+1,
							MPI_COMM_WORLD,&status);
					if(receiveTemp < temp){
						num[i] = receiveTemp;
						swap = 1;
					}
				}
			}
			else
				cmp_swap(&num[i], &num[i+1]);
		}
		// do nothing with rank 0 num[0]
		if(rank!=0 && odd_i && base){
			temp = num[0];
			MPI_Sendrecv(&temp, 1, MPI_FLOAT, rank-1, firstIndex,
					&receiveTemp, 1, MPI_FLOAT, rank-1, firstIndex-1,
					MPI_COMM_WORLD,&status);

			if(receiveTemp > temp){
				num[0] = receiveTemp;
				swap = 1;
			}
		}


		MPI_Allreduce(&swap, &swap_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

		if(swap_zero !=2){
			if(swap_sum == 0)
				swap_zero++;
			else
				swap_zero = 0;
		}
	}

	/****	Write Output File ****/
	MPI_File fout;
	MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_RDWR| MPI_MODE_CREATE, MPI_INFO_NULL, &fout);
	MPI_File_seek(fout, firstIndex*4, MPI_SEEK_SET);
	MPI_File_write(fout, num, base, MPI_FLOAT, &status);


	MPI_Finalize();
	//free(num);
}

