#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

int swap;

int compare(const void *a, const void *b)
{
	float c = *(float *)a;
	float d = *(float *)b;
	if(c < d) {return -1;}              
	else if (c == d) {return 0;}      
	else 
	{
		swap = 1;
		return 1;
	}
}


int isOdd(int a)
{
	if(a%2 != 0)
		return 1;
	return 0;
}


void merge(float *left, float *right, int lbase, int rbase, int isleft, float *tmp)
{

	int j=0, lptr=0, rptr=0, lbound, rbount;
	int copy_start, end, base;
	end = lbase + rbase;
	copy_start = isleft? 0:lbase;
	base = isleft? lbase: rbase;	

	int i;

	while(j<end){
	//while(lptr <lbase || rptr <rbase){
		//printf("lptr%d  rptr%d\n", lptr, rptr);
		//printf("left[lptr] %f right[rptr] %f\n",left[lptr], right[rptr]);
		if(lptr >= lbase)
			tmp[j++] = right[rptr++];
		else if (rptr >= rbase)
			tmp[j++] = left[lptr++];
		else
			tmp[j++] = (left[lptr] <= right[rptr])? left[lptr++]: right[rptr++];
		//printf("lptr%d  rptr%d\n", lptr, rptr); 
	}

	if(memcmp( (isleft?left:right), tmp+copy_start, base*4) !=0)
		swap=1;
	memcpy( (isleft?left:right), tmp+copy_start, base*4);
}


int main(int argc, char **argv){

	float temp, receiveTemp, *num;
	int i,total,all_swap=0, no_need_to_swap=0;
	int rank, world_size;
	int base, basePlus, firstIndex;
	int pre_base, next_base;
	float *tmp, *recvTmp, *sendTmp;

	double round=0;	

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

				if(rank == (basePlus-1))
					next_base = base-1;
				else
					next_base = base;
				pre_base = base;
			}
			else{
				firstIndex = basePlus*(base+1) + base*(rank-basePlus);
				next_base = base;
				if(rank == basePlus)
					pre_base = base +1;
				else
					pre_base = base;
			}
		}

		MPI_File_seek(fh, firstIndex*4, MPI_SEEK_SET);
		num = (float*) malloc(base* sizeof(float));
		tmp = (float*) malloc((base+1)*2* sizeof(float));
		recvTmp = (float*) malloc((base+1)* sizeof(float));
		sendTmp = (float*) malloc(base* sizeof(float));
		MPI_File_read(fh, num, base, MPI_FLOAT, &status);
	}
	// when (number of process > size of array) && (rank > size of array)
	else
		base = 0;



	/**** 	Sort ****/
	while(no_need_to_swap==0)	
	{
		qsort((void *)num, base, sizeof(num[0]), compare);	
		// even phase
		if(isOdd(rank) && rank!=0 && rank<total)
		{
			memcpy( sendTmp, num, base*4);
			MPI_Sendrecv(sendTmp, base, MPI_FLOAT, rank-1, rank,
					recvTmp, pre_base, MPI_FLOAT, rank-1, rank-1,
					MPI_COMM_WORLD, &status);
			merge(recvTmp, num, pre_base, base, 0, tmp); 
			//printf("***left rank:%d cpy_start %d base %d\n",rank, pre_base, base);
		}
		else if( !isOdd(rank) && rank< (world_size-1) && rank<(total-1))
		{
			memcpy( sendTmp, num, base*4);
			MPI_Sendrecv(sendTmp, base, MPI_FLOAT, rank+1, rank,
					recvTmp, next_base, MPI_FLOAT, rank+1, rank+1,
					MPI_COMM_WORLD, &status);
			merge(num, recvTmp, base, next_base, 1, tmp);
			//printf("---right rank:%d cpy_start %d base %d next_base%d\n",rank, 0, base, next_base);

		}
		//printf("even phase rank $d\n",rank);
		MPI_Barrier(MPI_COMM_WORLD);		


		// odd phase
		if( isOdd(rank) && rank< (world_size-1) && rank<(total-1))
		{
			memcpy( sendTmp, num, base*4);
			MPI_Sendrecv(sendTmp, base, MPI_FLOAT, rank+1, rank,
					recvTmp, next_base, MPI_FLOAT, rank+1, rank+1,
					MPI_COMM_WORLD, &status);
			merge(num, recvTmp, base, next_base, 1, tmp);
		}
		else if(!isOdd(rank) && rank!=0 && rank<total)
		{
			memcpy( sendTmp, num, base*4);
			MPI_Sendrecv(sendTmp, base, MPI_FLOAT, rank-1, rank,
					recvTmp, pre_base, MPI_FLOAT, rank-1, rank-1,
					MPI_COMM_WORLD, &status);

			merge(recvTmp, num, pre_base, base, 0, tmp);

		}

		//printf("odd phase rank%d\n",rank);
		// Use sum of swap (all_swap) to decide when to stop
		MPI_Reduce(&swap, &all_swap, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if(rank == 0){

			if(all_swap == 0)
				no_need_to_swap = 1;
			all_swap = 0;
		}

		swap = 0;

		// Broadcast the decition of end or not
		MPI_Bcast(&no_need_to_swap, 1, MPI_INT, 0, MPI_COMM_WORLD);

	}

	/****	Write Output File ****/
	MPI_File fout;
	MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_RDWR| MPI_MODE_CREATE, MPI_INFO_NULL, &fout);

	MPI_File_seek(fout, firstIndex*4, MPI_SEEK_SET);
	MPI_File_write(fout, num, base, MPI_FLOAT, &status);

	MPI_Finalize();
}

