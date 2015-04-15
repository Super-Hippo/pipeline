/* ### this is a .cu file### */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>

#include "header_def.cuh"
#include "data_partition.cuh"		/* o	nly for multi-GPU */

extern "C" bool* mymain(int sizeout, const char** s, const char *Chmm_path )
{
	const char *MODEL_PATH = Chmm_path;
	//Temporary code here..
	int ACTIVE_BLOCKS = 30;
	int ROW_IN_BLOCK = 20;

	/* Once we have multiple devices, display them and set them as 'current' for ready to run */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);

	int device;
	for (device = 0; device < deviceCount; ++device)
	{
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, device);
		//printf("\nDevice %d has compute capability %d.%d.\n", device, deviceProp.major, deviceProp.minor);
	}

	HMMER_PROFILE *hmm = NULL;						/* claim hmm         */
	hmm = hmmer_profile_Create(MODEL_PATH);					/* alloc mem for hmm */

	/* ************************************************ */
	/* 		1. Parameters, match, insert emission 		*/
	/* ************************************************ */

	/* get parameters (M, MU[], LAMBDA[]) */
	if(get_Parameters(hmm, MODEL_PATH) != 1) printf("error!\n");

	/* Get raw match & insert emission */
	if(get_Emission(hmm, MODEL_PATH) != 1) printf("error!\n");

	/* Transfer to prob from raw data (need do it before get degen)*/
	if(log_Odd_score(hmm) != 1) printf("error!\n");

	/* Get degenerate residues */
	if(mat_degen(hmm) != 1) printf("error in getiing 'degenerate residue'! \n");

	/* Rescaling and rounding */
	if(mf_conversion(hmm) != 1) printf("error!\n");

	//Until now, all preparation works of MSV has been done!


	/* ************************************************************************ */
	/* 		2. Transition probability only for Viterbi and Forward filter 		*/
	/* ************************************************************************ */

#if 1
	/**/
	if(get_transition(hmm, MODEL_PATH) != 1) printf("error!\n");

	/**/
	if(get_entryScore(hmm) != 1) printf("error!\n");

	/**/
	if(log_Trans(hmm) != 1) printf("error!\n");

	/**/
	if(xTrans(hmm) != 1) printf("error!\n");

	/**/
	if(vf_conversion(hmm) != 1) printf("error!\n");

	//Until now, all preparation works of Viterbi has been done!

#endif	/* Configure for Viterbi */

	/* **************************************************************************** */
	/* 		3. Protein sequence database read and pre-process               		*/
	/* **************************************************************************** */

	int number = 0;		/* the total number of seq in this database */

	/**/
	number = sizeout; //get_Seqnumber(DATABASE_PATH);
	//printf("\n In this database, we have %d seqs.\n", number);

	/**/
	const char** seq = (const char**)malloc(number * sizeof(char*));     					//dynamic memory for address value of each sequence: seq[i]->an address value
	int* seq_len = (int*)malloc(number * sizeof(int));       					//for cache length of each sequence
	//if(alloc_Eachseq(seq, seq_len, number, DATABASE_PATH) != 1) printf("error!\n");
	//if(fill_Eachseq(seq, seq_len, number, DATABASE_PATH) != 1) printf("error!\n");

	seq = s;

	for(int i =0;i< number;i++)//sets seq_length
	{
		seq_len[i] = strlen(seq[i]);
	}

#if 0
	/* verify reading data */
	FILE* vrd = NULL;
	vrd = fopen("E:\\verify_database", "w");
	int avg_len = 0;
	for(int i = 0; i < number; i++)
	{
		fprintf(vrd, "\n#%d===\n", i);
		avg_len = avg_len + seq_len[i];

		for(int j = 0; j < seq_len[i]; j++)
		{
			fprintf(vrd, "%c", seq[i][j]);
		}
	}
	avg_len = avg_len/number;
	fprintf(vrd, "\n==========\n, the avg_len = %d", avg_len);
	fclose(vrd);
#endif

	/* transfer seq into int values */
	char** iSeq = (char**)malloc(number * sizeof(char*));
	for(int i = 0; i < number; i++)
		iSeq[i] = (char*)malloc(seq_len[i] * sizeof(char));
	if(seq_To_char(iSeq, seq, seq_len, number) != 1) printf("error in transfer to short from char !\n");


	/* ***************************************** */
	/* 		6. Parameters for kernel launch	 	 */
	/* ***************************************** */

	/* gird & block */
	dim3 gridDim = dim3(1, ACTIVE_BLOCKS, 1);		//This is better since it is able to allocated paras dynamically.
	dim3 blockDim = dim3(32, ROW_IN_BLOCK, 1);

	//revolution begin
	int SPAN = number / deviceCount;
	int REMAIN = number % deviceCount;

	GPUs *kernel = new GPUs[deviceCount];	//dynamic allocation..

	//can try to use multiple threads...
	for(int i = 0; i < deviceCount; i++)
	{
		cudaSetDevice(i);
		kernel[i].get_duty(SPAN, REMAIN, i, deviceCount);	/* distribution */
		kernel[i].compress(iSeq, seq_len);					/* residues pack */
		kernel[i].get_Real_Length(seq_len);					/* real Length of everyone */
		kernel[i].trans_to_1D();							/* 2D array to 1D array */
		kernel[i].host_to_device(hmm);
	}

	//launch..
	for (int i = 0; i < deviceCount; i++)				 //Now we havent consider the time cost of data copy yet
	{
		cudaSetDevice(i);
		kernel[i].launch_MSV(gridDim, blockDim, hmm, i, ROW_IN_BLOCK);
	}
	cudaDeviceSynchronize();

	//get back results and counter..

	bool* passed_msv;
	passed_msv = (bool*)malloc(number*sizeof(bool));


	int counter = 0;
	for (int i = 0; i < deviceCount; i++)
	{
		cudaSetDevice(i);
		kernel[i].device_to_host(hmm,passed_msv);
	}

	//printf("\nTotally %d sequences pass through the MSV filter", counter);

	for(int i = 0; i < number ; i++)
	{
		//printf("\nsequence %i passed : %d",i,passed_msv[i]);
		if(passed_msv[i]==1)
		{
			counter++;
		}
	}

	/* ***************************************** */
	/* 		  6. Host & Device mem free   	 	 */
	/* ***************************************** */

	/**/
	delete[] kernel;


	/* Free */
	for(int i = 0; i < number; i++) {
		free(iSeq[i]);
		//free(seq[i]);
	}
	free(iSeq);				/* PLEASE REMIND: free child first, then parent */
	//free(seq);

	/* HMM model */
	freeHMM(hmm);

	cudaDeviceReset();		/* should put here since we have many cuda function before */

	return passed_msv;

}
