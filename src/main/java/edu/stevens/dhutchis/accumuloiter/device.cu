#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "header_def.cuh"

__global__ void MSV_TEST(int* seq, int total, unsigned int* offset, float* sc, int* L, int* L_6r, int* mat, int base, int bias, int tbm, int tec, float scale, int hmm_size, int rib)			/* total is the total # of residues */
{
	/*--------------------------------------------------------*/
	/*					  Shared memory 					  */
	/*--------------------------------------------------------*/
	extern __shared__ unsigned char sm[];
	unsigned int*  cache = (unsigned int*)sm;
	unsigned char* MMX   = (unsigned char*)&cache[rib];
	unsigned char* MAT   = (unsigned char*)&MMX[rib * (hmm_size + 1)];
	unsigned char* RED   = (unsigned char*)&MAT[hmm_size * PROTEIN_TYPE];

	/*--------------------------------------------------------*/
	/*					    registers 		   			      */
	/*--------------------------------------------------------*/
	const int row = blockIdx.y * blockDim.y + threadIdx.y;			/* individual index for each warp across whole grid */

	int xE, xJ, xB;
	int sv;
	int mmx;
	int i, j, p;
	int res;
	int count = 0;													/* used to times of cycle for each block */
	int tjb;

	unsigned int off;
	int Len;

	/*--------------------------------------------------------*/
	/*	    Global to Shared (ONLY FOR SHARED VERSION) 		  */
	/*--------------------------------------------------------*/
	for(p = threadIdx.y * 32 + threadIdx.x; p < hmm_size * PROTEIN_TYPE; p += rib * 32)	/* 这里 +=32 是因为一行有32个threads */
	{																								/* CAUTION：since each block needs MAT[] separately */
		MAT[p] = mat[p];																			/* So we dont use 'row' to do this job */
	}
	__syncthreads();


	/*--------------------------------------------------------*/
	/*				     OUTER LOOP BEGIN !!! 		   		  */
	/*--------------------------------------------------------*/
	while(row + rib * gridDim.y * count < total)			/* row + ROW_IN_BLOCK * gridDim.y * count 代表本row下一个目标seq的index */
	{
		Len = L_6r[row + rib * gridDim.y * count];			/* for reuse, so cache into "Len" of register */
		off = offset[row + rib * gridDim.y * count];		/* offset (beginning point) for each seq */

		i = 0;														/* must be refresh for each new seq come in */
		mmx = 0;
		xJ = 0;
		tjb = -1.0f * roundf(scale * logf(3.0f / (float) (L[row + rib * gridDim.y * count] + 3)));		/* EACH seq, we need recalculate tjbv */
		tjb = (tjb > 255.) ? 255 : (unsigned char) tjb;
		xB = subm(base, tjb);

		for(p = threadIdx.x; p < hmm_size + 1; p += 32)
		{						
			MMX[threadIdx.y * (hmm_size+1) + p] = 0;			/* 0 is -INFINITY, So here is initial for new seq come in */
		}

		/*--------------------------------------------------------*/
		/*				     MIDDLE LOOP BEGIN !!! 		   		  */
		/*--------------------------------------------------------*/
		while(i < Len)
		{
			cache[threadIdx.y]  = seq[off + i];						/* GLOBAL ACCESS: for read a compressed cell in */
			if( ((cache[threadIdx.y] >> 25) & 0x0000001F) == 31)
				break;												/* Immediately check the end of seq */

			#pragma unroll
			for(j = 0; j < 6 ; j++)									/* REPEAT 6 TIMES! SINCE this is a compressed cell */						
			{
				xE = 0;
				xB = subm(xB, tbm);
				res = ( cache[threadIdx.y] >> 25 - 5 * j ) & 0x0000001F;

				if(res == 31){
					break;
				}
					
				mmx = (int)MMX[threadIdx.y * (hmm_size+1) + threadIdx.x];		/* IMPORTANT: pull mmx back to the head of segemnt */

				/*--------------------------------------------------------*/
				/*				     INNER LOOP BEGIN !!! 		   		  */
				/*--------------------------------------------------------*/
				for(p = threadIdx.x; p < hmm_size; p += 32)						/* bank-conflict free ADDRESSING */
				{
					sv  = max(mmx, xB);
					sv  = addm(sv, bias);                     				

					sv  = subm(sv, MAT[res * hmm_size + p]);					/* SHARED VERSION */
					//sv  = subm(sv, __ldg(&mat[res * HMM_SIZE + p]));			/* READ-ONLY CACHE VERSION */


					xE  = max(xE, sv);

					if(p + 32 < hmm_size)
						mmx = (int)MMX[threadIdx.y * (hmm_size+1) + p + 32];	/* PROTECT FROM overflow of memory */              	

					MMX[threadIdx.y * (hmm_size+1) + p + 1] = sv;

				}	/* end inner alignment */

				RED[threadIdx.y * 32 + threadIdx.x] = (unsigned char)xE;

				/* thread 0 - 15 for each row */
				if(threadIdx.x < 16)
				{
					//problem under RELEASE running...
					RED[threadIdx.y * 32 + threadIdx.x] = max(RED[threadIdx.y * 32 + threadIdx.x], RED[threadIdx.y * 32 + threadIdx.x + 16]);
					RED[threadIdx.y * 32 + threadIdx.x] = max(RED[threadIdx.y * 32 + threadIdx.x], RED[threadIdx.y * 32 + threadIdx.x + 8]);
					RED[threadIdx.y * 32 + threadIdx.x] = max(RED[threadIdx.y * 32 + threadIdx.x], RED[threadIdx.y * 32 + threadIdx.x + 4]);
					RED[threadIdx.y * 32 + threadIdx.x] = max(RED[threadIdx.y * 32 + threadIdx.x], RED[threadIdx.y * 32 + threadIdx.x + 2]);
					RED[threadIdx.y * 32 + threadIdx.x] = max(RED[threadIdx.y * 32 + threadIdx.x], RED[threadIdx.y * 32 + threadIdx.x + 1]);
				}

				xE = (int)RED[threadIdx.y * 32];

				/* Imrediately check high score sequence */
				if( addm(xE, bias) == 255 ) 
				{
					sc[row + rib * gridDim.y * count] = 999999.0f;
					break;						/* break out MIDDLE LOOP, go next seq */
				}

				/* get rest of parameters */
				xE = subm(xE, tec);      		/* EC = EJ in MSV */
				xJ = max(xJ, xE);        		/* xJ = max (xJ, xE - tEJ) */

				xB = max(base, xJ);
				xB = subm(xB, tjb);     		/* xB = max (base, xJ) - tJB */				
			}

			i++;	/* IMPORTANT: index++ for next MIDDLE CYCLE (next compressed cell of this seq) */

		}	/* end loop over sequence residues 1..L */

		/* finally C->T, and add our missing precision on the NN,CC,JJ back */
		if( abs(sc[row + rib * gridDim.y * count] - 999999.0f) < 1e-6 )
		{
			/* do nothing */
		} else {
			sc[row + rib * gridDim.y * count] = ((float) (xJ - tjb) - (float) base);
			sc[row + rib * gridDim.y * count] /= scale;
			sc[row + rib * gridDim.y * count] -= 3.0;       /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */       
		}

		count++;	/* IMPORTANT: index++, for next OUTER CYCLE (index next new seq) */

	}	/* end loop over sequence database */
}

/* --------------------------------------------------
 * kernel warper
 * --------------------------------------------------
 *
 *
 * hs: input variable -> HMM_SIZE
 * rib: input variable -> ROW_IN_BLOCK
 */
void MSV_warp(dim3 grid, dim3 block, int* SEQ, int total, unsigned int* OFFSET, float* SC, int* LEN, int* LEN_6R, int* MAT, HMMER_PROFILE* HMM, int device, int rib)
{
	size_t shmem = rib * sizeof(unsigned int) + rib * (HMM->M + 1) * sizeof(unsigned char) + HMM->M * PROTEIN_TYPE * sizeof(unsigned char) + 32 * rib * sizeof(unsigned char);

	MSV_TEST<<<grid, block, shmem>>>(SEQ, total, OFFSET, SC, LEN, LEN_6R, MAT, HMM->base_b, HMM->bias_b, HMM->tbm_b, HMM->tec_b, HMM->scale_b, HMM->M, rib);

	//printf("GPU %d launched: %s\n", device, cudaGetErrorString(cudaGetLastError()));
}