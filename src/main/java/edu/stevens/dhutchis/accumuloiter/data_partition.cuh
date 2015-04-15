/* #### this is a .cuh file */


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "header_def.cuh"

/* This is a header file
 * 
 * I just want to simplify my multi-GPUs code through using C++
 * since here I can utilize "CLASS"
 */

#ifndef GPUS_PARTITION
#define GPUS_PARTITION

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

//define a class for duty allocation to each GPU hardware..
//each instance of this class will be loaded onto an individual GPU hardware..
//
class GPUs
{
  private:

    /**/
    unsigned int b_point;   /* BEGIN point of data partition */
    unsigned int e_point;   /* END point of ~ */

    int** seq_6r;
    int*  seq_len_6r;

    unsigned long acc;      /* largest data type here since this is sum of all seqs */

    int* real_L;            /* real length, not compressed length */

    int* seq_1D;            /* 1D array that caches all digitized seq */

    unsigned int* offset;   /* indexing of that 1D array */

    /**/
    int* d_seq;
    unsigned int* d_offset;
    int* d_len;
    int* d_len_6r;
    int* mat_m;
    float* score;

	public:

    GPUs()        /* Constructor: can be rewrite; can NOT be called explicitly */
    {
        /**/
        b_point = 0;
        e_point = 0;
        seq_6r = NULL;
        seq_len_6r = NULL;
        acc = 0;
        real_L = NULL;
        seq_1D = NULL;
        offset = NULL;
        /**/
        d_seq = NULL;
        d_offset = NULL;
        d_len = NULL;
        d_len_6r = NULL;
        mat_m = NULL;
        score = NULL;
    }

    ~GPUs()       /* Destructor: can NOT be rewrite; can be called explicitly */
    {
        int i;
        /**/
        free(seq_1D);
        free(offset);
        free(real_L);
        for(i = 0; i < (e_point - b_point); i++)
        {
            free(seq_6r[i]);
        }
        free(seq_6r);
        free(seq_len_6r);
        /**/
        cudaFree(d_seq);
        cudaFree(d_offset);
        cudaFree(d_len);
        cudaFree(d_len_6r);
        cudaFree(mat_m);
        cudaFree(score);
    }

	  void get_duty(unsigned int span, int remain, int i, int count);         /* STEP 1 */
	  void compress(char** seq, int* L);                                      /* STEP 2 */
    void get_Real_Length(int* L);                                           /* STEP 3 */
    void trans_to_1D();                                                     /* STEP 4 */
    void host_to_device(HMMER_PROFILE* HMM);                                /* STEP 5 */
    void launch_MSV(dim3 grid, dim3 block, HMMER_PROFILE* HMM, int i, int RIB);     /* STEP 6 */
    void device_to_host(HMMER_PROFILE* hmm, bool* passed_msv);                                 /* STEP 7 */
};

/* ------------------------------------------------------
 * Step 1: get_duty
 * ------------------------------------------------------
 * span:    the total number of seq this GPU need process
 * remain:  for the last GPU, it need process extra seqs
 * i:       the index of GPUs
 * count:   the total number of GPUs
 * ------------------------------------------------------
 */
void GPUs::get_duty(unsigned int span, int remain, int i, int count)
{
  if(i != (count - 1))  /* not for the last GPU */
  {
      b_point = i * span;
      e_point = i * span + span;
  }
  else           /* for the last GPU */
  {
      b_point = i * span;
      e_point = i * span + span + remain;
  }
}

/* -------------------------------------------------------------
 * Step 2: compress
 * -------------------------------------------------------------
 * Funciton:  Comress 6 chars into 1 int data type
 *            ** need release 'seq_6r', 'seq_len_6r', 'seq_6r[]'
 *            , 'seq_len_6r[]'.
 * -------------------------------------------------------------
 * seq:   2D array that cache raw seqs with 'char' data type
 * L:     1D array that cache raw length of each seq
 * -------------------------------------------------------------
 */
void GPUs::compress(char** seq, int* L)
{
  int a, b, c;      /* "c" is the length of each seq after compression */

  int s[6] = {0}; /* workers */

  int q = 0;    /* local index for this GPU */
  int i = 0;    /* global index for whole db */

  seq_6r = (int**)malloc((e_point - b_point) * sizeof(int*));
  seq_len_6r = (int*)malloc((e_point - b_point) * sizeof(int));

  for(i = b_point; i < e_point; i++)
  {
      a = (int)L[i] % 6;        /* # of rest residues in last segment */
      b = (int)L[i] / 6;        /* # of completed seqment */

      if(a != 0) {              /* has remainder */
        c = b + 1;
      }else {                   /* divisibility  */
        c = b;
      }

      seq_len_6r[q] = c;                /* length after compression */
      seq_6r[q] = (int*)malloc(c * sizeof(int));      /* allocate for compression */
      memset(seq_6r[q], 0, c * sizeof(int));          /* initial */

      /* refresh workers */
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;
      s[3] = 0;
      s[4] = 0;
      s[5] = 0;

      /* process completed elements */
      for(int j = 0; j < b; j++) {
        s[0] = (int)seq[i][6 * j + 0] << 25;
        s[1] = (int)seq[i][6 * j + 1] << 20;
        s[2] = (int)seq[i][6 * j + 2] << 15;
        s[3] = (int)seq[i][6 * j + 3] << 10;
        s[4] = (int)seq[i][6 * j + 4] << 5;
        s[5] = (int)seq[i][6 * j + 5];

        seq_6r[q][j] = s[0] | s[1] | s[2] | s[3] | s[4] | s[5];   /* j up to 'b-1'*/

      } /* end if divisibility */

      /* process the last element */
      if(a != 0) {                          /* has remainder */

        /* real residue compression */
        for(int j = 0; j < a; j++) {
          s[j] = seq[i][6 * b + j] << (25 - 5 * j);
          seq_6r[q][b] = seq_6r[q][b] | s[j];
        }

        /* wasteful residue compression */
        for(int j = 0; j < 6 - a; j++) {
          s[j] = 0x0000001F << (25 - 5 * a - 5 * j);
          seq_6r[q][b] = seq_6r[q][b] | s[j];
        }
      }

      acc += c;   /* accumulate for sum */
      q++;        /* local index + 1 */

  } /* all sequences have been compressed */

}

/* ------------------------------------------------
 * Step 3: get_Real_Length
 * ------------------------------------------------
 * L:   1D array that cache raw length of each seq
 * ------------------------------------------------
 */
void GPUs::get_Real_Length(int* L)
{
    int q = 0;  /* local index for this GPU */
    int i = 0;  /* global index for whole db */

    real_L = (int*)malloc((e_point - b_point) * sizeof(int));
    
    for(i = b_point; i < e_point; i++, q++)
    {
        real_L[q] = L[i];
    }
}

/* --------------------------------------------------
 * Step 4: trans_to_1D
 * --------------------------------------------------
 * Function:  transfer 2D to 1D array and give each
 *            seq a offset index for this 1D array.
 *            **need release 'seq_1D', 'offset'.
 * --------------------------------------------------
 */
void GPUs::trans_to_1D()
{
    unsigned int addr = 0; /* local index for this 1D array */
    int i, j;

    seq_1D = (int*)malloc(acc * sizeof(int));                                     /* ALLOC */
    offset = (unsigned int*)malloc((e_point - b_point) * sizeof(unsigned int));   /* ALLOC */

    for(i = 0; i < (e_point - b_point); i++)
    {
        for(j = 0; j < seq_len_6r[i]; j++)
        {
            seq_1D[addr] = seq_6r[i][j];
            addr += 1;
        }

        offset[i] = addr;   /* indexing within this 1D array */
    }
}

/* ------------------------------------------------
 * Step 5: host_to_device
 * ------------------------------------------------
 * Function: transfer data from host to device
 * ------------------------------------------------
 * 
 */
void GPUs::host_to_device(HMMER_PROFILE* HMM)
{
    cudaError_t r;
    /* copy 1D database to device */
    r = cudaMalloc((void**)&d_seq, acc * sizeof(int));                           CHECK(r);                                                        
    r = cudaMemcpy(d_seq, seq_1D, acc * sizeof(int), cudaMemcpyHostToDevice);    CHECK(r);

    /* copy offset of each seq to device */
    r = cudaMalloc((void**)&d_offset, (e_point - b_point) * sizeof(unsigned int));                         CHECK(r);                                                         
    r = cudaMemcpy(d_offset, offset, (e_point - b_point) * sizeof(unsigned int), cudaMemcpyHostToDevice);  CHECK(r);

    /* copy length of each seq to device */
    r = cudaMalloc((void**)&d_len, (e_point - b_point) * sizeof(int));                             CHECK(r);                                                          
    r = cudaMemcpy(d_len, real_L, (e_point - b_point) * sizeof(int), cudaMemcpyHostToDevice);      CHECK(r);

    /* for indexing each seq */
    r = cudaMalloc((void**)&d_len_6r, (e_point - b_point) * sizeof(int));                              CHECK(r);                                                          
    r = cudaMemcpy(d_len_6r, seq_len_6r, (e_point - b_point) * sizeof(int), cudaMemcpyHostToDevice);   CHECK(r);

    /* match */
    //r = cudaMalloc((void**)&mat_m, sizeof(HMM->mat_8bits));                                  CHECK(r);
    //r = cudaMemset(mat_m, 0, sizeof(HMM->mat_8bits));                                        CHECK(r);
    //r = cudaMemcpy(mat_m, HMM->mat_8bits, sizeof(HMM->mat_8bits), cudaMemcpyHostToDevice);   CHECK(r);
	r = cudaMalloc((void**)&mat_m, HMM->M * PROTEIN_TYPE * sizeof(int));                                  CHECK(r);
	r = cudaMemset(mat_m, 0, HMM->M * PROTEIN_TYPE * sizeof(int));                                        CHECK(r);
	r = cudaMemcpy(mat_m, HMM->mat_8bits, HMM->M * PROTEIN_TYPE * sizeof(int), cudaMemcpyHostToDevice);   CHECK(r);

    /**/
    r = cudaMalloc((void**)&score, (e_point - b_point) * sizeof(float));  CHECK(r);
    r = cudaMemset(score, 0, (e_point - b_point) * sizeof(float));        CHECK(r);   
}

/* -----------------------------------------------
 * Step 6: launch_MSV
 * -----------------------------------------------
 * Function:  launch MSV kernel
 * -----------------------------------------------
 * The last para is ROW_IN_BLOCK
 *
 */
void GPUs::launch_MSV(dim3 grid, dim3 block, HMMER_PROFILE* HMM, int i, int RIB)
{
    MSV_warp(grid, block, d_seq, (e_point - b_point), d_offset, score, d_len, d_len_6r, mat_m, HMM, i, RIB);
}

/* ------------------------------------------------------
 * Step 7: device_to_host
 * ------------------------------------------------------
 * Function:  copy the final results from device to host
 *            AND fprint into file..
 *            Get the number of seq can pass through
 * ------------------------------------------------------
 *
 */
void GPUs::device_to_host(HMMER_PROFILE* hmm, bool* passed_msv)
{
    float* sc = NULL;
    cudaError_t r;

    sc = (float*)malloc((e_point - b_point) * sizeof(float));
    r  = cudaMemcpy(sc, score, (e_point - b_point) * sizeof(float), cudaMemcpyDeviceToHost); CHECK(r);
    
    return seq_Pass(b_point, e_point, real_L, hmm, sc,passed_msv);
}

#ifdef __cplusplus
}
#endif

#endif /*HMMER_DEFS*/
