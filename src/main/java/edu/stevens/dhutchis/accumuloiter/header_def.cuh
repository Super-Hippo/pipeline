/* #### this is a .cuh file */


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
/* This is a header file
 * It includes all necessary marco value, enum numbers and function prototype
 * You can think it as a illustration file
 * 
 */

#ifndef HMMER_DEFS
#define HMMER_DEFS

#ifdef __cplusplus
extern "C" {
#endif


/* --------------------- EASY change for dead data ---------------------------*/
/**/  //#define HMM_SIZE 95                   /* hmm size (need change here for different hmm model) */
/**/  //static char *MODEL_PATH = "//Data//Hjiang//hmmer//100.hmm";
/**/  static char *DATABASE_PATH = "//home//echerin//tut//test_small.fasta";
/**/  //static char *DATABASE_PATH = "E:\\matrix\\Pfam-B.fasta\\Pfam-B.fasta";
/**/
/**/  //#define ROW_IN_BLOCK 32
/**/  //#define ACTIVE_BLOCKS 10
/**/  //#define PADD_SIZE 256
/**/  //#define NUM_GPUS 4
/**/  //#define ROW_IN_BLOCK_VIT 32     /*  For Viterbi only */
/**/ // #define ACTIVE_BLOCKS_VIT 30   
/* --------------------- EASY change for dead data ---------------------------*/


/* Define */
#define PROTEIN_TYPE 29           /* different types of protein, 20 basic + 1 degenerate (need change if there are more than 1 degen)*/
#define XTRANS_NODE 4             /* N,E,C,J*/
#define XTRANS_TYPE 2             /* Move, Loop*/

#define CHECK(res) if(res != cudaSuccess) {printf("ERROR in cuda! Now, ready to exit!"); getchar(); exit(1);}

#define fileERROR 0
#define fileOK 1

#define max(a,b) (((a) >= (b)) ? (a) : (b))
#define min(a,b) (((a) <= (b)) ? (a) : (b))

#define addm(a,b) (((a + b) > 255) ? 255 : (a + b))         /* saturated add operation for MSV */
#define subm(a,b) (((a - b) < 0) ? 0 : (a - b))             /* saturated sub operation for MSV */
//#define addv(a,b) (((a + b) > 32767) ? 32767 : (a + b))     /* saturated add operation for VIT */
//#define subv(a,b) (((a - b) < -32768) ? -32768 : (a - b))   /* saturated sub operation for VIT */
#define addv(a,b) (((a + b) < -32768) ? -32768 : (a + b))
#define clamp(a) (a > 32767) ? 32767 : a
 
#define INFINITY       999999999                    //���������
#define minusInfinity  logf(0.0f)                   //-infinity

#define eslCONST_LOG2  0.69314718055994529          //���Ӧ���� ln2

#define F1             0.02                         //���defaultֵû����Դ�����ҵ������ĵط�
#define F2             0.001
#define F3             1e-5

/* Define crossovers for numerical approximations.
 */
/* log(1+x) ~ x and  1-e^x = -x approximation.
 * Same threshold appears to be optimal for float or double x. xref STL9/138.
 */
#define eslSMALLX1    5e-9

/* INDEX ONLY for reading .hmm file */
enum r_trans {
  MM = 0, 
  MI = 1, 
  MD = 2, 
  IM = 3, 
  II = 4, 
  DM = 5, 
  DD = 6
};
#define p7P_NTRANS 7    /* MM,MI,MD,IM,II,DM,DD (order)*/

/* INDEX for configuring transiton matrix */
enum c_trans {
  B_M = 0,
  M_M = 1, 
  I_M = 2, 
  D_M = 3, 
  M_D = 4, 
  M_I = 5, 
  I_I = 6, 
  D_D = 7
};
#define TRANS_TYPE 8     /* BM,MM,DM,IM,MD,MI,II,DD (order)*/

/* Indices for special state types in the length model, gm->xsc[x][]
 */
enum p7p_xstates_e { 
  E = 0,
  N = 1,
  J = 2,
  C = 3
};
#define p7P_NXSTATES 4

/* Indices for transitions from the length modeling scores gm->xsc[][x]
 */
enum p7p_xtransitions_e {
  LOOP = 0,
  MOVE = 1
};
#define p7P_NXTRANS 2

//�����Ҫ����Ӷ�����ȥ...
typedef struct hmmer_profile_s {

  /* ============================= */
  /* Raw data read from .hmm files */
  /* ============================= */
  // float   mat_32bits[95  + 1][PROTEIN_TYPE];   /* match emission score */
  //float   ins_32bits[95  + 1][PROTEIN_TYPE];   /* insert emisssion score */
  
  //float   tran_32bits[95 + 1][7];              /* main mode transition. WHY '7'? Since there is not 'BM' */
  //float   log_tran_32bits[95][TRANS_TYPE];     /* i.e 1901*8 matrix (index manually), includes 'BM' */

  float** mat_32bits;
  float** ins_32bits;
  float** tran_32bits;
  float** log_tran_32bits;

  float   Xtran_32bits[XTRANS_NODE * XTRANS_TYPE];   /* special node transition: 4 x 2 = 8 */

  /* =============================================================== */
  /* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors */
  /* =============================================================== */

  //int   mat_8bits[95 * PROTEIN_TYPE];          /* rescaling rounding MATCH score. why change to "HMM_SIZE"��because we dont need Node 0 anymore */
  int* mat_8bits;
  //int  FUCK[95 * PROTEIN_TYPE];

  int   tbm_b;                                       /* constant B->Mk cost:    scaled log 2/M(M+1)       */
  int   tec_b;                                       /* constant E->C  cost:    scaled log 0.5            */
  int   tjb_b;                                       /* constant NCJ move cost: scaled log 3/(L+3)        */

  float scale_b;                                     /* typically 3 / log2: scores scale to 1/3 bits      */
  int   base_b;                                      /* typically +190: offset of uchar scores            */
  int   bias_b;                                      /* positive bias to emission scores, make them >=0   */

  /* ================================================================== */
  /* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors */
  /* ================================================================== */

  //int     mat_16bits[95 * PROTEIN_TYPE];     /* ����һάarray�������Kernel��Ѱַ����           */
  //int     tran_16bits[95 * TRANS_TYPE];      /* adopt 1 dimension addressing: BM,MM,IM,DM,MD,MI,II,DD */
  int*    mat_16bits;
  int*    tran_16bits;

  int     xw[XTRANS_NODE * XTRANS_TYPE];           /* NECJ state transition costs: N_move,N_loop; E;C;J */

  float   scale_w;                                 /* score units: typically 500 / log(2), 1/500 bits   */
  int     base_w;                                  /* offset of sword scores: typically +12000          */

  int     ddbound_w;                               /* threshold precalculated for lazy DD evaluation    */
  
  float   ncj_roundoff;                            /* missing precision on NN,CC,JJ after rounding      */

  /* ========================================================= */
  /* Information about current configuration, size, allocation */
  /* ========================================================= */

  int    M;            /* model length                                      */
  float  *MU;          /* parameters for getting P-value                    */
  float  *LAMBDA;      /* parameters for getting P-value                    */
  float  *f;           /* amino acid background frequencies                 */

  float  nj;           /* expected # of J's: 0 or 1, uni vs. multihit       */

} HMMER_PROFILE;


/* ------------------------------------------ Function Prototypes -------------------------------------------------------- */

/* degenerate.cpp */
int mat_degen(HMMER_PROFILE*);
void degen_enable_matrix(int**);
int Amino_Offset(char);

/* other_functions.cpp */
int NullOne(int, float *);
double esl_gumbel_surv(double, double, double);
int p7_AminoFrequencies(float*);
void freeHMM(HMMER_PROFILE*);
int round_DIY(float);
void seq_Pass(int, int, int*, HMMER_PROFILE*, float*,bool*);
void MSV_baseline(HMMER_PROFILE*, int, int*, char**);
void VIT_baseline(HMMER_PROFILE*, int, int*, char**);

/* read_files.cpp */
int get_hmm_size(const char*);
int get_Parameters(HMMER_PROFILE*, const char* );
void nextLine(FILE*, int);
void moveCursor(FILE*, int);
int get_Emission(HMMER_PROFILE*, const char*);
int get_transition(HMMER_PROFILE*, const char*);
int get_Seqnumber(const char*);                      
int get_seqID(const char*, int*);   /* Not always suitable in different database */
int alloc_Eachseq(char**, int*, int, const char*);
int fill_Eachseq(char**, int*, int, const char*);

/* model_config.cpp */
HMMER_PROFILE* hmmer_profile_Create(const char *MODEL_PATH);
int mf_conversion(HMMER_PROFILE*);
int unbiased_byteify(HMMER_PROFILE*, float);
int biased_byteify(HMMER_PROFILE*, float);
int seq_To_char(char**, const char**, int*, int);
int log_Odd_score(HMMER_PROFILE*);
int p7_hmm_CalculateOccupancy(HMMER_PROFILE*, float*);
int get_entryScore(HMMER_PROFILE*);
int log_Trans(HMMER_PROFILE*);
int xTrans(HMMER_PROFILE*);
int wordify(HMMER_PROFILE*, float);
int vf_conversion(HMMER_PROFILE*);
//void opt_seq_db(char** seq, int* L, int num, int number, int** seq_6r, int* seq_len_6r);
//unsigned long opt_seq_db(char** seq, int* L, int number, int** seq_6r, int* seq_len_6r);
unsigned long opt_seq_db(char**, int*, int**, int*, int, int);

/* device.cu */
extern void MSV_warp(dim3, dim3, int*, int, unsigned int*, float*, int*, int*, int*, HMMER_PROFILE*, int, int);

bool* mymain(int sizeout, const char** s,const char *Chmm_path );
//void MSV_warp(dim3 grid, dim3 block, int* SEQ, int total, unsigned int* OFFSET, float* SC, int* LEN, int* LEN_6R, int* MAT, HMMER_PROFILE* HMM, int device);
//void kernelTest(dim3 grid, dim3 block, int* SEQ, int total, unsigned int* OFFSET, float* SC, int* LEN, int* LEN_6R, int* MAT, HMMER_PROFILE* HMM, int* sc);
//void kernelTest_VIT(dim3 grid, dim3 block, int* SEQ, int total, float* SC, int* LEN, int* MAT, int* TRAN, HMMER_PROFILE* HMM);

//template <int HMM_SIZE, int ROW_IN_BLOCK>
//extern __global__ void MSV_TEST(int*, int, unsigned int*, float*, int*, int*, int*, int, int, int, int, float);

//void VIT_tail(dim3 grid, dim3 block, int* SEQ, int total, unsigned int* OFFSET, float* SC, int* LEN, int* MAT, int* TRAN, HMMER_PROFILE* HMM);

#ifdef __cplusplus
}
#endif

#endif /*HMMER_DEFS*/
