/* Configure hmm model
 * includes special nodes transition
 * expon to prob from score
 * entry score
 * log-odd to score from prob
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "header_def.cuh"


/*****************************************************************
 * 1. The HMMER_PROFILE structure: a score profile.
 *****************************************************************/

/* Function:  hmmer_profile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 *
 * Purpose:   Allocate for profiles of up to <allocM> nodes for digital alphabet <abc>.
 *
 */
HMMER_PROFILE* hmmer_profile_Create(const char *MODEL_PATH)
{
  int i,j = 0;
  HMMER_PROFILE *om  = NULL;                            /* claim here for return at last */
  om = (HMMER_PROFILE*)malloc(sizeof(HMMER_PROFILE));   /* this is an important step: after this, 
                                                         * we can do any assignment works for all
                                                         * elements of om. Otherwise, we cannot do
                                                         * none of below assignment works
                                                         */
  om->M = get_hmm_size(MODEL_PATH);   //temporary used

  /* Initialization */
  //memset(om->mat_32bits, 0, sizeof(om->mat_32bits));       /* initialize */
  //memset(om->ins_32bits, 0, sizeof(om->ins_32bits));       /* we read ins but never optimize and use them */

  /* Dynamic allocation 2D array (temperary) */
  om->mat_32bits = (float**)malloc((om->M+1) * sizeof(float*));
  for(i = 0; i < (om->M + 1); i++)
  {
	  om->mat_32bits[i] = (float*)malloc(PROTEIN_TYPE * sizeof(float));
	  for(j = 0; j < PROTEIN_TYPE; j++)
		  om->mat_32bits[i][j] = 0.0f;
  }


  om->ins_32bits = (float**)malloc((om->M+1) * sizeof(float*));
  for(i = 0; i < (om->M + 1); i++)
  {
	  om->ins_32bits[i] = (float*)malloc(PROTEIN_TYPE * sizeof(float));
	  for(j = 0; j < PROTEIN_TYPE; j++)
		  om->ins_32bits[i][j] = 0.0f;    
  }


  om->tran_32bits = (float**)malloc((om->M+1) * sizeof(float*));
  for(i = 0; i < (om->M + 1); i++)
  {
	  om->tran_32bits[i] = (float*)malloc(7 * sizeof(float));
	  for(j = 0; j < 7; j++)
		  om->tran_32bits[i][j] = 0.0f;
  }


  om->log_tran_32bits = (float**)malloc((om->M) * sizeof(float*));
  for(i = 0; i < (om->M); i++)
  {
	  om->log_tran_32bits[i] = (float*)malloc(TRANS_TYPE * sizeof(float));
	  for(j = 0; j < TRANS_TYPE; j++)
		  om->log_tran_32bits[i][j] = 0.0f; 
  }

  //
  om->MU         = (float*)malloc(3  * sizeof(float));     /* respectively, for MSV, VIT and FWD */
  om->LAMBDA     = (float*)malloc(3  * sizeof(float));

  om->f          = (float*)malloc(20 * sizeof(float));     /* above 3 pointer get their own space which is not included in om's space */
  if(p7_AminoFrequencies(om->f) != 1) printf("error\n");

  memset(om->Xtran_32bits, 0, sizeof(om->Xtran_32bits)); /* initialize */

  //om->M          = 0;                                      /* number of nodes in the model */

  //memset(om->tran_32bits, 0, sizeof(om->tran_32bits));            /* Initialize */
  //memset(om->log_tran_32bits, 0, sizeof(om->log_tran_32bits));
  
  /* =============================================================== */
  /* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors */
  /* =============================================================== */

  //memset(om->mat_8bits, 0, sizeof(om->mat_8bits));         //om->mat_8bits  = (short**)malloc(21 * sizeof(short*));
  om->mat_8bits = (int*)malloc(om->M * PROTEIN_TYPE * sizeof(int));
  for(j = 0; j < om->M * PROTEIN_TYPE; j++)
	 om->mat_8bits[j] = 0; 

  /*FUCKING TEST*/
  //printf("dynamic wrong size: %d\n", sizeof(om->mat_8bits));
  //printf("dynamic correct size: %d\n", om->M * PROTEIN_TYPE * sizeof(int));
  //printf("static correct size: %d\n", sizeof(om->FUCK));
  //getchar();



  om->tbm_b      = 0;
  om->tec_b      = 0;
  om->tjb_b      = 0;

  om->scale_b    = 0;
  om->base_b     = 0;
  om->bias_b     = 0;

  om->nj         = 1.0f;          /* H3 only uses multihit local alignment ( --fs mode) */

  /* ================================================================== */
  /* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors */
  /* ================================================================== */

   /* dynamic allocation (temprary solution) */
   om->mat_16bits = (int*)malloc(om->M * PROTEIN_TYPE * sizeof(int));
   for(j = 0; j < om->M * PROTEIN_TYPE; j++)
       om->mat_16bits[j] = 0; 

   om->tran_16bits = (int*)malloc(om->M * TRANS_TYPE * sizeof(int));
   for(j = 0; j < om->M * TRANS_TYPE; j++)
       om->tran_16bits[j] = 0; 


  //memset(om->mat_16bits,  0, sizeof(om->mat_16bits));     /* Initialize */
  //memset(om->tran_16bits, 0, sizeof(om->tran_16bits));

  memset(om->xw,          0, sizeof(om->xw));
  
  om->scale_w           =   0.0f;                  /* score units: typically 500 / log(2), 1/500 bits   */
  om->base_w            =   0;                     /* offset of sword scores: typically +12000          */
  om->ddbound_w         =   0;                     /* threshold precalculated for lazy DD evaluation    */
  om->ncj_roundoff      =   0.0f;                  /* missing precision on NN,CC,JJ after rounding      */

  return om;
}


/* mf_conversion()
 * 
 * This builds the MSVFilter() parts of the profile <om>, scores
 * in lspace uchars (16-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * (A cost of +255 is our -infinity "prohibited event")
 *
 * Returns <eslOK> on success;
 * throws <eslEINVAL> if <om> hasn't been allocated properly.
 */
int mf_conversion(HMMER_PROFILE *hmm)
{
  float   max_rs = 0.0;                           /* maximum residue score: used for unsigned emission score bias */
  int i, j;                                       /* index */

  /* Since we need get maximum residue score for emission score bias,
   * all log-odd 'match' score should be prepared before here. For the 
   * 'insert' score, all of them are zero and -INFINITY, so that we dont
   * need consider them because we are finding the maxium one. It must
   * be in the range of 'match' scores. 
   */

  //�����������Ҫ����Ҳ�����ռ� rsc[]�����ǰ� match �� insert ��������ֵ�ҳ����������Ǵ������ֵ
  //ʵ��������Ϊ��ֻҪ match �Ϳ����ˣ���Ϊ insert ���� 0 ���� ���������
  //������ֵ��ʵ����ʵ�� p7_ProfileConfig �ﴦ��ģ��������ȫ�ֱ��� rsc[][][]
  //�� p7_oprofile_Convert ��ֱ����������

  /* Here, we get the 'max_rs' values */
  for(i = 0; i <= hmm->M; i++)                /* 0 to M, total M+1 nodes that include node 0, [HMM_SIZE + 1] */
  {
    for(j = 0; j < 20; j++)                   /*   for (x = 0; x < gm->abc->K; x++) { max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (M+1)*2)); } */
    {
      max_rs = max(max_rs, hmm->mat_32bits[i][j]);
    }
  }

  /* First we determine the basis for the limited-precision MSVFilter scoring system. 
   * Default: 1/3 bit units, base offset 190:  range 0..255 => -190..65 => -63.3..21.7 bits
   * See J2/66, J4/138 for analysis.
   */
  
  hmm->scale_b = 3.0 / eslCONST_LOG2;                               /* scores in units of third-bits */
  hmm->base_b  = 190;
  hmm->bias_b  = unbiased_byteify(hmm, -1.0 * max_rs);              /* here, we use 'max_rs' */

  /* Match emission cost */
  //for(i = 1; i <= hmm->M; i++)                                    /* Row: 1-M */
  //{
  //  for(j = 0; j < PROTEIN_TYPE; j++)                             /* Column: 0-20 */
  //  {
  //    hmm->mat_8bits[(i - 1) * PROTEIN_TYPE + j] = biased_byteify(hmm, hmm->mat_32bits[i][j]);     /* Convert to 1-D from 2-D, and biased_byteify */
  //  }
  //}

  /* V2: Match emission cost */
  for(i = 0; i < PROTEIN_TYPE; i++)                  
  {
    for(j = 1; j <= hmm->M; j++)
    {
      hmm->mat_8bits[i * hmm->M + (j - 1)] = biased_byteify(hmm, hmm->mat_32bits[j][i]);    /* Convert to 1-D from 2-D, and biased_byteify */
    }
  }

  /* transition costs */
  hmm->tbm_b = unbiased_byteify(hmm, logf(2.0f / ((float) hmm->M * (float) (hmm->M+1))));   /* constant B->Mk penalty        */   //�� L û��ϵ�������� kernel �м���
  hmm->tec_b = unbiased_byteify(hmm, logf(0.5f));                                           /* constant multihit E->C = E->J */   //�� L û��ϵ�������� kernel �м���

  //hmm->tjb_b = unbiased_byteify(hmm, logf(3.0f / (float) (gm->L+3))); /* this adopts the L setting of the parent profile */   //�� L �йأ���Ҫ�� kernel �� reConfig

  return 1;
}

/* unbiased_byteify()
 * Convert original transition score to a rounded uchar cost
 *-----------------------------------------------------*
 * Transition scores for MSVFilter get this treatment. *
 *-----------------------------------------------------*
 *
 * e.g. a score of -2.1, with scale 3.0, becomes a cost of 6.
 * (A cost of +255 is our -infinity "prohibited event")
 */
int unbiased_byteify(HMMER_PROFILE *hmm, float sc)
{
  int b;

  sc  = -1.0f * round_DIY(hmm->scale_b * sc);       /* ugh. sc is now an integer cost represented in a float...    */
  b   = (sc > 255.) ? 255 : (int) sc;          /* and now we cast and saturate it to an unsigned char cost... */
  return b;
}

/* biased_byteify()
 * Converts original log-odds residue score to a rounded biased uchar cost.
 *--------------------------------------------------------*
 * Match emission scores for MSVFilter get this treatment.*
 *--------------------------------------------------------*
 * e.g. a score of +3.2, with scale 3.0 and bias 12, becomes 2.
 *    3.2*3 = 9.6; rounded = 10; bias-10 = 2.
 *
 * When used, we add the bias, then subtract this cost.
 * (A cost of +255 is our -infinity "prohibited event")
 */
int biased_byteify(HMMER_PROFILE *hmm, float sc)
{
  int b;

  sc  = -1.0f * round_DIY(hmm->scale_b * sc);                            /* ugh. sc is now an integer cost represented in a float...           */
  b   = (sc > 255 - hmm->bias_b) ? 255 : (int) sc + hmm->bias_b;    /* and now we cast, saturate, and bias it to an unsigned char cost... */
  return b;
  /* there is not roundf() function in standard C library 
   * it is under amp_math.h which is a parallel library for
   * C++. So we need implement roundf() by ourselves.
   */
}


/* Transfer seq into int value for optimization in Kernel */
int seq_To_char(char** iSEQ, const char** SEQ, int* SEQ_L, int num)
{
  int i,j;

  for(i = 0; i < num; i++) {
    for(j = 0; j < SEQ_L[i]; j++) {
      iSEQ[i][j] = Amino_Offset(SEQ[i][j]);
    }
  }

  return 1;
}


/* Get the log-odd score */
int log_Odd_score(HMMER_PROFILE *hmm) 
{
  int i,j;

  /* 1. Match Emission */
  for(j = 0; j < PROTEIN_TYPE; j++)
    hmm->mat_32bits[0][j] = minusInfinity;      /* For initialize Edge values. what ever the degen is, we set the Edge -infinity */

  for(i = 1; i <= hmm->M; i++) {                         
    for(j = 0; j < 20; j++)
      hmm->mat_32bits[i][j] = logf(hmm->mat_32bits[i][j] / hmm->f[j]);
  }

  return fileOK;
} 

 /*****************************************************************
 * 5. 
 *****************************************************************/

/* Function:  p7_hmm_CalculateOccupancy()
 *
 * Purpose:   Calculate a vector <mocc[1..M]> containing probability
 *            that each match state is used in a sampled path through
 *            the model. Caller provides allocated space (<M+1> floats)
 *            for <mocc>.
 *            
 *            Caller may optionally provide an array <iocc[0..M]> as
 *            well, which (if provided) will be set to contain the
 *            expected number of times that a sampled path would contain
 *            each insert state.
 *
 * Returns:   <eslOK> on success.
 */
int p7_hmm_CalculateOccupancy(HMMER_PROFILE *hmm, float *mocc)
{
  int k;

  mocc[0] = 0.;                                                       /* no M_0 state */
  mocc[1] = hmm->tran_32bits[0][MI] + hmm->tran_32bits[0][MM];        /* initialize w/ 1 - B->D_1 */
  for (k = 2; k <= hmm->M; k++) {
    mocc[k] = mocc[k-1] * (hmm->tran_32bits[k-1][MM] + hmm->tran_32bits[k-1][MI]) + (1.0 - mocc[k-1]) * hmm->tran_32bits[k-1][DM]; 
  }
 
  return fileOK;
}

/*
* CAUTION here: after this function the BMk values are already 'log'
* They are stored in TRAN[p7P_BM][0 to M-1], totaly M values represent entry score to M different nodes
* Although in TRAN[p7P_BM][] there are M+1 space, we dont use the last space which is TRAN[p7P_BM][M]
*/
int get_entryScore(HMMER_PROFILE *hmm) 
{
  /* Local mode entry:  occ[k] /( \sum_i occ[i] * (M-i+1))
   * (Reduces to uniform 2/(M(M+1)) for occupancies of 1.0)  */
  int k;
  float Z = 0.;
  float* occ = NULL;

  occ = (float*)malloc((hmm->M + 1) * sizeof(float));
  if (p7_hmm_CalculateOccupancy(hmm, occ) != fileOK) printf("error in get Occupancy!\n");

  for (k = 1; k <= hmm->M; k++){              /* since occ[0] is 0, so we start from k = 1 */
    Z += occ[k] * (float) (hmm->M - k + 1);
  }

  for (k = 1; k <= hmm->M; k++)
  {
    hmm->log_tran_32bits[k-1][0] = logf(occ[k] / Z);      /* note off-by-one: entry at Mk stored as [k-1][BM] */
  }                                                                                      

  free(occ);
  return fileOK;  
}

/* Get the log of main transition */
int log_Trans(HMMER_PROFILE *hmm)
{
  int i;

  /* only node 1 to M-1 have meaningful values */
  for(i = 1; i < hmm->M; i++) {                                     
    hmm->log_tran_32bits[i][M_M] = logf(hmm->tran_32bits[i][MM]);  // WHY we use [XX + 1]? Because log_tran_32bits[][0] is 'BM', and [][1] is 'MM' 
    hmm->log_tran_32bits[i][I_M] = logf(hmm->tran_32bits[i][IM]);
    hmm->log_tran_32bits[i][D_M] = logf(hmm->tran_32bits[i][DM]);
  }

  /* off-by-one */
  for(i = 0; i < (hmm->M - 1); i++) {
    hmm->log_tran_32bits[i][M_D] = logf(hmm->tran_32bits[i + 1][MD]);
    hmm->log_tran_32bits[i][M_I] = logf(hmm->tran_32bits[i + 1][MI]);
    hmm->log_tran_32bits[i][D_D] = logf(hmm->tran_32bits[i + 1][DD]);
    hmm->log_tran_32bits[i][I_I] = logf(hmm->tran_32bits[i + 1][II]);    
  }

  /* SET "-infinity" to Node 0 for MM,IM,DM */
  hmm->log_tran_32bits[0][M_M] = minusInfinity;
  hmm->log_tran_32bits[0][I_M] = minusInfinity; 
  hmm->log_tran_32bits[0][D_M] = minusInfinity;
  
  /* MD, MI, II, DD */
  hmm->log_tran_32bits[hmm->M - 1][M_D] = minusInfinity;
  hmm->log_tran_32bits[hmm->M - 1][M_I] = minusInfinity;
  hmm->log_tran_32bits[hmm->M - 1][I_I] = minusInfinity;
  hmm->log_tran_32bits[hmm->M - 1][D_D] = minusInfinity;

  /* Until now, we have finished all 'log_trans' of main Nodes */
  return 1;
}

/* Get special node transition (1-D arrary with 8 values)*/
int xTrans(HMMER_PROFILE *hmm)
{
  hmm->Xtran_32bits[E * XTRANS_TYPE + LOOP] = -eslCONST_LOG2;      /* H3 is --fs mode, which is multihit local alignment */
  hmm->Xtran_32bits[E * XTRANS_TYPE + MOVE] = -eslCONST_LOG2;      /* too */

  hmm->Xtran_32bits[N * XTRANS_TYPE + LOOP] = 0;                   /* Others need reconfigure based on different 'L', here just simply initialize to 0 */
  hmm->Xtran_32bits[N * XTRANS_TYPE + MOVE] = 0;
  hmm->Xtran_32bits[J * XTRANS_TYPE + LOOP] = 0;
  hmm->Xtran_32bits[J * XTRANS_TYPE + MOVE] = 0;
  hmm->Xtran_32bits[C * XTRANS_TYPE + LOOP] = 0;
  hmm->Xtran_32bits[C * XTRANS_TYPE + MOVE] = 0;

  return 0;
}


//====================================================== Configuration for Viterbi =======================================================================//

/* wordify()
 * Converts log probability score to a rounded signed 16-bit integer cost.      //���"log probability score" �� ����� "original transition score" ���� "original log-odds residue score"
 *--------------------------------------------------------------------------    //��ʲô������ϵ�أ�����������������������
 * Both emissions and transitions for ViterbiFilter get this treatment.
 *--------------------------------------------------------------------------
 * No bias term needed, because we use signed words. 
 *   e.g. a score of +3.2, with scale 500.0, becomes +1600.
 */
int wordify(HMMER_PROFILE* hmm, float sc)
{
  sc  = round_DIY(hmm->scale_w * sc);
  if      (sc >=  32767.0) return  32767;
  else if (sc <= -32768.0) return -32768;
  else return (int) sc;
}

/* vf_conversion(): 
 * 
 * This builds the ViterbiFilter() parts of the profile <om>, scores
 * in lspace swords (8-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 */
int vf_conversion(HMMER_PROFILE* hmm)
{
  int     ddtmp;    /* used in finding worst DD transition bound                    */

  /* First set the basis for the limited-precision scoring system. 
   * Default: 1/500 bit units, base offset 12000:  range -32768..32767 => -44768..20767 => -89.54..41.53 bits
   * See J4/138 for analysis.
   */
  hmm->scale_w = 500.0 / eslCONST_LOG2;
  hmm->base_w  = 12000;

  /* Using 1-D array: Match emission scores */
  //for(int i = 1; i <= hmm->M; i++) {                                                      /* Row: 1-M */
  //  for(int j = 0; j < PROTEIN_TYPE; j++) {                                               /* Column: 0-20 */
  //    hmm->mat_16bits[(i-1) * PROTEIN_TYPE + j] = wordify(hmm, hmm->mat_32bits[i][j]);    /* Convert to 1-D from 2-D, and wordify */
  //  }
  //}

  /* V2: Match emission cost */
  for(int i = 0; i < PROTEIN_TYPE; i++)                  
  {
    for(int j = 1; j <= hmm->M; j++)
    {
      hmm->mat_16bits[i * hmm->M + (j - 1)] = wordify(hmm, hmm->mat_32bits[j][i]);    /* Convert to 1-D from 2-D, and wordify */
    }
  }

  /* Using 1-D array, we begin Wordify log-transition score */
  //for(int i = 0; i < hmm->M; i++) {
  //  /* Convert to 1-D from 2-D, and wordify */
  //  hmm->tran_16bits[i * TRANS_TYPE + 0     ] = min(wordify(hmm, hmm->log_tran_32bits[i][0]),      0);   /* BM */
  //  hmm->tran_16bits[i * TRANS_TYPE + MM + 1] = min(wordify(hmm, hmm->log_tran_32bits[i][MM + 1]), 0);   /* MM */
  //  hmm->tran_16bits[i * TRANS_TYPE + IM + 1] = min(wordify(hmm, hmm->log_tran_32bits[i][IM + 1]), 0);   /* IM */
  //  hmm->tran_16bits[i * TRANS_TYPE + DM + 1] = min(wordify(hmm, hmm->log_tran_32bits[i][DM + 1]), 0);   /* DM */
  //  hmm->tran_16bits[i * TRANS_TYPE + MD + 1] = min(wordify(hmm, hmm->log_tran_32bits[i][MD + 1]), 0);   /* MD */
  //  hmm->tran_16bits[i * TRANS_TYPE + MI + 1] = min(wordify(hmm, hmm->log_tran_32bits[i][MI + 1]), 0);   /* MI */
  //  hmm->tran_16bits[i * TRANS_TYPE + II + 1] = min(wordify(hmm, hmm->log_tran_32bits[i][II + 1]), -1);  /* II */
  //  hmm->tran_16bits[i * TRANS_TYPE + DD + 1] = wordify(hmm, hmm->log_tran_32bits[i][DD + 1]);           /* DD */
  //}

  /* V2: exchange X-Y axis, that is better for addressing in kernel */
  /* KAO, �����index��Ȼ����ˣ����꣡�װ��˷��˺þã�����������ᣬindex�Ͳ�������"B_M * TRANS_TYPE"�� */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[B_M * hmm->M + i] = min(wordify(hmm, hmm->log_tran_32bits[i][B_M]), 0);   /* BM */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[M_M * hmm->M + i] = min(wordify(hmm, hmm->log_tran_32bits[i][M_M]), 0);   /* MM */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[I_M * hmm->M + i] = min(wordify(hmm, hmm->log_tran_32bits[i][I_M]), 0);   /* IM */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[D_M * hmm->M + i] = min(wordify(hmm, hmm->log_tran_32bits[i][D_M]), 0);   /* DM */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[M_D * hmm->M + i] = min(wordify(hmm, hmm->log_tran_32bits[i][M_D]), 0);   /* MD */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[M_I * hmm->M + i] = min(wordify(hmm, hmm->log_tran_32bits[i][M_I]), 0);   /* MI */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[I_I * hmm->M + i] = min(wordify(hmm, hmm->log_tran_32bits[i][I_I]), -1);  /* II */
  for(int i = 0; i < hmm->M; i++) hmm->tran_16bits[D_D * hmm->M + i] = wordify(hmm, hmm->log_tran_32bits[i][D_D]);           /* DD */

  /* Specials. (Actually in same order in om and gm, but we copy in general form anyway.)  */
  /* VF CC,NN,JJ transitions hardcoded zero; -3.0 nat approximation used instead; this papers
   * over a length independence problem, where the approximation weirdly outperforms the
   * exact solution, probably indicating that the model's Pascal distribution is problematic,
   * and the "approximation" is in fact closer to the One True Model, the mythic H4 supermodel.
   * [xref J5/36] 
   */
  hmm->xw[E * XTRANS_TYPE + LOOP] = wordify(hmm, hmm->Xtran_32bits[E * XTRANS_TYPE + LOOP]);  /* equal to MOVE, only need transfer one of them to Kernel */

  hmm->xw[E * XTRANS_TYPE + MOVE] = wordify(hmm, hmm->Xtran_32bits[E * XTRANS_TYPE + MOVE]);

  hmm->xw[N * XTRANS_TYPE + LOOP] = 0;                                 /* was wordify(om, gm->xsc[p7P_N][p7P_LOOP]); */
  
  hmm->xw[C * XTRANS_TYPE + LOOP] = 0;                                 /* was wordify(om, gm->xsc[p7P_C][p7P_LOOP]); */

  hmm->xw[J * XTRANS_TYPE + LOOP] = 0;                                 /* was wordify(om, gm->xsc[p7P_J][p7P_LOOP]); */

  hmm->ncj_roundoff = 0.0;                                             /* goes along with NN=CC=JJ=0, -3.0 nat approximation */
                                                                       /* otherwise, would be = om->scale_w * gm->xsc[p7P_N][p7P_LOOP] -  om->xw[p7O_N][p7O_LOOP];   */
                                                                       /* see J4/150 for discussion of VF error suppression, superceded by the -3.0 nat approximation */

  /* **************************************************************************************************** */
  /* Until now, we have finished MATCH emission, main node transition, and partly special node transition */ 
  /* **************************************************************************************************** */

  /* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
  hmm->ddbound_w = -32768;

  for(int i = 1; i < (hmm->M - 1) - 1; i++)                               // �����Ҹ���log_tran_32bits�Ĵ洢�ṹ����������Ҫ�Ĳ���i=2����i=1, ������ǰŲ1
  {                                                                       // ԭ������ i = 2 to M-1, ��ΪԴ����log_tran_32bits��Ӧ�ľ�����δ����ģ�Ҳ����i = 0ʱ��DD = �����
    ddtmp  = (int) wordify(hmm, hmm->log_tran_32bits[i][D_D]);            /* DD */ //����������ˣ�i = 0 ʱ����D1->D2, (����i = 1ʱ�����D2->D3)
    ddtmp += (int) wordify(hmm, hmm->log_tran_32bits[i + 2][D_M]);        /* DM */ //�����i+1��Ϊi+2  (i=1, ȡ���� D3->M4)
    ddtmp -= (int) wordify(hmm, hmm->log_tran_32bits[i + 2][B_M]);        /* BM */ //�����i+1��Ϊi+2  (i=1, ȡ���� B->M4)
    hmm->ddbound_w = max(hmm->ddbound_w, ddtmp);
  }

  return 1;
}
