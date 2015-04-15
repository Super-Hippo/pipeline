#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header_def.cuh"

/*****************************************************************
 * 2. Standard iid null model ("null1")
 *****************************************************************/

/* Function:  NullOne()
 *
 * Purpose:   Calculate the null1 lod score, for sequence <dsq>
 *            of length <L> "aligned" to the base null model <bg>. 
 * 
 * Note:      Because the residue composition in null1 <bg> is the
 *            same as the background used to calculate residue
 *            scores in profiles and null models, all we have to
 *            do here is score null model transitions.
 */
int NullOne(int L, float *ret_sc)
{
  float p1 = 0.0f;
  
  p1 = (float) L / (float) (L+1);

  *ret_sc = (float) L * logf(p1) + logf(1.-p1);

  return 1;
}

/* Function:  esl_gumbel_surv()
 * Synopsis:  Returns right tail mass above $x$, $P(S > x)$.
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ for a Gumbel 
 *            (that is, 1-cdf), the right tail's probability mass.
 * 
 *            Let $y=\lambda(x-\mu)$; for 64-bit doubles, 
 *            useful dynamic range for $y$ is $-3.6 <= y <= 708$.
 *            Returns 1.0 for $y$ below lower limit, and 0.0
 *            for $y$ above upper limit.
 */
double esl_gumbel_surv(double x, double mu, double lambda)         //��Ȼ����� x �� mu �� lambda �� double�ģ� ���Ǵ�������ֵ�� float �ġ�������֪���᲻����Ӱ��
{
  double y  = lambda*(x-mu);
  double ey = -exp(-y);

  /* Use 1-e^x ~ -x approximation here when e^-y is small. */
  if (fabs(ey) < eslSMALLX1) return -ey;                          // fabs(x) �� x �ľ��ֵ
  else                       return 1 - exp(ey);
}


/* For protein models, default iid background frequencies */
/* Function:  p7_AminoFrequencies()
 *
 * Purpose:   Fills a vector <f> with amino acid background frequencies,
 *            in [A..Y] alphabetic order, same order that Easel digital
 *            alphabet uses. Caller must provide <f> allocated for at
 *            least 20 floats.
 *            
 *            These were updated 4 Sept 2007, from SwissProt 50.8,
 *            (Oct 2006), counting over 85956127 (86.0M) residues.
 *
 * Returns:   <eslOK> on success.
 */
int p7_AminoFrequencies(float *f)
{
  f[0] = 0.0787945;   /* A */
  f[1] = 0.0151600;   /* C */
  f[2] = 0.0535222;   /* D */
  f[3] = 0.0668298;   /* E */
  f[4] = 0.0397062;   /* F */
  f[5] = 0.0695071;   /* G */
  f[6] = 0.0229198;   /* H */
  f[7] = 0.0590092;   /* I */
  f[8] = 0.0594422;   /* K */
  f[9] = 0.0963728;   /* L */
  f[10]= 0.0237718;   /* M */
  f[11]= 0.0414386;   /* N */
  f[12]= 0.0482904;   /* P */
  f[13]= 0.0395639;   /* Q */
  f[14]= 0.0540978;   /* R */
  f[15]= 0.0683364;   /* S */
  f[16]= 0.0540687;   /* T */
  f[17]= 0.0673417;   /* V */
  f[18]= 0.0114135;   /* W */
  f[19]= 0.0304133;   /* Y */
  return 1;
}

/* FREE mem of HMM model */
void freeHMM(HMMER_PROFILE* hmm)
{
  free(hmm->f);         /* 20 standard residues */
  free(hmm->MU);        /* 3 floats */
  free(hmm->LAMBDA);    /* 3 floats */

   //dynamic allocation free (temprary solution)
   free(hmm->tran_16bits);
   free(hmm->mat_16bits);
  free(hmm->mat_8bits);

  int i = 0;
  for(i = 0; i < hmm->M; i++)
    free(hmm->log_tran_32bits[i]);
  for(i = 0; i < (hmm->M + 1); i++)
  {
    free(hmm->tran_32bits[i]);
    free(hmm->ins_32bits[i]);
    free(hmm->mat_32bits[i]);
  }
  free(hmm->log_tran_32bits);
  free(hmm->tran_32bits);
  free(hmm->ins_32bits);
  free(hmm->mat_32bits);

  free(hmm);            /* free 'om', not include above extra space */
                        /* only include their pointer address: 3 float pointes: MU,LAMBDA,f */
}

/* Implement roundf() function
 * Example: roundf(2.1) = 2; roundf(2.6) = 3
 *          roundf(-2.1) = -2; roundf(-2.6) = -3
 *
 * Using existing functions: floor(); ceil();
 *
 * Example: floor(2.6) = 2; ceil(2.1) = 3
 */
int round_DIY(float input)
{ 
  int r;
  r = (input - floor(input) < 0.5 ) ? floor(input) : ceil(input);
  return r;
}

/* Get P-values */
void seq_Pass(int begin, int end, int* len, HMMER_PROFILE* hmm, float* sc,bool* passed_msv)
{
  float* nullsc = NULL;
  float  seq_score = 0.0f;         /* Final score of each seq for obtain P-values   */
  double P_value = 0.0f;           /* Used to filter seqs                           */
  long int n_past_msv = 0;         /* # of seqs pass MSV filter     */
  int span = end - begin;          /* null score for each seq with different length */

  nullsc = (float*)malloc(span * sizeof(float));     /* get Null score for each seq */
  for(int i = 0; i < span; i++)
    NullOne(len[i], &nullsc[i]);              /* len has been measured in void GPUs::get_Real_Length(int* L), it is real_L of private member */

  /* get # of seq pass through filter */
  for(int i = 0; i < span; i++)
  {
    seq_score = (sc[i] - nullsc[i]) / eslCONST_LOG2;
    P_value = esl_gumbel_surv(seq_score, hmm->MU[0], hmm->LAMBDA[0]);
    if(P_value <= F1)
       {
       	passed_msv[begin + i] = 1;
       }
       else
       {
       	passed_msv[begin + i] = 0;
       }
  }
  
  free(sc);
  free (nullsc);
  nullsc = NULL;
  sc = NULL;


}

#if 0
/* MSV: final score baseline calculation (and fprintf) */
void MSV_baseline(HMMER_PROFILE* hmm, int number, int* seq_len, char** iSeq)
{
  int L = 0;                    /* cache Length for each sequence */
  int tjb = 0;                  /**/
  float scale = hmm->scale_b;   /**/
  int base = hmm->base_b;       /**/
  int bias = hmm->bias_b;       /**/
  int tbm = hmm->tbm_b;
  int tec = hmm->tec_b;
  int xB, xJ, xE;
  int ceiling = 255;            /*�����д�ȷ��*/
  int SEQ;                      /**/

  int sv;
  int mmx;

  int MMX[HMM_SIZE + 1] = {0};                /* Why '+1'? Since we need use a wasteful node 0 */

  float* SC = (float*)malloc(number * sizeof(float));
  float* NULLSC = (float*)malloc(number * sizeof(float));

  for(int b = 0; b < number; b++)         /* Outer loop: reading each sequence */
  {
    L = seq_len[b];

      tjb = unbiased_byteify(hmm, logf(3.0f / (float) (L + 3)));
      xB  = subm(base, tjb);
      xJ  = 0;                   /* re-intialize */

      for(int z = 0; z < hmm->M + 1; z++)
      {
        MMX[z] = 0;               /* Need re-intialize */
      }

      for(int i = 1; i <= L; i++)         /* medium loop: reading each residue within this sequence */
      {
        xE = 0;
        xB = subm(xB, tbm);
        SEQ = iSeq[b][i-1];           /**/

        mmx = 0;                                                /* ���mmx����SSE�е�mpv�����ÿ����0ֻ��Ϊ�˱�ʾnode0 û��ֵ */

        for(int p = 0; p < hmm->M; p++)                         /* inner loop: here is M rather than M+1, since index p of mat_8bits[] from 0 to M-1 */
        {
          sv  = max(mmx, xB);                                   //sv = _mm_max_epu8(mpv, xBv);
          sv  = addm(sv, bias);                                 //sv = _mm_adds_epu8(sv, biasv); 
          sv  = subm(sv, hmm->mat_8bits[SEQ * hmm->M + p]);     //sv = _mm_subs_epu8(sv, *rsc);   rsc++;
          xE  = max(xE, sv);                                    //xEv  = _mm_max_epu8(xEv, sv);
          mmx = MMX[p + 1];                                     //���ﵱp=1900ʱ�������һ��node����:MMX[1901]                
          MMX[p + 1] = sv;
        }

        /* test for the overflow condition */
        if( addm(xE, bias) == 255 ) 
        {
            SC[b] = 999999.0f;
            break;
        }
            
        //������special state��dummy�治ͬ����Ϊ����û���� MOVE �� LOOP �Ƚϴ�С������ֻ�� MOVE
        xE = subm(xE, tec);      //EC = EJ
        xJ = max(xJ, xE);        //�����д�� xJ = max (xJ, xE - tEJ)
        
        xB = max(base, xJ);
        xB = subm(xB, tjb);     //xB = max (base, xJ) - tJB
      } /* end loop over sequence residues 1..L */

      /* finally C->T, and add our missing precision on the NN,CC,JJ back */
      if( abs(SC[b] - 999999.0f) < 1e-6 )
      {
        /* do nothing */
      } 
      else
      {
        SC[b] = ((float) (xJ - tjb) - (float) base);
        SC[b] /= scale;
        SC[b] -= 3.0;       /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */    
      } 
  } /* end loop over sequence database */

  /**/
  int msv_counter = 0;
  float SEQ_SCORE;
  double P_VALUE;

  for(int i = 0; i < number; i++)
  {
      NullOne(seq_len[i], &NULLSC[i]);
  }

  for(int i = 0; i < number; i++)
  {
     SEQ_SCORE = (SC[i] - NULLSC[i]) / eslCONST_LOG2;
     P_VALUE = esl_gumbel_surv(SEQ_SCORE, hmm->MU[0], hmm->LAMBDA[0]);

     if(P_VALUE <= F1) {
       msv_counter++;
     }
  }

  printf("\nTotally %d sequences pass through the MSV filter", msv_counter);

  free(SC);
  free(NULLSC);
}

/* VIT: final score baseline calculation (and fprintf) */
/*
 * here, -INFINITY is -32768
 * here, OVERFLOW is 32767
 *
 *
*/
void VIT_baseline(HMMER_PROFILE* hmm, int number, int* seq_len, char** iSeq)
{
  int L = 0;                        /* cache Length for each sequence */
  float scale = hmm->scale_w;       /**/
  int base = hmm->base_w;           /**/
  int xB, xJ, xE, xC, xN;
  int Dmax;                         /* maximum D cell score on row */
  int SEQ;                          /**/
  int ceiling = 32767;              /*�����д�ȷ��*/
  int NCJ_MOVE;                     /* Varying with L */
  float pmove;                      /* temp-value */

  /* For process final scores */
  float* SC = (float*)malloc(number * sizeof(float));
  float* NULLSC = (float*)malloc(number * sizeof(float));

  int sv;
  int iv;
  int dcv;

  int mmx, imx, dmx;

  int MMX[HMM_SIZE + 1] = {-32768};                /* Why '+1'? Since we need use a wasteful node 0 */
  int IMX[HMM_SIZE + 1] = {-32768};                /* ����ǲ��Ƕ�����Ϊ M + 1 �д���             */
  int DMX[HMM_SIZE + 1] = {-32768};

  /* -------------------------------------- */
  /* outer LOOP: read whole seqs one by one */
  /* -------------------------------------- */
  for(int b = 0; b < number; b++)
  {
    L = seq_len[b];     /* length of sequence */

    /* re-calculate NCJ_MOVE for different L (seq) */
    pmove = (2.0f + hmm->nj) / ((float) L + 2.0f + hmm->nj);
    NCJ_MOVE = wordify(hmm, logf(pmove));                      /* NCJ��MOVE��һ��ֻ��һ���͹��� */

    /* Initialization... */
    for(int z = 0; z < hmm->M + 1; z++)
    {
      MMX[z] = -32768;               /* Need re-intialize for every new seq */
      IMX[z] = -32768;
      DMX[z] = -32768;
    }

    xN = base;
    xB = xN + NCJ_MOVE;
    xJ = -32768;
    xC = -32768;
    xE = -32768;

    /* --------------------------------------------- */
    /* middle LOOP: read each residue of current seq */
    /* --------------------------------------------- */

    for(int i = 1; i <= L; i++)
    {
      SEQ = iSeq[b][i-1];     /**/
      xE = -32768;            /* need refresh for next residue (row-refresh)*/
      Dmax = -32768;          /* too */
      dcv = -32768;           /* for set DMX[node 1] to -32768 */

                              /* ����MSV�����ﲻ�ö�xB��ƫ�ƴ��� */

      mmx = -32768;           /* ���mmx����SSE�е�mpv�����ÿ����0ֻ��Ϊ�˱�ʾnode0 û��ֵ */
      imx = -32768;
      dmx = -32768;

      /* -------------------------------------------------- */
      /* inner LOOP: align whole model with current residue */
      /* -------------------------------------------------- */
      for(int p = 0; p < hmm->M; p++)
      {

        /* ------------------- M state -------------------------*/
        /* previous row -> current row (state transition)*/
        /* at the beginning, mmx,imx,dmx represent M0, which are all -32768 */
        sv = clamp(addv(xB, hmm->tran_16bits[B_M * hmm->M + p]));            //sv = _mm_adds_epi16(xBv, *tsc);  tsc++;����p=0��B->M1һ����ʤ����
        sv = max(sv, clamp(addv(mmx, hmm->tran_16bits[M_M * hmm->M + p])));  //sv   = _mm_max_epi16 (sv, _mm_adds_epi16(mpv, *tsc)); tsc++;����p=0��������ʵ��wasteful����Ϊû��M0->M1��= -32768
        sv = max(sv, clamp(addv(imx, hmm->tran_16bits[I_M * hmm->M + p])));  //sv   = _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;����p=0��������ʵ��wasteful����Ϊû��I0->M1��= -32768
        sv = max(sv, clamp(addv(dmx, hmm->tran_16bits[D_M * hmm->M + p])));  //sv   = _mm_max_epi16 (sv, _mm_adds_epi16(dpv, *tsc)); tsc++;����p=0��������ʵ��wasteful����Ϊû��D0->M1��= -32768

        /* state emission */
        sv = clamp(addv(sv, hmm->mat_16bits[SEQ * hmm->M + p]));             //sv   = _mm_adds_epi16(sv, *rsc); rsc++;  rsc��ָ��rwv�ģ�rwv�����match emission score, sv is Cur-row now

        /* get xE for this row */
        xE = max(xE, sv);                                                    //xEv  = _mm_max_epi16(xEv, sv);

        /* cache pre-row + right-col for calculation of next loop (p++) */
        mmx = MMX[p + 1];                                             /* ���ﵱp=1900ʱ�������һ��node����:MMX[1901] when p = 0, mmx will be M1 (pre-row) */
        imx = IMX[p + 1];                                             /* when p = 0, imx wiil be the I1 (Pre-row) */
        dmx = DMX[p + 1];                                             /* when p = 0, dmx will be caching D1 (Pre-row), that will be used to contribute for M2 in the next loop (p = 1) */

        //mmx = MMX[p]; //MMX[0] is M1
        //imx = IMX[p]; //I1
        //dmx = DMX[p]; //D1

        /* delay update cur-row for M,D */
        MMX[p + 1] = sv;           //when p=0,this is M1           //when p=0, this MMX[1] is the node 1, (store current), MMX will be Cur-row
        //MMX[p] = sv;  //new M1

        DMX[p + 1] = dcv;          //when p=0,this is D1, D1 is not existed, so it is -inf       //DMXo(q) = dcv = -32768 when p = 0, which means DMX[node 1] = -32768; why???????wing retract?? DMX will be Cur-row
        /* NOT RIGHT! here, dcv always be -inf when p = 0 in this loop */
        //DMX[p] = dcv; //new D1

        /* ------------------- D state -------------------------*/

        /* ONLY for M->D case here */
        /* when p = 0, we need calculate D2 by using above sv (M1) since we need M1->D2 */
        dcv = clamp(addv(sv, hmm->tran_16bits[M_D * hmm->M + p]));    /* when p = 0, here, sv = M1 (Cur-row), dcv = D2, so we need address M1->D2 */
        //if(p + 2 < HMM_SIZE)
        //  DMX[p + 2] = dcv;                       /* so, begin from DMX[node 2] is no longer -32768 (�Ҳ�֪�����㲻��wing retraction) */

        /* get the maximum DMX in cur-row */
        Dmax = max(Dmax, dcv);                                        /* why do this? for LAZY F loop below */

        /* ------------------- I state -------------------------*/
        /* update Cur-row for I */
        iv = clamp(addv(mmx, hmm->tran_16bits[M_I * hmm->M + p]));            /* mmx, here maintains pre-row's M1, so sv will be I1, come from M1->I1 that is cell of p = 0 for M_I*/
        iv = max(iv, clamp(addv(imx, hmm->tran_16bits[I_I * hmm->M + p])));   /* imx, here is the pre-row's I1, so here is I1->I1 */
        IMX[p + 1] = iv;

        /* Until now, for this p, all MID[p+1] has been updated to Cur-Node.
         * [i.e. Suppose p = 0, then cur-Node 1 has been done.
         * And, for next p cycle, mmx cache pre-Node 1 to calculate cur-Node 2
         * imx is the same; dmx is the same.
         * How about dcv ? dcv now is cur_D2, and it is not used to contribute
         * to Match state in the next cycle, which is the duty of dmx. so dcv is
         * only used to update DMX to Cur_Node and iterately get Dmax for this
         *  residue(row)]
         *
         */
      }

      /* After above Inner LOOP, we only need one thing: xE.
       * xE is only for current residue, for next residue
       * we will get another new xE.
       * In a word, one residue = one row alignment = one xE = one Dmax
       */

      /* immediately detect overflow */
      if(xE >= 32767 ) {
          SC[b] = 999999.0f;
          break;
      }

      /* Special states */
      xN = xN + hmm->xw[N * XTRANS_TYPE + LOOP];
      xC = max(xC + hmm->xw[C * XTRANS_TYPE + LOOP], xE + hmm->xw[E * XTRANS_TYPE + MOVE]);
      xJ = max(xJ + hmm->xw[J * XTRANS_TYPE + LOOP], xE + hmm->xw[E * XTRANS_TYPE + LOOP]);
      xB = max(xJ + NCJ_MOVE, xN + NCJ_MOVE);

      /* --------------------------------------------------------------------------- */
      /* Lazy F: need update DMX by using dcv sometimes BEFORE next i (middle loop)  */
      /* --------------------------------------------------------------------------- */
      /* Finally the "lazy F" loop (sensu [Farrar07]). We can often
       * prove that we don't need to evaluate any D->D paths at all.
       *
       * The observation is that if we can show that on the next row,
       * B->M(i+1,k) paths always dominate M->D->...->D->M(i+1,k) paths
       * for all k, then we don't need any D->D calculations.
       *
       * The test condition is:
       *      max_k D(i,k) + max_k ( TDD(k-2) + TDM(k-1) - TBM(k) ) < xB(i)
       * So:
       *   max_k (TDD(k-2) + TDM(k-1) - TBM(k)) is precalc'ed in om->dd_bound;
       *   max_k D(i,k) is why we tracked Dmaxv;
       *   xB(i) was just calculated above.
       */

      if(Dmax + hmm->ddbound_w > xB)
      {
        sv = -32768; //for max() below;
        for(int p = 0; p < hmm->M; p++)
        {
          DMX[p + 1] = max(sv, DMX[p + 1]);   //when p = 0, this step is wasteful since sv and DMX[p+1] are both -32768
          sv = clamp(addv(DMX[p + 1], hmm->tran_16bits[D_D * hmm->M + p]));
        }
       }

    } /* end loop over sequence residues 1..L */

    /* finally C->T */
    if( abs(SC[b] - 999999.0f) < 1e-6 ) {
        /* do nothing */
    }
    else if(xC > -32768) {
        SC[b] = (float) xC + (float) NCJ_MOVE - (float) base;   /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
        SC[b] /= scale;
        SC[b] -= 3.0;                                           /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }else {
        SC[b] = -999999.0f;
    }
  } /* end loop over sequence database */

  /* ------------------------ */
  /* Final score Calculation  */
  /* ------------------------ */
  int vit_counter = 0;
  float SEQ_SCORE;
  double P_VALUE;

  for(int i = 0; i < number; i++)
  {
      NullOne(seq_len[i], &NULLSC[i]);
  }


  for(int i = 0; i < number; i++)
  {
      SEQ_SCORE = (SC[i] - NULLSC[i]) / eslCONST_LOG2;
      P_VALUE = esl_gumbel_surv(SEQ_SCORE, hmm->MU[1], hmm->LAMBDA[1]);

      if(P_VALUE <= F2) {
        vit_counter++;
      }
  }

  printf("\nTotally %d sequences pass through the VIT filter", vit_counter);

  free(SC);
  free(NULLSC);
}
#endif
