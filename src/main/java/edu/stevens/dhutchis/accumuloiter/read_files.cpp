/* For reading useful data from .hmm files 
 * Now only focus on the MSV filter
 *
 * Simplest way to obtain all necessary data
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "header_def.cuh"

FILE* FilePointer = NULL;

/* read length (temprary solution) */
int get_hmm_size(const char* hmm_Path)
{
  char* LEN = "LENG";
  FilePointer = fopen(hmm_Path, "r");
  char* tag_len = (char*)malloc(5 * sizeof(char));    //read width of 5 bytes
  char* lenTok = (char*)malloc(100 * sizeof(char));   //enough width of 100 bytes to ensure get accurate 'length'
  int length = 0;

  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.hmm> file\n");
    exit(-1);
  }
  else {
    //printf("Got hmm length\n");
  }

  /* get the size of hmm */
  do{
    fgets(tag_len, 5, FilePointer);
    if(strcmp(tag_len, LEN))
      nextLine(FilePointer, 1);
  }while(strcmp(tag_len, LEN) && !feof(FilePointer));

  fgets(lenTok, 100, FilePointer);
  length = (int)atof(lenTok);       //reading raw .hmm file and get the length of it

  free(tag_len);
  free(lenTok);
  fclose(FilePointer);

 // printf("/n%d/n", length);
  return length;
}

/*****************************************************************
 * 1. 
 *****************************************************************/
 
/* Function:  get_Parameters()
 * Synopsis:  Read .hmm file and get specified data.
 *
 * Purpose:  We got several important parameters here:
 *           2. hmm->MU[0,1,2]
 *           3. hmm->LAMBDA[0,1,2]      
 * Returns:   
 *
 */

int get_Parameters(HMMER_PROFILE *hmm, const char* hmm_Path)
{
  char* LEN  = "LENG";
  char* PARA_MSV = "STATS LOCAL MSV";
  char* PARA_VIT = "STATS LOCAL VITERBI";   //没做
  char* PARA_FWD = "STATS LOCAL FORWARD";   

  FilePointer = fopen(hmm_Path, "r");
  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.hmm> file\n");
    getchar();
    return fileERROR;
  } else {
    //printf(".hmm file open well\n");
  }

  /* 2. get MSV parameters */
  char* tag_MSV = (char*)malloc(16 * sizeof(char)); //CAUTION: 16 is a dead data
  char* MSVTok = (char*)malloc(100 * sizeof(char));
  char *p_msv;
  do{
    fgets(tag_MSV, 16, FilePointer);
    if(strcmp(tag_MSV, PARA_MSV))
      nextLine(FilePointer, 1);
  }while(strcmp(tag_MSV, PARA_MSV) && !feof(FilePointer));
  fgets(MSVTok, 100, FilePointer);

  p_msv = strtok(MSVTok, "  ");
  hmm->MU[0] = (float)atof(p_msv);
  p_msv = strtok(NULL, "  "); 
  hmm->LAMBDA[0] = (float)atof(p_msv);

  /* demo for extracting values from a string
   * ==============================================
   * char* test1 = (char*)malloc(100*sizeof(char));
   * fgets(test1, 100, FilePointer);
   * char *p;
   * p = strtok(test1, "  ");
   * while(p)
   * {   
   *  printf("%s\n", p);   
   *  p = strtok(NULL, "  ");   
   * }
   * ==============================================
   */

  /* 3. get Viterbi parameters */
  char* tag_VIT = (char*)malloc(20 * sizeof(char));
  char* VITTok = (char*)malloc(100 * sizeof(char));
  char *p_vit;
  do{
    fgets(tag_VIT, 20, FilePointer);
    if(strcmp(tag_VIT, PARA_VIT))
      nextLine(FilePointer, 1);
  }while(strcmp(tag_VIT, PARA_VIT) && !feof(FilePointer));
  fgets(VITTok, 100, FilePointer);

  p_vit = strtok(VITTok, "  ");
  hmm->MU[1] = (float)atof(p_vit);
  p_vit = strtok(NULL, "  "); 
  hmm->LAMBDA[1] = (float)atof(p_vit);

  /* 4. get Forward parameters */
  char* tag_FWD = (char*)malloc(20 * sizeof(char));
  char* FWDTok = (char*)malloc(100 * sizeof(char));
  char *p_fwd;
  do{
    fgets(tag_FWD, 20, FilePointer);
    if(strcmp(tag_FWD, PARA_FWD))
      nextLine(FilePointer, 1);
  }while(strcmp(tag_FWD, PARA_FWD) && !feof(FilePointer));
  fgets(FWDTok, 100, FilePointer);

  p_fwd = strtok(FWDTok, "  ");
  hmm->MU[2] = (float)atof(p_fwd);
  p_fwd = strtok(NULL, "  "); 
  hmm->LAMBDA[2] = (float)atof(p_fwd);

  //printf("MSV MU: %f, LAMBDA: %f\n", hmm->MU[0], hmm->LAMBDA[0]);
  //printf("viterbi MU: %f, LAMBDA: %f\n", hmm->MU[1], hmm->LAMBDA[1]);
  //printf("Forward MU: %f, LAMBDA: %f\n", hmm->MU[2], hmm->LAMBDA[2]);

  /* 5. free */
  fclose(FilePointer);
  free(tag_MSV);
  free(MSVTok);
  free(tag_FWD);
  free(FWDTok);

  return fileOK;
}


/* always move build-in pointer before FIRST character in next 'times'line */
void nextLine(FILE* pointer_2, int times)
{
  int moveChar;
  for(int i = 0; i < times; i++)
  {
    moveChar = 0;
    while((char)moveChar != /*CAUTION*/ '\n')       //condition is based on Windows. It will be different in Linux and Mac OS
    {     
      moveChar = fgetc(pointer_2);
    }
  }
}

/* move file curser back with 'times' chars */
void moveCursor(FILE* pointer_2, int times)
{
  for(int i = 0; i < times; i++)
  { 
    fgetc(pointer_2);
  }
}

/*****************************************************************
 * 2. 
 *****************************************************************/
 
/* Function:  get_MatchEmission()
 * Synopsis:  Read .hmm file and get specified data.
 *
 */

int get_Emission(HMMER_PROFILE *hmm, const char* hmm_Path)
{
  char* BEGIN = "  COMPO   ";
  int i = 0;
  int j = 0;
  int q = 0;

  FilePointer = fopen(hmm_Path, "r");
  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.hmm> file\n");
    getchar();
    return fileERROR;
  } else {
    printf(".hmm file open well\n");
  }

  /* 1. locate the 'COMPO' */
  char* locate = (char*)malloc(11 * sizeof(char));  //why 11? since we need use 'fgets' to search BEGIN
  do{
    fgets(locate, 11, FilePointer);
    if(strcmp(locate, BEGIN))
      nextLine(FilePointer, 1);
  }while(strcmp(locate, BEGIN) && !feof(FilePointer));

  /* 2. Process node 0 for I (no M here) */
  char* m_Temp = (char*)malloc(179 * sizeof(char)); //why 179? since each line has 178 chars
  char* match;

  char* i_Temp = (char*)malloc(179 * sizeof(char));
  char* insert;

  nextLine(FilePointer, 1);
  moveCursor(FilePointer, 10);

  fgets(i_Temp, 179, FilePointer);
  insert = strtok(i_Temp, "  ");
  while(insert) {   
    //printf("%s\n", insert);
    hmm->ins_32bits[0][i] = expf(-1.0 * (float)atof(insert));       //这个要放在这里，因为你放在p赋值语句后，最后一次循环p会被赋值NULL，这里就错了
    i++;   
    insert = strtok(NULL, "  "); 
  }

  nextLine(FilePointer, 2);
  moveCursor(FilePointer, 10);

  /* 3. Process node 1 to M (both I and M) */
  
  j = j + 1;  //above, we read Node 0 of INSERT, but we wont use them since we will set all Node 0 of MATCH and INSERT into -infinity
              //Thus, here we j++, for both INSERT and MATCH, to make start from Node 1;  

  do{
    //Match Emission
    i = 0;
    fgets(m_Temp, 179, FilePointer);          //记住用 fgets的时候，要多出一个长度 num,因为他只截取 num-1个字符（这里情况特殊多出两个无所谓）
    match = strtok(m_Temp, "  ");
    //mat[i][j] = (float)atof(p);
    while(match) {   
      //printf("%s\n", match);
      hmm->mat_32bits[j][i] = expf(-1.0 * (float)atof(match)); //这个要放在这里，因为你放在p赋值语句后，最后一次循环p会被赋值NULL，这里就错了
      i++;   
      match = strtok(NULL, "  ");   
    }

    //Insert Emission
    i = 0;
    nextLine(FilePointer, 1);
    moveCursor(FilePointer, 10);

    fgets(i_Temp, 179, FilePointer);
    insert = strtok(i_Temp, "  ");
    while(insert) {   
      //printf("%s\n", insert);
      hmm->ins_32bits[j][i] = expf(-1.0 * (float)atof(insert)); //这个要放在这里，因为你放在p赋值语句后，最后一次循环p会被赋值NULL，这里就错了
      i++;   
      insert = strtok(NULL, "  "); 
    }

    //Finish one node..
    nextLine(FilePointer, 2);
    moveCursor(FilePointer, 10);

    j++;

  }while(j <= hmm->M);      //ensure that we can fill-up Mat/Ins[HMM_SIZE + 1][PROTEIN_TYPE]

  /* 4. free */
  fclose(FilePointer);
  free(locate);
  free(m_Temp);
  free(i_Temp);

  return fileOK;
}

/*==============================================================
 *
 *==============================================================
 * This function will return transition probability instead of
 * raw data score.
 * 
 * In other words, any value will be "expf(-1.0 * raw value)"
 * to get the PROBABILTIY.
 *
 */

int get_transition(HMMER_PROFILE *hmm, const char* hmm_Path) 
{
  int i;

  char* BEGIN = "  COMPO   ";
  FilePointer = fopen(hmm_Path, "r");
  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.hmm> file\n");
    getchar();
    return fileERROR;
  } else {
    printf(".hmm file open well_trans\n");
  }

  /* 1. locate the line of 1st stage */
  char* locate = (char*)malloc(11 * sizeof(char));  //why 11? since we need use 'fgets' to search BEGIN
  do{
    fgets(locate, 11, FilePointer);
    if(strcmp(locate, BEGIN))
      nextLine(FilePointer, 1);
  }while(strcmp(locate, BEGIN) && !feof(FilePointer));

  nextLine(FilePointer, 2);                         //move to the 'transition' line...

  char* t_Temp = (char*)malloc(72 * sizeof(char));  //why 72? since we need use 'fgets' to get TRAN line
  char* p;

  /* 2. Process node 0 */
  fgets(t_Temp, 72, FilePointer);

  p = strtok(t_Temp, "  ");
  hmm->tran_32bits[0][MM] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");
  hmm->tran_32bits[0][MI] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");
  hmm->tran_32bits[0][MD] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");
  hmm->tran_32bits[0][IM] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");
  hmm->tran_32bits[0][II] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");
  hmm->tran_32bits[0][DM] = 1.0f;

  hmm->tran_32bits[0][DD] = 0;

  nextLine(FilePointer, 3);             //move to next 'transition' line... 

  /* 3. Process node 1 to M-1 */
  for(i = 1; i < hmm->M; i++) {         //不去用 i = 0 (B) 和 i = length + 1 (E), 其中 E 的空间应该是永远都用不上的 

    fgets(t_Temp, 72, FilePointer);     //需要测试一下：p 设置在rolling外是不是可以？ 答案是肯定的，可以！！下一次用fgets时会覆盖掉原来的 72 个空间

    p = strtok(t_Temp, "  ");
    hmm->tran_32bits[i][MM] = expf(-1.0 * (float)atof(p));
    p = strtok(NULL, "  ");
    hmm->tran_32bits[i][MI] = expf(-1.0 * (float)atof(p));
    p = strtok(NULL, "  ");
    hmm->tran_32bits[i][MD] = expf(-1.0 * (float)atof(p));
    p = strtok(NULL, "  ");
    hmm->tran_32bits[i][IM] = expf(-1.0 * (float)atof(p));
    p = strtok(NULL, "  ");
    hmm->tran_32bits[i][II] = expf(-1.0 * (float)atof(p));
    p = strtok(NULL, "  ");
    hmm->tran_32bits[i][DM] = expf(-1.0 * (float)atof(p));
    p = strtok(NULL, "  ");
    hmm->tran_32bits[i][DD] = expf(-1.0 * (float)atof(p));

    nextLine(FilePointer, 3); //move to next 'transition' line...
  }

  /* 4. Process node M */
  fgets(t_Temp, 72, FilePointer);

  p = strtok(t_Temp, "  ");
  hmm->tran_32bits[hmm->M][MM] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");
  hmm->tran_32bits[hmm->M][MI] = expf(-1.0 * (float)atof(p));

  hmm->tran_32bits[hmm->M][MD] = expf(-1.0 * (float)atof(p));

  p = strtok(NULL, "        *  ");
  hmm->tran_32bits[hmm->M][IM] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");
  hmm->tran_32bits[hmm->M][II] = expf(-1.0 * (float)atof(p));
  p = strtok(NULL, "  ");

  hmm->tran_32bits[hmm->M][DM] = 1.0f;
    
  hmm->tran_32bits[hmm->M][DD] = 0;

  /* 3. free memory */
  fclose(FilePointer);
  free(locate);
  free(t_Temp);

  return fileOK;
}

/*****************************************************************
 * 3. 
 *****************************************************************/
 
/* Function:  get_Seqnumber()
 */

int get_Seqnumber(char* seq_Path)
{
  int fileChar = 0;
  int times = 0;

  FilePointer = fopen(seq_Path, "r");
  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.fsa> file\n");
    getchar();
    return fileERROR;
  } else {
    printf(".fsa file open well\n");
  }

  fileChar = fgetc(FilePointer);
  while(!feof(FilePointer)) {             //SAME AS： while(fileChar != EOF)
    if((char)fileChar == '>') {
      times++;
      nextLine(FilePointer, 1);
    }
    fileChar = fgetc(FilePointer);          //this is important because u need 'fgetc()' make build-in index move to the next
  }

  fclose(FilePointer);

  return times;
}

/* get the ID for each target sequence */
int get_seqID(char* seq_Path, int* ID_)
{
  int fileChar = 0;
  int count = 0;
  char* id_Temp = (char*)malloc(19 * sizeof(char)); //这里的数字 19 是根据 .fsa 文件中的格式得出的，如果 ID 长度有变，这个数值也需要改变
  char* id;

  FilePointer = fopen(seq_Path, "r");
  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.fsa> file\n");
    getchar();
    return fileERROR;
  } else {
    printf(".fsa file open well\n");
  }

  fileChar = fgetc(FilePointer);
  while(!feof(FilePointer)) {             //SAME AS： while(fileChar != EOF)
    if((char)fileChar == '>') {

      fgets(id_Temp, 19, FilePointer);
      /**********************************/    //这里重复两次是因为 .fsa 格式的原因
      /**/ id = strtok(id_Temp, "|"); /**/
      /**/ id = strtok(NULL, "|");    /**/
      /**/ id = strtok(NULL, "|");    /**/
      /**********************************/
      ID_[count] = (int)atof(id);
      //printf("\n%d", ID_[count]);
      count++;
    }
    fileChar = fgetc(FilePointer);          //this is important because u need 'fgetc()' make build-in index move to the next
  }

  fclose(FilePointer);
  free(id_Temp);

  return fileOK;
}

/* dynamic memory allocation for each sequence on host */
int alloc_Eachseq(char** addr, int* len_, int num, char* seq_Path)
{
  int fileChar = 0;
  int index = 0;

  FilePointer = fopen(seq_Path, "r");
  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.fsa> file\n");
    getchar();
    return fileERROR;
  } else {
    printf(".fsa file open well\n");
  }

  fileChar = fgetc(FilePointer);            //为什么fileChat要在while循环做一次，因为他是判断 == '>'，fileChat_in不做是因为他是判断 != '>'
  while(!feof(FilePointer)) {
    if((char)fileChar == '>') {           //fileChat的值从第一个‘>’进来后，就没有变过，多以每次都是满足条件进入if体
      nextLine(FilePointer, 1);
      int count = 0;
      int fileChar_in = 0;

      while((char)fileChar_in != '>' && fileChar_in != EOF) { //  here,it has already ensured security on the end of file
        fileChar_in = fgetc(FilePointer);
        if(fileChar_in >= 65 && fileChar_in <= 90) {
          count++;
        }
      }

      if(index < num) {                   //'if' condition can be removed actually.array cannot overflow.
        addr[index] = (char*)malloc(count*sizeof(char));  //allocate memory for this sequence.
        len_[index] = count;                //length of this sequence.

        index++;
      }
    }
  }

  fclose(FilePointer);

  return fileOK;
}

/* cache each sequence into corresponding position */
int fill_Eachseq(char** addr, int* len_, int num, char* seq_Path)
{
  int fileChar = 0;
  int index = 0;

  FilePointer = fopen(seq_Path, "r");
  if(FilePointer == NULL) {
    printf("Fatal error: Cannot open or find <.fsa> file\n");
    getchar();
    return fileERROR;
  } else {
    printf(".fsa file open well\n");
  }

  fileChar = fgetc(FilePointer);              //为什么fileChat要在while循环做一次，因为他是判断 == '>'，fileChat_in不做是因为他是判断 != '>'
  while(!feof(FilePointer)) {
    if((char)fileChar == '>') {               //fileChat的值从第一个‘>’进来后，就没有变过，多以每次都是满足条件进入if体
      nextLine(FilePointer, 1);
      int i = 0;
      int fileChar_in = 0;

      while((char)fileChar_in != '>' && fileChar_in != EOF) {         //  here,it has already ensured security on the end of file
        fileChar_in = fgetc(FilePointer);
        if(fileChar_in >= 65 && fileChar_in <= 90 && i < len_[index]) {
          addr[index][i] = (char)fileChar_in;
          i++;
        }
      }

      if(index < num) {                           //'if' condition can be removed actually.array cannot overflow.
        index++;
      }
    }
  }

  fclose(FilePointer);
  
  return fileOK;
}