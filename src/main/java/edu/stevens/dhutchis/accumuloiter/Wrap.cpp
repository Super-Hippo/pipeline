#include "Wrap.h"
#include <stdlib.h>
#include <string.h>
#include "header_def.cuh"

#ifdef __cplusplus
extern "C" {
#endif

JNIEXPORT jbooleanArray JNICALL Java_edu_stevens_dhutchis_accumuloiter_Wrap_seqpass(JNIEnv *env, jobject obj,jobjectArray array , jstring hmm_path)
{
	
	jint sizeout = env->GetArrayLength(array);
	jboolean carray[sizeout];
	bool* narray;
	const char *stringholder[sizeout];
	
	jint i;
	i=0;
	for(i = 0; i < sizeout;i++)//there is no standard function to convert java array of strings(which are seen as objects) to a native array of strings
	{

		 jstring str = (jstring) env->GetObjectArrayElement( array, i);//convert a jobject in an array to a jstring
		 stringholder[i] = env->GetStringUTFChars( str, 0);		
		 printf("\nresult %i is: %s\n", i,   stringholder[i]);
	}	
	const char *Chmm_path = env->GetStringUTFChars( hmm_path, 0);		
	
	printf("					I am in jni and about to go to cuda, path is: %s", Chmm_path);
	
	narray = mymain((int)sizeout, stringholder, Chmm_path); // DH->ERIC: pass sizeout as an int, stringholder as a char**
	
	//Release the hmm paths, we are done with them
	env->ReleaseStringUTFChars( hmm_path, Chmm_path);		
	
	i=0;
	for(i = 0; i < sizeout;i++)//cast all native booleans to jboolean
	{
		carray[i] = (jboolean)narray[i];
		
		//printf("\nresult %i is: %d", i,  narray[i]);
	}
	
	//each element in the array is true if the sequence made it through filter
	jbooleanArray boolholder = env->NewBooleanArray(sizeout);  
	env->SetBooleanArrayRegion( boolholder, 0, sizeout, carray);
	
	// DH->ERIC: Free the manually allocated memory once finished.
	//free(carray);
	//free(i);
	printf("\n leaving cpp");
	return boolholder;
}

JNIEXPORT jbooleanArray JNICALL Java_Wrap_seqpass(JNIEnv *env, jobject obj,jobjectArray array , jstring hmm_path)
{
	
	jint sizeout = env->GetArrayLength(array);
	jboolean carray[sizeout];
	bool* narray;
	const char *stringholder[sizeout];
	
	jint i;
	for(i = 0; i < sizeout;i++)//there is no standard function to convert java array of strings(which are seen as objects) to a native array of strings
	{
		 jstring str = (jstring) env->GetObjectArrayElement( array, i);//convert a jobject in an array to a jstring
		 stringholder[i] = env->GetStringUTFChars( str, 0);		
	}	
	const char *Chmm_path = env->GetStringUTFChars( hmm_path, 0);		
	
	//printf("I am in jni and about to go to cuda, path is: %s", Chmm_path);
	
	narray = mymain((int)sizeout, stringholder, Chmm_path); // DH->ERIC: pass sizeout as an int, stringholder as a char**
	
	//Release the hmm paths, we are done with them
	env->ReleaseStringUTFChars( hmm_path, Chmm_path);		
	
	for(i = 0; i < sizeout;i++)//cast all native booleans to jboolean
	{
		carray[i] = (jboolean)narray[i];
		
		//printf("\nresult %i is: %d", i,  narray[i]);
	}
	
	//each element in the array is true if the sequence made it through filter
	jbooleanArray boolholder = env->NewBooleanArray(sizeout);  
	env->SetBooleanArrayRegion( boolholder, 0, sizeout, carray);
	
	// DH->ERIC: Free the manually allocated memory once finished.
	//free(carray);
//printf("\n leaving cpp");
	return boolholder;
}

#ifdef __cplusplus
}
#endif

