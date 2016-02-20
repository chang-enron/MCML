/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

/****
 *	THINKCPROFILER is defined to generate profiler calls in 
 *	Think C. If 1, remember to turn on "Generate profiler 
 *	calls" in the options menu. 
 ****/
#define THINKCPROFILER 0	

/* GNU cc does not support difftime() and CLOCKS_PER_SEC.*/
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif

#include "mcml.h"

/*	Declare before they are used in main(). */
FILE *GetFile(char *);
short ReadNumRuns(FILE* );
void ReadParm(FILE* , InputStruct * );
void CheckParm(FILE* , InputStruct * );
void InitOutputData(InputStruct, OutStruct *);
void FreeData(InputStruct, OutStruct *);
REAL Rspecular(LayerStruct * );
void LaunchPhoton(REAL, LayerStruct *, PhotonStruct *);
void HopDropSpin(InputStruct  *,PhotonStruct *,OutStruct *);
void SumScaleResult(InputStruct, OutStruct *);
void WriteResult(InputStruct, OutStruct, char *);
void ShowVersion(char *version);
void cl_CheckParm(FILE* File_Ptr, cl_InputStruct * In_Ptr);
void cl_ReadParm(FILE* File_Ptr, cl_InputStruct * In_Ptr);
void cl_DoOneRun(short NumRuns, cl_InputStruct *In_Ptr,
					  int set_groupitems,
					  int set_batchitems);
#define ACL_ALIGNMENT 64

void* acl_aligned_malloc (size_t size) {
    void *result = NULL;
    posix_memalign (&result, ACL_ALIGNMENT, size);
    return result;
}
void acl_aligned_free (void *ptr) {
    free (ptr);
}


time_t PunchTime(char F, char *Msg)
{
#if GNUCC
  return(0);
#else
  static clock_t ut0;	/* user time reference. */
  static time_t  rt0;	/* real time reference. */
  REAL secs;
  char s[STRLEN];
  
  if(F==0) {
    ut0 = clock();
    rt0 = time(NULL);
    return(0);
  }
  else if(F==1)  {
    secs = (clock() - ut0)/(REAL)CLOCKS_PER_SEC;
    if (secs<0) secs=0;	/* clock() can overflow. */
    sprintf(s, "User time: %8.2lf sec = %8.2lf hr.  %s\n", 
	    secs, secs/3600.0, Msg);
    puts(s);
    strcpy(Msg, s);
    return(difftime(time(NULL), rt0));
//    return(time(NULL)-rt0);
  }
  else if(F==2) return(difftime(time(NULL), rt0));
  else return(0);
#endif
}

void PredictDoneTime(long P1, long Pt)	
{
  time_t now, done_time;
  struct tm *date;
  char s[80];
  
  now = time(NULL);
  date = localtime(&now);
  strftime(s, 80, "%H:%M %x", date);
  printf("Now %s, ", s);
  
  done_time = now + 
			(time_t) (PunchTime(2,"")/(REAL)P1*(Pt-P1));
  date = localtime(&done_time);
  strftime(s, 80, "%H:%M %x", date);
  printf("End %s\n", s);
}

void ReportResult(InputStruct In_Parm, OutStruct Out_Parm)
{
  char time_report[STRLEN];
  
  strcpy(time_report, " Simulation time of this run.");
  // PunchTime(1, time_report);
//  printf ("Time duration : %8.2lf sec\n\n", PunchTime(1, time_report)/3600.0);
  printf ("Time duration : %8.2lf sec\n\n", (REAL)PunchTime(1, time_report));

  SumScaleResult(In_Parm, &Out_Parm);
  WriteResult(In_Parm, Out_Parm, time_report);
}

void GetFnameFromArgv(int argc,
					  char * argv[],
					  int *set_groupitems,
					  int *set_batchitems,
					  char * input_filename)
					  
{
    if(argc>=2) {			/* filename in command line */
  	    *set_groupitems = atoi(argv[1]);
  	    *set_batchitems = atoi(argv[2]);
        strcpy(input_filename, argv[3]);
  	}
  	else
  	{  	
        printf("command line error \n");
        exit(1);
    }
} 

int main(int argc, char *argv[]) 
{
    char input_filename[STRLEN];
    FILE *input_file_ptr;
    short num_runs;	/* number of independent runs. */
    int set_groupitems=0;
    int set_batchitems=0;


    cl_InputStruct *cl_in_parm = (cl_InputStruct *)acl_aligned_malloc(sizeof(cl_InputStruct ));      
    ShowVersion("Version 1.2.2, 2000");
    GetFnameFromArgv(argc, argv, &set_groupitems, &set_batchitems, input_filename);
    printf("set_grouptimes = %d , set_batchitems = %d \n" ,set_groupitems,set_batchitems);
    input_file_ptr = GetFile(input_filename);
    cl_CheckParm(input_file_ptr, cl_in_parm);	
    num_runs = ReadNumRuns(input_file_ptr);

    while(num_runs--){
        cl_ReadParm(input_file_ptr, cl_in_parm);
        cl_DoOneRun(num_runs, cl_in_parm, set_groupitems, set_batchitems);
    }    
    return(0);
}
