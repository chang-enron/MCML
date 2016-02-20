#include "mcml.h"
#include "CL/opencl.h"


// ACL runtime configuration
static cl_platform_id platform;
static cl_device_id device;
static cl_context context;
static cl_command_queue queue , queue_read ;
static cl_kernel kernel;
static cl_program program;
static cl_int status;
cl_event event ;
int cnt = 0 ;
void* acl_aligned_malloc (size_t size) ;
void acl_aligned_free (void *ptr) ;  

static void dump_error(const char *str, cl_int status) {
  printf("%s\n", str);
  printf("Error code: %d\n", status);
}

// free the resources allocated during initialization
static void freeResources() {
  if(kernel) 
    clReleaseKernel(kernel);  
  if(program) 
    clReleaseProgram(program);
  if(queue) 
    clReleaseCommandQueue(queue);
  if(context) 
    clReleaseContext(context);
}

int RSEED=1245;

REAL Rspecular(LayerStruct * );
time_t PunchTime(char F, char *Msg);
void cl_SumScaleResult(cl_InputStruct, OutStruct *);
void cl_WriteResult(cl_InputStruct, OutStruct, char *);
void cl_InitOutputData(cl_InputStruct, OutStruct *);
void cl_FreeData(cl_InputStruct, OutStruct *);

void collect(OutStruct *Out_Ptr, tmpOutStruct *cl_OUTstruct, int item)      
{
    int i;
    if (cl_OUTstruct->Rd_valid){
        Out_Ptr->Rd_ra[cl_OUTstruct->Rdra.x][cl_OUTstruct->Rdra.y]+=cl_OUTstruct->Rdra.w;
    }
    if (cl_OUTstruct->Tt_valid){
        Out_Ptr->Tt_ra[cl_OUTstruct->Ttra.x][cl_OUTstruct->Ttra.y]+=cl_OUTstruct->Ttra.w;        
    }
    cnt += cl_OUTstruct->total_steps;
    for (i=0; i<cl_OUTstruct->total_steps; i++)
    {
        Out_Ptr->A_rz[cl_OUTstruct->data[i].x][cl_OUTstruct->data[i].y]+=cl_OUTstruct->data[i].w;
    }
}
void cl_ReportResult(cl_InputStruct In_Parm, OutStruct Out_Parm)
{
    char time_report[STRLEN];

    strcpy(time_report, " Simulation time of this run.");
    printf ("Duration : %8.2lf sec\n\n", (REAL)PunchTime(1, time_report));
    cl_SumScaleResult(In_Parm, &Out_Parm);
    cl_WriteResult(In_Parm, Out_Parm, time_report);
}



void cl_DoOneRun(    short NumRuns, 
                    cl_InputStruct *In_Ptr, 
                    int set_groupitems,
                    int set_batchitems) 
{

    cl_uint num_platforms;
    cl_uint num_devices;
    size_t vectorSize = set_batchitems;
    size_t workSize = set_groupitems;
    cl_event rdEvent_out_parm, MkernelEvent;
    cl_ulong eventStart, eventEnd; //, totalTime=0;
    cl_ulong kerneltime=0, rdbacktime=0;
    long i_photon;    
    OutStruct out_parm;       /* distribution of photons.*/
    long num_photons = In_Ptr->num_photons; //, photon_rep=10;
    int i ;
    size_t maxCmptUnits;
    size_t Batches;
    size_t endbatch_items;

    cl_InitOutputData(*In_Ptr, &out_parm);
    out_parm.Rsp = Rspecular(In_Ptr->layerspecs); 
    printf("out_parm.Rsp = %f \n" ,out_parm.Rsp);
    
    i_photon = num_photons;
    
    Batches = num_photons/vectorSize;
    endbatch_items = num_photons - Batches * vectorSize;
    if (0 != endbatch_items)
        Batches += 1;
    
    printf("Batches = %d \n" , Batches);

    printf("photon = %d ,Wth =  %f , dz = %f , dr = %f , da = %f , nz = %d ,nr = %d , na = %d ,layer = %d \n" , 
            In_Ptr->num_photons , In_Ptr->Wth, In_Ptr->dz, In_Ptr->dr, In_Ptr->da, In_Ptr->nz, In_Ptr->nr, In_Ptr->na, In_Ptr->num_layers ) ;
    for(int i = 0 ; i < 7 ; i++)   
        printf("layers = %d , n = %f , mua = %f , mus = %f , g = %f \n" , i , In_Ptr->layerspecs[i].n, In_Ptr->layerspecs[i].mua, In_Ptr->layerspecs[i].mus, In_Ptr->layerspecs[i].g);
    
    PunchTime(0, "");
    // get the platform ID
    status = clGetPlatformIDs(1, &platform, &num_platforms);
    if(status != CL_SUCCESS) {
        dump_error("Failed clGetPlatformIDs.", status);
        freeResources();
        return;
    }
    if(num_platforms != 1) {
        printf("Found %d platforms!\n", num_platforms);
        freeResources();
        return;
    }

    // get the device ID
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, &num_devices);
    if(status != CL_SUCCESS) {
        dump_error("Failed clGetDeviceIDs.", status);
        freeResources();
        return;
    }
    printf("num_devices = %d \n" , num_devices);
    if(num_devices != 1) {
        printf("Found %d devices!\n", num_devices);
        freeResources();
        return;
    }

    // create a context
    context = clCreateContext(0, 1, &device, NULL, NULL, &status);
    if(status != CL_SUCCESS) {
        dump_error("Failed clCreateContext.", status);
        freeResources();
        return;
    }

    // create a command queue
    queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
    if(status != CL_SUCCESS) {
        dump_error("Failed clCreateCommandQueue.", status);
        freeResources();
        return;
    }
    queue_read = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
    if(status != CL_SUCCESS) {
        dump_error("Failed clCreateCommandQueue.", status);
        freeResources();
        return;
    }
    // create the kernel
    const char *kernel_name = "photon";
    cl_int kernel_status;

    // create the program using binary already compiled offline using aoc (i.e. the .aocx file)
    FILE* fp = fopen("photon.aocx", "rb");
    if (fp == NULL) {
        printf("Failed to open photon.aocx file (fopen).\n");
        exit(1) ;
    }
    fseek(fp, 0, SEEK_END);
    size_t binary_length = ftell(fp);
    unsigned char*binary = (unsigned char*) malloc(sizeof(unsigned char) * binary_length);
    assert(binary && "Malloc failed");
    rewind(fp);
    if (fread((void*)binary, binary_length, 1, fp) == 0) {
        printf("Failed to read from photon.aocx file (fread).\n");
        exit(1) ;
    }
    fclose(fp);
    program = clCreateProgramWithBinary(context, 1, &device, &binary_length, (const unsigned char**)&binary, &kernel_status, &status);
    if(status != CL_SUCCESS || kernel_status != CL_SUCCESS) {
        dump_error("Failed clCreateProgramWithBinary.", status);
        freeResources();
        return;
    }

    // build the program
    status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
    if(status != CL_SUCCESS) {
        dump_error("Failed clBuildProgram.", status);
        freeResources();
        return;
    }

    // create the kernel
    kernel = clCreateKernel(program, kernel_name, &status);
    if(status != CL_SUCCESS) {
        dump_error("Failed clCreateKernel.", status);
        freeResources();
        return;
    }
    
    int ctrlflag = 66;


    tmpOutStruct     *out       = (tmpOutStruct    *)acl_aligned_malloc(sizeof(tmpOutStruct     ) * vectorSize);
    tmpOutStruct     *outA      = (tmpOutStruct    *)acl_aligned_malloc(sizeof(tmpOutStruct     ) * vectorSize);
    tmpOutStruct     *outB      = (tmpOutStruct    *)acl_aligned_malloc(sizeof(tmpOutStruct     ) * vectorSize);
    cl_InputStruct   *in        = (cl_InputStruct  *)acl_aligned_malloc(sizeof(cl_InputStruct   ) * vectorSize);
    unsigned int     *ran_seed  = (unsigned int *)   acl_aligned_malloc(sizeof(unsigned int     ) * vectorSize); 
    // create the input buffer
    
    //cl_mem cl_flag      = clCreateBuffer(context, CL_MEM_READ_ONLY , sizeof(short         )            , NULL, &status);
    cl_mem cl_in_ptr    = clCreateBuffer(context, CL_MEM_READ_ONLY , sizeof(cl_InputStruct)            , NULL, &status);
    cl_mem cl_seed      = clCreateBuffer(context, CL_MEM_READ_ONLY , sizeof(unsigned int  )* vectorSize, NULL, &status);
    cl_mem cl_outA      = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(tmpOutStruct  )* vectorSize, NULL, &status);
    cl_mem cl_outB      = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(tmpOutStruct  )* vectorSize, NULL, &status);
    
    if(status != CL_SUCCESS) {
        dump_error("Failed CreateBuffer .", status);
        return;
    }
    
    int a = 0;       
    // set the arguments
    status = clSetKernelArg(kernel, a++, sizeof(cl_mem), (void*)&cl_in_ptr);
    if(status != CL_SUCCESS) {
        dump_error("Failed set arg 0 .", status);
        return;
    }
    status = clSetKernelArg(kernel, a++, sizeof(cl_mem), (void*)&cl_seed);
    if(status != CL_SUCCESS) {
        dump_error("Failed set arg 1.", status);
        return;
    }
    status = clSetKernelArg(kernel, a++, sizeof(REAL)  , (void *)&out_parm.Rsp);
    if(status != CL_SUCCESS) {
        dump_error("Failed set arg 2.", status);
        return;
    }
    status = clSetKernelArg(kernel, a++, sizeof(cl_mem), (void*)&cl_outA);
    if(status != CL_SUCCESS) {
        dump_error("Failed set arg 3.", status);
        return;
    }
    status = clSetKernelArg(kernel, a++, sizeof(cl_mem), (void*)&cl_outB);
    if(status != CL_SUCCESS) {
        dump_error("Failed set arg 4.", status);
        return;
    }
    status = clSetKernelArg(kernel, a++, sizeof(int), (void*)&ctrlflag);
    if(status != CL_SUCCESS) {
        dump_error("Failed set arg 5.", status);
        return;
    }
    status = clEnqueueWriteBuffer(queue, cl_in_ptr  , CL_TRUE, 0, sizeof(cl_InputStruct), In_Ptr  , 0, NULL, NULL);
    if(status != CL_SUCCESS) {
        dump_error("Failed WriteBuffer .", status);
        return;
    }
       
    for(i = 0 ; i < vectorSize ; i++)
        ran_seed[i] = lrand48();
    status = clEnqueueWriteBuffer(queue, cl_seed  , CL_TRUE, 0, sizeof(unsigned int) * vectorSize , ran_seed, 0, NULL, NULL);

    // launch kernel
    status = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &vectorSize, &workSize, 0, NULL, &MkernelEvent);
    if (status != CL_SUCCESS) {
        dump_error("Failed to launch kernel.", status);
        freeResources();
        return;
    }
    clFlush(queue);
    clFinish(queue);
    status = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_START, sizeof(cl_ulong),&eventStart,NULL);
    status = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_END,   sizeof(cl_ulong),&eventEnd,NULL);
    if(status != CL_SUCCESS) {
        dump_error("Failed to Profile Kernel Time", status);
        freeResources();
        return;
    }
    kerneltime += (eventEnd-eventStart);
    

    int BatCnt;
    for( BatCnt = 1 ; BatCnt < Batches ; BatCnt++){
        
        printf("BatCnt %d \n" , BatCnt);
        //  printf("Read Buffer \n") ;
        
        if (ctrlflag == 66) //read CL data done in flag 66
        {
            status = clEnqueueReadBuffer(queue_read, cl_outA, CL_TRUE , 0, vectorSize * sizeof(tmpOutStruct), out, 0, NULL, &rdEvent_out_parm);   

        }
        if (ctrlflag == 67) //read CL data done in flag 67
        {
            status = clEnqueueReadBuffer(queue_read, cl_outB, CL_TRUE , 0, vectorSize * sizeof(tmpOutStruct), out, 0, NULL, &rdEvent_out_parm);
        }
        

        ctrlflag = ctrlflag ^ 1;

        status = clSetKernelArg(kernel, 5, sizeof(int), (void*)&ctrlflag);
        if(status != CL_SUCCESS) {
            dump_error("Failed set arg 5.", status);
            return;
        }
        for(i = 0 ; i < vectorSize ; i++)
            ran_seed[i] = lrand48();
        status = clEnqueueWriteBuffer(queue, cl_seed  , CL_TRUE, 0, sizeof(unsigned int) * vectorSize , ran_seed, 0, NULL, NULL);
        status = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&cl_seed);
        if(status != CL_SUCCESS) {
            dump_error("Failed set arg 1.", status);
            return;
        }

        // launch kernel
        status = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &vectorSize, &workSize, 0, NULL, &MkernelEvent);
        if (status != CL_SUCCESS) {
            dump_error("Failed to launch kernel.", status);
            freeResources();
            return;
        }
        clFlush(queue);

        clFinish(queue_read);
        if (ctrlflag == 67) //read CL data done in flag 66
        {
            for (i = 0; i < vectorSize; i++)
            {
                collect(&out_parm, &out[i], i);
            }

        }
        if (ctrlflag == 66) //read CL data done in flag 67
        {
            for (i = 0; i < vectorSize; i++)
            {
                collect(&out_parm, &out[i], i);
            }

        }
        clFinish(queue);
        status = clGetEventProfilingInfo(rdEvent_out_parm,CL_PROFILING_COMMAND_START, sizeof(cl_ulong),&eventStart,NULL);
        status = clGetEventProfilingInfo(rdEvent_out_parm,CL_PROFILING_COMMAND_END,   sizeof(cl_ulong),&eventEnd,NULL);
        if(status != CL_SUCCESS) {
            dump_error("Failed to Profile Read Back Time", status);
            freeResources();
            return;
        }
        rdbacktime += (eventEnd-eventStart);  
        status = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_START , sizeof(cl_ulong),&eventStart,NULL);
        status = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_END   , sizeof(cl_ulong),&eventEnd  ,NULL);
        if(status != CL_SUCCESS) {
            dump_error("Failed to Profile Kernel Time", status);
            freeResources();
            return;
        }
        kerneltime += (eventEnd-eventStart);
    }      

    printf("====== Last Collect Data =====\n") ;
    if (ctrlflag == 66) //read CL data done in flag 66
    {
        status = clEnqueueReadBuffer(queue_read, cl_outA, CL_TRUE , 0, vectorSize * sizeof(tmpOutStruct), out, 0, NULL, &rdEvent_out_parm);   
    }

    if (ctrlflag == 67) //read CL data done in flag 67
    {
        status = clEnqueueReadBuffer(queue_read, cl_outB, CL_TRUE , 0, vectorSize * sizeof(tmpOutStruct), out, 0, NULL, &rdEvent_out_parm);
    }
    clFinish(queue_read);
    if (ctrlflag == 67) //read CL data done in flag 66
    {
        for (i = 0; i < vectorSize; i++)
        {
            collect(&out_parm, &out[i], i);
        }

    }
    if (ctrlflag == 66) //read CL data done in flag 67
    {
        for (i = 0; i < vectorSize; i++)
        {
            collect(&out_parm, &out[i], i);
        }
    }
    status = clGetEventProfilingInfo(rdEvent_out_parm,CL_PROFILING_COMMAND_START,
            sizeof(cl_ulong),&eventStart,NULL);
    status = clGetEventProfilingInfo(rdEvent_out_parm,CL_PROFILING_COMMAND_END,
            sizeof(cl_ulong),&eventEnd,NULL);
    if(status != CL_SUCCESS) {
        dump_error("Failed to Profile Read Back Time", status);
        freeResources();
        return;
    }
    printf("BatCnt = %d , ReadBack Time = %6f \n" , BatCnt - 1 , (float)(eventEnd-eventStart)/1.e9);
    rdbacktime += (eventEnd-eventStart);  

    printf("Monte Carlo time = %f sec\n", (float)kerneltime / 1.e9);
    printf("Read back time = %f sec\n", (float)rdbacktime / 1.e9);
    printf("Total Steps = %d \n" , cnt);
    cl_ReportResult(*In_Ptr, out_parm);
    cl_FreeData(*In_Ptr, &out_parm);

    freeResources();

    return ;
}

