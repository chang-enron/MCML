#include <stdio.h>
#include <stdlib.h>
#include "CL/cl.h"
//#include "cl_trace.h"
//#include "oclUtils.h"
//#include "/root/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc/oclUtils.h"

#define MAX_SOURCE_SIZE (0x100000)

cl_uint numPlatforms;
cl_uint numDevices;
cl_platform_id *platforms = NULL;
cl_device_id* Devices = NULL;

int cl_errChk(cl_int ret, char *errinfo)
{
     if (ret != CL_SUCCESS) 
    {
        printf("%s\n", errinfo);
        printf("Error code = %d\n", ret);
        exit(ret);
    }
    return 0;
}   

int GetHW()
{
    char local_plat_buf[100];
    int i;
    
  // Get & Set OpenCL Platforms
    // get Platform numbers
    cl_int ret = clGetPlatformIDs(1, NULL, &numPlatforms);
    cl_errChk(ret,"Error 0>> clGetPlatformIDs");
    
    // get memory to store platform IDs
    platforms = (cl_platform_id*)malloc(numPlatforms *
                                    sizeof(cl_platform_id));
    // store IDs into memory
    ret = clGetPlatformIDs(numPlatforms, platforms, NULL);
    cl_errChk(ret,"Error 1>> clGetPlatformIDs");
	
// Get OpenCL Platforms & Devices Info.
    for (i = 0; i < numPlatforms; i++)
    {
    // Get Platform Info.
        ret = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
              sizeof(local_plat_buf), local_plat_buf, NULL);
        cl_errChk(ret,"Error >> clGetPlatformInfo");
                
        // get Devices numbers
        ret = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 
          0, NULL, &numDevices);
        cl_errChk(ret,"Error >> clGetDeviceIDs");
        
        // get memory to store device IDs
        Devices = (cl_device_id*)malloc(sizeof(cl_device_id)* numDevices);
        if (numDevices == 0)
        {
            printf("!! There is no device in platform #%d\n", i);
            exit(0);
        }
        else
        {
            ret = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 
                  numDevices, Devices, NULL);
        }
    }
    return 0;	
}

int QueryHWinfo(size_t *maxCmptUnits)
{
    cl_ulong globalmemSize, localmemSize, maxConstBufSize;
    size_t maxWGroupSize;
    size_t maxWIdims;
    size_t maxWItemSize3D[3];
    char device_str[100];
    char local_plat_buf[100];
    char local_dev_buf[100];
    int i;
     
  // Get & Set OpenCL Platforms
    // get Platform numbers
    cl_int ret = clGetPlatformIDs(1, NULL, &numPlatforms);
    cl_errChk(ret,"Error 0>> clGetPlatformIDs");
    
    printf(">> Get Platform num = %d\n\n", numPlatforms);
    
    // get memory to store platform IDs
    platforms = (cl_platform_id*)malloc(numPlatforms * sizeof(cl_platform_id));
    // store IDs into memory
    ret = clGetPlatformIDs(numPlatforms, platforms, NULL);
    cl_errChk(ret,"Error 1>> clGetPlatformIDs");
	
// Get OpenCL Platforms & Devices Info.
    for (i = 0; i < numPlatforms; i++)
    {
    // Get Platform Info.
        ret = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
              sizeof(local_plat_buf), local_plat_buf, NULL);
        cl_errChk(ret,"Error >> clGetPlatformInfo");
        // Vendor Info.
        printf(">> Platform #%d: Vendor => %s\n", i, local_plat_buf);
                
        // get Devices numbers
        ret = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 
          0, NULL, &numDevices);
        cl_errChk(ret,"Error >> clGetDeviceIDs");
        
        // get memory to store device IDs
        Devices = (cl_device_id*)malloc(sizeof(cl_device_id)* numDevices);
        if (numDevices == 0)
        {
            printf("!! There is no device in platform #%d\n", i);
            exit(0);
        }
        else
        {
            ret = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 
                  numDevices, Devices, NULL);
            printf(">> %d Device(s) in platform #%d\n", numDevices, i);
        }

        // Get Devices info.
        int j = 0;
       // cl_device_svm_capabilities caps;
        
        for (j=0; j< numDevices; j++)
        {
            printf("\n>> [ Device: %d ]\n", j);

         /*   // Get SVM support
            ret = clGetDeviceInfo(Devices[j], CL_DEVICE_SVM_CAPABILITIES,
                     sizeof(caps), &caps, 0);
            cl_errChk(ret,"Error >> clGetDeviceInfo_dev_svm");
              printf("\t>> SVM Capabilities:\n");
            if (ret == CL_SUCCESS){
              if (caps & CL_DEVICE_SVM_COARSE_GRAIN_BUFFER)
                 printf("\t\t>> CL_DEVICE_SVM_COARSE_GRAIN_BUFFER\n");
              if (caps & CL_DEVICE_SVM_FINE_GRAIN_BUFFER)
                 printf("\t\t>> CL_DEVICE_SVM_FINE_GRAIN_BUFFER\n");
              if (caps & CL_DEVICE_SVM_FINE_GRAIN_SYSTEM)
                 printf("\t\t>> CL_DEVICE_SVM_FINE_GRAIN_SYSTEM\n");
              if (caps & CL_DEVICE_SVM_ATOMICS)
                 printf("\t\t>> CL_DEVICE_SVM_ATOMICS\n");
             }
         */
            // Get Vendor info.
            ret = clGetDeviceInfo(Devices[j], CL_DEVICE_VENDOR, 
                    sizeof(device_str), device_str, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_dev_vendor");
            printf("\t>> Vendor: %s\n", device_str);

            // Get Name info.
            ret = clGetDeviceInfo(Devices[j], CL_DEVICE_NAME, 
                    sizeof(local_dev_buf), local_dev_buf, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_dev_name");
            printf("\t>> Model: %s\n", local_dev_buf);

            // Get Max Work Group Size
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_MAX_WORK_GROUP_SIZE, 
                    sizeof(maxWGroupSize), &maxWGroupSize, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_maxWGroupSize");
            printf("\t>> CL_DEVICE_MAX_WORK_GROUP_SIZE (WIs/WG): %d\n", (int)maxWGroupSize);

            // Get Max Compute Units Size
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_MAX_COMPUTE_UNITS, 
                    sizeof(*maxCmptUnits), maxCmptUnits, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_maxCmptUnits");
            printf("\t>> CL_DEVICE_MAX_COMPUTE_UNITS : %d\n", (int)*maxCmptUnits);

            // Get Max WORK_ITEM_DIMENSIONS
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, 
                    sizeof(maxWIdims), &maxWIdims, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_maxWorkItemD");
            printf("\t>> CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: %d\n", (int)maxWIdims);

            // Get Max WORK_ITEM_SIZES
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_MAX_WORK_ITEM_SIZES, 
                    sizeof(maxWItemSize3D), &maxWItemSize3D, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_maxWItemSize3D");
            printf("\t>> CL_DEVICE_MAX_WORK_ITEM_SIZES: %d, %d, %d\n", 
            (int)maxWItemSize3D[0], (int)maxWItemSize3D[1], (int)maxWItemSize3D[2]);

            // Get GLOBAL_MEM_SIZE
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_GLOBAL_MEM_SIZE, 
                    sizeof(globalmemSize), &globalmemSize, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_globalmemSize");
            printf("\t>> CL_DEVICE_GLOBAL_MEM_SIZE(B): %.1f\n", 
                (float)globalmemSize);

            // Get MAX_CONSTANT_BUFFER_SIZE
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, 
                    sizeof(maxConstBufSize), &maxConstBufSize, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_maxConstBufSize");
            printf("\t>> CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE(B): %.1f\n", 
                (float)maxConstBufSize);

            // Get LOCAL_MEM_SIZE
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_LOCAL_MEM_SIZE, 
                    sizeof(localmemSize), &localmemSize, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_localmemSize");
            printf("\t>> CL_DEVICE_LOCAL_MEM_SIZE(B): %.1f\n", 
                (float)localmemSize);
                
            // Get CL_DEVICE_MAX_CLOCK_FREQUENCY 
            ret = clGetDeviceInfo(Devices[j], 
                    CL_DEVICE_MAX_CLOCK_FREQUENCY , 
                    sizeof(localmemSize), &localmemSize, NULL);
            cl_errChk(ret,"Error >> clGetDeviceInfo_MAX_CLOCK_FREQUENCY");
            printf("\t>> CL_DEVICE_MAX_CLOCK_FREQUENCY (MHz): %lu\n", 
                localmemSize);
        }
    }
    return 0;	
}

int input_judge(int *inNum, int MaxNum, char *SelMsg)
{
    int sel_true = 0;
    while(sel_true != 1)
    {
        printf("%s\n", SelMsg);
	      scanf("%d", inNum);
        printf("\n");
        if (*inNum < MaxNum)
        {
            sel_true = 1;
        }
        else
        {
            printf("Wrong selection, please select again. \n");
            sel_true = 0;
        }
    }
    return 0;
}
	
int cl_OperateSelect(int *plt_sel, int *dev_sel, int Ntraces,
                     float mean, float stdv,
                     size_t *local_items_inAwgroup, size_t *Kernel_WorkNum)
{

    input_judge(plt_sel, numPlatforms, "Select Target Platform : ");
    input_judge(dev_sel, numDevices, "Select Target Device : ");
    
    printf("The TOTAL Traces Number is ");
    printf("%d", Ntraces);
    printf("\n");
    
    printf("The Iccident Engergy (MeV) is ");
    printf("%f\n", mean);
    printf("Deviation (MeV) is ");
    printf("%f\n\n", stdv);
//    input_judge(mev, MaxEnergy, "Input Iccident Engergy (MeV) (Max. 200): ");
    
//    printf("Input Random Seed : ");
//    scanf("%d", seed);
//    printf("\n");
    
    printf("Input (local items/work group) : ");
    scanf("%zd", local_items_inAwgroup);
    printf("\n");
    
    printf("Input Traces/Batch (Work Items/Kernel) : ");
    scanf("%zd", Kernel_WorkNum);
    printf("\n");	
    
    return 0;
}

int SetCont(int plt_sel, int dev_sel, cl_context *cont)
{
    // check if the environment is valid
    // gain the device numbers
    cl_int ret = clGetDeviceIDs(platforms[plt_sel], CL_DEVICE_TYPE_ALL, 
            0, NULL, &numDevices);
    cl_errChk(ret,"Error >> clGetDeviceIDs");
    if (dev_sel > numDevices)
    {
        printf("Invalid Device Number\n");
        exit(1);
    }
    // get the select platform & device
    ret = clGetDeviceIDs(platforms[plt_sel], CL_DEVICE_TYPE_ALL,
             numDevices, Devices, NULL);

    // check the device is a CPU or GPU
    cl_device_type DeviceTyep;
    //cl_device_id DeviceID;
    //DeviceID = Devices[dev_sel];
    ret = clGetDeviceInfo(Devices[dev_sel], CL_DEVICE_TYPE,
            sizeof(DeviceTyep),(void *)&DeviceTyep,NULL);
    cl_errChk(ret,"Error >> clGetDeviceInfo_DeviceTyep");

    if(DeviceTyep == CL_DEVICE_TYPE_GPU) 
        printf("Creating GPU Context\n");
    else if (DeviceTyep == CL_DEVICE_TYPE_CPU) 
        printf("Creating CPU Context\n");
    else
        printf("This Context Type not Supported.\n");

    // Create a context
    cl_context_properties prop[3] = {CL_CONTEXT_PLATFORM, 
                    (cl_context_properties) platforms[plt_sel], 0 };
    *cont = clCreateContextFromType(prop, 
                    (cl_device_type)DeviceTyep, NULL, NULL, &ret);

    if (*cont == 0)
    {
        printf("Cannot create OpenCL context\n");
        return 0;
    }    
    
	  return 0;
}

int LoadKernel(char *cl_filename, char **source_str, size_t *source_size)
{
    FILE *fp;
    fp = fopen(cl_filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }

    *source_str = (char*)malloc(MAX_SOURCE_SIZE);
    *source_size = fread(*source_str, 1, MAX_SOURCE_SIZE, fp);	
    fclose(fp);
    
    return 0;
}

int SetProgKernel(cl_program *prog, cl_kernel *ker, cl_context context, char *source_str, size_t source_size,
                  cl_device_id* Devices, int dev_sel, char *kername)
{
	  cl_int ret;
    // Create a program form the Kernel source (string from .cl)
    *prog = clCreateProgramWithSource(context, 1,
            (const char **)&source_str, (const size_t *)&source_size, &ret);
    cl_errChk(ret,"Error >> clCreateProgramWithSource");

    // Build Program
    // ret = clBuildProgram(*prog, dev_sel, Devices, NULL, NULL, NULL);
    ret = clBuildProgram(*prog, 0, NULL, NULL, NULL, NULL);
    //ret = clBuildProgram(*prog, 0, NULL, "-cl-std=CL2.0", NULL, NULL);
    if (ret != CL_SUCCESS) 
    {
        size_t len;
        char buffer[2048];
        clGetProgramBuildInfo(*prog, Devices[dev_sel], CL_PROGRAM_BUILD_LOG,
                              sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        printf("Error Build Prog.\n");
        printf("Error code = %d\n", ret);
        exit(ret);
    }
    // Create the OpenCL kernel
    *ker = clCreateKernel(*prog, kername, &ret);
    cl_errChk(ret,"Error >> clCreateKernel");
	  return 0;
}

