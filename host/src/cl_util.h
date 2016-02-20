#include <stdio.h>
#include <stdlib.h>
#include "CL/cl.h"


extern int cl_errChk(cl_int ret, char *errinfo);
extern int QueryHWinfo(size_t *maxCmptUnits);
extern int GetHW();
extern int input_judge(int *inNum, int MaxNum, char *SelMsg);
extern int cl_OperateSelect(int *plt_sel, int *dev_sel, int Ntraces,
                     float mean, float stdv,
                     size_t *local_items_inAwgroup, size_t *Kernel_WorkNum);

extern int SetCont(int plt_sel, int dev_sel, cl_context *cont);

extern int LoadKernel(char *cl_filename, char **source_str, size_t *source_size);
extern int SetProgKernel(cl_program *prog, cl_kernel *ker, cl_context context, char *source_str, size_t source_size,
                  cl_device_id* Devices, int dev_sel, char *kername);

extern cl_uint numPlatforms;
extern cl_uint numDevices;
extern cl_platform_id *platforms;
extern cl_device_id* Devices;
