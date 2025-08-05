/* 046267 Computer Architecture - HW #3 */
/* Implementation (skeleton)  for the dataflow statistics calculator */

#include "dflow_calc.h"

#define MAX_REGS 32

typedef struct
{
    int *instDepth, *instLat, *instDepSrc1, *instDepSrc2;
    unsigned int instCount;
} CtxObj;

ProgCtx analyzeProg(const unsigned int opsLatency[], const InstInfo progTrace[], unsigned int numOfInsts)
{
    if ((opsLatency == NULL || progTrace == NULL) && numOfInsts > 0)
    {
        return PROG_CTX_NULL; // Return null context if inputs are invalid
    }
    CtxObj *ctxObj = malloc(sizeof(CtxObj));                // Allocate memory for context object
    ctxObj->instDepth = malloc(numOfInsts * sizeof(int));   // Allocate memory for instruction depths
    ctxObj->instLat = malloc(numOfInsts * sizeof(int));     // Allocate memory for instruction latencies
    ctxObj->instDepSrc1 = malloc(numOfInsts * sizeof(int)); // Allocate memory for source 1 dependencies
    ctxObj->instDepSrc2 = malloc(numOfInsts * sizeof(int)); // Allocate memory for source 2 dependencies
    if (!ctxObj || !ctxObj->instDepth || !ctxObj->instLat || !ctxObj->instDepSrc1 || !ctxObj->instDepSrc2)
    {
        // Free allocated memory if any allocation fails
        free(ctxObj->instDepth);
        free(ctxObj->instLat);
        free(ctxObj->instDepSrc1);
        free(ctxObj->instDepSrc2);
        free(ctxObj);
        return PROG_CTX_NULL;
    }

    int regs[MAX_REGS];
    for (int i = 0; i < MAX_REGS; regs[i++] = -1) // Initialize register dependencies to -1

        ctxObj->instCount = numOfInsts;
    // Pre-fetching and analyzing instructions
    for (int i = 0; i < numOfInsts; ++i)
    {
        int src1Reg = progTrace[i].src1Idx, src2Reg = progTrace[i].src2Idx, dstReg = progTrace[i].dstIdx,
            currInstLat = opsLatency[progTrace[i].opcode], dep1InstIdx = regs[src1Reg], dep2InstIdx = regs[src2Reg],
            dep1InstDepthWithLat =
                dep1InstIdx != -1 ? ctxObj->instDepth[dep1InstIdx] + ctxObj->instLat[dep1InstIdx] : 0,
            dep2InstDepthWithLat =
                dep2InstIdx != -1 ? ctxObj->instDepth[dep2InstIdx] + ctxObj->instLat[dep2InstIdx] : 0;

        ctxObj->instLat[i] = currInstLat;     // Set current instruction latency
        ctxObj->instDepSrc1[i] = dep1InstIdx; // Set source 1 dependency
        ctxObj->instDepSrc2[i] = dep2InstIdx; // Set source 2 dependency
        ctxObj->instDepth[i] = dep1InstDepthWithLat > dep2InstDepthWithLat
                                   ? dep1InstDepthWithLat
                                   : dep2InstDepthWithLat; // Calculate instruction depth
        regs[dstReg] = i;                                  // Update the register written to by this instruction
    }
    ProgCtx handle = (ProgCtx)ctxObj; // Cast context object to ProgCtx
    return handle;
}

void freeProgCtx(ProgCtx ctx)
{
    if (!ctx)
        return;
    CtxObj *ctxObj = (CtxObj *)ctx;
    // Free allocated memory for context object
    free(ctxObj->instDepth);
    free(ctxObj->instLat);
    free(ctxObj->instDepSrc1);
    free(ctxObj->instDepSrc2);
    free(ctxObj);
}

int getInstDepth(ProgCtx ctx, unsigned int theInst)
{
    CtxObj *ctxObj = (CtxObj *)ctx;
    if (!ctxObj || theInst >= ctxObj->instCount)
        return -1;                     // Return -1 if context is invalid or instruction index is out of bounds
    return ctxObj->instDepth[theInst]; // Return the depth of the specified instruction
}

int getInstDeps(ProgCtx ctx, unsigned int theInst, int *src1DepInst, int *src2DepInst)
{
    CtxObj *ctxObj = (CtxObj *)ctx;
    if (!ctxObj || theInst >= ctxObj->instCount)
        return -1; // Return -1 if context is invalid or instruction index is out of bounds
    *src1DepInst = ctxObj->instDepSrc1[theInst]; // Set source 1 dependency
    *src2DepInst = ctxObj->instDepSrc2[theInst]; // Set source 2 dependency
    return 0;                                    // Return 0 for success
}

int getProgDepth(ProgCtx ctx)
{
    CtxObj *ctxObj = (CtxObj *)ctx;
    int instCount = ctxObj->instCount, progDepth = 0;
    // Calculate the longest execution path
    for (int i = 0; i < instCount; ++i)
    {
        if (ctxObj->instDepth[i] + ctxObj->instLat[i] > progDepth)
            progDepth = ctxObj->instDepth[i] + ctxObj->instLat[i];
    }
    return progDepth; // Return the longest execution path duration
}