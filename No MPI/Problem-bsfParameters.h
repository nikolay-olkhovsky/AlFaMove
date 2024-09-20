/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (No MPI)
Module: Problem-bsfParameters.h (BSF-skeleton Parameters)
Prefix: PP_BSF
Authors: Nikolay A. Olkhovsky & Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/

//=========================== BSF Skeleton Parameters =========================
#define PP_BSF_PRECISION (PP_SETW/2)// Decimal precision on output
#define PP_BSF_MAX_MPI_SIZE 400		// Maximal MPI Size
//#define PP_BSF_ITER_OUTPUT			// If PP_BSF_ITER_OUTPUT is defined then Iteration Output is performed
#define PP_BSF_TRACE_COUNT	1		// Each PP_BSF_TRACE_COUNT-th iteration to be outputted
#define PP_BSF_MAX_JOB_CASE 0

//--------------------------- OpenMP Parameters ---------------------------
// OpenMP mode is impossible here!

//--------------- BSF Lists parameters (For "No MPI" only) ----------------
#include "Problem-Parameters.h"
#define PP_BSF_MAP_LIST_LENGTH		PP_KK
#define PP_BSF_REDUCE_LIST_LENGTH	PP_KK
#define PP_BSF_REDUCE_LIST_1_LENGTH	1
#define PP_BSF_REDUCE_LIST_2_LENGTH	1
#define PP_BSF_REDUCE_LIST_3_LENGTH	1