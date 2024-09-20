/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (No MPI)
Module: Problem-bsfTypes.h (Predefined BSF Problem Types)
Prefix: PT_bsf
Authors: Nikolay A. Olkhovsky & Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {		// Type of Parameter for workers (current approximation)
	PT_vector_T x;				// Current surface point
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	int* faceCode;
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T d;		// d = PP_PROBE_LENGTH*(w-u)/||w-u||
	double objF_p;	// F(p), where p=u+d
	int faceCode;
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Not used
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	// Not used
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
// Not used
};