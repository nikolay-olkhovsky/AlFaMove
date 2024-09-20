/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (MPI)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Authors: Nikolay A. Olkhovsky & Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;
//========================== Problem variables ====================================
static int PD_m;					// Current number of inequalities
static int PD_n;					// Space dimension
static int PD_mh;					// Number of hyperplanes that include surface point
static int PD_ma;					// Number of hyperplanes used for pseudoprojection
static int PD_K;					// Number of faces of all possible dimensions
static int PD_iterNo;				// Number of iterations
//========================== Problem structures ====================================
static PT_matrix_T PD_A;			// Matrix of coefficients of inequalities
static PT_column_T PD_b;			// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;			// Gradient of Objective Function
static PT_vector_T PD_u;			// Current surface point
static PT_vector_T PD_hi;			// Higher bound
static PT_vector_T PD_lo;			// Lower bound
static PT_column_T PD_norm_a;		// Column of norms of matrix rows
static PT_vector_T PD_objVector;	// Used for pseudoprojecting
static int PD_pointHyperplanes[PP_MM];	// Index of hyperplanes that include surface point u
static int PD_faceCodes[PP_KK];			// Face codes
static int PD_faceHyperplanes[PP_MM];	// Index of hyperplanes used for pseudoprojection
//========================== Input/Output ====================================
static string PD_problemName;