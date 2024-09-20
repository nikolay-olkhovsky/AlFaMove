/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (No MPI)
Module: _Problems-Klee-Minty-Cube.h (Problems from the LP-Klee-Minty-Cube Set)
Prefix: PP
Authors: Nikolay A. Olkhovsky & Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
==============================================================================*/
#pragma once

//============================== Problem Parameters ======================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
// PP_EPS_PROJECTION_ROUND - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
//-----------------------------------------------------------------------

#define PP_PROBE_LENGTH				1		// Length of probe shift

//=============================================================================

/*============================== Klee-Minty5 LP problem =======================*/
// Starting point:	0 ... 0
// Exact solution:	0 ... 0	3125
#define PP_PROBLEM_NAME	"Klee-Minty5"
#define PP_KK 32		// Maximal number of faces that include surface point 2^m
#define PP_D 5			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 3125
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)// Precision for moving on polytope
#define PP_EPS_PROJECTION_ROUND		1E-7			// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+10			// Length of Objective Vector
//-----------------------------------------------------------------------------
// Elapsed time: 0.2348548
// Number of iterations: 9
// Computed objective value: 3124.999999985131
// Maximal objective value:  3125
// Relative error = 4.76e-12
//==========================================================================

/*============================== Klee-Minty6 LP problem =======================*
// Exact solution: 0 ... 0	15625
#define PP_PROBLEM_NAME	"Klee-Minty6"
#define PP_KK 64		// Maximal number of faces that include surface point 2^m
#define PP_D 6			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 15625
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)// Precision for moving on polytope
#define PP_EPS_PROJECTION_ROUND		1E-6			// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+13			// Length of Objective Vector
//-----------------------------------------------------------------------------
// Elapsed time: 1.6791758
// Number of iterations: 11
// Computed objective value: 15624.99999997016
// Maximal objective value:  15625
// Relative error = 1.91e-12
//==========================================================================

/*============================== Klee-Minty7 LP problem =======================*
// Start point:	   0 ... 0
// Exact solution: 0 ... 0	78125
#define PP_PROBLEM_NAME	"Klee-Minty7"
#define PP_KK 128		// Maximal number of faces that include surface point 2^m
#define PP_D 7			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 78125
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-8			// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)// Precision for moving on polytope
#define PP_EPS_PROJECTION_ROUND		1E-4			// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+14			// Length of Objective Vector
//-----------------------------------------------------------------------------
// Elapsed time: 6.1532229
// Number of iterations: 13
// Computed objective value: 78124.9999991611
// Maximal objective value:  78125
// Relative error = 1.07e-11
//=============================================================================

/*============================== Klee-Minty8 LP problem =======================*
// Start point:	   0 ... 0
// Exact solution: 0 ... 0	390625
#define PP_PROBLEM_NAME	"Klee-Minty8"
#define PP_KK 256		// Maximal number of faces that include surface point 2^m
#define PP_D 8			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 390625
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-7			// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)// Precision for moving on polytope
#define PP_EPS_PROJECTION_ROUND		1E-2			// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+17			// Length of Objective Vector
//-----------------------------------------------------------------------------
// Elapsed time: 310.2361
// Number of iterations: 15
// Computed objective value: 390624.999989941
// Maximal objective value:  390625
// Relative error = 2.58e-11
//=============================================================================

/*============================== Klee-Minty9 LP problem =======================*
// Start point:	   0 ... 0
// Exact solution: 0 ... 0	1953125
#define PP_PROBLEM_NAME	"Klee-Minty9"
#define PP_KK 512		// Maximal number of faces that include surface point 2^m
#define PP_D 9			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 1953125
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-6			// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)// Precision for moving on polytope
#define PP_EPS_PROJECTION_ROUND		1E-2			// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+19			// Length of Objective Vector
//-----------------------------------------------------------------------------
// Elapsed time: 655.64446
// Number of iterations: 17
// Computed objective value: 1953124.999629887
// Maximal objective value:  1953125
// Relative error = 1.89e-10
//=============================================================================

/*=============================================================================*/