/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (MPI)
Module: Problems15-1.h (LP problems of dimension 15 with 1 randome inequality)
Prefix: PP
Authors: Nikolay A. Olkhovsky & Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
Initial surface points for these problems were calculated using BSF-Apex-Quest.
==============================================================================*/
#pragma once

//============================== Problem Parameters =============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
// PP_EPS_PROJECTION_ROUND - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9				// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO			// Precision for point to be in halfspace
#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)	// Precision for moving on polytope (affects Shift = 0)
#define PP_EPS_PROJECTION_ROUND		1E-8				// Precision of rounding pseudoprojecting vectors
//-------------------------------------------------------------------------------
#define PP_MAX_PROJECTING_ITER	1E+7	// Maximum acceptable number of iterations in PseudoprojectionOnFace()
//=============================================================================

/*============================== rnd16-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd16-0"
#define PP_KK	65536	// Maximal number of faces that include surface point 2^16 (compilator limit: 2^24 = 16 777 216)
#define PP_M	17		// Number of equations (number of rows in *.mtx)
#define PP_N	33		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 27100
//-----------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+2				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.1					// Length of probe shift
//-----------------------------------------------------------------------------

/*============================== rnd17-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd17-0"
#define PP_KK	131072	// Maximal number of faces that include surface point 2^17 (compilator limit: 2^24 = 16 777 216)
#define PP_M	18		// Number of equations (number of rows in *.mtx)
#define PP_N	35		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 30500
//-----------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+2				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.1					// Length of probe shift
//-----------------------------------------------------------------------------

/*============================== rnd18-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd18-0"
#define PP_KK	262144	// Maximal number of faces that include surface point 2^18 (compilator limit: 2^24 = 16 777 216)
#define PP_M	19		// Number of equations (number of rows in *.mtx)
#define PP_N	37		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 34100
//-----------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+2				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.1					// Length of probe shift
//-----------------------------------------------------------------------------

/*============================== rnd19-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd19-0"
#define PP_KK	524288	// Maximal number of faces that include surface point 2^19 (compilator limit: 2^24 = 16 777 216)
#define PP_M	20		// Number of equations (number of rows in *.mtx)
#define PP_N	39		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 37900
//-----------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+2				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.1					// Length of probe shift
//-----------------------------------------------------------------------------

/*============================== rnd20-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd20-0"
#define PP_KK	1048576	// Maximal number of faces that include surface point 2^20 (compilator limit: 2^24 = 16 777 216)
#define PP_M	21		// Number of equations (number of rows in *.mtx)
#define PP_N	41		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 41900
//-----------------------------------------------------------------------------
// Elapsed time: 140.7833
// Number of iterations: 9
// Computed objective value: 41899.9999999701
// Maximal objective value:  41900
// Relative error = 7.14e-13
//-----------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+3				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.01					// Length of probe shift
//-----------------------------------------------------------------------------

/*============================== rnd22-0 LP problem ===========================*/
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd22-0"
#define PP_KK	4194304	// Maximal number of faces that include surface point 2^22 (compilator limit: 2^24 = 16 777 216)
#define PP_M	23		// Number of equations (number of rows in *.mtx)
#define PP_N	45		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 50500
//-----------------------------------------------------------------------------
// Elapsed time: 835.20186
// Number of iterations: 10
// Computed objective value: 50499.99999997041
// Maximal objective value:  50500
// Relative error = 5.86e-13//-----------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+4				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001				// Length of probe shift
//-----------------------------------------------------------------------------

/*============================== rnd24-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd24-0"
#define PP_KK	16777216// Maximal number of faces that include surface point 2^24 (compilator limit: 2^24 = 16 777 216) 
#define PP_M	25		// Number of equations (number of rows in *.mtx)
#define PP_N	49		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 59900
//-----------------------------------------------------------------------------
// Elapsed time: 1471.8844
// Number of iterations: 11
// Computed objective value: 59899.99999996762
// Maximal objective value:  59900
// Relative error = 5.41e-13
//-----------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+4				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001					// Length of probe shift
//-----------------------------------------------------------------------------

/*=============================================================================*/