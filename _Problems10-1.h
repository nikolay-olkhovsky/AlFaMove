/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (MPI)
Module: Problems10-1.h (LP problems of dimension 10 with 1 randome inequality)
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
#define PP_EPS_ZERO					1E-9				// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO			// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)	// Precision for moving on polytope (affects Shift = 0)
#define PP_EPS_PROJECTION_ROUND		1E-8				// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+9				// Length of Objective Vector
//-------------------------------------------------------------------------------
#define PP_KK	1024	// Maximal number of faces that include surface point (compilator limit: 2^24 = 16 777 216)
#define PP_M	11		// Number of equations (number of rows in *.mtx)
#define PP_N	21		// Number of variables (number of cols in *.mtx)

#define PP_PROBE_LENGTH			0.1		// Length of probe shift

//-------------------------------------------------------------------------------

/*============================== rnd10-0 LP problem ==============================*
// Exact solution:	100  200  200  200  200  200  200  200  200  200
#define PP_PROBLEM_NAME	"rnd10-0"
#define PP_MAX_OBJ_VALUE		10900			
//-------------------------------------------------------------------------------

/*============================== rnd10-1-1 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-1"
#define PP_MAX_OBJ_VALUE 9551.382889057777
//-------------------------------------------------------------------------------

/*============================== rnd10-1-2 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-2"
#define PP_MAX_OBJ_VALUE 9123.80496737513
//-------------------------------------------------------------------------------

/*============================== rnd10-1-3 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-3"
#define PP_MAX_OBJ_VALUE 9960.789316531023
//-------------------------------------------------------------------------------

/*============================== rnd10-1-4 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-4"
#define PP_MAX_OBJ_VALUE 9440.134567875428
//-------------------------------------------------------------------------------

/*============================== rnd10-1-5 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-5"
#define PP_MAX_OBJ_VALUE 10248.35536348364
//-------------------------------------------------------------------------------

/*============================== rnd10-1-6 LP problem ============================*/
#define PP_PROBLEM_NAME	"rnd10-1-6"
#define PP_MAX_OBJ_VALUE 10022.93578020061
//-------------------------------------------------------------------------------

/*===============================================================================*/