/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (MPI)
Module: Problems20-1.h (LP problems of dimension 20 with 1 randome inequality)
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
#define PP_KK	1048576	// Maximal number of faces that include surface point 2^20 (compilator limit: 2 147 483 647) 
#define PP_M	21		// Number of equations (number of rows in *.mtx)
#define PP_N	41		// Number of variables (number of cols in *.mtx)

#define PP_PROBE_LENGTH			0.1		// Length of probe shift
//=============================================================================

/*============================== rnd20-0 LP problem ===========================*/
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd20-0"
#define PP_MAX_OBJ_VALUE 31000
//-----------------------------------------------------------------------------

/*=============================================================================*/