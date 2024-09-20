/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (No MPI)
Module: Problem-Forwards.h (Problem Function Forwards)
Authors: Nikolay A. Olkhovsky & Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-bsfTypes.h"
#include "Problem-Types.h"
//====================== Private Functions ===========================
namespace PF {
	void	CodeToSubset(int code, int subset[PP_MM], int* ma);
//
	void	MakeFaceList(int* faceCodeList, int K);
	void	PreparationForIteration(PT_vector_T u);
	void	Print_Number_of_faces(PT_vector_T x);
//
// 
// 
}
//====================== Shared Functions ===========================
namespace SF {
	double	Distance_PointToHalfspace_i(PT_vector_T x, int i);
	double	Distance_PointToHyperplane_i(PT_vector_T x, int i);
	double	Distance_PointToPoint(PT_vector_T x, PT_vector_T y);
	double	Distance_PointToPolytope(PT_vector_T x);
	double	DistanceSQR_PointToPoint(PT_vector_T x, PT_vector_T y);
	void	JumpingOnPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double eps);
	void	MakeColumnOfNorms(PT_matrix_T A, PT_column_T norm_a);
	void	MakeListOfNotIncludingHalfspaces(PT_vector_T x, int* notIncludingHalfspacesList, double eps);
	void	MakePointHyperplaneList(PT_vector_T u, int* pointHyperplaneList, int* mh, double eps);
	void	MovingOnPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving);
	void	MovingToPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving);
	void	MTX_Conversion();
	void	MTX_ConversionSimple();
	bool	MTX_Load__Problem();
	bool	MTX_Load_A();
	bool	MTX_Load_b();
	bool	MTX_Load_c();
	bool	MTX_Load_hi();
	bool	MTX_Load_lo();
	bool	MTX_LoadPoint(PT_vector_T x, string postfix);
	void	MTX_RemoveFreeVariables(int m_equation, int m_inequality, int m_lowerBound, int m_higherBound);
	bool	MTX_SavePoint(PT_vector_T x, string postfix);
	void	MTX_SkipComments(FILE* stream);
	double	ObjF(PT_vector_T x);
	void	ObliqueProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T g, PT_vector_T o, double eps, int* exitCode);
	void	OrthogonalProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T r, double eps, int* exitcode);
	void	OrthogonalProjectingVectorOntoHyperplane_i(PT_vector_T x, int i, PT_vector_T p);
	bool	PointBelongsHalfspace_i(PT_vector_T point, int i, double eps);
	bool	PointBelongsHyperplane_i(PT_vector_T z, int i, double eps);
	bool	PointBelongsOuterCone(PT_vector_T x, int* notIncludingHalfspacesList, double eps);
	bool	PointBelongsPolytope(PT_vector_T x, double eps);
	void	PointHomothety(PT_vector_T x, PT_vector_T center, double ratio);
	bool	PointInsideHalfspace_i(PT_vector_T x, int i, double eps);
	int		PointLocation_i(PT_vector_T x, int i, double eps, double* a_DoT_x_MinuS_b);
	void	PolytopeHomothety(PT_vector_T center, double ratio);
	void	Print_Inequalities();
	void	Print_HalfspacesIncludingPoint(PT_vector_T x, double eps);
	void	Print_HalfspacesOutOfPoint(PT_vector_T x, double eps);
	void	Print_HyperplanesIncludingPoint(PT_vector_T x, double eps);
	void	Print_Vector(PT_vector_T x);
	void	PseudoprojectionOnFlat(int* flatHyperplanes, int m_flat, PT_vector_T v, double eps, int maxProjectingIter, PT_vector_T w, int* success);
	double	RelativeError(double trueValue, double calculatedValue);
	void	Shift(PT_vector_T point, PT_vector_T shiftVector, double factor, PT_vector_T shiftedPoint);
	void	Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z);
	void	Vector_Copy(PT_vector_T x, PT_vector_T y);
	void	Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y);
	void	Vector_DivideEquals(PT_vector_T x, double r);
	double	Vector_DotProduct(PT_vector_T x, PT_vector_T y);
	bool	Vector_Is_Tiny(PT_vector_T x, double eps);
	void	Vector_MakeLike(PT_vector_T x, double lengthOfLikeVector, PT_vector_T likeVector);
	void	Vector_MakeMinus_e(PT_vector_T minus_e);
	void	Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector);
	void	Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y);
	void	Vector_MultiplyEquals(PT_vector_T x, double r);
	double	Vector_Norm(PT_vector_T x);
	double	Vector_NormSquare(PT_vector_T x);
	void	Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector);
	void	Vector_Round(PT_vector_T x, double eps);
	void	Vector_SetValue(PT_vector_T x, double v);
	void	Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
	void	Vector_Zeroing(PT_vector_T x);
}
//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)
#define PF_MAP_LIST_INDEX (BSF_sv_addressOffset + BSF_sv_numberInSublist)