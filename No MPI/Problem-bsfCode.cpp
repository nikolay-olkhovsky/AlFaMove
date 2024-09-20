/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlFaMove - Along Faces Movement Method (No MPI)
Module: Problem-bsfCode.cpp (Implementation of Problem Code)
Prefix:	PC_bsf	- BSF Predefined Problem Functions
		SF		- Shared Functionc
		PF		- Private Functions
Authors: Nikolay A. Olkhovsky & Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;
using namespace SF;
using namespace PF;

//---------------------------------- BSF Predefined Problem Functions -----------------------------

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	Vector_Copy(parameterIn.x, parameterOutP->x);
}

void PC_bsf_Init(bool* success) {
	PD_problemName = PP_PROBLEM_NAME;

	*success = MTX_Load__Problem();
	if (*success == false)
		return;

	*success = MTX_LoadPoint(PD_u, PP_MTX_POSTFIX_U0);
	if (*success == false)
		return;

	MakeColumnOfNorms(PD_A, PD_norm_a);

	if (!PointBelongsPolytope(PD_u, PP_EPS_ZERO)) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Starting point does not belong to the feasible polytope with precision PP_EPS_ZERO = "
			<< PP_EPS_ZERO << "!!!\n";
		*success = false;
		return;
	}
	PD_iterNo = 0;
	Vector_MakeLike(PD_c, PP_OBJECTIVE_VECTOR_LENGTH, PD_objVector);
	PreparationForIteration(PD_u);
}

void PC_bsf_IterInit(PT_bsf_parameter_T parameter) {
	PreparationForIteration(parameter.x);
}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {

	cout << "# " << BSF_sv_iterCounter << "\tTime " << round(elapsedTime);
	cout << "\tx =";
	Print_Vector(parameter.x);
	cout << "\tF(x) = " << setw(PP_SETW) << ObjF(parameter.x);
}

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_JobDispatcher(PT_bsf_parameter_T* parameter, int* job, bool* toExit, double t) {
	// Not used
}

void PC_bsf_MainArguments(int argc, char* argv[]) {
	// Not used
}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success) {
	int faceCode = *mapElem->faceCode;
	PT_vector_T u;		// current surface point
	PT_vector_T v;		// v = u + PD_objVector (objVector = PP_OBJECTIVE_VECTOR_LENGTH*e_c)
	PT_vector_T w;		// pseudiprojection of v

	if (faceCode == 0) {
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = -PP_DBL_MAX;
		reduceElem->faceCode = 0;
		return;
	}

	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
	cout << "------------------------------------ Map(" << PF_MAP_LIST_INDEX << ") ------------------------------------" << endl;
#endif // PP_DEBUG /**/

	// Condition for breakpoint: PD_iterNo == 2 && (BSF_sv_addressOffset + BSF_sv_numberInSublist == 2)
	Vector_Copy(BSF_sv_parameter.x, u);
	double objF_u = ObjF(u);
	reduceElem->faceCode = faceCode;

	CodeToSubset(faceCode, PD_faceHyperplanes, &PD_ma);

	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
	cout << "Face hyperplanes: {";
	for (int i = 0; i < PD_ma - 1; i++) {
		cout << PD_faceHyperplanes[i] << ", ";
	}
	cout << PD_faceHyperplanes[PD_ma - 1] << "}.\n";
#endif // PP_DEBUG /**/

	Vector_Addition(u, PD_objVector, v);
	PseudoprojectionOnFlat(PD_faceHyperplanes, PD_ma, v, PP_EPS_PROJECTION_ROUND, PP_MAX_PSEUDOPROJECTING_ITER, w, success);

	if (!*success) {
		cout << "\n\nProcess " << BSF_sv_mpiRank
			<< ". Error in PC_bsf_MapF: Exceeded the maximum number of iterations when calculating pseudoprojection (PP_MAX_PSEUDOPROJECTING_ITER = "
			<< PP_MAX_PSEUDOPROJECTING_ITER << "). It is impossible to calculate Map function for element "
			<< BSF_sv_addressOffset + BSF_sv_numberInSublist << "!\n Perhaps you should decrease parameter PP_EPS_PROJECTION_ROUND.";
		return;
	}

	Vector_Round(w, PP_EPS_PROJECTION_ROUND * 10);
	Vector_Subtraction(w, u, (*reduceElem).d);

	double norm_d = Vector_Norm((*reduceElem).d);
	if (norm_d < PP_EPS_ZERO) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t||d|| = ||w-u|| < PP_EPS_ZERO ===>>> movement is impossible.\n";
#endif // PP_DEBUG /**/
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
	}

	Vector_MultiplyEquals((*reduceElem).d, PP_PROBE_LENGTH / norm_d);
	Vector_Round((*reduceElem).d, PP_EPS_PROJECTION_ROUND);

	PT_vector_T p;
	Vector_Addition(u, (*reduceElem).d, p);

	if (!PointBelongsPolytope(p, PP_EPS_POINT_IN_HALFSPACE)) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "Shifted point p = ";
		Print_Vector(p);
		cout << "\tnot in feasible polytope ===>>> movement is impossible." << endl;
#endif // PP_DEBUG /**/
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
	}

	reduceElem->objF_p = ObjF(p);

	if (RelativeError(objF_u, reduceElem->objF_p) < PP_EPS_ZERO) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "u =\t    ";
		Print_Vector(u);
		cout << "\tF(u) = " << setw(PP_SETW) << objF_u << endl;
		cout << "|F(u1)-F(u2)|/|F(u1)| = " << RelativeError(objF_u, reduceElem->objF_p) << " < PP_EPS_ZERO = " << PP_EPS_ZERO << " ===>>> movement is impossible.\n";
#endif // PP_DEBUG /**/

		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
	}

	if (reduceElem->objF_p < objF_u) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "u =\t    ";
		Print_Vector(u);
		cout << "\tF(w) = " << setw(PP_SETW) << reduceElem->objF_p << " < F(u) = " << objF_u << " ===>>> movement is impossible.\n";
#endif // PP_DEBUG /**/
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
	}

	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
	cout << "d = ";
	Print_Vector((*reduceElem).d);
	cout << endl;
	cout << "Shifted point p = ";
	Print_Vector(p);
	cout << "\tF(p) ="
		<< setw(PP_SETW) << reduceElem->objF_p << "\t\t---> Movement is possible." << endl;
#endif // PP_DEBUG /**/
} // end PC_bsf_MapF

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// Not used
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== " << PP_METHOD_NAME << " ====================================================" << endl;
	cout << "Problem name: " << PD_problemName << endl;

#ifdef PP_MPI
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#else
	cout << "No MPI" << endl;
#endif // PP_MPI

#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP

#ifdef PP_BSF_FRAGMENTED_MAP_LIST
	cout << "Map List is Fragmented" << endl;
#else
	cout << "Map List is not Fragmented" << endl;
#endif

#ifdef PP_SIMPLE_CONVERSION
	cout << "Conversion mode: simple (with preservation of free variables)" << endl;
#else
	cout << "Conversion mode: full (with elimination of free variables)" << endl;
#endif

	cout << "Before conversion: m =\t" << PP_M << "\tn = " << PP_N << endl;
	cout << "After conversion:  m =\t" << PD_m << "\tn = " << PD_n << endl;
	cout << "PP_EPS_ZERO\t\t\t" << PP_EPS_ZERO << endl;
	cout << "PP_EPS_POINT_IN_HALFSPACE\t" << PP_EPS_POINT_IN_HALFSPACE << endl;
	cout << "PP_EPS_PROJECTION_ROUND\t\t" << PP_EPS_PROJECTION_ROUND << endl;
	cout << "PP_OBJECTIVE_VECTOR_LENGTH\t" << PP_OBJECTIVE_VECTOR_LENGTH << endl;
	cout << "PP_PROBE_LENGTH\t\t\t" << PP_PROBE_LENGTH << endl;
	cout << "--------------- Data ---------------\n";
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	Print_Inequalities();
#endif // PP_MATRIX_OUTPUT
	cout << "Obj Function:\t";
	Print_Vector(PD_c);
	cout << endl;
	cout << "u0 =\t\t";
	Print_Vector(PD_u); cout << "\tF(x) = " << setw(PP_SETW) << ObjF(PD_u) << endl;

#ifdef PP_DEBUG
	if (!PointBelongsPolytope(PD_u, PP_EPS_POINT_IN_HALFSPACE))
		cout << "u0 is outside feasible polytope!!!\n";
	else
		cout << "u0 is belongs feasible polytope.\n";
	cout << "Including hyperplanes:\t"; Print_HyperplanesIncludingPoint(PD_u, PP_EPS_POINT_IN_HALFSPACE); cout << endl;
	cout << "Including faces:\t"; Print_Number_of_faces(PD_u);
#endif // PP_DEBUG
}

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	cout << setprecision(PP_SETW / 2);

	cout << "================================================" << endl;
	cout << "// Elapsed time: " << t << endl;
	cout << "// Number of iterations: " << PD_iterNo << endl;
	cout << "// Computed objective value: " << setprecision(16) << ObjF(PD_u) << endl;
	cout << "// Maximal objective value:  " << PP_MAX_OBJ_VALUE << endl;
	cout << "// Relative error = " << setprecision(3) << RelativeError(PP_MAX_OBJ_VALUE, ObjF(PD_u)) << setprecision(PP_SETW / 2) << endl;
	cout << "================================================" << endl;

#ifdef PP_SAVE_RESULT
	if (MTX_SavePoint(PD_u, PP_MTX_POSTFIX_SO))
		cout << "Calculated solution point is saved into file *.so" << endl;
#endif // PP_SAVE_RESULT

	cout << "Solution point:\t";
	Print_Vector(PD_u);	cout << endl;
#ifdef PP_DEBUG
	cout << "Distance to polytope: " << Distance_PointToPolytope(PD_u) << endl;
#endif // PP_DEBUG

} // end PC_bsf_ProblemOutput

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProcessResults(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	PT_vector_T u_moved;
	double shiftLength;

	if (Vector_Norm(reduceResult->d) < PP_EPS_ZERO) {
		*toExit = true;
		return;
	}

	//MovingOnPolytope(PD_u, reduceResult->d, u_moved, PP_EPS_MOVING_ON_POLYTOPE);
	JumpingOnPolytope(PD_u, reduceResult->d, u_moved, PP_EPS_ZERO);
	shiftLength = Distance_PointToPoint(PD_u, u_moved);

#ifdef PP_DEBUG
	if (shiftLength < PP_EPS_ZERO) {
		cout << "\nError in MovingOnPolytope(): Shift length < PP_EPS_ZERO! Possible you should increase PP_EPS_MOVING_ON_POLYTOPE or decrease PP_EPS_ZERO.\n";
		exit(1);
	}
#endif // PP_DEBUG

#ifdef PP_DEBUG
	CodeToSubset(reduceResult->faceCode, PD_faceHyperplanes, &PD_ma);
	cout << "_________________________________________________ " << PD_iterNo << " _____________________________________________________" << endl;
	cout << "Face dimension = " << PD_n - PD_ma << endl;
	cout << "Point:\t";
	Print_Vector(PD_u); cout << "\tF(x) = " << ObjF(PD_u) << endl;
	cout << "Point hyperplanes:\t"; Print_HyperplanesIncludingPoint(PD_u, PP_EPS_POINT_IN_HALFSPACE); cout << endl;
	cout << "Face Hyperplanes:\t{";
	for (int i = 0; i < PD_ma - 1; i++)
		cout << PD_faceHyperplanes[i] << ", ";
	cout << PD_faceHyperplanes[PD_ma - 1]
		<< "}\tShift = " << shiftLength << endl;
#endif // PP_DEBUG

	Vector_Copy(u_moved, PD_u);
	Vector_Copy(u_moved, parameter->x);
	PD_iterNo++;
}

void PC_bsf_ProcessResults_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	// Not used
}

void PC_bsf_ProcessResults_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	// Not used
}

void PC_bsf_ProcessResults_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	// Not used
}

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	if (x->objF_p > y->objF_p) {
		z->objF_p = x->objF_p;
		for (int j = 0; j < PD_n; j++)
			(*z).d[j] = (*x).d[j];
		z->faceCode = x->faceCode;
	}
	else {
		z->objF_p = y->objF_p;
		for (int j = 0; j < PD_n; j++)
			(*z).d[j] = (*y).d[j];
		z->faceCode = y->faceCode;
	}
}

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	// Not used
}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	// Not used
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// Not used
}

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	Vector_Copy(PD_u, parameter->x);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_KK;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->faceCode = &(PD_faceCodes[i]);
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; }
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; }
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; }
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; }
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; }
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; }
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; }
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; }

//---------------------------------- Shared Functions -------------------------
namespace SF {

	static inline double Distance_PointToHalfspace_i(PT_vector_T x, int i) {
		double a_DoT_z_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];

		if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
			return 0;

		if (a_DoT_z_MinuS_b < 0) // Point belongs to halfspace
			return 0;

		return a_DoT_z_MinuS_b / PD_norm_a[i];
	}

	static inline double Distance_PointToHyperplane_i(PT_vector_T x, int i) {
		if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
			return 0;
		else
			return fabs(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / PD_norm_a[i];
	}

	static inline double Distance_PointToPoint(PT_vector_T x, PT_vector_T y) {
		PT_vector_T z;
		Vector_Subtraction(x, y, z);
		return Vector_Norm(z);
	}

	static inline double Distance_PointToPolytope(PT_vector_T x) { // Measure of distance from point to polytope
		double maxDistance = 0;
		double distance;

		for (int i = 0; i < PD_m; i++) {
			distance = Distance_PointToHalfspace_i(x, i);
			if (distance > 0)
				maxDistance = PF_MAX(maxDistance, distance);
		}
		return maxDistance;
	}

	static inline double DistanceSQR_PointToPoint(PT_vector_T x, PT_vector_T y) {
		PT_vector_T z;
		Vector_Subtraction(x, y, z);
		return Vector_NormSquare(z);
	}

	static inline void JumpingOnPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double eps) {
		PT_vector_T o; // Oblique projection vector
		PT_vector_T o_min; // Oblique projection vector with minimum length
		double lengthSQR_o;
		double* z = startPoint;
		double* d = directionVector;
		double a_DoT_d;
		int location_z;
		double a_DoT_z_MinuS_b;
		double minLengthSQR = PP_INFINITY;

		Vector_Zeroing(o_min);

		for (int i = 0; i < PD_m; i++) {
			location_z = PointLocation_i(z, i, eps, &a_DoT_z_MinuS_b);
			assert(location_z != PP_DEGENERATE_INEQUALITY);

			switch (location_z) {
			case PP_ON_HYPERPLANE:
				continue;
			case PP_OUTSIDE_HALFSPACE:
				continue;
			case PP_INSIDE_HALFSPACE:
				a_DoT_d = Vector_DotProduct(PD_A[i], d); // <a,d>

				if (a_DoT_d < PP_EPS_ZERO)   // <a,d> <= 0
					continue;

				// Oblique projection vector: o = -(<a,z> - b)d/<a, d>
				Vector_MultiplyByNumber(d, -a_DoT_z_MinuS_b / a_DoT_d, o);
				lengthSQR_o = Vector_NormSquare(o);
				if (minLengthSQR > lengthSQR_o) {
					minLengthSQR = lengthSQR_o;
					Vector_Copy(o, o_min);
				}
				break;
			default:
				assert(false);
			}
		}
		Vector_Addition(startPoint, o_min, finishPoint);
	}

	static inline void MakeColumnOfNorms(PT_matrix_T A, PT_column_T norm_a) {
		for (int i = 0; i < PD_m; i++)
			norm_a[i] = Vector_Norm(A[i]);
	}

	static inline void MakeListOfNotIncludingHalfspaces(PT_vector_T x, int* notIncludingHalfspacesList, double eps) {
		int mo = 0;
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsHalfspace_i(x, i, eps)) {
				notIncludingHalfspacesList[mo] = i;
				mo++;
			}
		if (mo < PD_m)
			notIncludingHalfspacesList[mo] = -1;
	}

	static inline void MakePointHyperplaneList(PT_vector_T u, int* pointHyperplaneList, int* mh, double eps) {
		*mh = 0;
		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHyperplane_i(u, i, eps)) {
				pointHyperplaneList[*mh] = i;
				(*mh)++;
			}
		}
	}

	static inline void MovingOnPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving) {
		double leftBound = 0;
		double rightBound = PP_DBL_MAX;
		double factor = 1;
		double delta;

		assert(Vector_Norm(directionVector) >= PP_EPS_ZERO);

		delta = factor / 2;

		while (rightBound - leftBound >= PP_EPS_ZERO && delta > 0) {
			Shift(startPoint, directionVector, factor, finishPoint);
			if (PointBelongsPolytope(finishPoint, PP_EPS_POINT_IN_HALFSPACE)) {
				leftBound = factor;
				delta *= 2;
				factor += delta;
			}
			else {
				rightBound = factor;
				delta /= 2;
				factor -= delta;
			}
		}

		Shift(startPoint, directionVector, factor, finishPoint);
		delta = epsMoving;
		while (!PointBelongsPolytope(finishPoint, epsMoving) && delta > 0) {
			factor -= delta;
			delta *= 2;
			Shift(startPoint, directionVector, factor, finishPoint);
		}
	}

	static inline void MovingToPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving) {
		double leftBound = 0;
		double rightBound = PP_DBL_MAX;
		double factor = 1;
		double delta;
		static int outerHalspace_i[PP_MM];	// Index of out half-spaces
		int mo;								// Number of out half-spaces
		bool pointInsideCone;

		assert(Vector_Norm(directionVector) >= PP_EPS_ZERO);

		mo = 0;
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsHalfspace_i(startPoint, i, PP_EPS_POINT_IN_HALFSPACE)) {
				outerHalspace_i[mo] = i;
				mo++;
			}

		delta = factor / 2;

		while (rightBound - leftBound >= PP_EPS_ZERO && delta > 0) {
			Shift(startPoint, directionVector, factor, finishPoint);

			pointInsideCone = true;
			for (int i = 0; i < mo; i++)
				if (PointBelongsHalfspace_i(finishPoint, outerHalspace_i[i], PP_EPS_POINT_IN_HALFSPACE)) {
					pointInsideCone = false;
					break;
				}
			if (pointInsideCone) {
				leftBound = factor;
				delta *= 2;
				factor += delta;
			}
			else {
				rightBound = factor;
				delta /= 2;
				factor -= delta;
				assert(factor > 0);
			}
		}

		Shift(startPoint, directionVector, factor, finishPoint);
		delta = epsMoving;
		do {
			pointInsideCone = false;
			for (int i = 0; i < mo; i++)
				if (!PointBelongsHalfspace_i(finishPoint, outerHalspace_i[i], epsMoving)) {
					pointInsideCone = true;
					factor -= delta;
					delta *= 2;
					assert(factor > 0);
					Shift(startPoint, directionVector, factor, finishPoint);
					break;
				}
		} while (pointInsideCone && delta > 0);
	}

	static void MTX_Conversion() { // Transformation to inequalities & dimensionality reduction
		int m_equation = PD_m;
		int m_inequality;
		int m_lowerBound;
		int m_higherBound;

		for (int i = 0; i < m_equation; i++) { // Conversion to inequalities
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = -PD_A[i][j];
			PD_b[PD_m] = -PD_b[i];
			PD_m++; assert(PD_m <= PP_MM);
		}

		for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
			for (int j = 0; j < PD_n; j++)
				if (PD_A[i][j] == 0)
					PD_A[i][j] = 0;
		m_inequality = PD_m;

		for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
			for (int j = 0; j < PD_n; j++)
				PD_A[i + PD_m][j] = 0;
			PD_A[i + PD_m][i] = -1;
			if (PD_lo[i] == 0)
				PD_b[i + PD_m] = 0;
			else
				PD_b[i + PD_m] = -PD_lo[i];
		}
		PD_m += PD_n; assert(PD_m <= PP_MM);
		m_lowerBound = PD_m;

		for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
			if (PD_hi[i] != PP_INFINITY) {
				for (int j = 0; j < PD_n; j++)
					PD_A[PD_m][j] = 0;
				PD_A[PD_m][i] = 1;
				PD_b[PD_m] = PD_hi[i];
				PD_m++; assert(PD_m <= PP_MM);
			}
		}
		m_higherBound = PD_m;

		/**
		cout << "-----------------------------------------------------\n";
		Print_Inequalities();
		cout << "-----------------------------------------------------\n";
		cout << "PD_c: "; Print_Vector(PD_c); cout << endl;/**/

		MTX_RemoveFreeVariables(m_equation, m_inequality, m_lowerBound, m_higherBound);

		/**
		cout << "-----------------------------------------------------\n";
		Print_Inequalities();
		cout << "-----------------------------------------------------\n";/**/
	}

	static inline void MTX_ConversionSimple() { // Transformation to inequalities
		int m_equation = PD_m;

		for (int i = 0; i < m_equation; i++) { // Conversion to inequalities
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = -PD_A[i][j];
			PD_b[PD_m] = -PD_b[i];
			PD_m++; assert(PD_m <= PP_MM);
		}

		for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
			for (int j = 0; j < PD_n; j++)
				if (PD_A[i][j] == 0)
					PD_A[i][j] = 0;

		for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
			for (int j = 0; j < PD_n; j++)
				PD_A[i + PD_m][j] = 0;
			PD_A[i + PD_m][i] = -1;
			if (PD_lo[i] == 0)
				PD_b[i + PD_m] = 0;
			else
				PD_b[i + PD_m] = -PD_lo[i];
		}
		PD_m += PD_n; assert(PD_m <= PP_MM);

		for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
			if (PD_hi[i] != PP_INFINITY) {
				for (int j = 0; j < PD_n; j++)
					PD_A[PD_m][j] = 0;
				PD_A[PD_m][i] = 1;
				PD_b[PD_m] = PD_hi[i];
				PD_m++; assert(PD_m <= PP_MM);
			}
		}
	}

	static bool MTX_Load__Problem() {

		//--------------- Reading A ------------------
		if (!MTX_Load_A())
			return false;

		//--------------- Reading b ------------------
		if (!MTX_Load_b())
			return false;

		//--------------- Reading lo ------------------
		if (!MTX_Load_lo())
			return false;

		//--------------- Reading c ------------------
		if (!MTX_Load_c())
			return false;

		//--------------- Reading hi ------------------
		if (!MTX_Load_hi())
			return false;

		//---------- Conversion to inequalities -----------
#ifdef PP_SIMPLE_CONVERSION
		MTX_ConversionSimple();
#else
		MTX_Conversion();
#endif

		/**
		cout << "-----------------------------------------------------\n";
		Print_Inequalities();
		cout << "-----------------------------------------------------\n";
		cout << "PD_c: "; Print_Vector(PD_c); cout << endl;/**/

		return true;
	}

	static inline bool MTX_Load_A() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		int non;	// Number of non-zero elements
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string PD_MTX_File;

		PD_MTX_File = PP_PATH;
		PD_MTX_File += PP_MTX_PREFIX;
		PD_MTX_File += PD_problemName;
		PD_MTX_File += PP_MTX_POSTFIX_A;
		mtxFile = PD_MTX_File.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d%d", &nor, &noc, &non) < 3) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file " << mtxFile << endl;
			return false;
		}

		if (nor >= noc) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Number of rows m = " << nor << " must be < " << "Number of columns n = " << noc << "\n";
			return false;
		}

		if (noc != PP_N) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Invalid input data: PP_N must be = " << noc << "\n";
			return false;
		}

		if (nor != PP_M) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Invalid input data: PP_M must be = " << nor << "\n";
			return false;
		}

		PD_m = nor;
		PD_n = noc;

		if (2 * nor + noc > PP_MM) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Invalid input data: number of inequalities m = " << 2 * nor + noc
				<< " must be < PP_MM + 1 =" << PP_MM + 1 << "\n";
			return false;
		}

		for (int k = 0; k < non; k++) {
			int i, j;

			if (fscanf(stream, "%d%d%s", &i, &j, str) < 3) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file'" << mtxFile << "'." << endl;
				return false;
			}

			i -= 1;
			j -= 1;
			if (i < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Negative row index in'" << mtxFile << "'.\n" << endl;
				return false;
			}
			if (j < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Negative column index in'" << mtxFile << "'.\n" << endl;
				return false;
			}
			PD_A[i][j] = strtod(str, &chr);
		}

		fclose(stream);

		return true;
	}

	static inline bool MTX_Load_b() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string PD_MTX_File;

		PD_MTX_File = PP_PATH;
		PD_MTX_File += PP_MTX_PREFIX;
		PD_MTX_File += PD_problemName;
		PD_MTX_File += PP_MTX_POSTFIX_B;
		mtxFile = PD_MTX_File.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (PD_m != nor) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
			return false;
		}

		for (int i = 0; i < PD_m; i++) {
			if (fscanf(stream, "%s", str) < 1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file '" << mtxFile << "'." << endl;
				return false;
			}
			PD_b[i] = strtod(str, &chr);
		}
		fclose(stream);

		return true;
	}

	static inline bool MTX_Load_c() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string PD_MTX_File;

		PD_MTX_File = PP_PATH;
		PD_MTX_File += PP_MTX_PREFIX;
		PD_MTX_File += PD_problemName;
		PD_MTX_File += PP_MTX_POSTFIX_C;
		mtxFile = PD_MTX_File.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file" << endl;
				return false;
			}
			PD_c[j] = -strtod(str, &chr);
		}
		fclose(stream);

		return true;
	}

	static inline bool MTX_Load_hi() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string PD_MTX_File;

		PD_MTX_File = PP_PATH;
		PD_MTX_File += PP_MTX_PREFIX;
		PD_MTX_File += PD_problemName;
		PD_MTX_File += PP_MTX_POSTFIX_HI;
		mtxFile = PD_MTX_File.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "Unexpected end of file '" << mtxFile << "'." << endl;
				return false;
			}
			PD_hi[j] = strtod(str, &chr);
		}
		fclose(stream);
		return true;
	}

	static inline bool MTX_Load_lo() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string PD_MTX_File;

		PD_MTX_File = PP_PATH;
		PD_MTX_File += PP_MTX_PREFIX;
		PD_MTX_File += PD_problemName;
		PD_MTX_File += PP_MTX_POSTFIX_LO;
		mtxFile = PD_MTX_File.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file '" << mtxFile << "'." << endl;
				return false;
			}
			PD_lo[j] = strtod(str, &chr);
		}

		fclose(stream);

		return true;
	}

	static inline bool MTX_LoadPoint(PT_vector_T x, string postfix) {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string PD_MTX_File;

		PD_MTX_File = PP_PATH;
		PD_MTX_File += PP_MTX_PREFIX;
		PD_MTX_File += PD_problemName;
		PD_MTX_File += postfix;
		mtxFile = PD_MTX_File.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file" << endl;
				return false;
			}
			x[j] = strtod(str, &chr);
		}
		fclose(stream);

		return true;
	}

	static inline void MTX_RemoveFreeVariables(int m_equation, int m_inequality, int m_lowerBound, int m_higherBound) {
		int freeVariable_j_col[PP_M];
		int m_freeVariable_j_col = 0;
		int rowToDelete_i[PP_MM];
		int m_rowToDelete = 0;
		bool unique;
		bool ok;

		// Find free variables
		for (int j_col = 0; j_col < PD_n; j_col++) { // Find zero element in PD_c
			if (PD_c[j_col] != 0)
				continue;
			for (int i_row = 0; i_row < m_equation; i_row++) { // Find PD_A i_row with nonzero element in column j_col
				if (PD_A[i_row][j_col] == 0)
					continue;

				// Check uniqueness in column j_col
				unique = true;
				for (int i = 0; i < m_equation; i++) {
					if (i == i_row)
						continue;
					if (PD_A[i][j_col] != 0) {
						unique = false;
						break;
					}
				}
				if (!unique)
					continue;

				// Check lower bound
				ok = true;
				for (int i_lowerBound = m_inequality; i_lowerBound < m_lowerBound; i_lowerBound++) {
					if (PD_A[i_lowerBound][j_col] == 0)
						continue;
					if (PD_A[i_lowerBound][j_col] != -1) {
						ok = false;
						break;
					}
					if (PD_b[i_lowerBound] != 0) {
						ok = false;
						break;
					}
					rowToDelete_i[m_rowToDelete] = i_lowerBound;
					m_rowToDelete++; assert(m_rowToDelete <= PP_MM);
					break;
				}
				if (!ok)
					continue;

				// Check higher bound
				ok = true;
				for (int i_higherBound = m_lowerBound; i_higherBound < m_higherBound; i_higherBound++) {
					if (PD_A[i_higherBound][j_col] != 0) {
						ok = false;
						m_rowToDelete--;
						break;
					}
				}
				if (!ok)
					continue;

				freeVariable_j_col[m_freeVariable_j_col] = j_col;
				m_freeVariable_j_col++; assert(m_freeVariable_j_col <= PP_M);
				rowToDelete_i[m_rowToDelete] = m_equation + i_row;
				m_rowToDelete++; assert(m_rowToDelete <= PP_MM);
				// Check sign of free variable
				if (PD_A[i_row][j_col] < 0) {
					// Change sign of inequality
					for (int j = 0; j < PD_n; j++)
						PD_A[i_row][j] = -PD_A[i_row][j];
					PD_b[i_row] = -PD_b[i_row];
				}
				break;
			}
		}

		{// Eliminate columns with free variables
			static bool colToDeleteLable[PP_N];
			for (int j = 0; j < m_freeVariable_j_col; j++) {
				assert(freeVariable_j_col[j] < PP_N);
				colToDeleteLable[freeVariable_j_col[j]] = true;
			}

			for (int j = 0; j < PD_n; j++) {
				if (colToDeleteLable[j]) {
					for (int i = 0; i < PD_m; i++)
						PD_A[i][j] = PD_A[i][PD_n - 1];
					PD_c[j] = PD_c[PD_n - 1];
					colToDeleteLable[j] = colToDeleteLable[PD_n - 1];
					j--; assert(j >= 0);
					PD_n--; assert(PD_n >= 0);

					/**
					cout << "-----------------------------------------------------\n";
					Print_Inequalities();
					cout << "-----------------------------------------------------\n";/**/

				}
			}
		}

		{// Eliminate rows corresponding to free variables
			static bool rowToDeleteLable[PP_MM];
			for (int i = 0; i < m_rowToDelete; i++) {
				rowToDeleteLable[rowToDelete_i[i]] = true;
			}
			for (int i = 0; i < PD_m; i++) {
				if (rowToDeleteLable[i]) {
					for (int j = 0; j < PD_n; j++)
						PD_A[i][j] = PD_A[PD_m - 1][j];
					PD_b[i] = PD_b[PD_m - 1];
					rowToDeleteLable[i] = rowToDeleteLable[PD_m - 1];
					i--; assert(i >= 0);
					PD_m--; assert(PD_m >= 0);

					/**
					cout << "-----------------------------------------------------\n";
					Print_Inequalities();
					cout << "-----------------------------------------------------\n";/**/
				}
			}
		}
	}

	static bool MTX_SavePoint(PT_vector_T x, string postfix) {
		const char* mtxFile;
		FILE* stream;// Input stream
		string PD_MTX_File;

		PD_MTX_File = PP_PATH;
		PD_MTX_File += PP_MTX_PREFIX;
		PD_MTX_File += PD_problemName;
		PD_MTX_File += postfix;
		mtxFile = PD_MTX_File.c_str();
		stream = fopen(mtxFile, "w");
		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		fprintf(stream, "%d %d\n", PD_n, 1);

		for (int j = 0; j < PD_n; j++)
			fprintf(stream, "%.16f\n", x[j]);

		fclose(stream);
		return true;
	}

	static inline void MTX_SkipComments(FILE* stream) {
		fpos_t pos;	// Position in the input stream
		int res;
		res = fscanf(stream, "\n");
		fgetpos(stream, &pos);
		while (getc(stream) == '%') {
			while (getc(stream) != 10);
			res = fscanf(stream, "\n");
			fgetpos(stream, &pos);
		}
		fsetpos(stream, &pos);
	}

	static inline double ObjF(PT_vector_T x) {
		double s = 0;
		for (int j = 0; j < PD_n; j++)
			s += PD_c[j] * x[j];
		return s;
	}

	static inline void	ObliqueProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T d, PT_vector_T o, double eps, int* exitCode) {
		// Oblique projecting vector o of point z onto Half-space H_i with respect to vector d
		double a_DoT_g;	// <a,d>
		double a_DoT_z_MinuS_b;	// <a,z> - b
		double factor;	// (b - <a,z>) / <a,d>

		if (PD_norm_a[i] < PP_EPS_ZERO) {
			Vector_Zeroing(o);
			*exitCode = PP_DEGENERATE_INEQUALITY;
			return;
		}

		a_DoT_z_MinuS_b = Vector_DotProduct(PD_A[i], z) - PD_b[i]; // <a,z> - b

		if (fabs(a_DoT_z_MinuS_b) / PD_norm_a[i] < PP_EPS_ZERO) { // |<a,z> - b|/||a|| = 0
			*exitCode = PP_ON_HYPERPLANE;
			Vector_Zeroing(o);
			return;
		}

		if (a_DoT_z_MinuS_b < 0) { // <a,z> - b < 0
			*exitCode = PP_INSIDE_HALFSPACE; 
			Vector_Zeroing(o);
			return;
		}

		a_DoT_g = Vector_DotProduct(PD_A[i], d); // <a,d>


		if (fabs(a_DoT_g) < PP_EPS_ZERO) {
			*exitCode = PP_PARALLEL; 
			Vector_Zeroing(o);
			return;
		}

		if (a_DoT_g >= PP_EPS_ZERO) {
			*exitCode = PP_RECESSIVE; 
			Vector_Zeroing(o);
			return;
		}

		factor = a_DoT_z_MinuS_b / a_DoT_g; // (<a,z> - b) / <a,d>

		// Oblique projection vector: o = -(<a,z> - b)d/<a, d> = -factor * d
		Vector_MultiplyByNumber(d, -factor, o);

		*exitCode = PP_NONDEGENERATE_PROJECTING;
		return;
	}

	static inline void OrthogonalProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T r, double eps, int* exitCode) {
		double factor;
		double a_DoT_z_MinuS_b = Vector_DotProduct(PD_A[i], z) - PD_b[i]; // <a,z> - b
		double distance = fabs(a_DoT_z_MinuS_b) / PD_norm_a[i];

		if (PD_norm_a[i] < PP_EPS_ZERO) {
			Vector_Zeroing(r);
			*exitCode = PP_DEGENERATE_INEQUALITY; 
			return;
		}

		if (distance < eps) {
			Vector_Zeroing(r);
			*exitCode = PP_ON_HYPERPLANE; 
			return;
		}

		if (a_DoT_z_MinuS_b < 0) { // <a,z> - b < 0
			Vector_Zeroing(r);
			*exitCode = PP_INSIDE_HALFSPACE; 
			return;
		}

		factor = -a_DoT_z_MinuS_b / (PD_norm_a[i] * PD_norm_a[i]); // (b - <z,a>) / ||a||^2
		Vector_MultiplyByNumber(PD_A[i], factor, r); // r = a(b - <z,a>) / ||a||^2
		*exitCode = PP_NONDEGENERATE_PROJECTING;
	}

	static inline void OrthogonalProjectingVectorOntoHyperplane_i(PT_vector_T x, int i, PT_vector_T p) {
		assert(Vector_NormSquare(PD_A[i]));
		Vector_MultiplyByNumber(PD_A[i], -(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / Vector_NormSquare(PD_A[i]), p);
	}

	static inline bool PointBelongsHalfspace_i(PT_vector_T x, int i, double eps) {
		if (PD_norm_a[i] < eps) //Degenerate equation
			return true;
		double a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_DoT_x_MinuS_b) / PD_norm_a[i];
		if (distanceToHyperplane < eps)
			return true;
		if (a_DoT_x_MinuS_b < 0)
			return true;
		return false;
	}

	static inline bool PointBelongsHyperplane_i(PT_vector_T x, int i, double eps) {
		if (Distance_PointToHyperplane_i(x, i) < eps)
			return true;
		else
			return false;
	}

	static inline bool PointBelongsPolytope(PT_vector_T x, double eps) { // If the point belongs to the polytope with prescigion of eps
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsHalfspace_i(x, i, eps))
				return false;
		return true;
	}

	static inline bool PointBelongsOuterCone(PT_vector_T x, int* notIncludingHalfspacesList, double eps) { // If the point belongs to the outer cone with prescigion of eps
		for (int i = 0; i < PD_m && notIncludingHalfspacesList[i] >= 0; i++)
			if (PointBelongsHalfspace_i(x, i, eps))
				return false;
		return true;
	}

	static inline void PointHomothety(PT_vector_T x, PT_vector_T center, double ratio) { // https://en.wikipedia.org/wiki/Homothety
		if (ratio == 1)
			return;
		assert(ratio > 0);
		for (int j = 0; j < PD_n; j++)
			x[j] = ratio * x[j] - (ratio - 1) * center[j];
	}

	static inline bool PointInsideHalfspace_i(PT_vector_T x, int i, double eps) {
			if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
			return true;
		double a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_DoT_x_MinuS_b) / PD_norm_a[i];
		if (distanceToHyperplane < eps)
			return false;
		if (a_DoT_x_MinuS_b < 0)
			return true;
		return false;
	}

	static inline int PointLocation_i(PT_vector_T x, int i, double eps, double* a_DoT_x_MinuS_b) {

		if (PD_norm_a[i] < PP_EPS_ZERO) 
			return PP_DEGENERATE_INEQUALITY;

		*a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];

		if (fabs(*a_DoT_x_MinuS_b) / PD_norm_a[i] < PP_EPS_ZERO)// <a,x> = b
			return PP_ON_HYPERPLANE;

		if (*a_DoT_x_MinuS_b < 0)								// <a,x> < b
			return PP_INSIDE_HALFSPACE;	

		return PP_OUTSIDE_HALFSPACE;							// <a,x> > b
	}

	static inline void PolytopeHomothety(PT_vector_T center, double ratio) { // https://en.wikipedia.org/wiki/Homothety
		if (ratio == 1)
			return;
		assert(ratio > 0);

		for (int i = 0; i < PD_m; i++) {
			PD_b[i] = ratio * PD_b[i] - (ratio - 1) * Vector_DotProduct(PD_A[i], center);
		}
	}

	static inline void Print_Inequalities() {
		for (int i = 0; i < PD_m; i++) {
			cout << i << ")";
			for (int j = 0; j < PD_n; j++)
				cout << setw(PP_SETW) << PD_A[i][j];
			cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
		}
	}

	static inline void Print_HalfspacesIncludingPoint(PT_vector_T x, double eps) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHalfspace_i(x, i, eps)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_HalfspacesOutOfPoint(PT_vector_T x, double eps) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (!PointBelongsHalfspace_i(x, i, eps)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_HyperplanesIncludingPoint(PT_vector_T x, double eps) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHyperplane_i(x, i, eps)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_Vector(PT_vector_T x) {
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << x[j];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	}

	static inline void PseudoprojectionOnFlat(int* flatHyperplanes, int m_flat, PT_vector_T v, double eps, int maxProjectingIter, PT_vector_T w, int* success) {
		PT_vector_T r;
		PT_vector_T w_previous;
		double distSQR;
		int iterCount = 0;
		double eps_distSQR = (eps * eps) / 100;

		Vector_Copy(v, w);

		do {
			Vector_Zeroing(r);
			Vector_Copy(w, w_previous);

			for (int i = 0; i < m_flat; i++) {
				PT_vector_T p;
				OrthogonalProjectingVectorOntoHyperplane_i(w, flatHyperplanes[i], p);
				Vector_PlusEquals(r, p);
			}

			Vector_DivideEquals(r, m_flat);
			Vector_Round(r, eps);
			Vector_PlusEquals(w, r);

			distSQR = DistanceSQR_PointToPoint(w, w_previous);
			iterCount++;
			if (iterCount > maxProjectingIter) {
				*success = false;
				break;
			}
		} while (distSQR >= eps_distSQR);
	}

	static inline double RelativeError(double trueValue, double calculatedValue) {
		if (fabs(trueValue) >= PP_EPS_ZERO)
			return fabs(calculatedValue - trueValue) / fabs(trueValue);
		else
			return fabs(calculatedValue - trueValue);
	}

	static inline void Shift(PT_vector_T point, PT_vector_T shiftVector, double factor, PT_vector_T shiftedPoint) {
		for (int j = 0; j < PD_n; j++)
			shiftedPoint[j] = point[j] + shiftVector[j] * factor;
	}

	static inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] + y[j];
	}

	static inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
		for (int j = 0; j < PD_n; j++)
			toPoint[j] = fromPoint[j];
	}

	static inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = x/r
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] / r;
	}

	static inline void Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
		for (int j = 0; j < PD_n; j++)
			x[j] /= r;
	}

	static inline double Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
		double sum = 0;
		for (int j = 0; j < PD_n; j++)
			sum += x[j] * y[j];
		return sum;
	}

	static inline bool Vector_Is_Tiny(PT_vector_T x, double eps) {
		return Vector_Norm(x) < eps;
	}

	static inline void Vector_MakeLike(PT_vector_T x, double lengthOfLikeVector, PT_vector_T likeVector) {
		double norm_x = Vector_Norm(x);
		if (norm_x == 0)
			Vector_Zeroing(likeVector);
		else
			Vector_MultiplyByNumber(x, lengthOfLikeVector / norm_x, likeVector);
	}

	static inline void Vector_MakeMinus_e(PT_vector_T minus_e) {
		for (int j = 0; j < PD_n; j++)
			minus_e[j] = -1;
	}

	static inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
		for (int j = 0; j < PD_n; j++)
			equalPoint[j] -= minusVector[j];
	}

	static inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] * r;
	}

	static inline void Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
		for (int j = 0; j < PD_n; j++)
			x[j] *= r;
	}

	static inline double Vector_Norm(PT_vector_T x) {
		double norm_x = sqrt(Vector_NormSquare(x));
		return norm_x;
	}

	static inline double Vector_NormSquare(PT_vector_T x) {
		double sum = 0;

		for (int j = 0; j < PD_n; j++) {
			sum += x[j] * x[j];
		}
		return sum;
	}

	static inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
		for (int j = 0; j < PD_n; j++)
			equalVector[j] += plusVector[j];
	}

	static inline void Vector_Round(PT_vector_T x, double eps) {
		double floorValue;
		double fractionalPart;
		double sign;
		double absValue;

		for (int j = 0; j < PD_n; j++) {
			if (fabs(x[j]) < eps) {
				x[j] = 0;
				continue;
			}
			absValue = fabs(x[j]);
			sign = x[j] > 0 ? 1 : -1;
			floorValue = floor(absValue);
			fractionalPart = absValue - floorValue;
			if (1 - fractionalPart < eps) {
				x[j] = sign * (floorValue + 1);
				continue;
			}
			if (fractionalPart < eps)
				x[j] = sign * floorValue;
		}
	}

	static inline void Vector_SetValue(PT_vector_T x, double v) {  // x = (v,...,v)
		for (int j = 0; j < PD_n; j++) x[j] = v;
	}

	static inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] - y[j];
	}

	static inline void Vector_Zeroing(PT_vector_T x) {  // x = 0
		for (int j = 0; j < PD_n; j++) x[j] = 0;
	}
}

//---------------------------------- Private Functions -------------------------
namespace PF {
	using namespace SF;

	static inline void CodeToSubset(int code, int subset[PP_MM], int* ma) {
		*ma = 0;
		for (int i = 0; i < PD_mh; i++) {
			if (code % 2 == 1) {
				subset[*ma] = PD_pointHyperplanes[i];
				(*ma)++; assert(*ma <= PP_MM);
			}
			code /= 2;
			if (code == 0)
				break;
		}
	}

	static inline void MakeFaceList(int* faceCodeList, int K) {
		int index;

		for (int k = 0; k < PP_KK; k++) {
			faceCodeList[k] = 0;
		}

		if (PP_KK <= (PP_RND_MAX + 1)) {
			for (int k = 0; k < K; k++) {
				index = rand() % PP_KK;
				if (faceCodeList[index] == 0)
					faceCodeList[index] = k;
				else {
					for (int ki = 1; ki < K; ki++) {
						if (faceCodeList[(index + ki) % PP_KK] == 0) {
							faceCodeList[(index + ki) % PP_KK] = k;
							break;
						}
					}
				}
			}
			return;
		}

		assert(K <= PP_KK);

		int segmentCount = PP_KK / (PP_RND_MAX + 1);
		assert(PP_KK % (PP_RND_MAX + 1) == 0);

		if (K < segmentCount) {
			for (int k = 1; k < K; k++)
				faceCodeList[(k - 1) * (PP_RND_MAX + 1)] = k;
			return;
		}

		int segmentK = K / segmentCount;
		assert(K % segmentK == 0);

		for (int l = 0; l < segmentCount; l++) {
			for (int k = l * segmentK; k < (l + 1) * segmentK; k++) {
				index = rand() % (PP_RND_MAX + 1) + l * (PP_RND_MAX + 1);
				if (faceCodeList[index] == 0)
					faceCodeList[index] = k;
				else {
					for (int ki = 1; ki < (PP_RND_MAX + 1); ki++) {
						if (faceCodeList[l * (PP_RND_MAX + 1) + (index + ki) % (PP_RND_MAX + 1)] == 0) {
							faceCodeList[l * (PP_RND_MAX + 1) + (index + ki) % (PP_RND_MAX + 1)] = k;
							break;
						}
					}
				}
			}
		}
	}

	static inline void PreparationForIteration(PT_vector_T u) {
		MakePointHyperplaneList(u, PD_pointHyperplanes, &PD_mh, PP_EPS_POINT_IN_HALFSPACE);
		assert(PD_mh <= PP_MM);

		if (!PointBelongsPolytope(u, PP_EPS_POINT_IN_HALFSPACE)) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
				cout << "\nPoint does not belong to the feasible polytope with precision PP_EPS_ZERO = "
					<< PP_EPS_POINT_IN_HALFSPACE << "!!! You should decrease this parameter.\n";
				exit(1);
			}
		}

		PD_K = (int)pow(2, (double)PD_mh);
		if (PD_K > PP_KK) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Parameter PP_KK = " << PP_KK << " must be greater than or equal to " << PD_K << "\n";
			exit(1);
		}

		//
		MakeFaceList(PD_faceCodes, PD_K);
		//
		//
		//
		//
	}

	static inline void Print_Number_of_faces(PT_vector_T x) {
		int mh;
		int K;

		mh = 0;
		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHyperplane_i(x, i, PP_EPS_POINT_IN_HALFSPACE))
				mh++;
		}

		//
		//
		//
		//
		//
		//
		//
		//
		//
		//
		//

		K = (int)pow(2, (double)mh) - 1;
		cout << K << endl;
		//
	}
}