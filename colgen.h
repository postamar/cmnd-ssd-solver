/*
 * colgen.h, by Marius Posta, 2009
 *
 * Column-generation data structure for solving a Capacitated Multicommodity Network Flow problem associated to the CMND problem.
 *
 */


#ifndef COLGEN_H
#define COLGEN_H

#include "cmnd.h"
#include "shortestpath.h"
#include "grbenv.h"


typedef struct
{
	cmnd_t* instp; // CMND instance used
	simple_graph_t* graph; // graph used for pricing new columns

	GRBmodel* lp; // master problem
	int overflow; // non-zero if at least one commodity overflows on the artificial arc from its source to its sink (i.e. the design is infeasible)

	double z_ub; // upper bound on the optimal objective value (obtained by querying the LP solver of the master problem)
	double z_lb; // lower bound on the optimal objective value (obtained by solving a lagrangian relaxation)
	double* arc_cap; // current arc capacities
	double* arc_mul; // current arc multipliers (i.e. multipliers of the capacity constraints)
	double* overflow_mul; // current arc multipliers of the artificial arcs
	double* commod_mul; // current commodity multipliers (i.e. multipliers of the bundle constraints)
	double** commod_arc_lfactor; // linearization factors for each arc for each commodity, used by the slope-scaling descent algorithm

	double* col_val; // values in a column of the constraint coefficient matrix of the master problem (either 0.0 or 1.0)
	int* col_ind; // indices of the same

	int prev_alloc; // memory allocation counter for storing column (i.e. path) attributes
	int current_alloc; // idem
	int n_paths; // current number of columns (i.e. paths)

	double* path_ucost; // cost for a unit of flow on each path
	double* path_flow; // amount of flow on each path
	int* path_commod; // commodity which is carried on each path
	char** path_arcs; // 1 if a given arc is on a path, 0 if not

	double* arc_total_flow; // sum of all flows routed through each arc
	double* arc_total_ucost; // sum of all unit costs paid for each arc
	double** commod_arc_flow; // sum of all flows of each commodity routed through each arc
	int* arc_open; // 1 if a given arc carries any flow, 0 if not
	int* best_arc_open; // best arc_open configuration found during slope-scaling descent
} colgen_t;



colgen_t colgen_create (cmnd_t* instp);

void colgen_destroy (colgen_t* colgenp);

void colgen_change_arc_capacities (colgen_t* colgenp, int* arc_open); // modifies right-hand side of current master problem according to values in arc_open

double colgen_solve // solves a capacitated multicommodity network flow problem
	(colgen_t* colgenp, 
	double tolerance, // numerical error tolerance for comparing the lower and upper bounds on the optimum objective value
	int max_iter, // maximum number of column generation iterations 
	int max_paths); // maximum number of columns in the master problem

double colgen_slope_scaling // slope-scaling heuristic
	(colgen_t* colgenp, 
	 int* arc_open, // initial arc configuration: -1 if arc is free, 0 if not allowed open, 1 if forced open; if arc_open is NULL, then all arcs are free
	 double tolerance, // numerical error tolerance for each internal call to colgen_solve
	 int max_iter, // maximum number of column generation iterations for each internal call to colgen_solve
	 int max_paths, // maximum number of columns in the master problem for each internal call to colgen_solve
	 int max_ss_iter, // maximum number of slope-scaling descent iterations
	 int max_ss_iter_no_improv); // maximum number of consecutive non-improving slope-scaling descent iterations

#endif

