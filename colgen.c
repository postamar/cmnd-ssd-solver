/*
 * colgen.c, by Marius Posta, 2009
 * Check LICENSE.txt for the legal blah-blah.
 */


#include "colgen.h"

int grb_status; // status variable for Gurobi


void colgen_destroy (colgen_t* colgenp)
{
	int k;

	for (k = 0; k < colgenp->current_alloc; k++) 
		free(colgenp->path_arcs[k]);
	free(colgenp->path_arcs);

	for (k = 0; k < colgenp->instp->n_commods; k++) {
		free(colgenp->commod_arc_flow[k]);
		free(colgenp->commod_arc_lfactor[k]);
		simple_graph_destroy(&colgenp->graph[k]);
	}
	GRB(freemodel(colgenp->lp));
	free(colgenp->graph);
	free(colgenp->commod_arc_lfactor);
	free(colgenp->arc_cap);
	free(colgenp->arc_mul);
	free(colgenp->overflow_mul);
	free(colgenp->commod_mul);
	free(colgenp->col_val);
	free(colgenp->col_ind);
	free(colgenp->path_ucost);
	free(colgenp->path_flow);
	free(colgenp->path_commod);
	free(colgenp->arc_total_flow);
	free(colgenp->arc_total_ucost);
	free(colgenp->arc_open);
	free(colgenp->best_arc_open);
	free(colgenp->commod_arc_flow);
	
}


colgen_t colgen_create (cmnd_t* instp)
{
	colgen_t colgen;
	int i, k;
	char* sense = (char*) (malloc((instp->n_commods + instp->n_arcs) * sizeof(char)));

	colgen.instp = instp;
	colgen.graph = (malloc(instp->n_commods * sizeof(simple_graph_t)));
	colgen.overflow = 1;
	colgen.n_paths = colgen.prev_alloc = colgen.current_alloc = instp->n_commods;
	colgen.arc_cap = (malloc((instp->n_arcs + instp->n_commods) * sizeof(double)));
	colgen.arc_mul = (malloc(instp->n_arcs * sizeof(double)));
	colgen.overflow_mul = (malloc(instp->n_commods* sizeof(double)));
	colgen.commod_mul = (malloc(instp->n_commods * sizeof(double)));
	colgen.col_val = (malloc((2 + instp->n_arcs) * sizeof(double)));
	colgen.col_ind = (malloc((2 + instp->n_arcs) * sizeof(int)));
	colgen.path_ucost = (malloc(instp->n_commods * sizeof(double)));
	colgen.path_flow = (malloc(instp->n_commods * sizeof(double)));
	colgen.path_commod = (malloc(instp->n_commods * sizeof(int)));
	colgen.arc_total_flow = (malloc(instp->n_arcs * sizeof(double)));
	colgen.arc_total_ucost = (malloc(instp->n_arcs * sizeof(double)));
	colgen.arc_open = (malloc(instp->n_arcs * sizeof(int)));
	colgen.best_arc_open = (malloc(instp->n_arcs * sizeof(int)));
	colgen.commod_arc_flow = (malloc(instp->n_commods * sizeof(double*)));
	colgen.path_arcs = (malloc(colgen.current_alloc * sizeof(char*)));
	colgen.commod_arc_lfactor = (malloc(instp->n_commods * sizeof(double*)));

	for (i = 0; i < instp->n_arcs; i++) {
		sense[i] = GRB_LESS_EQUAL; // sense of capacity constraints
		colgen.arc_cap[i] = instp->arc_capacity[i]; // initial right-hand sides of capacity constraints
		colgen.col_val[i] = 1.0;
	}
	colgen.col_val[instp->n_arcs] = 1.0;
	for (k = 0; k < instp->n_commods; k++) {
		colgen.commod_arc_lfactor[k] = (malloc(instp->n_arcs * sizeof(double)));
		for (i = 0; i < instp->n_arcs; i++)
			colgen.commod_arc_lfactor[k][i] = 0.0; // default linearization factors
		colgen.path_arcs[k] = (calloc(instp->n_arcs, sizeof(char)));
		sense[instp->n_arcs + k] = GRB_EQUAL; // sense of bundle constraints
		colgen.arc_cap[instp->n_arcs + k] = instp->commod_supply[k]; // right-hand sides of bundle constraints
		colgen.graph[k] = simple_graph_create(instp); // pricing graph creation
		colgen.commod_arc_flow[k] = (malloc(instp->n_arcs * sizeof(double)));
	}

	// master problem creation
	GRB(newmodel(env, &colgen.lp, "CMNF path-based formulation", 0, NULL, NULL, NULL, NULL, NULL));
	// the first n_arcs constraints are the capacity constraints
	GRB(addconstrs(colgen.lp, instp->n_arcs, 0, NULL, NULL, NULL, sense, colgen.arc_cap, NULL)); 
	// followed by n_commods bundle constraints
	GRB(addconstrs(colgen.lp, instp->n_commods, 0, NULL, NULL, NULL, &sense[instp->n_arcs], instp->commod_supply, NULL));
	GRB(updatemodel(colgen.lp));

	// overflow paths creation (paths containing only the artificial arc from source to sink for each commodity)
	for (k = 0; k < instp->n_commods; k++) {
		colgen.path_ucost[k] = instp->commod_overflow_ucost[k];
		colgen.path_flow[k] = 0.0;
		colgen.path_commod[k] = k;
		colgen.col_ind[0] = instp->n_arcs + k;
		GRB(addvar(colgen.lp, 1, colgen.col_ind, colgen.col_val, colgen.path_ucost[k], 0.0, instp->commod_supply[k], GRB_CONTINUOUS, NULL));
	}
	GRB(updatemodel(colgen.lp));

	// solver parameter settings
	GRB(setintparam(env, GRB_INT_PAR_METHOD, 1));
	// optimization
	GRB(optimize(colgen.lp));
	GRB(getintattr(colgen.lp, GRB_INT_ATTR_STATUS, &grb_status));
	assert(grb_status == GRB_OPTIMAL);
	// initial primal and dual values
	GRB(getdblattrarray(colgen.lp, GRB_DBL_ATTR_X, 0, colgen.n_paths, colgen.path_flow));
	GRB(getdblattrarray(colgen.lp, GRB_DBL_ATTR_PI, 0, instp->n_arcs, colgen.arc_mul));
	GRB(getdblattrarray(colgen.lp, GRB_DBL_ATTR_PI, instp->n_arcs, instp->n_commods, colgen.commod_mul));

	// initial bounds on the optimal objective value
	colgen.z_lb = 0.0;
	colgen.z_ub = inf;

	return colgen;
}


int arc_cmp (const void* a, const void* b)
// comparison function for sorting arc indices in ascending order
{
	int na = *((int*) a);
	int nb = *((int*) b);

	if (na == nb) 
		return 0;
	else 
		return (na < nb) ? -1 : 1;
}


void colgen_compute_shortest_paths (colgen_t* colgenp)
// column pricing subproblem for each commodity
{
	int i, k;

	for (k = 0; k < colgenp->instp->n_commods; k++) {
		simple_graph_reset(&colgenp->graph[k], k);

		for (i = 0; i < colgenp->instp->n_arcs; i++) {
			// the arc multiplier is subtracted from the unit cost
			// in case of a slope-scaling descent, the linearization factor is added
			colgenp->graph[k].arc_cost[i] += colgenp->commod_arc_lfactor[k][i] - colgenp->arc_mul[i];
			assert(colgenp->graph[k].arc_cost[i] >= 0.0);
		}
		
		simple_graph_solve_shortest_path(&colgenp->graph[k]);
	}
}


int colgen_generate (colgen_t* colgenp, int commodity)
// attempts to generate a new column for the master problem, with a negative reduced cost, 
// returns 1 on success, 0 on failure.
{
	int i, c;
	int size, col_nzcnt = 0;
	double ucost = 0.0; // unit cost of the path

	// walk through the shortest path
	for (c = colgenp->instp->commod_dest_node[commodity]; c != colgenp->instp->commod_orig_node[commodity]; c = colgenp->graph[commodity].pred_node[c]) {
		i = colgenp->graph[commodity].pred_arc[c];
		if (i == colgenp->instp->n_arcs) // test if the shortest path has the artificial arc
			return 0;
		ucost += colgenp->instp->arc_commod_ucost[i][commodity] + colgenp->commod_arc_lfactor[commodity][i]; // update path unit cost
		colgenp->col_ind[col_nzcnt++] = i; // add index to constraint coefficient matrix column
	}

	qsort(colgenp->col_ind, col_nzcnt, sizeof(int), arc_cmp); // sort column indices in ascending order
	colgenp->col_ind[col_nzcnt++] = colgenp->instp->n_arcs + commodity; // add bundle constraint coefficient

	// reallocate array space for path attributes, using a fibonacci sequence: new_size = old_size + older_size; 
	if (colgenp->n_paths == colgenp->current_alloc) { 
		size = colgenp->current_alloc;
		colgenp->current_alloc += colgenp->prev_alloc;
		colgenp->prev_alloc = size;
		colgenp->path_ucost = (realloc(colgenp->path_ucost, colgenp->current_alloc * sizeof(double)));
		colgenp->path_flow = (realloc(colgenp->path_flow, colgenp->current_alloc * sizeof(double)));
		colgenp->path_commod = (realloc(colgenp->path_commod, colgenp->current_alloc * sizeof(int)));
		colgenp->path_arcs = (realloc(colgenp->path_arcs, colgenp->current_alloc * sizeof(char*)));
		for (c = size; c < colgenp->current_alloc; c++) 
			colgenp->path_arcs[c] = (calloc(colgenp->instp->n_arcs, sizeof(char)));
	}

	// store the arc composition of the new column 
	for (i = 0; i < colgenp->instp->n_arcs; i++)
		colgenp->path_arcs[colgenp->n_paths][i] = 0;
	for (c = 0; c < col_nzcnt - 1; c++) {
		i = colgenp->col_ind[c];
		colgenp->path_arcs[colgenp->n_paths][i] = 1;
	}
	// store its other attributes
	colgenp->path_ucost[colgenp->n_paths] = ucost;
	colgenp->path_flow[colgenp->n_paths] = 0.0;
	colgenp->path_commod[colgenp->n_paths] = commodity;

	++colgenp->n_paths;

	// add the column
	GRB(addvar(colgenp->lp, col_nzcnt, colgenp->col_ind, colgenp->col_val, ucost, 0.0, colgenp->instp->commod_supply[commodity], GRB_CONTINUOUS, NULL));
	return 1;
}


int colgen_iteration (colgen_t* colgenp)
// performs one iteration of the column-generation algorithm
{
	int i, k, flag;
	int newcols = 0;

	// solve the column pricing subproblem
	colgen_compute_shortest_paths(colgenp); 

	// generate new columns with non-positive reduced costs
	for (k = 0; k < colgenp->instp->n_commods; k++) 
		if (colgenp->graph[k].dist[colgenp->instp->commod_dest_node[k]] <= colgenp->commod_mul[k]) 
			newcols += colgen_generate(colgenp, k); 

	// compute the lower bound on the optimal objective value
	colgenp->z_lb = 0.0;
	for (i = 0; i < colgenp->instp->n_arcs; i++) 
		colgenp->z_lb += colgenp->arc_mul[i] * colgenp->arc_cap[i];
	for (k = 0; k < colgenp->instp->n_commods; k++) {
		colgenp->z_lb += colgenp->graph[k].dist[colgenp->instp->commod_dest_node[k]] * colgenp->instp->commod_supply[k];
	}

	if (newcols > 0) { 
		// reoptimize
		GRB(updatemodel(colgenp->lp));
		GRB(setintparam(env, GRB_INT_PAR_METHOD, 0));
		GRB(optimize(colgenp->lp));
		GRB(getintattr(colgenp->lp, GRB_INT_ATTR_STATUS, &grb_status));
		assert(grb_status == GRB_OPTIMAL);
		// retrieve simplex multipliers
		GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_PI, 0, colgenp->instp->n_arcs, colgenp->arc_mul));
		GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_PI, colgenp->instp->n_arcs, colgenp->instp->n_commods, colgenp->commod_mul));
		if (colgenp->overflow) { 
			// retrieve artificial arc flow values
			GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_X, 0, colgenp->instp->n_commods, colgenp->path_flow));
			flag = 1;
			for (k = 0; k < colgenp->instp->n_commods && flag; k++) 
				if (colgenp->path_flow[k] > 0.0) 
					flag = 0;
			// if all are zero, then the solution is primal-feasible
			if (flag) 
				colgenp->overflow = 0;
		}
	}

	// retrieve the upper bound on the optimal objective value
	GRB(getdblattr(colgenp->lp, GRB_DBL_ATTR_OBJVAL, &colgenp->z_ub));

	//printf("%.40f <= Z* <= %.40f\n", colgenp->z_lb, colgenp->z_ub);
	//fflush(stdout);

	return newcols;
}


double colgen_add_flows (colgen_t* colgenp)
// compute arc flow values from the optimally-solved path-based formulation
{
	int i, k, p;
	double z = 0.0;

	// retrieve all primal values, i.e. path flow values
	GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_X, 0, colgenp->n_paths, colgenp->path_flow));

	// reset arc flows and costs
	for (k = 0; k < colgenp->instp->n_commods; k++) 
		for (i = 0; i < colgenp->instp->n_arcs; i++) 
			colgenp->commod_arc_flow[k][i] = 0.0;
	for (i = 0; i < colgenp->instp->n_arcs; i++) 
		colgenp->arc_total_ucost[i] = colgenp->arc_total_flow[i] = 0.0;

	for (p = colgenp->instp->n_commods; p < colgenp->n_paths; p++)
		if (colgenp->path_flow[p] > 0.0) 
			for (i = 0; i < colgenp->instp->n_arcs; i++) 
				if (colgenp->path_arcs[p][i]) {
					// add flow values to arc
					colgenp->commod_arc_flow[colgenp->path_commod[p]][i] += colgenp->path_flow[p]; 
					colgenp->arc_total_flow[i] += colgenp->path_flow[p];
					// add flow cost to arc
					colgenp->arc_total_ucost[i] += colgenp->path_flow[p] * colgenp->instp->arc_commod_ucost[i][colgenp->path_commod[p]];
				}

	// deduce the state of the arcs and compute the total cost
	for (i = 0; i < colgenp->instp->n_arcs; i++) {
		colgenp->arc_open[i] = 0;
		if (colgenp->arc_total_flow[i] > 0.0) {
			colgenp->arc_open[i] = 1;
			z += colgenp->instp->arc_fcost[i] + colgenp->arc_total_ucost[i];
		}
	}
	
	return z;
}


void colgen_change_arc_capacities (colgen_t* colgenp, int* arc_open)
// modifies right-hand side of current master problem according to values in arc_open
{
	int i;
	int nchangec = 0; // amount of open arcs which we closed
	int nchangeo = 0; // amount of closed arcs which we opened

	// close open arcs
	for (i = 0; i < colgenp->instp->n_arcs; i++) {
		if (!arc_open[i] && colgenp->arc_cap[i] > 0.0) {
			colgenp->col_val[nchangec] = colgenp->arc_cap[i] = 0.0;
			colgenp->col_ind[nchangec++] = i;
		}
	}

	if (nchangec) {
		// assume the solution to be infeasible
		if (colgenp->overflow == 0) 
			colgenp->overflow = 1;
		// update right-hand sides of affected capacity constraints and reoptimize
		GRB(setdblattrlist(colgenp->lp, GRB_DBL_ATTR_RHS, nchangec, colgenp->col_ind, colgenp->col_val));
		GRB(updatemodel(colgenp->lp));
		GRB(setintparam(env, GRB_INT_PAR_METHOD, 1));
		GRB(optimize(colgenp->lp));
		GRB(getintattr(colgenp->lp, GRB_INT_ATTR_STATUS, &grb_status));
		assert(grb_status == GRB_OPTIMAL);
	}

	// open closed arcs
	for (i = 0; i < colgenp->instp->n_arcs; i++) {
		if (arc_open[i] && colgenp->arc_cap[i] == 0.0) {
			colgenp->col_val[nchangeo] = colgenp->arc_cap[i] = colgenp->instp->arc_capacity[i];
			colgenp->col_ind[nchangeo++] = i;
		}
	}

	if (nchangeo) {
		// update right-hand sides of affected capacity constraints and reoptimize
		GRB(setdblattrlist(colgenp->lp, GRB_DBL_ATTR_RHS, nchangeo, colgenp->col_ind, colgenp->col_val));
		GRB(updatemodel(colgenp->lp));
		GRB(setintparam(env, GRB_INT_PAR_METHOD, 0));
		GRB(optimize(colgenp->lp));
		GRB(getintattr(colgenp->lp, GRB_INT_ATTR_STATUS, &grb_status));
		assert(grb_status == GRB_OPTIMAL);
	}

	if (nchangec || nchangeo) { // if any change was made
		// undo changes to col_val
		for (i = 0; i <= nchangec || i <= nchangeo; i++) 
			colgenp->col_val[i] = 1.0;
		// retrieve simplex multipliers
		GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_PI, 0, colgenp->instp->n_arcs, colgenp->arc_mul));
		GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_PI, colgenp->instp->n_arcs, colgenp->instp->n_commods, colgenp->commod_mul));
		// retrieve artificial arc path flows
		GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_X, 0, colgenp->instp->n_commods, colgenp->path_flow));
		// reset bounds
		colgenp->z_lb = 0.0;
		colgenp->z_ub = inf;
	}
}



void colgen_trim (colgen_t* colgenp)
// removes all columns from the master problem except those corresponding to the artificial arc paths
{
	int i;
	
	// assume infeasibility
	colgenp->overflow = 1;
	// delete the added columns
	for (i = colgenp->instp->n_commods; i < colgenp->n_paths; i++)
		GRB(delvars(colgenp->lp, 1, &i));
	GRB(updatemodel(colgenp->lp));
	colgenp->n_paths = colgenp->instp->n_commods;
	// reoptimize
	GRB(setintparam(env, GRB_INT_PAR_METHOD, 1));
	GRB(optimize(colgenp->lp));
	GRB(getintattr(colgenp->lp, GRB_INT_ATTR_STATUS, &grb_status));
	assert(grb_status == GRB_OPTIMAL);
	// retrieve simplex multipliers
	GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_PI, 0, colgenp->instp->n_arcs, colgenp->arc_mul));
	GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_PI, colgenp->instp->n_arcs, colgenp->instp->n_commods, colgenp->commod_mul));
	// retrieve artificial arc path flows
	GRB(getdblattrarray(colgenp->lp, GRB_DBL_ATTR_X, 0, colgenp->n_paths, colgenp->path_flow));
	// reset bounds
	colgenp->z_lb = 0.0;
	colgenp->z_ub = inf;
}


double colgen_solve (colgen_t* colgenp, double tolerance, int max_iter, int max_paths)
// solves a capacitated multicommodity network flow problem
{
	int iter;

	// reset bounds
	colgenp->z_ub = inf;
	colgenp->z_lb = 0.0;

	// remove all added columns if the column limit is exceeded
	if (colgenp->n_paths > max_paths)
		colgen_trim(colgenp);

	// add columns, until the gap falls below the tolerance, or the iteration limit is reached
	for (iter = 0; iter < max_iter && colgenp->z_lb * (1.0 + tolerance) < colgenp->z_ub; iter++) 
		if (!colgen_iteration(colgenp)) 
			break;
	
	// recompute the node potential multipliers
	colgen_compute_shortest_paths(colgenp);

	// add all path flows on each arc
	return colgen_add_flows(colgenp);
}



int colgen_update_lfactors (colgen_t* colgenp)
// updates the linearization factors for the slope-scaling descent heuristic
{
	int i, k, nupd = 0;

	for (k = 0; k < colgenp->instp->n_commods; k++) 
		for (i = 0; i < colgenp->instp->n_arcs; i++)
			if (colgenp->commod_arc_flow[k][i] > 0.0) {
				colgenp->commod_arc_lfactor[k][i] = colgenp->instp->arc_fcost[i] / colgenp->arc_total_flow[i];
				++nupd;
			}

	return nupd;
}


double colgen_slope_scaling(colgen_t* colgenp, int* arc_open, double tolerance, int max_iter, int max_paths, int max_ss_iter, int max_ss_iter_no_improv)
// slope-scaling heuristic implementation
{
	double lfactor, z_ub, old_z_ub = inf, best_z_ub = inf;
	int i, k, iter;
	int no_improv = 0; // counter for consecutive iterations with no improvement

	// remove all added columns
	colgen_trim(colgenp);

	if (arc_open == NULL) {
		// set the linearization factors to their initial values
		for (i = 0; i < colgenp->instp->n_arcs; i++) {
			lfactor = (arc_open == NULL || arc_open[i] == -1) ? (colgenp->instp->arc_fcost[i] / colgenp->instp->arc_capacity[i]) : 0.0;
			for (k = 0; k < colgenp->instp->n_commods; k++)
				colgenp->commod_arc_lfactor[k][i] = lfactor;
		}
	}
	else {
		// set the linearization factors for open arcs to 0
		for (i = 0; i < colgenp->instp->n_arcs; i++) 
			if (arc_open[i] == 1)
				for (k = 0; k < colgenp->instp->n_commods; k++)
					colgenp->commod_arc_lfactor[k][i] = 0.0;
		// remove all forbidden arcs
		colgen_change_arc_capacities(colgenp, arc_open);
	}


	// perform slope-scaling descent, until iteration limits are reached or the linearization factors remain the same
	for (iter = 0; iter < max_ss_iter && no_improv < max_ss_iter_no_improv; iter++) {
		// solve the CMNF problem for the upper bound on the objective value
		z_ub = colgen_solve(colgenp, tolerance, max_iter, max_paths);
		// update the no-improvement counter
		no_improv = (z_ub == old_z_ub) ? (no_improv + 1) : 0;
		old_z_ub = z_ub;
		if (z_ub < best_z_ub) {
			// replace the best arc configuration
			best_z_ub = z_ub;
			memcpy(colgenp->best_arc_open, colgenp->arc_open, colgenp->instp->n_arcs * sizeof(int));
		}
		// remove all added columns
		colgen_trim(colgenp);
		// update the linearization factors
		if (!colgen_update_lfactors(colgenp))
			break;
	}

	// reset the linearization factors to zero
	for (k = 0; k < colgenp->instp->n_commods; k++)
		for (i = 0; i < colgenp->instp->n_arcs; i++) 
			colgenp->commod_arc_lfactor[k][i] = 0.0;

	// remove all added columns
	colgen_trim(colgenp);
	// open and close the arcs according to the best configuration
	colgen_change_arc_capacities(colgenp, colgenp->best_arc_open);
	// solve the CMNF
	return colgen_solve(colgenp, tolerance, max_iter, max_paths);
}



void colgen_check_flows (colgen_t* colgenp)
// debugging function, checks that the numbers add up
{
	int i, j, k, l;
	double f, z = 0.0;
	int nofail = 1;

	for (i = 0; i < colgenp->instp->n_arcs; i++) {
		f = 0.0;
		for (k = 0; k < colgenp->instp->n_commods; k++) 
			f += colgenp->commod_arc_flow[k][i];
		assert(f == colgenp->arc_total_flow[i]);
	}


	for (k = 0; k < colgenp->instp->n_commods; k++) {
		for (i = 0; i < colgenp->instp->n_arcs; i++) 
			if (colgenp->commod_arc_flow[k][i] > 0.0) 
				z += colgenp->commod_arc_flow[k][i] * colgenp->instp->arc_commod_ucost[i][k];

		for (j = 0; j < colgenp->instp->n_nodes; j++) {
			f = colgenp->instp->node_commod_supply[j][k];
			if (f < 0.0) 
				f += colgenp->path_flow[k];
			if (f > 0.0) 
				f -= colgenp->path_flow[k];
			for (l = 0; l < colgenp->instp->node_n_outgoing_arcs[j]; l++) {
				i = colgenp->instp->node_outgoing_arc[j][l];
				f -= colgenp->commod_arc_flow[k][i];
			}
			for (l = 0; l < colgenp->instp->node_n_ingoing_arcs[j]; l++) {
				i = colgenp->instp->node_ingoing_arc[j][l];
				f += colgenp->commod_arc_flow[k][i];
			}
			if (f != 0.0) {
				printf("Flow conservation violation of %f for commodity %i at node %i\n", f, k, j);
				printf("\tsupply = %f\n", colgenp->instp->node_commod_supply[j][k]);
				printf("\toverflow = %f\n", colgenp->path_flow[k]);
				for (l = 0; l < colgenp->instp->node_n_outgoing_arcs[j]; l++) {
					i = colgenp->instp->node_outgoing_arc[j][l];
					f = colgenp->commod_arc_flow[k][i];
					if (f != 0.0)
						printf("\toutgoing =(%i)=> %i  : %f\n", i, colgenp->instp->arc_dest_node[i], f);
				}
				for (l = 0; l < colgenp->instp->node_n_ingoing_arcs[j]; l++) {
					i = colgenp->instp->node_ingoing_arc[j][l];
					f = colgenp->commod_arc_flow[k][i];
					if (f != 0.0)
						printf("\tingoing =(%i)=> %i  : %f\n", i, colgenp->instp->arc_orig_node[i], f);
				}
				nofail = 0;
			}
		}
	}

	assert(z == colgenp->z_ub);
	assert(nofail);
}


