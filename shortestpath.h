/*
 * shortestpath.h, by Marius Posta, 2009.
 * Check LICENSE.txt for the legal blah-blah.
 *
 * simple directed graph and residual graph data structure definitions for the CMND problem.
 *
 */


#ifndef SHORTESTPATH_H
#define SHORTESTPATH_H

#include "cmnd.h"


typedef struct
{
	cmnd_t* instp;
	int commodity; // commodity used for the arc costs
	double* arc_cost;

	int* node_flag; // required by Dijkstra's algorithm

	double* dist; // distance from node predecessor
	int* pred_node; // node predecessor
	int* pred_arc; // arc from node predecessor
} simple_graph_t;


simple_graph_t simple_graph_create (cmnd_t* instp);
void simple_graph_destroy (simple_graph_t* graphp);
void simple_graph_reset (simple_graph_t* graphp, int commodity); // resets all arc costs according to commodity
double simple_graph_solve_shortest_path (simple_graph_t* graphp); // applies Dijkstra's algorithm



typedef struct // neighboring node data structure
{
	int node; // node index
	int arc; // arc index
	double* distp; // pointer to arc cost value
} resnode_t;


typedef struct
{
	cmnd_t* instp;
	double* arc_cost; // arc cost for adding flow
	double* arc_cost_inv; // arc cost for subtracting flow

	resnode_t** node_outward; // accessible neighboring nodes
	int* node_n_outward;

	int list_size; // required by Ford-Bellman algorithm
	int* list; // idem

	double* dist; // distance from node predecessor
	int* pred_node; // node predecessor
	int* pred_arc; // arc from node predecessor
} residual_graph_t;


residual_graph_t residual_graph_create (cmnd_t* instp);
void residual_graph_destroy (residual_graph_t* graphp);
void residual_graph_reset (residual_graph_t* graphp); // resets all arc costs and node neighbors
void residual_graph_prune (residual_graph_t* graphp); // finds all node neighbors linked without infinite cost arcs
void residual_graph_solve_shortest_path (residual_graph_t* graphp, int source_node); // finds shortest or almost-shortest paths from source_node (Ford-Bellman)

#endif


