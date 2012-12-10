/*
 * shortestpath.c, by Marius Posta, 2009
 *
 */

#include "shortestpath.h"



void simple_graph_destroy (simple_graph_t* graphp)
{
	free(graphp->arc_cost);
	free(graphp->dist);
	free(graphp->node_flag);
	free(graphp->pred_node);
	free(graphp->pred_arc);
}


simple_graph_t simple_graph_create (cmnd_t* instp)
{
	simple_graph_t graph;
	int n_nodes = instp->n_nodes;
	int n_arcs = instp->n_arcs;

	graph.commodity = -1;
	graph.instp = instp;
	graph.arc_cost = malloc(n_arcs * sizeof(double));
	graph.dist = malloc(n_nodes * sizeof(double));
	graph.node_flag = malloc(n_nodes * sizeof(int));
	graph.pred_node = (malloc(n_nodes * sizeof(double)));
	graph.pred_arc = (malloc(n_nodes * sizeof(double)));

	return graph;
}


void simple_graph_reset (simple_graph_t* graphp, int commodity)
// resets all arc costs according to commodity 
{
	int i, j, src, snk;

	graphp->commodity = commodity;

	for (i = 0; i < graphp->instp->n_arcs; i++) 
		graphp->arc_cost[i] = graphp->instp->arc_commod_ucost[i][commodity];

	for (j = 0; j < graphp->instp->n_nodes; j++) {
		graphp->node_flag[j] = 1;
		graphp->dist[j] = inf;
		graphp->pred_node[j] = -1;
		graphp->pred_arc[j] = -1;
	}

	// sets the artificial arc as linking the source and the sink
	// indentified by the value graphp->instp->n_arcs, which is not a valid arc index
	// this remains so until a better path is found
	src = graphp->instp->commod_orig_node[commodity];
	snk = graphp->instp->commod_dest_node[commodity];
	graphp->dist[src] = 0.0;
	graphp->dist[snk] = graphp->instp->commod_overflow_ucost[commodity];
	graphp->pred_node[snk] = src;
	graphp->pred_arc[snk] = graphp->instp->n_arcs;
}


double simple_graph_solve_shortest_path (simple_graph_t* graphp)
// applies Dijkstra's algorithm
{
	int i, j, k, l, c = -1;
	double cdist, alt;
	
	for (l = 0; l < graphp->instp->n_nodes; l++) {
		cdist = inf; 
		c = -1;
		for (j = 0; j < graphp->instp->n_nodes; j++)
			if (graphp->node_flag[j] && graphp->dist[j] < cdist) {
				c = j;
				cdist = graphp->dist[j];
			}

		if (c == -1)
			break;
			
		graphp->node_flag[c] = 0;
		for (k = 0; k < graphp->instp->node_n_outgoing_arcs[c]; k++) {
			i = graphp->instp->node_outgoing_arc[c][k];
			j = graphp->instp->arc_dest_node[i];
			alt = cdist + graphp->arc_cost[i];
			if (alt < graphp->dist[j]) {
				graphp->dist[j] = alt;
				graphp->pred_node[j] = c;
				graphp->pred_arc[j] = i;
			}
		}
	}

	return graphp->dist[graphp->instp->commod_dest_node[graphp->commodity]];
}





void residual_graph_destroy (residual_graph_t* graphp)
{
	int j;

	for (j = 0; j < graphp->instp->n_nodes; j++) 
		free(graphp->node_outward[j]);

	free(graphp->list);
	free(graphp->arc_cost);
	free(graphp->arc_cost_inv);
	free(graphp->dist);
	free(graphp->pred_node);
	free(graphp->pred_arc);
	free(graphp->node_n_outward);
	free(graphp->node_outward);
}


residual_graph_t residual_graph_create (cmnd_t* instp)
{
	residual_graph_t graph;
	int j;
	int n_nodes = instp->n_nodes;
	int n_arcs = instp->n_arcs;

	graph.list_size = 0;
	graph.instp = instp;
	graph.list = (malloc(n_nodes * sizeof(double)));
	graph.arc_cost = (malloc(n_arcs * sizeof(double)));
	graph.arc_cost_inv = (malloc(n_arcs * sizeof(double)));
	graph.dist = (malloc(n_nodes * sizeof(double)));
	graph.pred_node = (malloc(n_nodes * sizeof(double)));
	graph.pred_arc = (malloc(n_nodes * sizeof(double)));
	graph.node_n_outward = (malloc(n_nodes * sizeof(int)));
	graph.node_outward = (malloc(n_nodes * sizeof(resnode_t*)));

	for (j = 0; j < instp->n_nodes; j++) {
		graph.node_n_outward[j] = 0;
		graph.node_outward[j] = (malloc(2*n_nodes * sizeof(resnode_t)));
	}

	return graph;
}


int resnode_cmp (const void* a, const void* b) 
// comparison function for qsort, comparing neighbors, preferring those linked with the cheapest arcs
{
	resnode_t* rsap = (resnode_t*) a;
	resnode_t* rsbp = (resnode_t*) b;

	if (*rsap->distp == *rsbp->distp)
		return 0;
	else
		return ((*rsap->distp < *rsbp->distp) ? -1 : 1);
}


void residual_graph_reset (residual_graph_t* graphp)
// resets arc costs and node neighbors
{
	int i, j;

	for (i = 0; i < graphp->instp->n_arcs; i++) 
		graphp->arc_cost[i] = graphp->arc_cost_inv[i] = inf;

	for (j = 0; j < graphp->instp->n_nodes; j++) 
		graphp->node_n_outward[j] = 0;


	for (j = 0; j < graphp->instp->n_nodes; j++) {
		graphp->dist[j] = inf;
		graphp->pred_node[j] = -1;
		graphp->pred_arc[j] = -1;
	}

	graphp->list_size = 0;
}


void residual_graph_prune (residual_graph_t* graphp) 
// finds all node neighbors linked without infinite cost arcs
{
	int i, j, o, d;
	resnode_t resnode;

	for (i = 0; i < graphp->instp->n_arcs; i++) {
		o = graphp->instp->arc_orig_node[i];
		d = graphp->instp->arc_dest_node[i];
	
		if (graphp->arc_cost[i] != inf) {
			resnode.node = d;
			resnode.arc = i; 
			resnode.distp = &graphp->arc_cost[i];
			assert(graphp->node_n_outward[o] < 2*graphp->instp->n_nodes);
			graphp->node_outward[o][(graphp->node_n_outward[o])++] = resnode;
		}

		if (graphp->arc_cost_inv[i] != inf) {
			resnode.node = o;
			resnode.arc = i; 
			resnode.distp = &graphp->arc_cost_inv[i];
			assert(graphp->node_n_outward[d] < 2*graphp->instp->n_nodes);
			graphp->node_outward[d][(graphp->node_n_outward[d])++] = resnode;
		}
	}

	for (j = 0; j < graphp->instp->n_nodes; j++) 
		// for each node, sort the neighbors in order of ascending distance
		qsort(graphp->node_outward[j], graphp->node_n_outward[j], sizeof(resnode_t), resnode_cmp);
}


int is_in_list (residual_graph_t* graphp, int node_maybe_in_list) 
// returns true if node_maybe_in_list is in visited nodes list
{
	int i;
				
	for (i = 0; i < graphp->list_size; i++) 
		if (graphp->list[i] == node_maybe_in_list)
			return 1;

	return 0;
}


int is_in_path (residual_graph_t* graphp, int node_current, int node_maybe_in_path)
// returns true if node_maybe_in_path is in the path leading from the source to node_current
{
	int c;
	
	for (c = node_current; c != -1; c = graphp->pred_node[c])
		if (node_maybe_in_path == c)
			return 1;

	return 0;
}


void residual_graph_solve_shortest_path (residual_graph_t* graphp, int source_node) 
// finds shortest or almost-shortest paths from source_node
{
	int j, k, o;
	resnode_t resnode;

	// reset precedence information
	for (j = 0; j < graphp->instp->n_nodes; j++) {
		graphp->dist[j] = inf;
		graphp->pred_node[j] = -1;
		graphp->pred_arc[j] = -1;
	}
	graphp->list_size = 1;
	graphp->list[0] = source_node;
	graphp->dist[source_node] = 0.0;

	// apply Ford-Bellman
	while (graphp->list_size > 0) {
		o = graphp->list[--graphp->list_size];
		for (k = 0; k < graphp->node_n_outward[o]; k++) {
			resnode = graphp->node_outward[o][k];
			if (graphp->dist[resnode.node] > graphp->dist[o] + (*resnode.distp) && !is_in_path(graphp, o, resnode.node)) {
				graphp->dist[resnode.node] = graphp->dist[o] + (*resnode.distp);
				graphp->pred_node[resnode.node] = o;
				graphp->pred_arc[resnode.node] = resnode.arc;

				if (!is_in_list(graphp, resnode.node))
					graphp->list[graphp->list_size++] = resnode.node;
			}
		}
	}
}


