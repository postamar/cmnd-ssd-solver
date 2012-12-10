/*
* cmnd.c, by Marius Posta, 2009
* Check LICENSE.txt for the legal blah-blah.
*/

#include "cmnd.h"





void cmnd_destroy (cmnd_t* inst)
{
	int i;

	free(inst->arc_orig_node);
	free(inst->arc_dest_node);
	free(inst->arc_capacity);
	free(inst->arc_fcost);
	for (i = 0; i < inst->n_commods; i++) {
		free(inst->arc_commod_ucost[i]);
		free(inst->arc_commod_bound[i]);
	}
	free(inst->arc_commod_ucost);
	free(inst->arc_commod_bound);

	free(inst->commod_orig_node);
	free(inst->commod_dest_node);
	free(inst->commod_supply);
	free(inst->commod_overflow_ucost);

	free(inst->node_n_outgoing_arcs);
	free(inst->node_n_ingoing_arcs);
	for(i = 0; i < inst->n_nodes; i++) {
		free(inst->node_outgoing_arc[i]);
		free(inst->node_ingoing_arc[i]);
		free(inst->node_commod_supply[i]);
	}
	free(inst->node_outgoing_arc);
	free(inst->node_ingoing_arc);
	free(inst->node_commod_supply);
}


cmnd_t cmnd_alloc (int n_nodes, int n_arcs, int n_commods)
{
	cmnd_t inst;
	int i, j;

	inst.n_arcs = n_arcs;
	inst.n_nodes = n_nodes;
	inst.n_commods = n_commods;
	inst.avg_arc_fcost = 0.0;
	inst.min_arc_fcost = inf;
	inst.max_arc_fcost = -inf;

	inst.arc_orig_node = ((int*) calloc(inst.n_arcs, sizeof(int)));
	inst.arc_dest_node = ((int*) calloc(inst.n_arcs, sizeof(int)));
	inst.arc_capacity = ((double*) calloc(inst.n_arcs, sizeof(double)));
	inst.arc_fcost = ((double*) calloc(inst.n_arcs, sizeof(double)));
	inst.arc_commod_ucost = ((double**) calloc(inst.n_arcs, sizeof(double*)));
	inst.arc_commod_bound = ((double**) calloc(inst.n_arcs, sizeof(double*)));
	for (i = 0; i < inst.n_arcs; i++) {
		inst.arc_commod_ucost[i] = ((double*) calloc(inst.n_commods, sizeof(double)));
		inst.arc_commod_bound[i] = ((double*) calloc(inst.n_commods, sizeof(double)));
	}

	inst.commod_orig_node = ((int*) calloc(inst.n_commods, sizeof(int)));
	inst.commod_dest_node = ((int*) calloc(inst.n_commods, sizeof(int)));
	inst.commod_supply = ((double*) calloc(inst.n_commods, sizeof(double)));
	inst.commod_overflow_ucost = ((double*) calloc(inst.n_commods, sizeof(double)));

	inst.node_n_outgoing_arcs = ((int*) calloc(inst.n_nodes, sizeof(int)));
	inst.node_n_ingoing_arcs = ((int*) calloc(inst.n_nodes, sizeof(int)));
	inst.node_outgoing_arc = ((int**) calloc(inst.n_nodes, sizeof(int*)));
	inst.node_ingoing_arc = ((int**) calloc(inst.n_nodes, sizeof(int*)));
	inst.node_commod_supply = ((double**) calloc(inst.n_nodes, sizeof(double*)));
	for (i = 0; i < inst.n_nodes; i++) {
		inst.node_n_outgoing_arcs[i] = inst.node_n_ingoing_arcs[i] = 0;
		inst.node_outgoing_arc[i] = ((int*) calloc(inst.n_nodes, sizeof(int)));
		inst.node_ingoing_arc[i] = ((int*) calloc(inst.n_nodes, sizeof(int)));
		inst.node_commod_supply[i] = ((double*) calloc(inst.n_commods, sizeof(double)));
		for (j = 0; j < inst.n_commods; j++) 
			inst.node_commod_supply[i][j] = 0.0;
	}


	return inst;
}



cmnd_t cmnd_create (char* filename)
{
	FILE* f = fopen(filename, "r");
	cmnd_t inst;
	int i, j, n;
	int n_nodes, n_arcs, n_commods;
	int orig, dest;
	double ucost;
	char buf1[1024], buf2[1024];

	/* read file contents */ 
	n = fscanf(f, "%s \n", buf1);
	assert(1 == n);
	n = fscanf(f, "%i %i %i \n", &n_nodes, &n_arcs, &n_commods);
	assert(3 == n);

	inst = cmnd_alloc(n_nodes, n_arcs, n_commods);

	for (i = strlen(filename); i > 0; i--) 
		if (filename[i - 1] == '/')
			break;
	strcpy(inst.name, &filename[i]);
	for (i = strlen(inst.name) - 1; i >= 0; --i)
		if (inst.name[i] == '.') {
			inst.name[i] = '\0';
			break;
		}

	for (i = 0; i < n_arcs; i++) {
		n = fscanf(f, "%i %i %lf %lf %lf %s %s \n", &orig, &dest, &ucost, &inst.arc_capacity[i], &inst.arc_fcost[i], buf1, buf2);
		assert(7 == n);
		inst.avg_arc_fcost += inst.arc_fcost[i] / ((double) n_arcs);
		if (inst.arc_fcost[i] < inst.min_arc_fcost)
			inst.min_arc_fcost =  inst.arc_fcost[i];
		if (inst.arc_fcost[i] > inst.max_arc_fcost)
			inst.max_arc_fcost =  inst.arc_fcost[i];
		inst.arc_orig_node[i] = orig - 1;
		inst.arc_dest_node[i] = dest - 1;
		inst.node_outgoing_arc[inst.arc_orig_node[i]][inst.node_n_outgoing_arcs[inst.arc_orig_node[i]]++] = i;
		inst.node_ingoing_arc[inst.arc_dest_node[i]][inst.node_n_ingoing_arcs[inst.arc_dest_node[i]]++] = i;
		for (j = 0; j < n_commods; j++) {
			inst.arc_commod_ucost[i][j] = ucost;
		}
	}

	for (i = 0; i < n_commods; i++) {
		inst.commod_overflow_ucost[i] = 0.0;
		n = fscanf(f, "%i %i %lf \n", &orig, &dest, &inst.commod_supply[i]);
		assert(3 == n);
		inst.commod_orig_node[i] = orig - 1;
		inst.commod_dest_node[i] = dest - 1;
		inst.node_commod_supply[inst.commod_orig_node[i]][i] = inst.commod_supply[i];
		inst.node_commod_supply[inst.commod_dest_node[i]][i] = -inst.commod_supply[i];
		for (j = 0; j < n_arcs; j++) {
			inst.commod_overflow_ucost[i] += 10.0 * inst.arc_commod_ucost[j][i]; // 10 times the total unit costs of all arcs for this commodity 
			inst.arc_commod_bound[j][i] = (inst.arc_capacity[j] < inst.commod_supply[i]) ? inst.arc_capacity[j] : inst.commod_supply[i];
		}
	}

	fclose(f);

	return inst;
}


