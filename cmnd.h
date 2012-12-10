/*
 * cmnd.h, by Marius Posta, 2009.
 *
 * Capacitated Multicommodity Network Design problem instance structure definition.
 */


#ifndef CMND_H
#define CMND_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

#ifndef INFINITY
	#define INFINITY DBL_MAX
#endif


static const double inf = INFINITY;

typedef uint8_t IPC_header_t;

static const IPC_header_t IPC_error = 0;
static const IPC_header_t IPC_nothing = 1;
static const IPC_header_t IPC_end_search = 2;
static const IPC_header_t IPC_new_z_ub = 3;
static const IPC_header_t IPC_design = 4;
static const IPC_header_t IPC_solution = 5;


typedef struct 
{
	char name[256];
	int n_nodes;
	int n_arcs;
	int n_commods;

	// arc attributes
	int* arc_orig_node;
	int* arc_dest_node;
	double* arc_capacity;
	double* arc_fcost; // fixed cost
	double** arc_commod_ucost; // unit cost
	double** arc_commod_bound; // min(arc capacity, commodity supply)

	// commodity attributes
	int* commod_orig_node; // source
	int* commod_dest_node; // sink
	double* commod_supply; 
	double* commod_overflow_ucost; // unit cost of artificial arc linking source to sink (unit cost = bigM, fixed cost = 0, capacity = inf)

	// node attributes
	int* node_n_outgoing_arcs; 
	int* node_n_ingoing_arcs;
	int** node_outgoing_arc;
	int** node_ingoing_arc;
	double** node_commod_supply; // +supply if source, -supply if sink, 0 otherwise

	// statistical info
	double avg_arc_fcost;
	double min_arc_fcost;
	double max_arc_fcost;
} cmnd_t;


cmnd_t cmnd_create (char* filename);
void cmnd_destroy (cmnd_t* inst);

#endif

