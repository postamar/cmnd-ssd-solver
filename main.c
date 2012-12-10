/*
 * main.c, by Marius Posta, 2009
 * check LICENSE.txt for legal blah-blah.
 */

#include "colgen.h"

GRBenv* env;

int main (int argc, char** argv)
{
	cmnd_t inst;
	colgen_t ssd;
	int grb_version_major, grb_version_minor, grb_version_technical;
	double best_z, z;
	int *design;

	// initialize Gurobi
	GRB(loadenv(&env, "rscmnd.log"));
	GRB(setintparam(env, "OutputFlag", 0));
	GRBversion(&grb_version_major, &grb_version_minor, &grb_version_technical);
    printf("Gurobi initialized, version %d.%d.%d.\n", grb_version_major, grb_version_minor, grb_version_technical);

	// load instance
	if (1 == argc) {
		fprintf(stderr, "ERROR: no instance filename passed as argument\n");
		return 1;
	}
	printf("Loading CMND instance %s\n", argv[1]);
	inst = cmnd_create(argv[1]);
	ssd = colgen_create(&inst);

	// perform slope-scaling descent
	printf("Performing slope-scaling descent...\n");
	fflush(stdout);
	best_z = z = colgen_slope_scaling(&ssd, NULL, 0.00001, 1000, 500, 1000, 5);
	printf("Best solution value: %.2f\n", z);

	design = calloc(inst.n_arcs, sizeof(int));
	while (1) {
		int i;

		for (i = 0; i < inst.n_arcs; i++)
			design[i] = (random() % 10) ? -1 : ((random() % 3) ? 0 : 1);
		z = colgen_slope_scaling(&ssd, design, 0.00001, 1000, 500, 1000, 5);
		printf("%.2f", z);
		if (best_z > z) {
			best_z = z;
			printf("*");
		}
		printf("\n");

	}

	// clean up
	colgen_destroy(&ssd);
	cmnd_destroy(&inst);
    GRBfreeenv(env);

	return 0;
}


