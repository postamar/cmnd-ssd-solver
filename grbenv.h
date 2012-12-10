/*
 * grbenv.h, by Marius Posta, 2009
 * Check LICENSE.txt for the legal blah-blah.
 * 
 * Gurobi LP solver header file.
 * 
 */

#ifndef GRBENV_H
#define GRBENV_H

#include "gurobi_c.h"
#include "assert.h"

/* error-checking macro borrowed from Paul Khuong (originally for CPLEX) */
#define GRB(CALL)             \
        do {                                                            \
                int status = (GRB ## CALL);                             \
                if (status) {                                           \
                        fprintf(stderr, "GRB error: %i\n", status);     \
                        assert(0);                                      \
                }                                                       \
        } while (0)


extern GRBenv* env; // our Gurobi environment, defined in main.c

#endif

