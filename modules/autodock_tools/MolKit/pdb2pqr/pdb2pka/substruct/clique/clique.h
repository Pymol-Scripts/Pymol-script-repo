#ifndef CLIQUE_CLIQUE_H
#define CLIQUE_CLIQUE_H

#include "set.h"

/* possible return values from the clique detection routines */
/* Must keep these < 0 to make life simpler with dfmax implementation */
/* and keep them != -1 so Python parsing is easier */
enum {
  CLIQUE_FOUND = -11,           /* clique found; success */
  CLIQUE_ABORT = -12,           /* system problem; no success */
  CLIQUE_CONTINUE = -13         /* nothing happened here */
};


/* Called whenever a clique is found */
/* must return one of the enums given earlier) */
typedef int (*clique_callback)(const _set *, void *context);

/* Generic form */
typedef int (*clique_solver)(int N, char **connection, void *params,
			     Set*, clique_callback, void *context);


#endif
