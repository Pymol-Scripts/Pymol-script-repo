#ifndef CLIQUE_BK_H
#define CLIQUE_BK_H

#include "set.h"
#include "clique.h"
#include <stdio.h>
#include <stdlib.h>
#include "bk.h"


int bron_kerbosch(int N, char **connected, void *, Set *best,
		  clique_callback cb, void *context);

static int bk_extend_v2(char **connected, int *old, int ne, int ce,
			Set *compsub, Set *best,
			clique_callback, void *context);


int bron_kerbosch(int N, char **connected, void *params, Set *best,
		  clique_callback cb, void *context)
{
  int *ALL = (int *) malloc(N*sizeof(int));
  Set compsub;
  int c;
  int result;
  
  if (!init_Set(&compsub, N)) {
    free(ALL);
    return CLIQUE_ABORT;
  }
  
  best -> size = 0;
  
  for (c=0; c<N; c++) {
    ALL[c] = c;
  }
  
  result = bk_extend_v2(connected, ALL, 0, N, &compsub, best, cb, context);
  if (result == CLIQUE_CONTINUE) {
    /* even if it wasn't found, we can return an empty set */
    result = CLIQUE_FOUND;
  }

  del_Set(&compsub);
  free(ALL);
  return result;
}

static int bk_extend_v2(char **connected, int *old, int ne, int ce,
			Set *compsub, Set *best,
			clique_callback cb, void *context)
{
  int *new_ = (int *)malloc(ce*sizeof(int));
  int nod, fixp;
  int newne, newce, i, j, count, pos, p, s, sel, minnod;
  int result = CLIQUE_CONTINUE;

  minnod = ce;
  nod = 0;

  /* Determine each counter value and look for minimum */

  for (i=0; i<ce && minnod != 0; i++) {
    p = old[i];
    count = 0;

    /* Count disconnections */
    for (j=ne; j<ce && count < minnod; j++) {
      if (!connected[p][old[j]]) {
	count++;

	/* Save position of potential candidate */
	pos = j;
      }
    }

    /* Test new minimum */
    if (count < minnod) {
      fixp = p;
      minnod = count;
      if (i<ne) {
	s = pos;
      } else {
	s = i;
	/* preincr */
	nod = 1;
      }
    }
  }

  /* If fixed point initially chosen from candidates then */
  /* number of diconnections will be preincreased by one */

  /* Backtrackcycle */
  for (nod=minnod+nod; nod>=1; nod--) {
    /* Interchange */
    p = old[s];
    old[s] = old[ne];
    sel = old[ne] = p;

    /* Fill new set "not" */
    newne = 0;
    for (i=0; i<ne; i++) {
      if (connected[sel][old[i]]) {
	new_[newne++] = old[i];
      }
    }

    /* Fill new set "cand" */
    newce = newne;
    for (i=ne+1; i<ce; i++) {
      if (connected[sel][old[i]]) {
	new_[newce++] = old[i];
      }
    }
    /* Add to compsub */
    compsub->vertex[compsub->size++] = sel;

    if (newce == 0) {
      if (best->size < compsub->size) {
	/* found a max clique */
	copy_Set(compsub, best);
      }
      if (cb) {
	result = cb(compsub, context);
	switch (result) {
	case CLIQUE_CONTINUE:
	  break;
	case CLIQUE_FOUND:
	  /* copy the current values */
	  copy_Set(compsub, best);
	case CLIQUE_ABORT:
	  goto END_CLIQUE;
	default:
	  /* out or range */
	  result = CLIQUE_ABORT;
	  goto END_CLIQUE;
	}
      }

    } else {
      if (newne < newce) {
	result = bk_extend_v2(connected, new_, newne, newce, compsub,
			      best, cb, context);
	if (result != CLIQUE_CONTINUE) {
	  goto END_CLIQUE;
	}
      }
    }
    
    /* Remove from compsub */
    compsub->size--;
    
    /* Add to "nod" */
    ne++;
    if (nod > 1) {
      /* Select a candidate disconnected to the fixed point */
      for (s=ne; connected[fixp][old[s]]; s++) {
	/* nothing */
      }

    } /* end selection */
  } /* Backtrackcycle */

 END_CLIQUE:
  free(new_);

  return result;
}



#endif
