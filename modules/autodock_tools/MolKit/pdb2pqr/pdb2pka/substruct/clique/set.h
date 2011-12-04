#ifndef CLIQUE_SET_H
#define CLIQUE_SET_H

#include <stdio.h>
#include "set.h"

/* Sets are really simple things and are really just used to implement
a fixed size stack.  They contain an allocated block of memory (which
isn't ever resized) and the number of elements in the set. */

typedef struct _set {
  int size;
  int *vertex;
} Set;

/* allocate the given amount of memory for the set. */
/* Returns 0 if failure (that is, out of memory), 1 if success */
int init_Set(Set *set, int size);

/* free memory allocated for this set */
void del_Set(Set *set);

void copy_Set(const Set *from, Set *to);

/* returns 0 if failure, !0 if success */
/* Creates a set of a predefined size */
int init_Set(Set *set, int size)
{
  set->size = 0;
  set->vertex = (int *) malloc(size * sizeof(int));
  return set->vertex != NULL;
}
  
void del_Set(Set *set)
{
  set->size = 0;
  free(set->vertex);
  set->vertex = NULL;
}

void copy_Set(const Set *from, Set *to)
{
  int loc = from->size;
  const int *v1 = from->vertex;
  int *v2 = to->vertex;
  to->size = loc;

  while (loc-- > 0) {
    *v2++=*v1++;
  }
}

#endif
