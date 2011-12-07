#ifndef CLIQUE_UTIL_H
#define CLIQUE_UTIL_H

int print_clique_callback(const _set *st, void *context)
{
  int i;
  //printf("Size = %2d found\n", st->size);
  for (i=0; i<st->size; i++) {
    //printf("%d ", st->vertex[i]);
  }
  //printf("\n");
  return CLIQUE_CONTINUE;
}

#endif
