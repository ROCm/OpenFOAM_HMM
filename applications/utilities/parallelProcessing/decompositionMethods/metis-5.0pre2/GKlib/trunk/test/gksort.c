/*!
\file  gksort.c
\brief Testing module for the various sorting routines in GKlib

\date   Started 4/4/2007
\author George
\version\verbatim $Id: gksort.c 1421 2007-04-06 14:37:41Z karypis $ \endverbatim
*/

#include <GKlib.h>

#define N       10000

/*************************************************************************/
/*! Testing module for gk_?isort() routine */
/*************************************************************************/
void test_isort()
{
  gk_idx_t i;
  int array[N];

  /* test the increasing sort */
  printf("Testing iisort...\n");
  for (i=0; i<N; i++)
    array[i] = RandomInRange(123432);

  gk_iisort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i] > array[i+1])
      printf("gk_iisort error at index %jd [%d %d]\n", (intmax_t)i, array[i], array[i+1]);
  }


  /* test the decreasing sort */
  printf("Testing disort...\n");
  for (i=0; i<N; i++)
    array[i] = RandomInRange(123432);

  gk_disort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i] < array[i+1])
      printf("gk_disort error at index %jd [%d %d]\n", (intmax_t)i, array[i], array[i+1]);
  }

}


/*************************************************************************/
/*! Testing module for gk_?fsort() routine */
/*************************************************************************/
void test_fsort()
{
  gk_idx_t i;
  float array[N];

  /* test the increasing sort */
  printf("Testing ifsort...\n");
  for (i=0; i<N; i++)
    array[i] = RandomInRange(123432)/(1.0+RandomInRange(645323));

  gk_ifsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i] > array[i+1])
      printf("gk_ifsort error at index %jd [%f %f]\n", (intmax_t)i, array[i], array[i+1]);
  }


  /* test the decreasing sort */
  printf("Testing dfsort...\n");
  for (i=0; i<N; i++)
    array[i] = RandomInRange(123432)/(1.0+RandomInRange(645323));

  gk_dfsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i] < array[i+1])
      printf("gk_dfsort error at index %jd [%f %f]\n", (intmax_t)i, array[i], array[i+1]);
  }

}


/*************************************************************************/
/*! Testing module for gk_?idxsort() routine */
/*************************************************************************/
void test_idxsort()
{
  gk_idx_t i;
  gk_idx_t array[N];

  /* test the increasing sort */
  printf("Testing iidxsort...\n");
  for (i=0; i<N; i++)
    array[i] = RandomInRange(123432);

  gk_iidxsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i] > array[i+1])
      printf("gk_iidxsort error at index %jd [%jd %jd]\n", (intmax_t)i, array[i], array[i+1]);
  }


  /* test the decreasing sort */
  printf("Testing didxsort...\n");
  for (i=0; i<N; i++)
    array[i] = RandomInRange(123432);

  gk_didxsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i] < array[i+1])
      printf("gk_didxsort error at index %jd [%jd %jd]\n", (intmax_t)i, array[i], array[i+1]);
  }

}



/*************************************************************************/
/*! Testing module for gk_?ikvsort() routine */
/*************************************************************************/
void test_ikvsort()
{
  gk_idx_t i;
  gk_ikv_t array[N];

  /* test the increasing sort */
  printf("Testing iikvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432);
    array[i].val = i;
  }

  gk_iikvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key > array[i+1].key)
      printf("gk_iikvsort error at index %jd [%d %d] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }


  /* test the decreasing sort */
  printf("Testing dikvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432);
    array[i].val = i;
  }

  gk_dikvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key < array[i+1].key)
      printf("gk_dikvsort error at index %jd [%d %d] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }

}



/*************************************************************************/
/*! Testing module for gk_?fkvsort() routine */
/*************************************************************************/
void test_fkvsort()
{
  gk_idx_t i;
  gk_fkv_t array[N];

  /* test the increasing sort */
  printf("Testing ifkvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432)/(1.0+RandomInRange(645323));
    array[i].val = i;
  }

  gk_ifkvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key > array[i+1].key)
      printf("gk_ifkvsort error at index %jd [%f %f] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }


  /* test the decreasing sort */
  printf("Testing dfkvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432)/(1.0+RandomInRange(645323));
    array[i].val = i;
  }

  gk_dfkvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key < array[i+1].key)
      printf("gk_dfkvsort error at index %jd [%f %f] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }

}


/*************************************************************************/
/*! Testing module for gk_?dkvsort() routine */
/*************************************************************************/
void test_dkvsort()
{
  gk_idx_t i;
  gk_dkv_t array[N];

  /* test the increasing sort */
  printf("Testing idkvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432)/(1.0+RandomInRange(645323));
    array[i].val = i;
  }

  gk_idkvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key > array[i+1].key)
      printf("gk_idkvsort error at index %jd [%lf %lf] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }


  /* test the decreasing sort */
  printf("Testing ddkvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432)/(1.0+RandomInRange(645323));
    array[i].val = i;
  }

  gk_ddkvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key < array[i+1].key)
      printf("gk_ddkvsort error at index %jd [%lf %lf] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }

}


/*************************************************************************/
/*! Testing module for gk_?skvsort() routine */
/*************************************************************************/
void test_skvsort()
{
  gk_idx_t i;
  gk_skv_t array[N];
  char line[256];

  /* test the increasing sort */
  printf("Testing iskvsort...\n");
  for (i=0; i<N; i++) {
    sprintf(line, "%d", RandomInRange(123432));
    array[i].key = gk_strdup(line);
    array[i].val = i;
  }

  gk_iskvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (strcmp(array[i].key, array[i+1].key) > 0)
      printf("gk_idkvsort error at index %jd [%s %s] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }


  /* test the decreasing sort */
  printf("Testing dskvsort...\n");
  for (i=0; i<N; i++) {
    sprintf(line, "%d", RandomInRange(123432));
    array[i].key = gk_strdup(line);
    array[i].val = i;
  }

  gk_dskvsort(N, array);

  for (i=0; i<N-1; i++) {
    /*printf("%s\n", array[i].key);*/
    if (strcmp(array[i].key, array[i+1].key) < 0)
      printf("gk_ddkvsort error at index %jd [%s %s] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }

}


/*************************************************************************/
/*! Testing module for gk_?idxkvsort() routine */
/*************************************************************************/
void test_idxkvsort()
{
  gk_idx_t i;
  gk_idxkv_t array[N];

  /* test the increasing sort */
  printf("Testing iidxkvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432);
    array[i].val = i;
  }

  gk_iidxkvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key > array[i+1].key)
      printf("gk_iidxkvsort error at index %jd [%jd %jd] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }


  /* test the decreasing sort */
  printf("Testing didxkvsort...\n");
  for (i=0; i<N; i++) {
    array[i].key = RandomInRange(123432);
    array[i].val = i;
  }

  gk_didxkvsort(N, array);

  for (i=0; i<N-1; i++) {
    if (array[i].key < array[i+1].key)
      printf("gk_didxkvsort error at index %jd [%jd %jd] [%jd %jd]\n", (intmax_t)i, array[i].key, array[i+1].key, (intmax_t)array[i].val, (intmax_t)array[i+1].val);
  }

}




int main()
{
  test_isort();
  test_fsort();
  test_idxsort();

  test_ikvsort();
  test_fkvsort();
  test_dkvsort();
  test_skvsort();
  test_idxkvsort();
}

