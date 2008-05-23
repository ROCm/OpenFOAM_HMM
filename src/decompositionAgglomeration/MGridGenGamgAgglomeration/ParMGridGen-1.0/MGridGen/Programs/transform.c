/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * transform.c
 *
 * This file is the driver for Transform
 *
 * Started 4/14/99
 * George
 *
 */

/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  char filename[256];

  if (argc != 2) {
    printf("Usage: %s <GraphFile>\n", argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);

  TransformGraph(filename);

  printf("Parameters: %s\n", filename);

  return(0);
}  
