#include "lib.h"
#include <stdio.h>
#include <stdlib.h>

struct Brain {};

int brain_ini(int nx, int ny, int nz, float *GM, float *WM, double Dw,
              double rho, struct Brain **pbrain) {
  struct Brain *brain;
  if ((brain = (struct Brain *)malloc(sizeof *brain)) == NULL) {
    fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
    return 1;
  }
  *pbrain = brain;
  return 0;
}

int brain_fin(struct Brain *brain) {
  free(brain);
  return 0;
}

float *brain_read(const char *path, int *nx, int *ny, int *nz) {
  float *d;
  FILE *file;
  int32_t header[6];
  int n;
  if ((file = fopen(path, "r")) == NULL) {
    fprintf(stderr, "%s:%d: error: fail to open '%s'\n", __FILE__, __LINE__,
            path);
    return NULL;
  }
  if (fread(header, sizeof header, 1, file) != 1) {
    fprintf(stderr, "%s:%d: error: fail to read '%s'\n", __FILE__, __LINE__,
            path);
    return NULL;
  }
  if (header[0] != 1234) {
    fprintf(stderr, "%s:%d: not a data file '%s'\n", __FILE__, __LINE__, path);
    return NULL;
  }
  if (header[1] != 3) {
    fprintf(stderr, "%s:%d: not a three dimensional file '%s'\n", __FILE__,
            __LINE__, path);
    return NULL;
  }
  if (header[5] != 1) {
    fprintf(stderr, "%s:%d: not a float32 file '%s'\n", __FILE__, __LINE__,
            path);
    return NULL;
  }
  n = header[2] * header[3] * header[4];
  *nx = header[2];
  *ny = header[3];
  *nz = header[4];
  if ((d = (float *)malloc(n * sizeof *d)) == NULL) {
    fprintf(stderr, "%s:%d: malloc failed for '%s'\n", __FILE__, __LINE__,
            path);
    return NULL;
  }
  if (fread(d, sizeof *d, n, file) != n) {
    fprintf(stderr, "%s:%d: error: fail to read '%s'\n", __FILE__, __LINE__,
            path);
    return NULL;
  }
  if (fclose(file) != 0) {
    fprintf(stderr, "%s:%d: error: fail to close '%s'\n", __FILE__, __LINE__,
            path);
    return NULL;
  }
  return d;
}