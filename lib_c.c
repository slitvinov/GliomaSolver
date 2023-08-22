#include <stdio.h>
#include <stdlib.h>
#include "lib.h"

float *brain_read(const char *path, int *n) {
  float *d;
  FILE *file;
  int32_t header[6];
  size_t m;
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
  n[0] = header[2];
  n[1] = header[3];
  n[2] = header[4];
  m = n[0] * n[1] * n[2];
  if ((d = (float *)malloc(m * sizeof *d)) == NULL) {
    fprintf(stderr, "%s:%d: malloc failed for '%s'\n", __FILE__, __LINE__,
	    path);
    return NULL;
  }
  if (fread(d, sizeof *d, m, file) != m) {
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
