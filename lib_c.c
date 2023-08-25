#include "lib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI (3.141592653589793)

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

static int sgn(double d) {
  double eps = 0.0;
  if (d < -eps) {
    return -1;
  } else {
    return d > eps;
  }
}

static long double TiLogLikelihood(int m, float *model, float *data,
                                   double slope, double uc) {
  int i;
  long double sum, diff, omega2, alpha;
  sum = 0.0;
  for (i = 0; i < m; i++) {
    diff = model[i] - uc;
    omega2 = (diff > 0.) ? 1. : diff * diff;
    alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp(-omega2 / slope));
    sum += data[i] == 1 ? log(alpha) : log(1. - alpha);
  }
  return sum;
}

int brain_likelihood(struct LikelihoodParams *params, float *model,
                     double *ans) {
  int N, m, i;
  long double sum, p1, p2, Lt1, Lt2;
  m = params->n[0] * params->n[1] * params->n[2];
  N = 0;
  sum = 0.;
  for (i = 0; i < m; i++)
    if (params->PET[i] > 0.) {
      sum += (model[i] - params->PETscale * params->PET[i]) *
             (model[i] - params->PETscale * params->PET[i]);
      N++;
    }
  p1 = -0.5 * N * log(2. * PI * params->PETsigma2);
  p2 = -0.5 * (1. / params->PETsigma2) * sum;
  Lt1 = TiLogLikelihood(m, model, params->T1c, params->slope, params->T1uc);
  Lt2 = TiLogLikelihood(m, model, params->FLAIR, params->slope, params->T2uc);
  *ans = -(p1 + p2 + Lt1 + Lt2);
  return 0;
}
