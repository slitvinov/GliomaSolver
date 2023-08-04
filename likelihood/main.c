#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double PETsigma2, PETscale, slope, T1uc, T2uc;
static size_t mNelements;
#define PI (3.141592653589793)

static int sgn(double d) {
  double eps = 0.0;
  if (d < -eps) {
    return -1;
  } else {
    return d > eps;
  }
}

float *D3D(const char *path) {
  FILE *file;
  int header[6];
  float *mData;
  file = fopen(path, "r");
  fread(header, sizeof header, 1, file);
  assert(header[0] == 1234);
  assert(header[1] == 3);
  assert(header[5] == 1);
  mNelements = header[2] * header[3] * header[4];
  mData = (float *)malloc(mNelements * sizeof *mData);
  fread(mData, mNelements, sizeof *mData, file);
  fclose(file);
  return mData;
}

long double LogBernoulli(double u, double y, int Ti) {
  double uc;
  if (Ti == 1) {
    uc = T1uc;
  } else {
    uc = T2uc;
  }
  double diff = u - uc;
  long double omega2 = (diff > 0.) ? 1. : diff * diff;
  long double alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp(-omega2 / slope));
  return (y == 1) ? log(alpha) : log(1. - alpha);
}
long double TiLogLikelihood(float *model, int Ti) {
  float *data;
  if (Ti == 1)
    data = D3D("tumT1c.dat");
  else
    data = D3D("tumFLAIR.dat");
  long double sum = 0.;
  for (int i = 0; i < mNelements; i++)
    sum += LogBernoulli(model[i], data[i], Ti);
  return sum;
}
int main(int argc, const char **argv) {
  float *model, *PETdata;
  int N, i;
  long double sum;

  PETsigma2 = 0.000361;
  PETscale = 0.85;
  T1uc = 0.7;
  T2uc = 0.25;
  slope = 2;
  model = D3D("HGG_data.dat");
  PETdata = D3D("tumPET.dat");
  N = 0;
  sum = 0.;
  for (i = 0; i < mNelements; i++)
    if (PETdata[i] > 0.) {
      sum += (model[i] - PETscale * PETdata[i]) *
             (model[i] - PETscale * PETdata[i]);
      N++;
    }
  long double p1 = -0.5 * N * log(2. * PI * PETsigma2);
  long double p2 = -0.5 * (1. / PETsigma2) * sum;
  long double Lpet = p1 + p2;
  long double Lt1 = TiLogLikelihood(model, 1);
  long double Lt2 = TiLogLikelihood(model, 2);
  printf("%.16Le\n", -(Lpet + Lt1 + Lt2));
};
