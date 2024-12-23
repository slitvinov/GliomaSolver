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

long double TiLogLikelihood(float *model, float *data, double uc) {
  int i;
  long double sum, diff, omega2, alpha;
  sum = 0.0;
  for (i = 0; i < mNelements; i++)  {
    diff = model[i] - uc;
    omega2 = (diff > 0.) ? 1. : diff * diff;
    alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp(-omega2 / slope));
    sum += data[i] == 1 ? log(alpha) : log(1. - alpha);
  }
  return sum;
}
int main(int argc, const char **argv) {
  float *model, *PETdata, *tumT1c, *tumFLAIR;
  int N, i;
  long double sum, p1, p2, Lt1, Lt2;
  double ans;

  PETsigma2 = 0.000361;
  PETscale = 0.85;
  T1uc = 0.7;
  T2uc = 0.25;
  slope = 2;
  model = D3D("HGG_data.dat");
  PETdata = D3D("tumPET.dat");
  tumT1c = D3D("tumT1c.dat");
  tumFLAIR = D3D("tumFLAIR.dat");
  N = 0;
  sum = 0.;
  for (i = 0; i < mNelements; i++)
    if (PETdata[i] > 0.) {
      sum += (model[i] - PETscale * PETdata[i]) *
             (model[i] - PETscale * PETdata[i]);
      N++;
    }
  p1 = -0.5 * N * log(2. * PI * PETsigma2);
  p2 = -0.5 * (1. / PETsigma2) * sum;
  Lt1 = TiLogLikelihood(model, tumT1c, T1uc);
  Lt2 = TiLogLikelihood(model, tumFLAIR, T2uc);
  ans = -(p1 + p2 + Lt1 + Lt2);
  printf("%.16le\n", ans);
};
