#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <string>
#include <vector>

static double PETsigma2, PETscale, slope, T1uc, T2uc;
static size_t mNx, mNy, mNz, mNelements;

static int sgn(double d) {
  double eps = 0.0;
  if (d < -eps) {
    return -1;
  } else {
    return d > eps;
  }
}

using namespace std;
class D3D {
public:
  D3D(const char *filename) {
    std::ifstream fin(filename, std::ios::binary);
    int size[3];
  int header[2];
  fin.read((char *)header, 2 * sizeof(int));
  fin.read((char *)size, 3 * sizeof(int));
    mNx = size[0];
    mNy = size[1];
    mNz = size[2];
    mNelements = mNx * mNy * mNz;
    mData = new float[mNelements];

  int data_type;
  fin.read((char *)&data_type, sizeof(int));
  fin.read((char *)mData, sizeof(float) * mNelements);
    fin.close();
  }
  ~D3D() { delete mData; }

public:
  float operator()(size_t i, size_t j, size_t k) const {
    assert(i < mNx && j < mNy && k < mNz);
    return mData[i + (j + k * mNy) * mNx];
  }
  float operator()(size_t i) const {
    return mData[i];
  }
private:
  float *mData;
};
long double PETLogLikelihood(D3D &model) {
  D3D PETdata("tumPET.dat");
  int N = 0;
  long double sum = 0.;
  for (int iz = 0; iz < mNz; iz++)
    for (int iy = 0; iy < mNy; iy++)
      for (int ix = 0; ix < mNx; ix++) {
        if (PETdata(ix, iy, iz) > 0.) {
          sum += (model(ix, iy, iz) - PETscale * PETdata(ix, iy, iz)) *
                 (model(ix, iy, iz) - PETscale * PETdata(ix, iy, iz));
          N++;
        }
      }
  long double p1 = -0.5 * N * log(2. * M_PI * PETsigma2);
  long double p2 = -0.5 * (1. / PETsigma2) * sum;
  return p1 + p2;
}
long double LogBernoulli(double u, double y, int Ti) {
  double uc, is2;
  if (Ti == 1) {
    uc = T1uc;
    is2 = 1. / slope;
  } else {
    uc = T2uc;
    is2 = 1. / slope;
  }
  double diff = u - uc;
  long double omega2 = (diff > 0.) ? 1. : diff * diff;
  long double alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp(-omega2 * is2));
  return (y == 1) ? log(alpha) : log(1. - alpha);
}
long double TiLogLikelihood(D3D &model, int Ti) {
  char filename[256];
  if (Ti == 1)
    sprintf(filename, "tumT1c.dat");
  else
    sprintf(filename, "tumFLAIR.dat");
  D3D data(filename);
  long double sum = 0.;
  for (int i = 0; i < mNelements; i++)
        sum += LogBernoulli(model(i), data(i), Ti);
  return sum;
}
int main(int argc, const char **argv) {
  ifstream mydata("LikelihoodInput.txt");
  if (mydata.is_open()) {
    mydata >> PETsigma2;
    mydata >> PETscale;
    mydata >> T1uc;
    mydata >> T2uc;
    mydata >> slope;
    mydata.close();
  } else {
    printf("Aborting: missing input file LikelihoodInput.txt \n");
    abort();
  }
  D3D model("HGG_data.dat");
  long double Lpet = PETLogLikelihood(model);
  long double Lt1 = TiLogLikelihood(model, 1);
  long double Lt2 = TiLogLikelihood(model, 2);
  printf("%.16Le\n", -(Lpet + Lt1 + Lt2));
};
