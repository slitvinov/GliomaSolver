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

static  double PETsigma2, PETscale, slope, T1uc, T2uc;
static int sgn(double d) {
  double eps = 0.0;
  if (d < -eps) {
    return -1;
  } else {
    return d > eps;
  }
}

using namespace std;
namespace Matrix {
enum TypeID { TID_DOUBLE = 0, TID_FLOAT = 1, TID_INVALID = -1 };
template <typename T> inline int typeId() { return TID_INVALID; }
template <> inline int typeId<double>() { return TID_DOUBLE; }
template <> inline int typeId<float>() { return TID_FLOAT; }
inline std::istream &deserializeHeader(std::istream &is, size_t dim,
                                       int *size) {
  int header[2];
  is.read((char *)header, 2 * sizeof(int));
  if (header[0] != 1234) {
    std::cout << "Magic number did not match. Aborting\n";
    is.clear(std::ios::badbit);
    return is;
  }
  if (header[1] != dim) {
    std::cout << "Dimensions of matrix do not match. Aborting\n";
    is.clear(std::ios::badbit);
    return is;
  }
  is.read((char *)size, dim * sizeof(int));
  return is;
}
inline void deserializeHeader(const char *filename, size_t dim, int *size) {
  std::ifstream fin(filename, std::ios::binary);
  if (!fin.is_open()) {
    std::cout << "ERROR while opening " << filename << std::endl;
    return;
  }
  if (!deserializeHeader(fin, dim, size)) {
    std::cout << "ERROR while reading " << filename << std::endl;
  }
  fin.close();
}
template <typename T2, typename T>
inline std::istream &deserializeConvert(std::istream &is, T *data, int n_elem) {
  T2 *tmp = new T2[n_elem];
  is.read((char *)tmp, sizeof(T2) * n_elem);
  for (int i = 0; i < n_elem; ++i) {
    data[i] = (T)tmp[i];
  }
  delete tmp;
  return is;
}
template <typename T>
inline std::istream &deserialize(std::istream &is, T *data, int n_elem) {
  int data_type;
  is.read((char *)&data_type, sizeof(int));
  if (!is || TypeID(data_type) == TID_INVALID) {
    is.clear(std::ios::badbit);
    return is;
  }
  if (data_type == typeId<T>()) {
    return is.read((char *)data, sizeof(T) * n_elem);
  } else {
    switch (TypeID(data_type)) {
    case TID_DOUBLE:
      return deserializeConvert<double>(is, data, n_elem);
    case TID_FLOAT:
      return deserializeConvert<float>(is, data, n_elem);
    default:
      is.clear(std::ios::badbit);
      return is;
    }
  }
}
template <typename T> class D3D {
public:
  typedef T ElementType;

public:
  D3D(const size_t nx, const size_t ny, const size_t nz) { init(nx, ny, nz); }
  D3D(const size_t nx, const size_t ny, const size_t nz, T *data) {
    init(nx, ny, nz);
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = data[i];
    }
  }
  D3D(const D3D &from) {
    init(from.mNx, from.mNy, from.mNz);
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = from.mData[i];
    }
  }
  D3D(const char *filename) : mNx(0), mData(NULL) { load(filename); }
  D3D(std::istream &is) : mNx(0), mData(NULL) { load(is); }
  ~D3D() { delete mData; }

private:
  void init(const size_t nx, const size_t ny, const size_t nz) {
    mNx = nx;
    mNy = ny;
    mNz = nz;
    mNelements = nx * ny * nz;
    mData = new T[mNelements];
  }

public:
  D3D &operator=(const D3D &from) {
    if (this == &from)
      return *this;
    if (from.mNx != mNx || from.mNy != mNy || from.mNz != mNz) {
      delete mData;
      init(from.mNx, from.mNy, from.mNz);
      assert(mNelements == from.mNelements);
    }
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = from.mData[i];
    }
    return *this;
  }

public:
  size_t getSizeX() const { return mNx; }
  size_t getSizeY() const { return mNy; }
  size_t getSizeZ() const { return mNz; }
  T operator()(size_t i, size_t j, size_t k) const {
    assert(i < mNx && j < mNy && k < mNz);
    return mData[i + (j + k * mNy) * mNx];
  }
  T &operator()(size_t i, size_t j, size_t k) {
    assert(i < mNx && j < mNy && k < mNz);
    return mData[i + (j + k * mNy) * mNx];
  }

public:
  void dump(const char *filename, TypeID tid = (TypeID)typeId<T>()) const {
    std::ofstream fout(filename, std::ios::binary);
    if (!fout.is_open()) {
      std::cout << "ERROR while opening " << filename << std::endl;
      return;
    }
    if (!dump(fout, tid)) {
      std::cout << "ERROR while writing to " << filename << std::endl;
    }
    fout.close();
  }
  std::ostream &dump(std::ostream &os, TypeID tid = (TypeID)typeId<T>()) const {
    int header[5] = {1234, 3, mNx, mNy, mNz};
    os.write((char *)header, 5 * sizeof(int));
    return serialize(os, mData, mNelements, tid);
  }
  void load(const char *filename) {
    std::ifstream fin(filename, std::ios::binary);
    if (!fin.is_open()) {
      std::cout << "ERROR while opening " << filename << std::endl;
      return;
    }
    if (!load(fin)) {
      std::cout << "ERROR while reading " << filename << std::endl;
    }
    fin.close();
  }
  std::istream &load(std::istream &is) {
    int size[3];
    deserializeHeader(is, 3, size);
    if (size[0] != mNx || size[1] != mNy || size[2] != mNz) {
      delete mData;
      init(size[0], size[1], size[2]);
    }
    return deserialize(is, mData, mNelements);
  }

private:
  size_t mNx;
  size_t mNy;
  size_t mNz;
  size_t mNelements;
  T *mData;
};
} // namespace Matrix
typedef Matrix::D3D<double> MatrixD3D;
class HGG_Likelihood {
private:
  long double _computePETLogLikelihood(MatrixD3D model);
  long double _computeTiLogLikelihood(MatrixD3D model, int Ti);
  long double _computeLogBernoulli(double y, double u, int Ti);
  void _writeToFile(long double output);
public:
  HGG_Likelihood(const int argc, const char **argv);
  ~HGG_Likelihood(){};
  void run();
};
HGG_Likelihood::HGG_Likelihood(const int argc, const char **argv) {
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
  printf("PET: PETsigma2=%f, PETscale=%f \n", PETsigma2, PETscale);
  printf("MRI: T1uc=%f, T2uc =%f, slope=%f \n", T1uc, T2uc, slope);
}
long double HGG_Likelihood::_computePETLogLikelihood(MatrixD3D model) {
  char filename[256];
  sprintf(filename, "tumPET.dat");
  MatrixD3D PETdata(filename);
  int dataX = PETdata.getSizeX();
  int dataY = PETdata.getSizeY();
  int dataZ = PETdata.getSizeZ();
  int N = 0;
  long double sum = 0.;
  for (int iz = 0; iz < dataZ; iz++)
    for (int iy = 0; iy < dataY; iy++)
      for (int ix = 0; ix < dataX; ix++) {
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
long double HGG_Likelihood::_computeTiLogLikelihood(MatrixD3D model, int Ti) {
  char filename[256];
  if (Ti == 1)
    sprintf(filename, "tumT1c.dat");
  else
    sprintf(filename, "tumFLAIR.dat");
  MatrixD3D data(filename);
  int dataX = data.getSizeX();
  int dataY = data.getSizeY();
  int dataZ = data.getSizeZ();
  long int N = dataX * dataY * dataZ;
  assert(N == model.getSizeX() * model.getSizeY() * model.getSizeZ());
  long double sum = 0.;
  for (int iz = 0; iz < dataZ; iz++)
    for (int iy = 0; iy < dataY; iy++)
      for (int ix = 0; ix < dataX; ix++)
        sum += _computeLogBernoulli(model(ix, iy, iz), data(ix, iy, iz), Ti);
  printf("LogLike of T%i = %Lf \n", Ti, sum);
  return sum;
}
long double HGG_Likelihood::_computeLogBernoulli(double u, double y, int Ti) {
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
void HGG_Likelihood::_writeToFile(long double output) {
  long double MinusLogLikelihood = -output;
  FILE *myfile = fopen("Likelihood.txt", "w");
  if (myfile != NULL)
    fprintf(myfile, "%Lf \n", MinusLogLikelihood);
  fclose(myfile);
}
void HGG_Likelihood::run() {
  char filename[256];
  sprintf(filename, "HGG_data.dat");
  MatrixD3D model(filename);
  long double Lpet = _computePETLogLikelihood(model);
  long double Lt1 = _computeTiLogLikelihood(model, 1);
  long double Lt2 = _computeTiLogLikelihood(model, 2);
  long double costFunction = Lpet + Lt1 + Lt2;
  printf("L_Pet=%Lf, L_T1=%Lf, L_T2 = %Lf \n", Lpet, Lt1, Lt2);
  printf("LogLike = %Lf \n", costFunction);
  _writeToFile(costFunction);
}

int main(int argc, const char **argv) {
  HGG_Likelihood l(argc, (const char **)argv);
  l.run();
  return 0;
};
