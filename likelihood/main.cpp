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
#include <sys/time.h>
#include <vector>

using namespace std;

class Value {
private:
  string content;

public:
  Value() : content("") {}

  Value(string content_) : content(content_) {}

  double asDouble(double def = 0) const {
    if (content == "")
      return def;
    return (double)atof(content.c_str());
  }

  int asInt(int def = 0) const {
    if (content == "")
      return def;
    return atoi(content.c_str());
  }

  bool asBool(bool def = false) const {
    if (content == "")
      return def;
    if (content == "0")
      return false;
    if (content == "false")
      return false;

    return true;
  }

  string asString(string def = "") const {
    if (content == "")
      return def;

    return content;
  }
};

class ArgumentParser {
private:
  map<string, Value> mapArguments;

  const int iArgC;
  const char **vArgV;

public:
  Value operator()(const string arg) {
    printf("%s is %s\n", arg.data(), mapArguments[arg].asString().data());
    return mapArguments[arg];
  }

  ArgumentParser(const int argc, const char **argv)
      : mapArguments(), iArgC(argc), vArgV(argv) {
    for (int i = 1; i < argc; i++)
      if (argv[i][0] == '-') {
        string values = "";
        int itemCount = 0;

        for (int j = i + 1; j < argc; j++)
          if (argv[j][0] == '-')
            break;
          else {
            if (strcmp(values.c_str(), ""))
              values += ' ';

            values += argv[j];
            itemCount++;
          }

        if (itemCount == 0)
          values += '1';
        mapArguments[argv[i]] = Value(values);
        i += itemCount;
      }
  }
};

/**\file Matrix.h
 * Collection of matrices (models of #CMatrix2D and #CMatrix3D).
 * Most of them also feature serialization and deserialization.
 *
 * #include "Matrix.h"
 *
 * @author Gerardo Tauriello
 * @date 7/19/10
 *
 * @see Matrix, CMatrix2D, CMatrix3D
 */


// system includes
#include <fstream>

/**
 * Matrices all put into namespace Matrix (also includes some helpers).
 *
 * Naming: XYD
 * - X: D = dynamic, S = static allocation
 * - Y: 2 or 3 dimensional matrix
 * Usually templatized over type.
 *
 * Dump layout:
 * - Header (see Matrix::deserializeHeader())
 *   - Magic number 1234 (as int) -> to check for Little vs BigEndian
 *   - DIM=2 or 3 (int)
 *   - nx,ny,(nz) (int)
 *   - data
 * - Generic way to store data (see Matrix::serialize(), Matrix::deserialize())
 *   - type-id (see Matrix::TypeID)
 *   - data[..] -> 1D array
 */
namespace Matrix {
// note: all global stuff must be inline or we get duplicate symbols!


/** Enumerator for supported data types. */
enum TypeID { TID_DOUBLE = 0, TID_FLOAT = 1, TID_INT = 2, TID_INVALID = -1 };

/**
 * Type id for data.
 * @return ID for given template-type or -1 if not supported.
 */
template <typename T> inline int typeId() { return TID_INVALID; }
template <> inline int typeId<double>() { return TID_DOUBLE; }
template <> inline int typeId<float>() { return TID_FLOAT; }
template <> inline int typeId<int>() { return TID_INT; }

/**
 * Read header from data stream.
 * @param is	Input stream from which to load data (open as binary).
 * @param dim   Number of dimensions (2 or 3).
 * @param size  Size of matrix (given array with dim elements, filled here).
 * @return Input stream after reading the data (check success with "if (is)
 * ...").
 */
inline std::istream &deserializeHeader(std::istream &is, size_t dim,
                                       int *size) {
  // HEADER
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
  // MATRIX SIZE
  is.read((char *)size, dim * sizeof(int));
  return is;
}

/** Read header from file (see deserializeHeader() for details). */
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

/**
 * Get and convert data from a stream.
 * Data is stored as n_elem elements of type T2 and should be stored in
 * data array which must have space for n_elem elements of type T.
 * Do not use this for T == T2 (use <code>is.read((char*)data,
 * sizeof(T)*n_elem)</code> instead).
 *
 * @param is		Input stream from which to load data (open as binary).
 * @param data		Pointer to where to write read data (please make sure to
 * have enough space!)
 * @param n_elem	Number of elements to be read.
 * @tparam T2		Source type for data (i.e. type in stream)
 * @tparam T		Target type for data (may differ from type in stream)
 * @return Input stream after reading the data (check success with "if (is)
 * ...").
 */
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

/**
 * Get data from a stream.
 * Data is stored as:
 * - type_id (int)
 * - n_elem elements of type given by type_id
 * The size of the data array must be n_elem * sizeof(T).
 *
 * @param is		Input stream from which to load data (open as binary).
 * @param data		Pointer to where to write read data (please make sure to
 * have enough space!)
 * @param n_elem	Number of elements to be read.
 * @tparam T		Target type for data (may differ from type in file)
 * @return Input stream after reading the data (check success with "if (is)
 * ...").
 */
template <typename T>
inline std::istream &deserialize(std::istream &is, T *data, int n_elem) {
  // read type of data
  int data_type;
  is.read((char *)&data_type, sizeof(int));
  if (!is || TypeID(data_type) == TID_INVALID) {
    // FAIL
    is.clear(std::ios::badbit);
    return is;
  }
  // compare with given type
  if (data_type == typeId<T>()) {
    // fast version
    return is.read((char *)data, sizeof(T) * n_elem);
  } else {
    // must convert
    switch (TypeID(data_type)) {
    case TID_DOUBLE:
      return deserializeConvert<double>(is, data, n_elem);
    case TID_FLOAT:
      return deserializeConvert<float>(is, data, n_elem);
    case TID_INT:
      return deserializeConvert<int>(is, data, n_elem);
    default:
      // FAIL
      is.clear(std::ios::badbit);
      return is;
    }
  }
}

/**
 * Get and convert data and write to stream.
 * Data is stored as n_elem elements of type T2 and should came from
 * data array which must contain n_elem elements of type T.
 * Do not use this for T == T2 (use <code>os.write((char*)data,
 * sizeof(T)*n_elem)</code> instead).
 *
 * @param os		Output stream to which to write data (open as binary).
 * @param data		Pointer to array with n_elem elements.
 * @param n_elem	Number of elements to be written.
 * @tparam T2		Target type for data (i.e. type in stream)
 * @tparam T		Source type for data (may differ from type in stream)
 * @return Output stream after writing the data (check success with "if (os)
 * ...").
 */
template <typename T2, typename T>
inline std::ostream &serializeConvert(std::ostream &os, const T *data,
                                      int n_elem) {
  T2 *tmp = new T2[n_elem];
  for (int i = 0; i < n_elem; ++i) {
    tmp[i] = (T2)data[i];
  }
  os.write((char *)tmp, sizeof(T2) * n_elem);
  delete tmp;
  return os;
}

/**
 * Write data to a stream.
 * Data is stored as:
 * - type_id (int)
 * - n_elem elements of type given by type_id
 * The size of the data array must be n_elem * sizeof(T).
 *
 * @param os		Output stream to which to write data (open as binary).
 * @param data		Pointer to array with n_elem elements.
 * @param n_elem	Number of elements to be written.
 * @param tid		(optional) Target type for data (default: type id for T)
 * @tparam T		Source type for data (may differ from type in stream)
 * @return Output stream after writing the data (check success with "if (os)
 * ...").
 */
template <typename T>
inline std::ostream &serialize(std::ostream &os, const T *data, int n_elem,
                               TypeID tid = (TypeID)typeId<T>()) {
  // type id check
  if (tid == TID_INVALID) {
    // FAIL
    os.clear(std::ios::badbit);
    return os;
  }

  // write type of data
  int data_type = tid;
  os.write((char *)&data_type, sizeof(int));

  // compare with given type
  if (tid == typeId<T>()) {
    // fast version
    return os.write((char *)data, sizeof(T) * n_elem);
  } else {
    // must convert
    switch (TypeID(data_type)) {
    case TID_DOUBLE:
      return serializeConvert<double>(os, data, n_elem);
    case TID_FLOAT:
      return serializeConvert<float>(os, data, n_elem);
    case TID_INT:
      return serializeConvert<int>(os, data, n_elem);
    default:
      // FAIL
      os.clear(std::ios::badbit);
      return os;
    }
  }
}

/**
 * 2D matrix with dynamic allocation.
 * Loops over elements should always go first in y and then in x.
 *
 * - Associated types:
 *   - <code>ElementType</code> - Type of elements
 * - Methods:
 *   - <code>size_t getSizeX() const</code> - Get number of elements in x
 *   - <code>size_t getSizeY() const</code> - Get number of elements in y
 *   - <code>ElementType operator()(size_t i, size_t j) const</code>\n
 *     Read only access to element at x = i (in [0, getSizeX()-1]), y = j (in
[0, getSizeY()-1])
 *   - <code>ElementType& operator()(size_t i, size_t j)</code>\n
 *     Access element at x = i (in [0, getSizeX()-1]), y = j (in [0,
getSizeY()-1])
 * .
 * Example using variable m of type X:
\code
for (int j = 0; j < m.getSizeY(); ++j) {
    for (int i = 0; i < m.getSizeX(); ++i) {
        m(i,j) = 0;
    }
}
\endcode
 */
template <typename T> class D2D {
  // TYPEDEFS
public:
  typedef T ElementType;

public:
  /**
   * Default constructor creating nx x ny matrix.
   */
  D2D(const size_t nx, const size_t ny) { init(nx, ny); }

  /**
   * Constructor filling matrix with data.
   * Data must have nx*ny elements.
   */
  D2D(const size_t nx, const size_t ny, T *data) {
    init(nx, ny);
    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = data[i];
    }
  }

  /**
   * Copy constructor.
   *
   * @param from The value to copy to this object.
   */
  D2D(const D2D &from) {
    init(from.mNx, from.mNy);
    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = from.mData[i];
    }
  }

  /**
   * Construct from file.
   * See #Matrix for dump layout.
   */
  D2D(const char *filename) : mNx(0), mData(NULL) { load(filename); }
  /**
   * Construct from binary stream.
   * See #Matrix for dump layout.
   */
  D2D(std::istream &is) : mNx(0), mData(NULL) { load(is); }

  /**
   * Destructor.
   */
  ~D2D() { delete mData; }

private:
  /** Initialize matrix (allocate and set stuff). */
  void init(const size_t nx, const size_t ny) {
    // use this only for primitive types (rest should be in init. list)
    mNx = nx;
    mNy = ny;
    mNelements = nx * ny;
    mData = new T[mNelements];
  }

public:
  /**
   * Assignment operator.
   *
   * @param from The value to assign to this object.
   *
   * @return A reference to this object.
   */
  D2D &operator=(const D2D &from) {
    // handle self-assignment
    if (this == &from)
      return *this;

    // resize needed?
    if (from.mNx != mNx || from.mNy != mNy) {
      delete mData;
      init(from.mNx, from.mNy);
      assert(mNelements == from.mNelements);
    }

    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = from.mData[i];
    }

    // don't forget to return
    return *this;
  }

public:
  /** Get size of matrix (X). */
  size_t getSizeX() const { return mNx; }
  /** Get size of matrix (Y). */
  size_t getSizeY() const { return mNy; }
  /** Access element (read-only). */
  T operator()(size_t i, size_t j) const {
    assert(i < mNx && j < mNy);
    return mData[i + j * mNx];
  }
  /** Access element. */
  T &operator()(size_t i, size_t j) {
    assert(i < mNx && j < mNy);
    return mData[i + j * mNx];
  }

public:
  /**
   * Dump matrix to file.
   * See #Matrix for dump layout.
   */
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

  /**
   * Dump matrix to binary stream.
   * See #Matrix for dump layout.
   */
  std::ostream &dump(std::ostream &os, TypeID tid = (TypeID)typeId<T>()) const {
    // HEADER
    int header[4] = {1234, 2, mNx, mNy};
    os.write((char *)header, 4 * sizeof(int));
    // DATA
    return serialize(os, mData, mNelements, tid);
  }

  /**
   * Load matrix from file.
   * See #Matrix for dump layout.
   */
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

  /**
   * Load matrix from binary stream.
   * See #Matrix for dump layout.
   */
  std::istream &load(std::istream &is) {
    // HEADER
    int size[2];
    deserializeHeader(is, 2, size);
    // resize needed?
    if (size[0] != mNx || size[1] != mNy) {
      delete mData;
      init(size[0], size[1]);
    }

    // LOAD DATA
    return deserialize(is, mData, mNelements);
  }

private:
  // mData is mNx x mNy array with mNelements elements
  // mNx = -1 && mData = NULL is for non-initialized matrix
  size_t mNx;
  size_t mNy;
  size_t mNelements;
  T *mData;
};

/**
 * 3D matrix with dynamic allocation.
 * Loops over elements should always go first in z, then in y and then in x.
 *
 * - Associated types:
 *   - <code>ElementType</code> - Type of elements
 * - Methods:
 *   - <code>size_t getSizeX() const</code> - Get number of elements in x
 *   - <code>size_t getSizeY() const</code> - Get number of elements in y
 *   - <code>size_t getSizeZ() const</code> - Get number of elements in z
 *   - <code>ElementType operator()(size_t i, size_t j, size_t k) const</code>\n
 *     Read only access to element at x = i (in [0, getSizeX()-1]), y = j (in
[0, getSizeY()-1]), z = k (in [0, getSizeZ()-1])
 *   - <code>ElementType& operator()(size_t i, size_t j, size_t k)</code>\n
 *     Access element at x = i (in [0, getSizeX()-1]), y = j (in [0,
getSizeY()-1]), z = k (in [0, getSizeZ()-1])
 * .
 * Example using variable m of type X:
\code
for (int k = 0; k < m.getSizeZ(); ++k) {
    for (int j = 0; j < m.getSizeY(); ++j) {
        for (int i = 0; i < m.getSizeX(); ++i) {
            m(i,j,k) = 0;
        }
    }
}
\endcode
 */
template <typename T> class D3D {
  // TYPEDEFS
public:
  typedef T ElementType;

public:
  /**
   * Default constructor creating nx x ny x nz matrix.
   */
  D3D(const size_t nx, const size_t ny, const size_t nz) { init(nx, ny, nz); }

  /**
   * Constructor filling matrix with data.
   * Data must have nx*ny*nz elements.
   */
  D3D(const size_t nx, const size_t ny, const size_t nz, T *data) {
    init(nx, ny, nz);
    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = data[i];
    }
  }

  /**
   * Copy constructor.
   *
   * @param from The value to copy to this object.
   */
  D3D(const D3D &from) {
    init(from.mNx, from.mNy, from.mNz);
    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = from.mData[i];
    }
  }

  /**
   * Construct from file.
   * See #Matrix for dump layout.
   */
  D3D(const char *filename) : mNx(0), mData(NULL) { load(filename); }
  /**
   * Construct from binary stream.
   * See #Matrix for dump layout.
   */
  D3D(std::istream &is) : mNx(0), mData(NULL) { load(is); }

  /**
   * Destructor.
   */
  ~D3D() { delete mData; }

private:
  /** Initialize matrix (allocate and set stuff). */
  void init(const size_t nx, const size_t ny, const size_t nz) {
    // use this only for primitive types (rest should be in init. list)
    mNx = nx;
    mNy = ny;
    mNz = nz;
    mNelements = nx * ny * nz;
    mData = new T[mNelements];
  }

public:
  /**
   * Assignment operator.
   *
   * @param from The value to assign to this object.
   *
   * @return A reference to this object.
   */
  D3D &operator=(const D3D &from) {
    // handle self-assignment
    if (this == &from)
      return *this;

    // resize needed?
    if (from.mNx != mNx || from.mNy != mNy || from.mNz != mNz) {
      delete mData;
      init(from.mNx, from.mNy, from.mNz);
      assert(mNelements == from.mNelements);
    }

    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = from.mData[i];
    }

    // don't forget to return
    return *this;
  }

public:
  /** Get size of matrix (X). */
  size_t getSizeX() const { return mNx; }
  /** Get size of matrix (Y). */
  size_t getSizeY() const { return mNy; }
  /** Get size of matrix (Z). */
  size_t getSizeZ() const { return mNz; }
  /** Access element (read-only). */
  T operator()(size_t i, size_t j, size_t k) const {
    assert(i < mNx && j < mNy && k < mNz);
    return mData[i + (j + k * mNy) * mNx];
  }
  /** Access element. */
  T &operator()(size_t i, size_t j, size_t k) {
    assert(i < mNx && j < mNy && k < mNz);
    return mData[i + (j + k * mNy) * mNx];
  }

public:
  /**
   * Dump matrix to file.
   * See #Matrix for dump layout.
   */
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

  /**
   * Dump matrix to binary stream.
   * See #Matrix for dump layout.
   */
  std::ostream &dump(std::ostream &os, TypeID tid = (TypeID)typeId<T>()) const {
    // HEADER
    int header[5] = {1234, 3, mNx, mNy, mNz};
    os.write((char *)header, 5 * sizeof(int));
    // DATA
    return serialize(os, mData, mNelements, tid);
  }

  /**
   * Load matrix from file.
   * See #Matrix for dump layout.
   */
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

  /**
   * Load matrix from binary stream.
   * See #Matrix for dump layout.
   */
  std::istream &load(std::istream &is) {
    // HEADER
    int size[3];
    deserializeHeader(is, 3, size);
    // resize needed?
    if (size[0] != mNx || size[1] != mNy || size[2] != mNz) {
      delete mData;
      init(size[0], size[1], size[2]);
    }

    // LOAD DATA
    return deserialize(is, mData, mNelements);
  }

private:
  // mData is mNx x mNy x mNz array with mNelements elements
  // mNx = -1 && mData = NULL is for non-initialized matrix
  size_t mNx;
  size_t mNy;
  size_t mNz;
  size_t mNelements;
  T *mData;
};

} // namespace Matrix

typedef Matrix::D3D<double> MatrixD3D;
typedef Matrix::D2D<double> MatrixD2D;

class HGG_Likelihood {
private:
  ArgumentParser parser;
  long double _computePETLogLikelihood(MatrixD3D model);
  long double _computeTiLogLikelihood(MatrixD3D model, int Ti);
  long double _computeLogBernoulli(double y, double u, int Ti);
  void _writeToFile(long double output);
  int sgn(double d);

  double PETsigma2, PETscale;
  double slope; // same for T1 & T2(sigma2 double sigmoid, k singel sigmoid)
  double T1uc, T2uc;
  bool bSelectivePoints;
  int stepPET, stepTi;

public:
  HGG_Likelihood(const int argc, const char **argv);
  ~HGG_Likelihood(){};
  void run();
};


HGG_Likelihood::HGG_Likelihood(const int argc, const char **argv)
    : parser(argc, argv) {
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

  bSelectivePoints = parser("-bROI").asBool(0);
  stepPET = parser("-stepPET").asInt(1);
  stepTi = parser("-stepMRI").asInt(1);

  printf("PET: PETsigma2=%f, PETscale=%f \n", PETsigma2, PETscale);
  printf("MRI: T1uc=%f, T2uc =%f, slope=%f \n", T1uc, T2uc, slope);
  printf("stepPET=%d, stepTi =%d \n", stepPET, stepTi);
}

long double HGG_Likelihood::_computePETLogLikelihood(MatrixD3D model) {
  /* 1) Read in PET data */
  char filename[256];
  sprintf(filename, "tumPET.dat");
  MatrixD3D PETdata(filename); // pet signal, u_real = PETscale*pet_signal

  int dataX = PETdata.getSizeX();
  int dataY = PETdata.getSizeY();
  int dataZ = PETdata.getSizeZ();

  int N = 0;
  long double sum = 0.;

  for (int iz = 0; iz < dataZ; iz++)
    for (int iy = 0; iy < dataY; iy++)
      for (int ix = 0; ix < dataX; ix++) {
        if ((ix % stepPET == 0) && (iy % stepPET == 0) && (iz % stepPET == 0)) {
          if (PETdata(ix, iy, iz) > 0.) {
            sum += (model(ix, iy, iz) - PETscale * PETdata(ix, iy, iz)) *
                   (model(ix, iy, iz) - PETscale * PETdata(ix, iy, iz));
            N++;
          }
        }
      }

  long double p1 = -0.5 * N * log(2. * M_PI * PETsigma2);
  long double p2 = -0.5 * (1. / PETsigma2) * sum;

  printf("Points=%i, stepPET=%i \n", N, stepPET);
  return (p1 + p2);
}

long double HGG_Likelihood::_computeTiLogLikelihood(MatrixD3D model, int Ti) {
  /* 1) Read in T1 data */
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

  if (!bSelectivePoints) {
    //#pragma omp parallel for reduction(+:sum)
    for (int iz = 0; iz < dataZ; iz++)
      for (int iy = 0; iy < dataY; iy++)
        for (int ix = 0; ix < dataX; ix++)
          sum += _computeLogBernoulli(model(ix, iy, iz), data(ix, iy, iz), Ti);
  } else {
    sprintf(filename, "ROI.dat");
    MatrixD2D Points(filename);
    int Npoints = Points.getSizeX();

    for (int i = 0; i < Npoints; i++) {
      int ix = Points(i, 0);
      int iy = Points(i, 1);
      int iz = Points(i, 2);

      if ((ix % stepTi == 0) && (iy % stepTi == 0) && (iz % stepTi == 0))
        sum += _computeLogBernoulli(model(ix, iy, iz), data(ix, iy, iz), Ti);
    }

    printf("Ti Npoints = %i \n", Npoints);
  }

  printf("LogLike of T%i = %Lf \n", Ti, sum);
  return sum;
}

/* Likelihood based on Bernoulli distr. of double logsitic sigmoid
 1) compute alpha
 - should be in (0,1)
 - rounding errros can make alpha = 0 or =1
 - if that is the case, correct it since it will be used in log() */
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

  // OLD version: long double alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp( -diff *
  // diff * is2));
  return (y == 1) ? log(alpha) : log(1. - alpha);
}

int HGG_Likelihood::sgn(double d) {
  double eps = 0.0;
  if (d < -eps) {
    return -1;
  } else {
    return d > eps;
  }
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

double my_gettime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1.0E-6;
}

int main(int argc, const char **argv) {
  printf("==========================================\n");
  printf("          Computing likelihood            \n");
  printf("==========================================\n");

  ArgumentParser parser(argc, argv);
  double t_begin = my_gettime();

  HGG_Likelihood *l = new HGG_Likelihood(argc, (const char **)argv);
  l->run();
  delete l;

  double t_end = my_gettime();
  printf("Likelihood: elapsed time=%f (sec) \n", t_end - t_begin);

  return 0;
};
