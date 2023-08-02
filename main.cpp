#include "../MRAG/MRAGHeaders.h"
#include <fstream>
namespace Matrix {
enum TypeID { TID_DOUBLE = 0, TID_FLOAT = 1, TID_INT = 2, TID_INVALID = -1 };
template <typename T> inline int typeId() { return TID_INVALID; }
template <> inline int typeId<double>() { return TID_DOUBLE; }
template <> inline int typeId<float>() { return TID_FLOAT; }
template <> inline int typeId<int>() { return TID_INT; }
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
template <typename T> class D3D {
  // TYPEDEFS
public:
  typedef T ElementType;

public:
  D3D(const size_t nx, const size_t ny, const size_t nz) { init(nx, ny, nz); }
  D3D(const size_t nx, const size_t ny, const size_t nz, T *data) {
    init(nx, ny, nz);
    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = data[i];
    }
  }
  D3D(const D3D &from) {
    init(from.mNx, from.mNy, from.mNz);
    // copy from data
    for (int i = 0; i < mNelements; ++i) {
      mData[i] = from.mData[i];
    }
  }
  D3D(const char *filename) : mNx(0), mData(NULL) { load(filename); }
  D3D(std::istream &is) : mNx(0), mData(NULL) { load(is); }
  ~D3D() { delete mData; }

private:
  void init(const size_t nx, const size_t ny, const size_t nz) {
    // use this only for primitive types (rest should be in init. list)
    mNx = nx;
    mNy = ny;
    mNz = nz;
    mNelements = nx * ny * nz;
    mData = new T[mNelements];
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
    // HEADER
    int header[5] = {1234, 3, mNx, mNy, mNz};
    os.write((char *)header, 5 * sizeof(int));
    // DATA
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
/*
 Computes rhs of the reaction-diffusio equation for tumor:
 //  ------------------------------------
 //    ∂φ / ∂t = ∇( D∇φ) + ρ * φ(1-φ)
 //     + no flux BC
 //
 //     φ - tumor density
 //     ρ - tumor proliferation rate
 //     D - diffusivity
 //  ------------------------------------
 */
struct ReactionDiffusionOperator {
  int stencil_start[3];
  int stencil_end[3];

  const Real Dw, Dg, rho;

  ReactionDiffusionOperator(const Real Dw_, const Real Dg_, const Real rho_)
      : Dw(Dw_), Dg(Dg_), rho(rho_) {
    stencil_start[0] = stencil_start[1] = -1;
    stencil_end[0] = stencil_end[1] = +2;
    stencil_start[2] = -1;
    stencil_end[2] = +2;
  }

  ReactionDiffusionOperator(const ReactionDiffusionOperator &copy)
      : Dw(copy.Dw), Dg(copy.Dg), rho(copy.rho) {
    stencil_start[0] = stencil_start[1] = -1;
    stencil_end[0] = stencil_end[1] = +2;
    stencil_start[2] = -1;
    stencil_end[2] = +2;
  }

  template <typename LabType, typename BlockType>
  inline void operator()(LabType &lab, const BlockInfo &info,
                         BlockType &o) const {
    double h = info.h[0];
    double ih2 = 1. / (h * h);
    Real df[6];  // diffusion coefficient
    Real chf[6]; // domain charact. func, chf=0 -> outside, chf=1 inside domain:
                 // use to apply BC
    Real df_loc; // diffusion at the current point (local)
    Real chf_loc; // diffusion at the current point (local)

    for (int iz = 0; iz < BlockType::sizeZ; iz++)
      for (int iy = 0; iy < BlockType::sizeY; iy++)
        for (int ix = 0; ix < BlockType::sizeX; ix++) {
          // check if we are in the brain domain
          if (lab(ix, iy, iz).p_w + lab(ix, iy, iz).p_g + lab(ix, iy, iz).phi >
              0.) {
            df_loc = lab(ix, iy, iz).p_w * Dw + lab(ix, iy, iz).p_g * Dg;
            df[0] = lab(ix - 1, iy, iz).p_w * Dw + lab(ix - 1, iy, iz).p_g * Dg;
            df[1] = lab(ix + 1, iy, iz).p_w * Dw + lab(ix + 1, iy, iz).p_g * Dg;
            df[2] = lab(ix, iy - 1, iz).p_w * Dw + lab(ix, iy - 1, iz).p_g * Dg;
            df[3] = lab(ix, iy + 1, iz).p_w * Dw + lab(ix, iy + 1, iz).p_g * Dg;
            df[4] = lab(ix, iy, iz - 1).p_w * Dw + lab(ix, iy, iz - 1).p_g * Dg;
            df[5] = lab(ix, iy, iz + 1).p_w * Dw + lab(ix, iy, iz + 1).p_g * Dg;

            _harmonic_mean(df, df_loc);

            chf[0] = lab(ix - 1, iy, iz).phi + lab(ix - 1, iy, iz).p_w +
                     lab(ix - 1, iy, iz).p_g;
            chf[1] = lab(ix + 1, iy, iz).phi + lab(ix + 1, iy, iz).p_w +
                     lab(ix + 1, iy, iz).p_g;
            chf[2] = lab(ix, iy - 1, iz).phi + lab(ix, iy - 1, iz).p_w +
                     lab(ix, iy - 1, iz).p_g;
            chf[3] = lab(ix, iy + 1, iz).phi + lab(ix, iy + 1, iz).p_w +
                     lab(ix, iy + 1, iz).p_g;
            chf[4] = lab(ix, iy, iz - 1).phi + lab(ix, iy, iz - 1).p_w +
                     lab(ix, iy, iz - 1).p_g;
            chf[5] = lab(ix, iy, iz + 1).phi + lab(ix, iy, iz + 1).p_w +
                     lab(ix, iy, iz + 1).p_g;

            _applyNoFluxBC(df, chf);

            // diffusion fluxes
            double diffusionFluxIn = ih2 * (df[0] * lab(ix - 1, iy, iz).phi +
                                            df[1] * lab(ix + 1, iy, iz).phi +
                                            df[2] * lab(ix, iy - 1, iz).phi +
                                            df[3] * lab(ix, iy + 1, iz).phi +
                                            df[4] * lab(ix, iy, iz - 1).phi +
                                            df[5] * lab(ix, iy, iz + 1).phi);

            double diffusionFluxOut =
                -((df[0] + df[1] + df[2] + df[3] + df[4] + df[5]) *
                  lab(ix, iy, iz).phi * ih2);
            double reactionFlux =
                rho * lab(ix, iy, iz).phi * (1. - lab(ix, iy, iz).phi);

            o(ix, iy, iz).dphidt =
                diffusionFluxOut + diffusionFluxIn + reactionFlux;
          } else
            o(ix, iy, iz).dphidt = 0.;
        }
  }

  // Di,j = 2 * (Di * Dj / (Di + Dj)
  // set Di,j to zero if (Di + Dj = 0) i.e. no update and avoid division by zero
  inline void _harmonic_mean(Real (&df)[6], Real df_loc) const {
    Real eps = 1.0e-08; // to avoid divisin by zero

    df[0] =
        (df[0] + df_loc < eps) ? 0. : 2. * df[0] * df_loc / (df[0] + df_loc);
    df[1] =
        (df[1] + df_loc < eps) ? 0. : 2. * df[1] * df_loc / (df[1] + df_loc);
    df[2] =
        (df[2] + df_loc < eps) ? 0. : 2. * df[2] * df_loc / (df[2] + df_loc);
    df[3] =
        (df[3] + df_loc < eps) ? 0. : 2. * df[3] * df_loc / (df[3] + df_loc);

    df[4] =
        (df[4] + df_loc < eps) ? 0. : 2. * df[4] * df_loc / (df[4] + df_loc);
    df[5] =
        (df[5] + df_loc < eps) ? 0. : 2. * df[5] * df_loc / (df[5] + df_loc);
  }

  inline void _applyNoFluxBC(Real (&df)[6], Real n[6]) const {
    // n is domain char. func, use to apply bc by modifying the df term by the
    // ghost point
    Real eps = 0.1;

    if (n[0] < eps) {
      df[1] *= 2.0;
    }
    if (n[1] < eps) {
      df[0] *= 2.0;
    }
    if (n[2] < eps) {
      df[3] *= 2.0;
    }
    if (n[3] < eps) {
      df[2] *= 2.0;
    }

    if (n[4] < eps) {
      df[5] *= 2.0;
    }
    if (n[5] < eps) {
      df[4] *= 2.0;
    }
  }
};

struct UpdateTumor {
  double dt;

  UpdateTumor(double dt_) : dt(dt_) {}

  UpdateTumor(const UpdateTumor &copy) : dt(copy.dt) {}

  template <typename BlockType>
  inline void operator()(const BlockInfo &info, BlockType &o) const {
    for (int iz = 0; iz < BlockType::sizeZ; iz++)
      for (int iy = 0; iy < BlockType::sizeY; iy++)
        for (int ix = 0; ix < BlockType::sizeX; ix++) {
          o(ix, iy, iz).phi += dt * o(ix, iy, iz).dphidt;
          o(ix, iy, iz).phi = max((Real)0., o(ix, iy, iz).phi);
          o(ix, iy, iz).phi = min((Real)1., o(ix, iy, iz).phi);
        }
  }
};

using namespace MRAG;
struct Cell {
  /* tumor */
  Real phi;
  Real dphidt;

  /* tissue percentage per voxel*/
  Real p_g, p_w, p_csf;

  /* tissue concentration */
  Real wm, gm, csf; //
  Real dwmdt, dgmdt, dcsfdt;

  // velocity field + helper fields for WENO
  Real ux, uy, uz;
  Real omega, domegadt;

  // pressure + auxiliary functions for pressure source, phase filed funtion,
  // charact. function
  Real p, dpdt;
  Real f;
  Real chi;         // domain char. funciton
  Real pff, dpffdt; // pahse field function of whole anatomy, of tissue

  // other helper fields
  Real exact;
  Real tmp;

  // Volume Perception
  Real vp;

  // Cahn-Hilliard
  Real mob, mu; // mobility and chemical potetnial mu

  Cell() {
    phi = 0.0;
    dphidt = 0.0;
    p_g = 0.0;
    p_w = 0.0;
    p_csf = 0.0;
    wm = gm = csf = 0.0;
    dwmdt = dgmdt = dcsfdt = 0.0;
    ux = uy = uz = 0.0;
    p = 0.0;
    dpdt = 0.0;
    omega = 0.0;
    domegadt = 0.0;
    exact = 0.0;
    tmp = 0.0;
    f = 0.0;
    chi = 0.0;
    pff = 0.0;
    dpffdt = 0.0;
    vp = 0.0;
    mob = 0.0;
    mu = 0.0;
  }

  Cell(Real phi_, Real dphidt_, Real p_g_, Real p_w_, Real p_csf_, Real wm_,
       Real gm_, Real csf_, Real dwmdt_, Real dgmdt_, Real dcsfdt_, Real ux_,
       Real uy_, Real uz_, Real p_, Real dpdt_, Real omega_, Real domegadt_,
       Real exact_, Real tmp_, Real f_, Real chi_, Real pff_, Real dpffdt_,
       Real vp_, Real mob_, Real mu_) {
    phi = phi_;
    dphidt = dphidt_;
    p_g = p_g_;
    p_w = p_w_;
    p_csf = p_csf_;
    wm = wm_;
    gm = gm_;
    csf = csf_;
    dwmdt = dwmdt_;
    dgmdt = dgmdt_;
    dcsfdt = dcsfdt_;
    ux = ux_;
    uy = uy_;
    uz = uz_;
    p = p_;
    dpdt = dpdt_;
    omega = omega_;
    domegadt = domegadt_;
    exact = exact_;
    tmp = tmp_;
    f = f_;
    chi = chi_;
    pff = pff_;
    dpffdt = dpffdt_;
    vp = vp_;
    mob = mob_;
    mu = mu_;
  }

  void operator+=(Cell t) {
    phi += t.phi;
    dphidt += t.dphidt;
    p_g += t.p_g;
    p_w += t.p_w;
    p_csf += t.p_csf;
    wm += t.wm;
    gm += t.gm;
    csf += t.csf;
    dwmdt += t.dwmdt;
    dgmdt += t.dgmdt;
    dcsfdt += t.dcsfdt;
    ux += t.ux;
    uy += t.uy;
    uz += t.uz;
    p += t.p;
    dpdt += t.dpdt;
    omega += t.omega;
    domegadt += t.domegadt;
    exact += t.exact;
    tmp += t.tmp;
    f += t.f;
    chi += t.chi;
    pff += t.pff;
    dpffdt += t.dpffdt;
    vp += t.vp;
    mob += t.mob;
    mu += t.mu;
  }

  operator Real() { return (Real)phi; }

  void integrate(float dt) {
    phi += dt * dphidt;
    dphidt = 0.0;
  }

  template <int i> Real evaluate_concentration(double dt) {
    return phi + dt * dphidt;
  }

  Real giveMe(int i, Real h = 0) {

    switch (i) {
    case 0:
      return phi;
    case 1:
      return phi + 0.1 * p_g + 0.2 * p_w;
    case 2:
      return p_w;
    case 3:
      return p_g;
    case 4:
      return p_csf;
    case 5:
      return wm;
    case 6:
      return gm;
    case 7:
      return csf;
    case 8:
      return ux;
    case 9:
      return uy;
    case 10:
      return uz;
    case 11:
      return p;
    case 12:
      return chi;
    case 13:
      return f;
    case 14:
      return pff;

    default:
      abort();
      return 0;
    }
  }
};

inline Cell operator*(const Cell &p, Real v) {
  Cell c;
  c.phi = p.phi * v;
  c.dphidt = p.dphidt * v;
  c.p_g = p.p_g * v;
  c.p_w = p.p_w * v;
  c.p_csf = p.p_csf * v;
  c.wm = p.wm * v;
  c.gm = p.gm * v;
  c.csf = p.csf * v;
  c.dwmdt = p.dwmdt * v;
  c.dgmdt = p.dgmdt * v;
  c.dcsfdt = p.dcsfdt * v;
  c.ux = p.ux * v;
  c.uy = p.uy * v;
  c.uz = p.uz * v;
  c.p = p.p * v;
  c.dpdt = p.dpdt * v;
  c.omega = p.omega * v;
  c.domegadt = p.domegadt * v;
  c.exact = p.exact * v;
  c.tmp = p.tmp * v;
  c.f = p.f * v;
  c.chi = p.chi * v;
  c.pff = p.pff * v;
  c.dpffdt = p.dpffdt * v;
  c.vp = p.vp * v;
  c.mob = p.mob * v;
  c.mu = p.mu * v;

  return c;
}

template <typename T, int i> inline Real RD_projector_impl_vtk(const T &t) {
  //	return (Real)(0.2*t.p_w + 0.1*t.p_g );
  return (Real)(0.1 * t.p_g + 0.2 * t.p_w + t.p_csf + t.phi);
}

template <typename T, int i> inline Real RD_projector_impl_wav(const T &t) {
  // return i==0 ? (Real)(t.phi) : (Real)(t.p_w);  // for refinment w.r.t 2
  // channels
  return (Real)(t.phi);
  //    return (Real)(t.pff) ;
}

make_projector(RD_Projector_VTK, RD_projector_impl_vtk)
    make_projector(RD_Projector_Wavelets, RD_projector_impl_wav)

#ifndef _BLOCKSIZE_
#define _BLOCKSIZE_ 16
#endif

#ifndef _BLOCKSIZE_Z_
#define _BLOCKSIZE_Z_ _BLOCKSIZE_
#endif

#ifndef _BPD_
#define _BPD_ 4
#endif

#ifndef _MAXLEVEL_
#define _MAXLEVEL_ 3
#endif

#ifndef _RESJUMP_
#define _RESJUMP_ 1
#endif

        static const int blockSize = _BLOCKSIZE_;
static const int blockSizeZ = _BLOCKSIZE_Z_;
static const int blocksPerDimension = _BPD_;

// Structural parameters
static const bool bIsCellCentered = true;
static const bool bVerbose = true;

// Multiresolution parameters
static const int maxLevel = _MAXLEVEL_; // 8bpd has maxLevel 3 since 2^3
static const int resJump =
    _RESJUMP_; // modulo(maxLevel,resJum) = 0, !!! and reJump < maxLevel
const double refinement_tolerance = 1e-4;
const double compression_tolerance = 1e-5;

typedef Block<Cell, blockSize, blockSize, blockSizeZ> B;
typedef _WAVELET_TYPE W;

typedef Matrix::D3D<double> MatrixD3D;

static const int nThreads = 1;
typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
class Glioma_ReactionDiffusion {
private:
  Grid<W, B> *grid;
  BlockProcessing blockProcessing;
  Refiner_SpaceExtension *refiner;
  Compressor *compressor;
  BlockFWT<W, B, RD_Projector_Wavelets>
      blockfwt; // refinment based on single channel
  SpaceTimeSorter stSorter;
  Profiler profiler;
  ArgumentParser parser;
  IO_VTK<W, B, RD_Projector_VTK> vtk;
  BlockLab<B> lab;
  int numberOfIterations;
  double whenToWrite;
  double whenToWriteOffset;
  bool isDone;
  bool bAdaptivity;
  bool bVerbose;
  bool bVTK;
  bool bUQ;
  bool bDumpIC;
  string PatientFileName;
  Real L;
  Real tumor_ic[3];

  static void _ic(Grid<W, B> &grid, string PatientFileName, Real &L,
                  Real tumor_ic[3]);
  void _reactionDiffusionStep(BoundaryInfo *boundaryInfo,
                              const int nParallelGranularity, const Real Dw,
                              const Real Dg, const Real rho, double dt);
  void _dump(int counter);
  void _dumpUQoutput();

public:
  Glioma_ReactionDiffusion(int argc, const char **argv);
  ~Glioma_ReactionDiffusion();
  void run();
};

static int maxStencil[2][3] = {-1, -1, -1, +2, +2, +2};

Glioma_ReactionDiffusion::Glioma_ReactionDiffusion(int argc, const char **argv)
    : parser(argc, argv) {
  bVerbose = parser("-verbose").asBool(1);
  bVTK = parser("-vtk").asBool(1);
  bUQ = parser("-UQ").asBool(0);
  bDumpIC = parser("-bDumpIC").asBool(1);
  bAdaptivity = parser("-adaptive").asBool(1);
  PatientFileName = parser("-PatFileName").asString();
  if (bVerbose)
    printf("Set up: blockSize=%d Wavelets=w%s (blocksPerDimension=%d, "
           "maxLevel=%d)\n",
           blockSize, "w", blocksPerDimension, maxLevel);

  refiner = new Refiner_SpaceExtension(resJump, maxLevel);
  compressor = new Compressor(resJump);
  grid = new Grid<W, B>(blocksPerDimension, blocksPerDimension,
                        blocksPerDimension, maxStencil);
  grid->setCompressor(compressor);
  grid->setRefiner(refiner);
  stSorter.connect(*grid);

  L = 1;

  if (bUQ) {
    ifstream mydata("TumorIC.txt");

    if (mydata.is_open()) {
      mydata >> tumor_ic[0];
      mydata >> tumor_ic[1];
      mydata >> tumor_ic[2];
      mydata.close();
    } else {
      printf("Aborting: missing input file TumorIC.txt \n");
      abort();
    }
  } else {
    tumor_ic[0] = parser("-icx").asDouble(0.28);
    tumor_ic[1] = parser("-icy").asDouble(0.67);
    tumor_ic[2] = parser("-icz").asDouble(0.35);
  }

  _ic(*grid, PatientFileName, L, tumor_ic);

  if (bDumpIC)
    _dump(0);

  isDone = false;
  whenToWriteOffset = parser("-dumpfreq").asDouble();
  whenToWrite = whenToWriteOffset;
  numberOfIterations = 0;
}

Glioma_ReactionDiffusion::~Glioma_ReactionDiffusion() {
  std::cout << "------Adios muchachos------" << std::endl;
}

// 1) Read in anatomies - rescaled to [0,1]^3
// 2) Initialize tumor
// 3) Set the characteristic length L as the length of the data
void Glioma_ReactionDiffusion::_ic(Grid<W, B> &grid, string PatientFileName,
                                   Real &L, Real tumor_ic[3]) {
  printf("Reading data from file: %s \n", PatientFileName.c_str());

  char anatomy[200];
  sprintf(anatomy, "%sGM.dat", PatientFileName.c_str());
  MatrixD3D GM(anatomy);
  sprintf(anatomy, "%sWM.dat", PatientFileName.c_str());
  MatrixD3D WM(anatomy);
  sprintf(anatomy, "%sCSF.dat", PatientFileName.c_str());
  MatrixD3D CSF(anatomy);

  int brainSizeX = (int)GM.getSizeX();
  int brainSizeY = (int)GM.getSizeY();
  int brainSizeZ = (int)GM.getSizeZ();
  printf("brainSizeX=%i, brainSizeY=%i, brainSizeZ=%i \n", brainSizeX,
         brainSizeY, brainSizeZ);

  int brainSizeMax = max(brainSizeX, max(brainSizeY, brainSizeZ));
  L = brainSizeMax *
      0.1; // voxel spacing 1mm, convert from mm to cm  // L = 25.6 cm
  printf("Characteristic Lenght L=%f \n", L);

  double brainHx = 1.0 / ((double)(brainSizeMax)); //  w.r.t. longest dimension
                                                   //  for correct aspect ratio
  double brainHy = 1.0 / ((double)(brainSizeMax)); //  w.r.t. longest dimension
                                                   //  for correct aspect ratio
  double brainHz = 1.0 / ((double)(brainSizeMax)); //  w.r.t. longest dimension
                                                   //  for correct aspect ratio

  // Tumor set up
  const Real tumorRadius = 0.005;
  const Real smooth_sup = 2.; // suppor of smoothening
  const Real h =
      1. / 128; // use fixed h, for same IC smoothening at all resolutions
  const Real iw = 1. / (smooth_sup * h); // widht of smoothening

  Real pGM, pWM, pCSF;

  vector<BlockInfo> vInfo = grid.getBlocksInfo();

  for (int i = 0; i < vInfo.size(); i++) {
    BlockInfo &info = vInfo[i];
    B &block = grid.getBlockCollection()[info.blockID];

    for (int iz = 0; iz < B::sizeZ; iz++)
      for (int iy = 0; iy < B::sizeY; iy++)
        for (int ix = 0; ix < B::sizeX; ix++) {
          Real x[3];
          info.pos(x, ix, iy, iz);

          int mappedBrainX = (int)floor(x[0] / brainHx);
          int mappedBrainY = (int)floor(x[1] / brainHy);
          int mappedBrainZ = (int)floor(x[2] / brainHz);

          //                    // aspect ratio correction
          //                    mappedBrainX -= (int) ( (brainSizeMax -
          //                    brainSizeX) * 0.5); mappedBrainY -= (int) (
          //                    (brainSizeMax - brainSizeY) * 0.5); mappedBrainZ
          //                    -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);

          if ((mappedBrainX >= 0 && mappedBrainX < brainSizeX) &
                  (mappedBrainY >= 0 && mappedBrainY < brainSizeY) &&
              (mappedBrainZ >= 0 && mappedBrainZ < brainSizeZ)) {
            // anatomy
            pGM = GM(mappedBrainX, mappedBrainY, mappedBrainZ);
            pWM = WM(mappedBrainX, mappedBrainY, mappedBrainZ);
            pCSF = CSF(mappedBrainX, mappedBrainY, mappedBrainZ);

            // separat tissue and fluid based on majority voting
            double tissue = pWM + pGM;
            pCSF = (pCSF > tissue) ? 1. : 0.;
            pWM = (pCSF > tissue) ? 0. : pWM;
            pGM = (pCSF > tissue) ? 0. : pGM;

            tissue = pWM + pGM;
            block(ix, iy, iz).p_w = (tissue > 0.) ? (pWM / tissue) : 0.;
            block(ix, iy, iz).p_g = (tissue > 0.) ? (pGM / tissue) : 0.;
            block(ix, iy, iz).p_csf = pCSF;

            // tumor
            const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1],
                               x[2] - tumor_ic[2]};
            const Real dist =
                sqrt(p[0] * p[0] + p[1] * p[1] +
                     p[2] * p[2]); // distance of curent voxel from tumor center
            const Real psi = (dist - tumorRadius) * iw;

            if ((psi < -1) && (pGM + pWM > 0.001)) // we are in tumor
              block(ix, iy, iz).phi = 1.0;
            else if (((-1 <= psi) && (psi <= 1)) && (pGM + pWM > 0))
              block(ix, iy, iz).phi =
                  1.0 * 0.5 * (1 - psi - sin(M_PI * psi) / (M_PI));
            else
              block(ix, iy, iz).phi = 0.0;
          }
        }

    grid.getBlockCollection().release(info.blockID);
  }
}

void Glioma_ReactionDiffusion::_reactionDiffusionStep(
    BoundaryInfo *boundaryInfo, const int nParallelGranularity, const Real Dw,
    const Real Dg, const Real rho, double dt) {

  vector<BlockInfo> vInfo = grid->getBlocksInfo();
  const BlockCollection<B> &collecton = grid->getBlockCollection();

  ReactionDiffusionOperator rhs(Dw, Dg, rho);
  UpdateTumor updateTumor(dt);

  blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
  BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}

void Glioma_ReactionDiffusion::_dump(int counter) {
  if (bVTK) {
    char filename[256];
    sprintf(filename, "Data_%04d", counter);

    IO_VTKNative3D<W, B, 5, 0> vtkdumper2;
    vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
  }
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_ReactionDiffusion::_dumpUQoutput() {
  int gpd = blocksPerDimension * blockSize;
  double hf = 1. / gpd;
  double eps = hf * 0.5;

  MatrixD3D tumor(gpd, gpd, gpd);
  vector<BlockInfo> vInfo = grid->getBlocksInfo();

  for (int i = 0; i < vInfo.size(); i++) {
    BlockInfo &info = vInfo[i];
    B &block = grid->getBlockCollection()[info.blockID];
    double h = info.h[0];

    for (int iz = 0; iz < B::sizeZ; iz++)
      for (int iy = 0; iy < B::sizeY; iy++)
        for (int ix = 0; ix < B::sizeX; ix++) {
          double x[3];
          info.pos(x, ix, iy, iz);

          // mapped coordinates
          int mx = (int)floor((x[0]) / hf);
          int my = (int)floor((x[1]) / hf);
          int mz = (int)floor((x[2]) / hf);

          if (h < hf + eps)
            tumor(mx, my, mz) = block(ix, iy, iz).phi;
          else if (h < 2. * hf + eps) {
            for (int cz = 0; cz < 2; cz++)
              for (int cy = 0; cy < 2; cy++)
                for (int cx = 0; cx < 2; cx++)
                  tumor(mx + cx, my + cy, mz + cz) = block(ix, iy, iz).phi;
          } else if (h < 3. * hf + eps) {
            for (int cz = 0; cz < 3; cz++)
              for (int cy = 0; cy < 3; cy++)
                for (int cx = 0; cx < 3; cx++)
                  tumor(mx + cx, my + cy, mz + cz) = block(ix, iy, iz).phi;
          } else {
            for (int cz = 0; cz < 4; cz++)
              for (int cy = 0; cy < 4; cy++)
                for (int cx = 0; cx < 4; cx++)
                  tumor(mx + cx, my + cy, mz + cz) = block(ix, iy, iz).phi;
          }
        }
  }

  char filename2[256];
  sprintf(filename2, "HGG_data.dat");
  tumor.dump(filename2);
}

void Glioma_ReactionDiffusion::run() {
  const int nParallelGranularity = (grid->getBlocksInfo().size() <= 8 ? 1 : 4);
  BoundaryInfo *boundaryInfo = &grid->getBoundaryInfo();

  /* Tumor growth parameters*/
  Real Dw, Dg, rho, tend;

  if (bUQ) {
    ifstream mydata("InputParameters.txt");
    if (mydata.is_open()) {
      mydata >> Dw;
      mydata >> rho;
      mydata >> tend;
      mydata.close();
    } else {
      printf("Aborting: missing input file InputParameters.txt \n");
      abort();
    }
  } else {
    Dw = (Real)parser("-Dw").asDouble(0.0013);
    rho = (Real)parser("-rho").asDouble(0.025);
    tend = (Real)parser("-Tend").asDouble(300);
  }

  /*rescale for correct space dimension*/
  Dw = Dw / (L * L);
  Dg = 0.1 * Dw;

  Real t = 0.0;
  Real h = 1. / (blockSize * blocksPerDimension);
  Real dt = 0.99 * h * h / (2. * 3 * max(Dw, Dg));
  int iCounter = 1;
  if (bVerbose)
    printf("Dg=%e, Dw=%e, dt= %f, rho=%f , h=%f\n", Dg, Dw, dt, rho, h);

  // Initial compression, later just refinment since tumor just grow
  Science::AutomaticRefinement<0, 0>(*grid, blockfwt, refinement_tolerance,
                                     maxLevel, 1, &profiler);
  Science::AutomaticCompression<0, 0>(*grid, blockfwt, compression_tolerance,
                                      -1, &profiler);

  while (t <= tend) {
    _reactionDiffusionStep(boundaryInfo, nParallelGranularity, Dw, Dg, rho, dt);
    t += dt;
    numberOfIterations++;

    if (t >= ((double)(whenToWrite))) {
      if (bAdaptivity) {
        Science::AutomaticRefinement<0, 0>(
            *grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
        // Science::AutomaticCompression	<0,0>(*grid, blockfwt,
        // compression_tolerance, -1, &profiler);
      }

      _dump(iCounter++);
      whenToWrite = whenToWrite + whenToWriteOffset;
      if ((bVerbose) && (bVTK))
        printf("Dumping data at time t=%f\n", t);
    }
  }

  // Refine + save the last one
  if (bAdaptivity)
    Science::AutomaticRefinement<0, 0>(*grid, blockfwt, refinement_tolerance,
                                       maxLevel, 1, &profiler);

  _dump(iCounter);
  if (bUQ)
    _dumpUQoutput();

  if (bVerbose)
    printf("**** Dumping done\n");
  if (bVerbose)
    printf("\n\n Run Finished \n\n");
}

using namespace MRAG;
int main(int argc, const char **argv) {
  ArgumentParser parser(argc, argv);
  Environment::setup(max(1, parser("-nthreads").asInt()));
  Glioma_ReactionDiffusion s(argc, (const char **)argv);
  s.run();
}
