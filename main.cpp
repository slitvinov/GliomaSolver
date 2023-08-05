#include "MRAGHeaders.h"
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGBlockCollection.h"
#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGGridNode.h"
#include "MRAGcore/MRAGrid.h"
#include <fstream>
using namespace MRAG;
static int mNx, mNy, mNz, mNelements;
static float* D3D(const char *path) {
    float *mData;
    FILE *file;
    int header[6];
    file = fopen(path, "r");
    fread(header, sizeof header, 1, file);
    assert(header[0] == 1234);
    assert(header[1] == 3);
    assert(header[5] == 1);
    mNelements = header[2] * header[3] * header[4];
    mNx = header[2];
    mNy = header[3];
    mNz = header[4];
    mData = (float *)malloc(mNelements * sizeof *mData);
    fread(mData, mNelements, sizeof *mData, file);
    fclose(file);
    return mData;
};
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

template <typename T, int i> inline Real RD_projector_impl_wav(const T &t) {
  // return i==0 ? (Real)(t.phi) : (Real)(t.p_w);  // for refinment w.r.t 2
  // channels
  return (Real)(t.phi);
  //    return (Real)(t.pff) ;
}

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
static const bool bIsCellCentered = true;
static const int maxLevel = _MAXLEVEL_;
static const int resJump = _RESJUMP_;
const double refinement_tolerance = 1e-4;
const double compression_tolerance = 1e-5;
typedef Block<Cell, blockSize, blockSize, blockSizeZ> B;
typedef _WAVELET_TYPE W;
static const int nThreads = 1;
typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
class Glioma_ReactionDiffusion {
private:
  Grid<W, B> *grid;
  BlockProcessing blockProcessing;
  Refiner_SpaceExtension *refiner;
  Compressor *compressor;
  BlockFWT<W, B, RD_Projector_Wavelets> blockfwt;
  SpaceTimeSorter stSorter;
  ArgumentParser parser;
  BlockLab<B> lab;
  int numberOfIterations;
  double whenToWrite;
  double whenToWriteOffset;
  bool isDone;
  bool bAdaptivity;
  bool bVerbose;
  bool bDumpIC;
  string PatientFileName;
  Real L;
  Real tumor_ic[3];

  static void _ic(Grid<W, B> &grid, string PatientFileName, Real &L,
                  Real tumor_ic[3]);
  void _reactionDiffusionStep(BoundaryInfo *boundaryInfo,
                              const int nParallelGranularity, const Real Dw,
                              const Real Dg, const Real rho, double dt);
  void _dumpUQoutput();

public:
  Glioma_ReactionDiffusion(int argc, const char **argv);
  void run();
};

static int maxStencil[2][3] = {-1, -1, -1, +2, +2, +2};

Glioma_ReactionDiffusion::Glioma_ReactionDiffusion(int argc, const char **argv)
    : parser(argc, argv) {
  bVerbose = parser("-verbose").asBool(1);
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

  _ic(*grid, PatientFileName, L, tumor_ic);

  isDone = false;
  whenToWriteOffset = parser("-dumpfreq").asDouble();
  whenToWrite = whenToWriteOffset;
  numberOfIterations = 0;
}
void Glioma_ReactionDiffusion::_ic(Grid<W, B> &grid, string PatientFileName,
                                   Real &L, Real tumor_ic[3]) {
  float *GM, *WM, *CSF;
  GM = D3D("GM.dat");
  WM = D3D("WM.dat");
  CSF = D3D("CSF.dat");
  int brainSizeX = mNx;
  int brainSizeY = mNy;
  int brainSizeZ = mNz;
  printf("brainSizeX=%i, brainSizeY=%i, brainSizeZ=%i \n", brainSizeX,
         brainSizeY, brainSizeZ);

  int brainSizeMax = max(brainSizeX, max(brainSizeY, brainSizeZ));
  L = brainSizeMax * 0.1;
  printf("Characteristic Lenght L=%f \n", L);
  double brainHx = 1.0 / ((double)(brainSizeMax));
  double brainHy = 1.0 / ((double)(brainSizeMax));
  double brainHz = 1.0 / ((double)(brainSizeMax));
  const Real tumorRadius = 0.005;
  const Real smooth_sup = 2.;
  const Real h = 1. / 128;
  const Real iw = 1. / (smooth_sup * h);
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
          if ((mappedBrainX >= 0 && mappedBrainX < brainSizeX) &
                  (mappedBrainY >= 0 && mappedBrainY < brainSizeY) &&
              (mappedBrainZ >= 0 && mappedBrainZ < brainSizeZ)) {
	    int index = mappedBrainX + (mappedBrainY + mappedBrainZ * mNy) * mNx;
            pGM = GM[index];
            pWM = WM[index];
            pCSF = CSF[index];
            double tissue = pWM + pGM;
            pCSF = (pCSF > tissue) ? 1. : 0.;
            pWM = (pCSF > tissue) ? 0. : pWM;
            pGM = (pCSF > tissue) ? 0. : pGM;
            tissue = pWM + pGM;
            block(ix, iy, iz).p_w = (tissue > 0.) ? (pWM / tissue) : 0.;
            block(ix, iy, iz).p_g = (tissue > 0.) ? (pGM / tissue) : 0.;
            block(ix, iy, iz).p_csf = pCSF;
            const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1],
                               x[2] - tumor_ic[2]};
            const Real dist = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            const Real psi = (dist - tumorRadius) * iw;
            if ((psi < -1) && (pGM + pWM > 0.001))
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

void Glioma_ReactionDiffusion::_dumpUQoutput() {
  float *d;
  FILE *file;
  int gpd = blocksPerDimension * blockSize;
  double hf = 1. / gpd;
  double eps = hf * 0.5;
  d = (float *)malloc(gpd * gpd * gpd * sizeof *d);
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
          int mx = (int)floor((x[0]) / hf);
          int my = (int)floor((x[1]) / hf);
          int mz = (int)floor((x[2]) / hf);
          if (h < hf + eps) {
            d[mx + (my + mz * gpd) * gpd] = block(ix, iy, iz).phi;
          } else if (h < 2. * hf + eps) {
            for (int cz = 0; cz < 2; cz++)
              for (int cy = 0; cy < 2; cy++)
                for (int cx = 0; cx < 2; cx++) {
                  d[mx + cx + (my + cy + (mz + cz) * gpd) * gpd] =
                      block(ix, iy, iz).phi;
                }
          } else if (h < 3. * hf + eps) {
            for (int cz = 0; cz < 3; cz++)
              for (int cy = 0; cy < 3; cy++)
                for (int cx = 0; cx < 3; cx++) {
                  d[mx + cx + (my + cy + (mz + cz) * gpd) * gpd] =
                      block(ix, iy, iz).phi;
                }
          } else {
            for (int cz = 0; cz < 4; cz++)
              for (int cy = 0; cy < 4; cy++)
                for (int cx = 0; cx < 4; cx++) {
                  d[mx + cx + (my + cy + (mz + cz) * gpd) * gpd] =
                      block(ix, iy, iz).phi;
                }
          }
        }
  }
  file = fopen("HGG_data.dat", "w");
  int header[6] = {1234, 3, gpd, gpd, gpd, 1};
  fwrite(header, sizeof header, 1, file);
  fwrite(d, gpd * gpd * gpd, sizeof *d, file);
  fclose(file);
  free(d);
}

void Glioma_ReactionDiffusion::run() {
  const int nParallelGranularity = (grid->getBlocksInfo().size() <= 8 ? 1 : 4);
  BoundaryInfo *boundaryInfo = &grid->getBoundaryInfo();
  Real Dw, Dg, rho, tend;

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
  Dw = Dw / (L * L);
  Dg = 0.1 * Dw;

  Real t = 0.0;
  Real h = 1. / (blockSize * blocksPerDimension);
  Real dt = 0.99 * h * h / (2. * 3 * max(Dw, Dg));
  int iCounter = 1;
  if (bVerbose)
    printf("Dg=%e, Dw=%e, dt= %f, rho=%f , h=%f\n", Dg, Dw, dt, rho, h);
  Science::AutomaticRefinement<0, 0>(*grid, blockfwt, refinement_tolerance,
                                     maxLevel, 1);
  Science::AutomaticCompression<0, 0>(*grid, blockfwt, compression_tolerance,
                                      -1);

  while (t <= tend) {
    _reactionDiffusionStep(boundaryInfo, nParallelGranularity, Dw, Dg, rho, dt);
    t += dt;
    numberOfIterations++;

    if (t >= ((double)(whenToWrite))) {
      if (bAdaptivity) {
        Science::AutomaticRefinement<0, 0>(*grid, blockfwt,
                                           refinement_tolerance, maxLevel, 1);
        // Science::AutomaticCompression	<0,0>(*grid, blockfwt,
        // compression_tolerance, -1, &profiler);
      }
      whenToWrite = whenToWrite + whenToWriteOffset;
    }
  }
  if (bAdaptivity)
    Science::AutomaticRefinement<0, 0>(*grid, blockfwt, refinement_tolerance,
                                       maxLevel, 1);
  _dumpUQoutput();
}

int main(int argc, const char **argv) {
  Glioma_ReactionDiffusion s(argc, (const char **)argv);
  s.run();
}
