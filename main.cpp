#include <math.h>
#include <assert.h>
#include <vector>
#include <set>
#include "MRAGcore/MRAGCommon.h"
#define _WAVELET_TYPE Wavelets_Interp2ndOrder
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGRefiner.h"
#include "MRAGcore/MRAGBlockLab.h"
#include "MRAGcore/MRAGBlockFWT.h"
#include "MRAGscience/MRAGScienceCore.h"
#include "MRAGcore/MRAGCommon.h"
#include "MRAGscience/MRAGAutomaticRefiner.h"
#include "MRAGscience/MRAGSpaceTimeSorter.h"
#include "MRAGscience/MRAGRefiner_SpaceExtension.h"
#include "MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
static int mNx, mNy, mNz, mNelements;
static float *D3D(const char *path) {
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
}
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
  inline void operator()(LabType &lab, const MRAG::BlockInfo &info,
                         BlockType &o) const {
    double h = info.h[0];
    double ih2 = 1. / (h * h);
    Real df[6];  // diffusion coefficient
    Real chf[6]; // domain charact. func, chf=0 -> outside, chf=1 inside domain:
                 // use to apply BC
    Real df_loc; // diffusion at the current point (local)
    Real chf_loc; // diffusion at the current point (local)
    Real eps;
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

            eps = 1.0e-08;
            /// Di,j = 2 * (Di * Dj / (Di + Dj)
            df[0] = (df[0] + df_loc < eps)
                        ? 0.
                        : 2. * df[0] * df_loc / (df[0] + df_loc);
            df[1] = (df[1] + df_loc < eps)
                        ? 0.
                        : 2. * df[1] * df_loc / (df[1] + df_loc);
            df[2] = (df[2] + df_loc < eps)
                        ? 0.
                        : 2. * df[2] * df_loc / (df[2] + df_loc);
            df[3] = (df[3] + df_loc < eps)
                        ? 0.
                        : 2. * df[3] * df_loc / (df[3] + df_loc);
            df[4] = (df[4] + df_loc < eps)
                        ? 0.
                        : 2. * df[4] * df_loc / (df[4] + df_loc);
            df[5] = (df[5] + df_loc < eps)
                        ? 0.
                        : 2. * df[5] * df_loc / (df[5] + df_loc);

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

            eps = 0.1;
            if (chf[0] < eps) {
              df[1] *= 2.0;
            }
            if (chf[1] < eps) {
              df[0] *= 2.0;
            }
            if (chf[2] < eps) {
              df[3] *= 2.0;
            }
            if (chf[3] < eps) {
              df[2] *= 2.0;
            }
            if (chf[4] < eps) {
              df[5] *= 2.0;
            }
            if (chf[5] < eps) {
              df[4] *= 2.0;
            }
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
};

struct UpdateTumor {
  double dt;

  UpdateTumor(double dt_) : dt(dt_) {}

  UpdateTumor(const UpdateTumor &copy) : dt(copy.dt) {}

  template <typename BlockType>
  inline void operator()(const MRAG::BlockInfo &info, BlockType &o) const {
    for (int iz = 0; iz < BlockType::sizeZ; iz++)
      for (int iy = 0; iy < BlockType::sizeY; iy++)
        for (int ix = 0; ix < BlockType::sizeX; ix++) {
          o(ix, iy, iz).phi += dt * o(ix, iy, iz).dphidt;
          o(ix, iy, iz).phi = max((Real)0., o(ix, iy, iz).phi);
          o(ix, iy, iz).phi = min((Real)1., o(ix, iy, iz).phi);
        }
  }
};

struct Cell {
  /* tumor */
  Real phi;
  Real dphidt;

  /* tissue percentage per voxel*/
  Real p_g, p_w, p_csf;

  /* tissue concentration */
  Real wm, gm, csf; //
  Real dwmdt, dgmdt, dcsfdt;

  // pressure + auxiliary functions for pressure source, phase filed funtion,
  // charact. function
  Real f;
  Real chi;         // domain char. funciton
  Real pff, dpffdt; // pahse field function of whole anatomy, of tissue

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
    f = 0.0;
    chi = 0.0;
    pff = 0.0;
    dpffdt = 0.0;
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
    f += t.f;
    chi += t.chi;
    pff += t.pff;
    dpffdt += t.dpffdt;
    mob += t.mob;
    mu += t.mu;
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
  c.f = p.f * v;
  c.chi = p.chi * v;
  c.pff = p.pff * v;
  c.dpffdt = p.dpffdt * v;
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

make_projector(RD_Projector_Wavelets, RD_projector_impl_wav);

static const int blockSize = _BLOCKSIZE_;
static const int blockSizeZ = _BLOCKSIZE_;
static const int blocksPerDimension = _BPD_;
static const int maxLevel = _MAXLEVEL_;
static const int resJump = _RESJUMP_;
const double refinement_tolerance = 1e-4;
const double compression_tolerance = 1e-5;
typedef MRAG::Block<Cell, blockSize, blockSize, blockSizeZ> B;
typedef MRAG::_WAVELET_TYPE W;
static const int nThreads = 1;
typedef MRAG::Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
static MRAG::Grid<W, B> *grid;
static BlockProcessing blockProcessing;
static MRAG::Refiner_SpaceExtension *refiner;
static MRAG::Compressor *compressor;
static MRAG::BlockFWT<W, B, RD_Projector_Wavelets> blockfwt;
static MRAG::SpaceTimeSorter stSorter;
static MRAG::BlockLab<B> lab;
static double whenToWrite;
static double whenToWriteOffset;
static bool isDone;
static Real L;
static Real tumor_ic[3];
static int maxStencil[2][3] = {-1, -1, -1, +2, +2, +2};
int main(int argc, const char **argv) {
  refiner = new MRAG::Refiner_SpaceExtension(resJump, maxLevel);
  compressor = new MRAG::Compressor(resJump);
  grid = new MRAG::Grid<W, B>(blocksPerDimension, blocksPerDimension,
                              blocksPerDimension, maxStencil);
  grid->setCompressor(compressor);
  grid->setRefiner(refiner);
  stSorter.connect(*grid);

  L = 1;
  tumor_ic[0] = 0.28;
  tumor_ic[1] = 0.75;
  tumor_ic[2] = 0.35;
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
  const Real h0 = 1. / 128;
  const Real iw = 1. / (smooth_sup * h0);
  Real pGM, pWM, pCSF;
  vector<MRAG::BlockInfo> vInfo = grid->getBlocksInfo();
  for (int i = 0; i < vInfo.size(); i++) {
    MRAG::BlockInfo &info = vInfo[i];
    B &block = grid->getBlockCollection()[info.blockID];

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
            int index =
                mappedBrainX + (mappedBrainY + mappedBrainZ * mNy) * mNx;
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

    grid->getBlockCollection().release(info.blockID);
  }

  isDone = false;
  whenToWriteOffset = 50;
  whenToWrite = whenToWriteOffset;

  const int nParallelGranularity = (grid->getBlocksInfo().size() <= 8 ? 1 : 4);
  MRAG::BoundaryInfo *boundaryInfo = &grid->getBoundaryInfo();
  Real Dw, Dg, rho, tend;
  Dw = 0.0013;
  rho = 0.025;
  tend = 300;
  Dw = Dw / (L * L);
  Dg = 0.1 * Dw;

  Real t = 0.0;
  Real h = 1. / (blockSize * blocksPerDimension);
  Real dt = 0.99 * h * h / (2. * 3 * max(Dw, Dg));
  int iCounter = 1;
  MRAG::Science::AutomaticRefinement<0, 0>(*grid, blockfwt,
                                           refinement_tolerance, maxLevel, 1);
  MRAG::Science::AutomaticCompression<0, 0>(*grid, blockfwt,
                                            compression_tolerance, -1);
  ReactionDiffusionOperator rhs(Dw, Dg, rho);
  UpdateTumor updateTumor(dt);
  const MRAG::BlockCollection<B> &collecton = grid->getBlockCollection();
  while (t <= tend) {
    vInfo = grid->getBlocksInfo();
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor,
                             nParallelGranularity);
    t += dt;
    if (t >= whenToWrite) {
      MRAG::Science::AutomaticRefinement<0, 0>(
          *grid, blockfwt, refinement_tolerance, maxLevel, 1);
      whenToWrite = whenToWrite + whenToWriteOffset;
    }
  }
  MRAG::Science::AutomaticRefinement<0, 0>(*grid, blockfwt,
                                           refinement_tolerance, maxLevel, 1);
  float *d;
  FILE *file;
  int gpd = blocksPerDimension * blockSize;
  double hf = 1. / gpd;
  double eps = hf * 0.5;
  d = (float *)malloc(gpd * gpd * gpd * sizeof *d);
  vInfo = grid->getBlocksInfo();
  for (int i = 0; i < vInfo.size(); i++) {
    MRAG::BlockInfo &info = vInfo[i];
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
