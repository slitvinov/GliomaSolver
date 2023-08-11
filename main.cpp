#include <math.h>
#include <assert.h>
#include <vector>
#include <set>
#include "MRAGcore/MRAGCommon.h"
#define _WAVELET_TYPE Wavelets_Interp2ndOrder
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
using namespace std;
#include "MRAGcore/MRAGCompressor.h"
#include "MRAGcore/MRAGBoundaryBlockInfo.h"
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
#include "write.h"
#include "lib.h"

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
    Real eps;
    int iz, iy, ix;
    double diffusionFluxIn, diffusionFluxOut, reactionFlux;
    for (iz = 0; iz < BlockType::sizeZ; iz++)
      for (iy = 0; iy < BlockType::sizeY; iy++)
        for (ix = 0; ix < BlockType::sizeX; ix++) {
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
            diffusionFluxIn = ih2 * (df[0] * lab(ix - 1, iy, iz).phi +
                                     df[1] * lab(ix + 1, iy, iz).phi +
                                     df[2] * lab(ix, iy - 1, iz).phi +
                                     df[3] * lab(ix, iy + 1, iz).phi +
                                     df[4] * lab(ix, iy, iz - 1).phi +
                                     df[5] * lab(ix, iy, iz + 1).phi);

            diffusionFluxOut =
                -((df[0] + df[1] + df[2] + df[3] + df[4] + df[5]) *
                  lab(ix, iy, iz).phi * ih2);
            reactionFlux =
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
  inline void operator()(const MRAG::BlockInfo &, BlockType &o) const {
    int ix, iy, iz;
    for (iz = 0; iz < BlockType::sizeZ; iz++)
      for (iy = 0; iy < BlockType::sizeY; iy++)
        for (ix = 0; ix < BlockType::sizeX; ix++) {
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
  Real p_g, p_w;

  Cell() {
    phi = 0.0;
    dphidt = 0.0;
    p_g = 0.0;
    p_w = 0.0;
  }

  void operator+=(Cell t) {
    phi += t.phi;
    dphidt += t.dphidt;
    p_g += t.p_g;
    p_w += t.p_w;
  }
};

inline Cell operator*(const Cell &p, Real v) {
  Cell c;
  c.phi = p.phi * v;
  c.dphidt = p.dphidt * v;
  c.p_g = p.p_g * v;
  c.p_w = p.p_w * v;
  return c;
}

template <typename T, int i> inline Real RD_projector_impl_wav(const T &t) {
  return (Real)(t.phi);
}
make_projector(RD_Projector_Wavelets, RD_projector_impl_wav);

int main(int, char **) {
  struct Brain *brain;
  const int blockSize = _BLOCKSIZE_;
  const int blockSizeZ = _BLOCKSIZE_;
  int blocksPerDimension = 16;
  int maxLevel = 4;
  int resJump = 1;;
  const double refinement_tolerance = 1e-4;
  const double compression_tolerance = 1e-5;
  typedef MRAG::Block<Cell, blockSize, blockSize, blockSizeZ> B;
  typedef MRAG::_WAVELET_TYPE W;
  typedef MRAG::Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
  BlockProcessing blockProcessing;
  MRAG::BlockFWT<W, B, RD_Projector_Wavelets> blockfwt;
  MRAG::SpaceTimeSorter stSorter;
  MRAG::BlockLab<B> lab;
  double whenToWrite;
  double whenToWriteOffset;
  Real L;
  double ic[3];
  int maxStencil[2][3] = {-1, -1, -1, +2, +2, +2};
  MRAG::Refiner_SpaceExtension refiner(resJump, maxLevel);
  MRAG::Compressor compressor(resJump);
  MRAG::Grid<W, B> grid(blocksPerDimension, blocksPerDimension,
                        blocksPerDimension, maxStencil);
  float *GM, *WM;
  int brainSizeMax;
  double brainHx, brainHy, brainHz;
  Real pGM, pWM;
  double tissue;
  int i, ix, iy, iz, cx, cy, cz;
  Real x[3], dist, psi;
  char path[FILENAME_MAX - 9];
  int step;
  int mappedBrainX, mappedBrainY, mappedBrainZ;
  int index;
  int mx, my, mz;
  int nx, ny, nz;
  Real rho, tend;
  double Dw, Dg;
  double tumorRadius, smooth_sup, h0, iw;

  grid.setCompressor(&compressor);
  grid.setRefiner(&refiner);
  stSorter.connect(grid);

  L = 1;
  ic[0] = 0.6497946102507519;
  ic[1] = 0.5908331665234543;
  ic[2] = 0.3715947899171972;
  GM = brain_read("GM.dat", &nx, &ny, &nz);
  WM = brain_read("WM.dat", &nx, &ny, &nz);
  Dw = 0.0013;
  rho = 0.025;
  tend = 300;
  brain_ini(nx, ny, nz, GM, WM, Dw, rho, &brain);
  brainSizeMax = max(nx, max(ny, nz));
  L = brainSizeMax * 0.1;
  printf("Characteristic Lenght L=%f \n", L);
  brainHx = 1.0 / ((double)(brainSizeMax));
  brainHy = 1.0 / ((double)(brainSizeMax));
  brainHz = 1.0 / ((double)(brainSizeMax));
  tumorRadius = 0.005;
  smooth_sup = 2.;
  h0 = 1. / 128;
  iw = 1. / (smooth_sup * h0);
  vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
  for (i = 0; i < vInfo.size(); i++) {
    MRAG::BlockInfo &info = vInfo[i];
    B &block = grid.getBlockCollection()[info.blockID];
    for (iz = 0; iz < B::sizeZ; iz++)
      for (iy = 0; iy < B::sizeY; iy++)
        for (ix = 0; ix < B::sizeX; ix++) {
          info.pos(x, ix, iy, iz);

          mappedBrainX = (int)floor(x[0] / brainHx);
          mappedBrainY = (int)floor(x[1] / brainHy);
          mappedBrainZ = (int)floor(x[2] / brainHz);
          if ((mappedBrainX >= 0 && mappedBrainX < nx) &
                  (mappedBrainY >= 0 && mappedBrainY < ny) &&
              (mappedBrainZ >= 0 && mappedBrainZ < nz)) {
            index = mappedBrainX + (mappedBrainY + mappedBrainZ * ny) * nx;
            pGM = GM[index];
            pWM = WM[index];
            tissue = pWM + pGM;
            tissue = pWM + pGM;
            block(ix, iy, iz).p_w = (tissue > 0.) ? (pWM / tissue) : 0.;
            block(ix, iy, iz).p_g = (tissue > 0.) ? (pGM / tissue) : 0.;
            const Real p[3] = {x[0] - (Real)ic[0], x[1] - (Real)ic[1],
                               x[2] - (Real)ic[2]};
            dist = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            psi = (dist - tumorRadius) * iw;
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
  free(GM);
  free(WM);

  whenToWriteOffset = 50;
  whenToWrite = whenToWriteOffset;
  MRAG::BoundaryInfo *boundaryInfo = &grid.getBoundaryInfo();
  Dw = Dw / (L * L);
  Dg = 0.1 * Dw;
  Real t = 0.0;
  Real h = 1. / (blockSize * blocksPerDimension);
  Real dt = 0.99 * h * h / (2. * 3 * max(Dw, Dg));
  MRAG::Science::AutomaticRefinement<0, 0>(grid, blockfwt, refinement_tolerance,
                                           maxLevel, 1);
  MRAG::Science::AutomaticCompression<0, 0>(grid, blockfwt,
                                            compression_tolerance, -1);
  ReactionDiffusionOperator rhs(Dw, Dg, rho);
  UpdateTumor updateTumor(dt);
  const MRAG::BlockCollection<B> &collecton = grid.getBlockCollection();
  step = 0;
  while (t <= tend) {
    vInfo = grid.getBlocksInfo();
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor);
    t += dt;
    step++;
    if (t >= whenToWrite) {
      MRAG::Science::AutomaticRefinement<0, 0>(
          grid, blockfwt, refinement_tolerance, maxLevel, 1);
      sprintf(path, "a.%09d", step);
      write<W, B, MRAG::BlockLab<B>>(&grid, boundaryInfo, path);
      whenToWrite = whenToWrite + whenToWriteOffset;
    }
  }
  MRAG::Science::AutomaticRefinement<0, 0>(grid, blockfwt, refinement_tolerance,
                                           maxLevel, 1);
  float *d;
  FILE *file;
  int gpd = blocksPerDimension * blockSize;
  double hf = 1. / gpd;
  double eps = hf * 0.5;
  d = (float *)malloc(gpd * gpd * gpd * sizeof *d);
  vInfo = grid.getBlocksInfo();
  for (int i = 0; i < vInfo.size(); i++) {
    MRAG::BlockInfo &info = vInfo[i];
    B &block = grid.getBlockCollection()[info.blockID];
    double h = info.h[0];
    for (iz = 0; iz < B::sizeZ; iz++)
      for (iy = 0; iy < B::sizeY; iy++)
        for (ix = 0; ix < B::sizeX; ix++) {
          info.pos(x, ix, iy, iz);
          mx = (int)floor((x[0]) / hf);
          my = (int)floor((x[1]) / hf);
          mz = (int)floor((x[2]) / hf);
          if (h < hf + eps) {
            d[mx + (my + mz * gpd) * gpd] = block(ix, iy, iz).phi;
          } else if (h < 2. * hf + eps) {
            for (cz = 0; cz < 2; cz++)
              for (cy = 0; cy < 2; cy++)
                for (cx = 0; cx < 2; cx++) {
                  d[mx + cx + (my + cy + (mz + cz) * gpd) * gpd] =
                      block(ix, iy, iz).phi;
                }
          } else if (h < 3. * hf + eps) {
            for (cz = 0; cz < 3; cz++)
              for (cy = 0; cy < 3; cy++)
                for (cx = 0; cx < 3; cx++) {
                  d[mx + cx + (my + cy + (mz + cz) * gpd) * gpd] =
                      block(ix, iy, iz).phi;
                }
          } else {
            for (cz = 0; cz < 4; cz++)
              for (cy = 0; cy < 4; cy++)
                for (cx = 0; cx < 4; cx++) {
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
  brain_fin(brain);
}
