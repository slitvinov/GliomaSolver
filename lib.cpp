#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <set>
#include "MRAGcore/MRAGCommon.h"
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
const int blockSize = _BLOCKSIZE_;
typedef MRAG::Block<struct Cell, blockSize, blockSize, blockSize> B;
typedef MRAG::Wavelets_Interp2ndOrder W;
typedef MRAG::Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;

struct Cell {
  Real phi, dphidt, p_g, p_w;
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

static Cell operator*(const Cell &p, Real v) {
  Cell c;
  c.phi = p.phi * v;
  c.dphidt = p.dphidt * v;
  c.p_g = p.p_g * v;
  c.p_w = p.p_w * v;
  return c;
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
    Real df[6];
    Real chf[6];
    Real df_loc;
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

template <typename T, int i> inline Real RD_projector_impl_wav(const T &t) {
  return (Real)(t.phi);
}
make_projector(RD_Projector_Wavelets, RD_projector_impl_wav);

struct Brain {
  MRAG::Refiner_SpaceExtension *refiner;
  MRAG::Compressor *compressor;
  MRAG::Grid<W, B> *grid;
  MRAG::BlockFWT<W, B, RD_Projector_Wavelets> *blockfwt;
  MRAG::SpaceTimeSorter *stSorter;
  BlockProcessing *blockProcessing;
  ReactionDiffusionOperator *rhs;
  UpdateTumor *updateTumor;
};

int brain_ini(struct BrainParams *params, struct Brain **pbrain) {
  struct Brain *brain;
  int maxLevel = 4;
  int resJump = 1;
  double refinement_tolerance = 1e-4;
  double compression_tolerance = 1e-5;
  int maxStencil[2][3] = {-1, -1, -1, +2, +2, +2};
  int brainSizeMax;
  double brainHx, brainHy, brainHz;
  Real pGM, pWM;
  double tissue;
  int i, ix, iy, iz;
  Real x[3], dist, psi;
  int mX, mY, mZ;
  int index;
  double tumorRadius, smooth_sup, h0, iw;
  if ((brain = (struct Brain *)malloc(sizeof *brain)) == NULL) {
    fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
    return 1;
  }
  brain = new Brain;
  brain->grid = new MRAG::Grid<W, B>(params->blocksPerDimension, params->blocksPerDimension,
                                     params->blocksPerDimension, maxStencil);
  brain->blockfwt = new MRAG::BlockFWT<W, B, RD_Projector_Wavelets>;
  brain->blockProcessing = new BlockProcessing;
  brain->refiner = new MRAG::Refiner_SpaceExtension(resJump, maxLevel);
  brain->compressor = new MRAG::Compressor(resJump);
  brain->stSorter = new MRAG::SpaceTimeSorter;

  brain->grid->setCompressor(brain->compressor);
  brain->grid->setRefiner(brain->refiner);
  brain->stSorter->connect(*brain->grid);

  brainSizeMax = std::max(params->n[0], std::max(params->n[1], params->n[2]));
  brainHx = 1.0 / ((double)(brainSizeMax));
  brainHy = 1.0 / ((double)(brainSizeMax));
  brainHz = 1.0 / ((double)(brainSizeMax));
  tumorRadius = 0.005;
  smooth_sup = 2.;
  h0 = 1. / 128;
  iw = 1. / (smooth_sup * h0);
  vector<MRAG::BlockInfo> vInfo = brain->grid->getBlocksInfo();
  for (i = 0; i < vInfo.size(); i++) {
    MRAG::BlockInfo &info = vInfo[i];
    B &block = brain->grid->getBlockCollection()[info.blockID];
    for (iz = 0; iz < B::sizeZ; iz++)
      for (iy = 0; iy < B::sizeY; iy++)
        for (ix = 0; ix < B::sizeX; ix++) {
          info.pos(x, ix, iy, iz);
          mX = (int)floor(x[0] / brainHx);
          mY = (int)floor(x[1] / brainHy);
          mZ = (int)floor(x[2] / brainHz);
          if (mX >= 0 && mX < params->n[0] && mY >= 0 && mY < params->n[1] &&
              mZ >= 0 && mZ < params->n[2]) {
            index = mX + (mY + mZ * params->n[1]) * params->n[0];
            pGM = params->GM[index];
            pWM = params->WM[index];
            tissue = pWM + pGM;
            block(ix, iy, iz).p_w = (tissue > 0.) ? (pWM / tissue) : 0.;
            block(ix, iy, iz).p_g = (tissue > 0.) ? (pGM / tissue) : 0.;
            const Real p[3] = {x[0] - (Real)params->ic[0],
                               x[1] - (Real)params->ic[1],
                               x[2] - (Real)params->ic[2]};
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

    brain->grid->getBlockCollection().release(info.blockID);
  }
  MRAG::Science::AutomaticRefinement<0, 0>(*brain->grid, *brain->blockfwt,
                                           refinement_tolerance, maxLevel, 1);
  MRAG::Science::AutomaticCompression<0, 0>(*brain->grid, *brain->blockfwt,
                                            compression_tolerance, -1);
  brain->rhs =
      new ReactionDiffusionOperator(params->Dw, 0.1 * params->Dw, params->rho);
  brain->updateTumor = new UpdateTumor(params->dt);
  *pbrain = brain;
  return 0;
}

int brain_fin(struct Brain *brain) {
  delete brain->refiner;
  delete brain->compressor;
  delete brain->grid;
  delete brain->blockfwt;
  delete brain->stSorter;
  delete brain->blockProcessing;
  delete brain->rhs;
  delete brain->updateTumor;
  free(brain);
  return 0;
}

int brain_step(struct Brain *brain) {
  int maxLevel = 4;
  double refinement_tolerance = 1e-4;
  const MRAG::BlockCollection<B> &collecton = brain->grid->getBlockCollection();
  vector<MRAG::BlockInfo> vInfo = brain->grid->getBlocksInfo();
  MRAG::BoundaryInfo *boundaryInfo = &brain->grid->getBoundaryInfo();
  brain->blockProcessing->pipeline_process(vInfo, collecton, *boundaryInfo,
                                           *brain->rhs);
  BlockProcessing::process(vInfo, collecton, *brain->updateTumor);
  MRAG::Science::AutomaticRefinement<0, 0>(*brain->grid, *brain->blockfwt,
                                           refinement_tolerance, maxLevel, 1);
  return 0;
}

int brain_dump(struct Brain *brain, const char *path) {
  MRAG::BoundaryInfo *boundaryInfo = &brain->grid->getBoundaryInfo();
  write<W, B, MRAG::BlockLab<B>>(brain->grid, boundaryInfo, path);
  return 0;
}

int brain_project(struct Brain *brain, float *d) {
  int blocksPerDimension = 16;
  int ix, iy, iz, cx, cy, cz, mx, my, mz, gpd = blocksPerDimension * blockSize;
  Real x[3];
  double h, hf = 1. / gpd, eps = hf * 0.5;
  vector<MRAG::BlockInfo> vInfo = brain->grid->getBlocksInfo();
  for (int i = 0; i < vInfo.size(); i++) {
    MRAG::BlockInfo &info = vInfo[i];
    B &block = brain->grid->getBlockCollection()[info.blockID];
    h = info.h[0];
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
  return 0;
}

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
