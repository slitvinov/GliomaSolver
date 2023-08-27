#include "lib.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
int main(void) {
  struct Brain *brain;
  double whenToWrite, whenToWriteOffset, tend, h, t, L, ans;
  float *model;
  char path[FILENAME_MAX - 9];
  int step, gpd, n[3];
  struct BrainParams params;
  struct LikelihoodParams lparams;

  params.blocksPerDimension = 32;
  params.ic[0] = 0.6497946102507519;
  params.ic[1] = 0.5908331665234543;
  params.ic[2] = 0.3715947899171972;
  params.GM = brain_read("GM.dat", params.n);
  params.WM = brain_read("WM.dat", params.n);
  gpd = (_BLOCKSIZE_)*params.blocksPerDimension;
  h = 1. / gpd;
  L = max(params.n[0], max(params.n[1], params.n[2])) * 0.1;
  params.Dw = 0.0013 / (L * L);
  params.rho = 0.025;
  tend = 300;
  params.dt = 0.99 * h * h / (2. * 3 * params.Dw);
  brain_ini(&params, &brain);
  free(params.GM);
  free(params.WM);
  whenToWriteOffset = 50;
  whenToWrite = whenToWriteOffset;
  t = 0.0;
  step = 0;
  while (t <= tend) {
    brain_step(brain);
    t += params.dt;
    step++;
    if (t >= whenToWrite) {
      sprintf(path, "a.%09d", step);
      brain_dump(brain, path);
      whenToWrite = whenToWrite + whenToWriteOffset;
    }
  }
  model = (float *)malloc(gpd * gpd * gpd * sizeof *model);
  if (model == NULL) {
    fprintf(stderr, "%s:%d: error: malloc failed\n", __FILE__, __LINE__);
    exit(1);
  }

  brain_project(brain, model);
  brain_fin(brain);

  lparams.PET = brain_read("tumPET.dat", n);
  lparams.T1c = brain_read("tumT1c.dat", n);
  lparams.FLAIR = brain_read("tumFLAIR.dat", n);
  assert(n[0] == gpd);
  assert(n[1] == gpd);
  assert(n[2] == gpd);
  lparams.n[0] = lparams.n[1] = lparams.n[2] = gpd;
  lparams.PETsigma2 = 0.000361;
  lparams.PETscale = 0.85;
  lparams.T1uc = 0.7;
  lparams.T2uc = 0.25;
  lparams.slope = 2;

  brain_likelihood(&lparams, model, &ans);
  free(model);
  free(lparams.PET);
  free(lparams.T1c);
  free(lparams.FLAIR);
  printf("%.16e\n", ans);
}
