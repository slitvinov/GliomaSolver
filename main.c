#include <stdio.h>
#include <stdlib.h>
#include "lib.h"

const int blockSize = _BLOCKSIZE_;
#define max(a, b) ((a) > (b) ? (a) : (b))

int main(void) {
  struct Brain *brain;
  double whenToWrite, whenToWriteOffset, tend, h, t;
  float *d;
  char path[FILENAME_MAX - 9];
  int nx, ny, nz, step, gpd;
  int32_t header[6];
  FILE *file;
  struct BrainParams params;

  params.blocksPerDimension = 16;
  params.ic[0] = 0.6497946102507519;
  params.ic[1] = 0.5908331665234543;
  params.ic[2] = 0.3715947899171972;
  params.GM = brain_read("GM.dat", params.n);
  params.WM = brain_read("WM.dat", params.n);
  params.Dw = 0.0013;
  params.rho = 0.025;
  tend = 300;
  gpd = blockSize * params.blocksPerDimension;
  h = 1. / gpd;
  params.dt = 0.99 * h * h / (2. * 3 * params.Dw);

  params.blocksPerDimension = 16;
  params.n[0] = nx;
  params.n[1] = ny;
  params.n[2] = nz;

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
  d = (float *)malloc(gpd * gpd * gpd * sizeof *d);
  brain_project(brain, d);
  file = fopen("HGG_data.dat", "w");
  header[0] = 1234;
  header[1] = 3;
  header[2] = gpd;
  header[3] = gpd;
  header[4] = gpd;
  header[5] = 1;
  fwrite(header, sizeof header, 1, file);
  fwrite(d, gpd * gpd * gpd, sizeof *d, file);
  fclose(file);
  free(d);
  brain_fin(brain);
}
