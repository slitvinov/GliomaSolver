#include <stdio.h>
#include <stdlib.h>
#include "lib.h"

const int blockSize = _BLOCKSIZE_;
#define max(a, b) (a) > (b) ? (a) : (b)

int main(void) {
  struct Brain *brain;
  int blocksPerDimension = 16;
  double whenToWrite;
  double whenToWriteOffset;
  double ic[3];
  float *GM, *WM, *d;
  char path[FILENAME_MAX - 9];
  int nx, ny, nz, step, gpd, brainSizeMax;
  int32_t header[6];
  double Dw, Dg, rho, tend, h, L, t, dt;
  FILE *file;

  ic[0] = 0.6497946102507519;
  ic[1] = 0.5908331665234543;
  ic[2] = 0.3715947899171972;
  GM = brain_read("GM.dat", &nx, &ny, &nz);
  WM = brain_read("WM.dat", &nx, &ny, &nz);
  Dw = 0.0013;
  rho = 0.025;
  tend = 300;
  L = 1;
  Dw = Dw / (L * L);
  Dg = 0.1 * Dw;
  h = 1. / (blockSize * blocksPerDimension);
  dt = 0.99 * h * h / (2. * 3 * max(Dw, Dg));
  brain_ini(nx, ny, nz, GM, WM, ic, Dw, Dg, rho, dt, &brain);
  free(GM);
  free(WM);
  whenToWriteOffset = 50;
  whenToWrite = whenToWriteOffset;
  t = 0.0;
  brainSizeMax = max(nx, max(ny, nz));
  L = brainSizeMax * 0.1;
  h = 1. / (blockSize * blocksPerDimension);
  Dw = Dw / (L * L);
  Dg = 0.1 * Dw;
  dt = 0.99 * h * h / (2. * 3 * max(Dw, Dg));
  step = 0;
  while (t <= tend) {
    brain_step(brain);
    t += dt;
    step++;
    if (t >= whenToWrite) {
      sprintf(path, "a.%09d", step);
      brain_dump(brain, path);
      whenToWrite = whenToWrite + whenToWriteOffset;
    }
  }
  gpd = blocksPerDimension * blockSize;
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
