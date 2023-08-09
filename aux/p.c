#include <stdio.h>

int main() {
  int i;
  float pGM, pWM, pCSF, tissue;
  float Q[8][3] = {
      {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
      {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1},
  };

  pGM = 1.0;
  pWM = 1.0;
  pCSF = 1.0;

  for (i = 0; i < 8; i++) {
    pWM = Q[i][0];
    pGM = Q[i][1];
    pCSF = Q[i][2];
    printf("%g %g %g : ", pGM, pWM, pCSF);
    tissue = pWM + pGM;
    pCSF = (pCSF > tissue) ? 1. : 0.;
    pWM = (pCSF > tissue) ? 0. : pWM;
    pGM = (pCSF > tissue) ? 0. : pGM;
    tissue = pWM + pGM;
    printf("%g %g %g\n", tissue, pGM, pWM);
  }
}
