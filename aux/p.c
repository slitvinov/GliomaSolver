#include <stdio.h>

int main() {
  int i;
  float pGM, pWM, pCSF;
  float Q[8][3] = {
      {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
      {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1},
  };
  for (i = 0; i < 8; i++) {
    pWM = Q[i][0];
    pGM = Q[i][1];
    pCSF = Q[i][2];
    printf("%g %g %g : ", pGM, pWM, pCSF);
    double tissue = pWM + pGM;
    pCSF = (pCSF > tissue) ? 1. : 0.;
    pWM = (pCSF > tissue) ? 0. : pWM;
    pGM = (pCSF > tissue) ? 0. : pGM;
    tissue = pWM + pGM;
    printf("%g %g %g\n", pGM, pWM, tissue);
  }
}
