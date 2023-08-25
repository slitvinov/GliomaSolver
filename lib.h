#ifdef __cplusplus
extern "C" {
#endif

struct Brain;
struct BrainParams {
  int blocksPerDimension;
  int n[3];
  float *GM, *WM;
  double ic[3];
  double Dw, rho, dt;
};
struct LikelihoodParams {
  int n[3];
  float *PET, *T1c, *FLAIR;
  double PETsigma2, PETscale, slope, T1uc, T2uc;
};
float *brain_read(const char *, int *);
int brain_ini(struct BrainParams *, struct Brain **);
int brain_step(struct Brain *);
int brain_dump(struct Brain *, const char *);
int brain_project(struct Brain *, float *);
int brain_fin(struct Brain *);
int brain_likelihood(struct LikelihoodParams *, float *, double *);

#ifdef __cplusplus
}
#endif
