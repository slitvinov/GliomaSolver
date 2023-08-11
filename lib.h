#ifdef __cplusplus
extern "C" {
#endif

struct Brain;
float *brain_read(const char *, int *, int *, int *);
  int brain_ini(int, int, int, const float *, const float *, const double *, double, double, double, double, struct Brain **);
int brain_step(struct Brain *);
int brain_dump(struct Brain *, const char *);
int brain_project(struct Brain *, float *);
int brain_fin(struct Brain *);

#ifdef __cplusplus
}
#endif
