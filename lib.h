#ifdef __cplusplus
extern "C" {
#endif

struct Brain;
float *brain_read(const char *, int *, int *, int *);
int brain_ini(int, int, int, float *, float *, double, double, struct Brain *);
int brain_step(struct Brain *, double);
int brain_dump(struct Brain *, char *path);
int brain_project(struct Brain *, float *);
int brain_fin(struct Brain*);

#ifdef __cplusplus
}
#endif
