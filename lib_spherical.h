#include <complex.h>
/*LEBEDEV GRID */
int get_lebedev_grid(double **theta, double **phi, double **wleb, int nleb, unsigned long *ngrid, char *path);
/* SPHERICAL HARMONICS */
double ylm_real(int l, int m, double th, double phi);
double *get_ylm_real_matrix(int lmax, double *th, double *phi, unsigned long nleb);
double *get_sph_real_coeff(double *v, double *wleb, double *y, unsigned long nleb, unsigned long nl);
double *get_ylm_real_orthonormality(double *wleb, double *y, unsigned long nl, unsigned long nleb);
double *get_rev_project_real(double *vll, double *y, unsigned long nl, unsigned long nleb);
double *get_sph_real_coeff_1d(double *vleb, double *wleb, double *y, unsigned long nl, unsigned long nleb);
double *get_ylm_real_vec(int lmax, double th, double phi);
double complex ylm(int l, int m, double th, double phi);
double complex *get_ylm_matrix(int lmax, double *th, double *phi, unsigned long nleb);
double complex *get_sph_coeff(double *v, double *wleb, double complex *y, unsigned long nleb, unsigned long nl);
double complex *get_ylm_orthonormality(double *wleb, double complex *y, unsigned long nl, unsigned long nleb);
double complex *get_rev_project(double complex *vll, double complex *y, unsigned long nl, unsigned long nleb);
/* ROATATION */
int rodrigues_rot_3d(double *v, double *z, double *vr);
/* TESTS */
double test_get_lebedev_grid(int leb, int tfun, char *leb_path);
double test_ylm_orthonormality(int leb, int lmax, char *leb_path);
double test_ylm_real_orthonormality(int leb, int lmax, char *leb_path);
double test_ylm_real_projection_1d(unsigned long np, int leb, int lmax, int tfun, char *path);
double test_ylm_real_projection_2d(unsigned long np, int leb, int lmax, int tfun, char *path);
