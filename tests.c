#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include "lib_spherical.h"
#include "../lib_arr/lib_arr.h"

#define PI (3.1415926535897)

double test_get_lebedev_grid(int leb, int tfun, char *leb_path)
{

	double *th, *phi, *wleb, x, y, z, x2, y2, z2, II[4], Icomp[4], Idiff[4];
	unsigned long ngrid, l;

	fprintf(stderr, "test_get_lebedev_grid() test #%d %s:%d\n", tfun, __FILE__, __LINE__);

	Icomp[0] = 19.388114662154152;
	Icomp[1] = 0;
	Icomp[2] = 12.566370614359172;
	Icomp[3] = 12.566370614359172;

	II[0] = 0;
	II[1] = 0;
	II[2] = 0;
	II[3] = 0;

	get_lebedev_grid(&th, &phi, &wleb, leb, &ngrid, leb_path);
	assert(th);
	assert(phi);
	assert(wleb);

	for (l = 0; l < ngrid; l++) {
		x = cos(phi[l]) * sin(th[l]);
		y = sin(phi[l]) * sin(th[l]);
		z = cos(th[l]);

		x2 = x * x;
		y2 = y * y;
		z2 = z * z;

		II[0] += wleb[l] * (1 + x + y2 + x2 * y + x2 * x2 + y2 * y2 * y + x2 * y2 * z2);
		II[1] += wleb[l] * x * y * z;
		II[2] += wleb[l] * (1 + tanh(z - x - y));
		II[3] += wleb[l];
	}

	for (l = 0; l < 4; l++) {
		Idiff[l] = fabs(II[l] - Icomp[l]);
	}

	free(th);
	free(phi);
	free(wleb);

	return Idiff[tfun];
}

double test_ylm_orthonormality(int leb, int lmax, char *leb_path)
{
	double complex *orth, *y;
	double *th, *phi, *wleb, ret;
	unsigned long nl, nleb;

	fprintf(stderr, "test_ylm_orthonormality() %s:%d\n", __FILE__, __LINE__);

	nl = (unsigned long)(lmax * lmax);

	get_lebedev_grid(&th, &phi, &wleb, leb, &nleb, leb_path);
	assert(th);
	assert(phi);
	assert(wleb);

	y = get_ylm_matrix(lmax, th, phi, nleb);
	assert(y);

	orth = get_ylm_orthonormality(wleb, y, nl, nleb);
	assert(orth);

	ret = test_identity_complex(orth, nl);

	free(th);
	free(phi);
	free(wleb);
	free(y);
	free(orth);

	return ret;
}

double test_ylm_real_orthonormality(int leb, int lmax, char *leb_path)
{
	double *orth, *y, *th, *phi, *wleb, ret;
	unsigned long nl, nleb;

	fprintf(stderr, "test_ylm_real_orthonormality() %s:%d\n", __FILE__, __LINE__);

	nl = (unsigned long)(lmax * lmax);

	get_lebedev_grid(&th, &phi, &wleb, leb, &nleb, leb_path);
	assert(th);
	assert(phi);
	assert(wleb);

	y = get_ylm_real_matrix(lmax, th, phi, nleb);
	assert(y);

	orth = get_ylm_real_orthonormality(wleb, y, nl, nleb);
	assert(orth);

	ret = test_identity_real(orth, nl);

	free(th);
	free(phi);
	free(wleb);
	free(y);
	free(orth);

	return ret;
}

static double fun(int tfun, double th, double phi)
{
	double ret[3], r, x, y, z, x2, y2, z2;

	r = 1.0;
	x = r * cos(phi) * sin(th);
	y = r * sin(phi) * sin(th);
	z = r * cos(th);

	x2 = x * x;
	y2 = y * y;
	z2 = z * z;

	ret[0] = 1 + x + y2 + x2 * y + x2 * x2 + y2 * y2 * y + x2 * y2 * z2;
	ret[1] = x * y * z;
	ret[2] = 1 + tanh(0.1 * (z - x - y));

	return ret[tfun];
}

double test_ylm_real_projection_1d(unsigned long np, int leb, int lmax, int tfun, char *path)
{
	double *th, *phi, *wleb, *vleb, *y, *yp, *vl, *vp, *va, *thp, *phip, dnorm, tmp;
	unsigned long nleb, l, i, nl;

	fprintf(stderr, "test_ylm_real_projection_1d() test #%d %s:%d\n", tfun, __FILE__, __LINE__);

	nl = (unsigned long)(lmax * lmax);

	get_lebedev_grid(&th, &phi, &wleb, leb, &nleb, path);
	assert(th);
	assert(phi);
	assert(wleb);

	vleb = malloc(nleb * sizeof(double));
	assert(vleb);

	for (i = 0; i < nleb; i++) {
		vleb[i] = fun(tfun, th[i], phi[i]);
	}

	y = get_ylm_real_matrix(lmax, th, phi, nleb);
	assert(y);

	vl = get_sph_real_coeff_1d(vleb, wleb, y, nl, nleb);
	assert(vl);

	thp = malloc(np * sizeof(double));
	assert(thp);
	phip = malloc(np * sizeof(double));
	assert(phip);
	vp = malloc(np * sizeof(double));
	assert(vp);
	va = malloc(np * sizeof(double));
	assert(va);

	for (i = 0; i < np; i++) {
		thp[i] = PI * rand() / RAND_MAX;
		phip[i] = 2 * PI * rand() / RAND_MAX;
		va[i] = fun(tfun, thp[i], phip[i]);
	}

	yp = get_ylm_real_matrix(lmax, thp, phip, np);
	assert(yp);

	for (i = 0; i < np; i++) {
		tmp = 0;
		for (l = 0; l < nl; l++) {
			tmp += vl[l] * yp[l * np + i];
		}
		vp[i] = tmp;
	}

	for (i = 0; i < np; i++) {
		vp[i] = fabs(vp[i] - va[i]);
	}

	dnorm = get_norm_real(vp, np);

	free(yp);
	free(va);
	free(vp);
	free(phip);
	free(thp);
	free(vl);
	free(y);
	free(vleb);
	free(wleb);
	free(phi);
	free(th);

	return dnorm;
}

static double fun2d(double k_th, double k_phi, double kp_th, double kp_phi, int tfun)
{
	double v[3], k, kp, lam, kx, ky, kz, kpx, kpy, kpz, kd2;

	k = 4.0;
	kp = 3.0;
	kx = k * cos(k_phi) * sin(k_th);
	ky = k * sin(k_phi) * sin(k_th);
	kz = k * cos(k_th);

	kpx = kp * cos(kp_phi) * sin(kp_th);
	kpy = kp * sin(kp_phi) * sin(kp_th);
	kpz = kp * cos(kp_th);

	kd2 = (kx - kpx) * (kx - kpx) + (ky - kpy) * (ky - kpy) + (kz - kpz) * (kz - kpz);

	lam = 5000.0;

	v[0] = lam / (938 * 938 + kd2 * 197.0);
	v[1] = kd2;
	v[2] = kx * ky * kx + kpx * kpy * kpz;

	return v[tfun];
}

double test_ylm_real_projection_2d(unsigned long np, int leb, int lmax, int tfun, char *path)
{
	double *vleb, *vll, *y, *yp, *th, *phi, *wleb, *thp, *phip, *va, *vp, dnorm;
	unsigned long i, j, nl, np2, nleb;

	fprintf(stderr, "test_ylm_real_projection_2d() test #%d %s:%d\n", tfun, __FILE__, __LINE__);

	nl = (unsigned long)(lmax * lmax);

	get_lebedev_grid(&th, &phi, &wleb, leb, &nleb, path);
	assert(th);
	assert(phi);
	assert(wleb);

	vleb = malloc(nleb * nleb * sizeof(double));
	assert(vleb);

	for (i = 0; i < nleb; i++) {
		for (j = 0; j < nleb; j++) {
			vleb[i * nleb + j] = fun2d(th[i], phi[i], th[j], phi[j], tfun);
		}
	}

	y = get_ylm_real_matrix(lmax, th, phi, nleb);
	assert(y);

	vll = get_sph_real_coeff(vleb, wleb, y, nleb, nl);
	assert(vll);

	thp = malloc(np * sizeof(double));
	assert(thp);
	phip = malloc(np * sizeof(double));
	assert(phip);
	va = malloc(np * np * sizeof(double));
	assert(va);

	for (i = 0; i < np; i++) {
		thp[i] = PI * rand() / RAND_MAX;
		phip[i] = 2 * PI * rand() / RAND_MAX;
	}

	for (i = 0; i < np; i++) {
		for (j = 0; j < np; j++) {
			va[i * np + j] = fun2d(thp[i], phip[i], thp[j], phip[j], tfun);
		}
	}

	yp = get_ylm_real_matrix(lmax, thp, phip, np);
	assert(yp);

	vp = get_rev_project_real(vll, yp, nl, np);
	assert(vp);

	np2 = np * np;
	for (i = 0; i < np2; i++) {
		vp[i] = fabs(va[i] - vp[i]);
	}

	dnorm = get_norm_real(vp, np2);

	free(yp);
	free(vp);
	free(va);
	free(phip);
	free(thp);
	free(vll);
	free(y);
	free(vleb);
	free(wleb);
	free(phi);
	free(th);

	return dnorm;
}
