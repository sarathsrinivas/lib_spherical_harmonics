#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include "lib_spherical.h"

double complex ylm(int l, int m, double th, double phi)
{
	double complex y, II, ret;
	double fac;
	II = (double complex)I;
	y = gsl_sf_legendre_sphPlm(l, abs(m), cos(th)) * cexp(m * phi * II);
	if (m >= 0)
		fac = 1.0;
	else
		fac = (m % 2) ? -1.0 : 1.0;

	ret = fac * y;
	return ret;
}

double ylm_real(int l, int m, double th, double phi)
{
	double complex y, II;
	double pf, yr[3], ret, sm[2], plm;
	int id;

	pf = sqrt(2);
	II = (double complex)I;

	plm = gsl_sf_legendre_sphPlm(l, abs(m), cos(th));
	y = plm * cexp(abs(m) * phi * II);

	yr[0] = cimag(y);
	yr[1] = creal(y);
	yr[2] = plm;

	sm[0] = 1.0;
	sm[1] = -1.0;

	if (m == 0) {
		id = 2;
		pf = 1.0;
	} else {
		id = (m + abs(m)) / (2 * m);
	}

	ret = pf * sm[abs(m % 2)] * yr[id];

	return ret;
}

double *get_ylm_real_vec(int lmax, double th, double phi)
{
	double *ylm_vec;
	unsigned long i, nl;
	int l, m;

	nl = (unsigned long)(lmax * lmax);
	ylm_vec = malloc(nl * sizeof(double));
	assert(ylm_vec);

	i = 0;
	for (l = 0; l < lmax; l++) {
		for (m = -1 * l; m <= l; m++) {
			ylm_vec[i] = ylm_real(l, m, th, phi);
			i++;
		}
	}

	return ylm_vec;
}

double *get_sph_real_coeff_1d(double *vleb, double *wleb, double *y, unsigned long nl, unsigned long nleb)
{
	double *vl, tmp;
	unsigned long i, l;

	vl = malloc(nl * sizeof(double));
	assert(vl);

	for (l = 0; l < nl; l++) {
		tmp = 0;
		for (i = 0; i < nleb; i++) {
			tmp += wleb[i] * vleb[i] * y[l * nleb + i];
		}
		vl[l] = tmp;
	}

	return vl;
}

double complex *get_ylm_matrix(int lmax, double *th, double *phi, unsigned long nleb)
{
	unsigned long nl, i, j;
	int l, m;
	double complex *y;

	nl = (unsigned long)(lmax * lmax);
	y = malloc(nl * nleb * sizeof(double complex));
	assert(y);

	i = 0;
	for (l = 0; l < lmax; l++) {
		for (m = -1 * l; m <= l; m++) {
			for (j = 0; j < nleb; j++) {
				y[i * nleb + j] = ylm(l, m, th[j], phi[j]);
			}
			i++;
		}
	}

	return y;
}

double *get_ylm_real_matrix(int lmax, double *th, double *phi, unsigned long nleb)
{
	unsigned long nl, i, j;
	int l, m;
	double *y;

	nl = (unsigned long)(lmax * lmax);
	y = malloc(nl * nleb * sizeof(double));
	assert(y);

	i = 0;
	for (l = 0; l < lmax; l++) {
		for (m = -1 * l; m <= l; m++) {
			for (j = 0; j < nleb; j++) {
				y[i * nleb + j] = ylm_real(l, m, th[j], phi[j]);
			}
			i++;
		}
	}

	return y;
}

double complex *get_sph_coeff(double *v, double *wleb, double complex *y, unsigned long nleb,
			      unsigned long nl)
{
	double complex *vll, yli, ylpj, vllp;
	double wi, wj, vij;
	unsigned long l, lp, i, j;

	vll = malloc(nl * nl * sizeof(double complex));
	assert(vll);

	for (l = 0; l < nl; l++) {
		for (lp = 0; lp < nl; lp++) {
			vllp = 0;
			for (i = 0; i < nleb; i++) {
				yli = y[l * nleb + i];
				wi = wleb[i];
				for (j = 0; j < nleb; j++) {
					ylpj = y[lp * nleb + j];
					wj = wleb[j];
					vij = v[i * nleb + j];

					vllp += wi * wj * vij * yli * ylpj;
				}
			}
			vll[l * nl + lp] = vllp;
		}
	}

	return vll;
}

double *get_sph_real_coeff(double *v, double *wleb, double *y, unsigned long nleb, unsigned long nl)
{
	double *vll, yli, ylpj, vllp;
	double wi, wj, vij;
	unsigned long l, lp, i, j;

	vll = malloc(nl * nl * sizeof(double));
	assert(vll);

	for (l = 0; l < nl; l++) {
		for (lp = 0; lp < nl; lp++) {
			vllp = 0;
			for (i = 0; i < nleb; i++) {
				yli = y[l * nleb + i];
				wi = wleb[i];
				for (j = 0; j < nleb; j++) {
					ylpj = y[lp * nleb + j];
					wj = wleb[j];
					vij = v[i * nleb + j];

					vllp += wi * wj * vij * yli * ylpj;
				}
			}
			vll[l * nl + lp] = vllp;
		}
	}

	return vll;
}

double complex *get_ylm_orthonormality(double *wleb, double complex *y, unsigned long nl, unsigned long nleb)
{
	double complex *orth, tmp, yli, ylpi;
	unsigned long i, l, lp;
	orth = malloc(nl * nl * sizeof(double complex));
	assert(orth);
	for (l = 0; l < nl; l++) {
		for (lp = 0; lp < nl; lp++) {
			tmp = 0;
			for (i = 0; i < nleb; i++) {
				yli = y[l * nleb + i];
				ylpi = y[lp * nleb + i];

				tmp += wleb[i] * yli * conj(ylpi);
			}
			orth[l * nl + lp] = tmp;
		}
	}

	return orth;
}

double *get_ylm_real_orthonormality(double *wleb, double *y, unsigned long nl, unsigned long nleb)
{
	double *orth, tmp, yli, ylpi;
	unsigned long i, l, lp;
	orth = malloc(nl * nl * sizeof(double));
	assert(orth);
	for (l = 0; l < nl; l++) {
		for (lp = 0; lp < nl; lp++) {
			tmp = 0;
			for (i = 0; i < nleb; i++) {
				yli = y[l * nleb + i];
				ylpi = y[lp * nleb + i];

				tmp += wleb[i] * yli * ylpi;
			}
			orth[l * nl + lp] = tmp;
		}
	}

	return orth;
}

double complex *get_rev_project(double complex *vll, double complex *y, unsigned long nl, unsigned long nleb)
{
	double complex *v, vij, vllp, yli, ylpj;
	unsigned long i, j, l, lp;

	v = malloc(nleb * nleb * sizeof(double complex));
	assert(v);

	for (i = 0; i < nleb; i++) {
		for (j = 0; j < nleb; j++) {
			vij = 0;
			for (l = 0; l < nl; l++) {
				yli = y[l * nleb + i];
				for (lp = 0; lp < nl; lp++) {
					ylpj = y[lp * nleb + j];
					vllp = vll[l * nl + lp];

					vij += vllp * yli * ylpj;
				}
			}
			v[i * nleb + j] = vij;
		}
	}
	return v;
}

double *get_rev_project_real(double *vll, double *y, unsigned long nl, unsigned long nleb)
{
	double *v, vij, vllp, yli, ylpj;
	unsigned long i, j, l, lp;

	v = malloc(nleb * nleb * sizeof(double));
	assert(v);

	for (i = 0; i < nleb; i++) {
		for (j = 0; j < nleb; j++) {
			vij = 0;
			for (l = 0; l < nl; l++) {
				yli = y[l * nleb + i];
				for (lp = 0; lp < nl; lp++) {
					ylpj = y[lp * nleb + j];
					vllp = vll[l * nl + lp];

					vij += vllp * yli * ylpj;
				}
			}
			v[i * nleb + j] = vij;
		}
	}
	return v;
}

int rodrigues_rot_3d(double *v, double *z, double *vr)
{
	double nx, ny, nz, vx, vy, vz, cos_th, sin_th, nxv_x, nxv_y, nxv_z, ndotv, vrx, vry, vrz, v_th, v_phi,
	    z_th, z_phi;

	v_th = v[1];
	v_phi = v[2];

	z_th = z[1];
	z_phi = z[2];

	vx = cos(v_phi) * sin(v_th);
	vy = sin(v_phi) * sin(v_th);
	vz = cos(v_th);

	cos_th = cos(z_th);
	sin_th = sin(z_th);

	nx = sin(z_phi);
	ny = -cos(z_phi);
	nz = 0;

	nxv_x = ny * vz;
	nxv_y = -nx * vz;
	nxv_z = nx * vy - ny * vx;

	ndotv = nx * vx + ny * vy + nz * vz;

	vrx = vx * cos_th + nxv_x * sin_th + nx * ndotv * (1 - cos_th);
	vry = vy * cos_th + nxv_y * sin_th + ny * ndotv * (1 - cos_th);
	vrz = vz * cos_th + nxv_z * sin_th + nz * ndotv * (1 - cos_th);

	vr[0] = v[0] * sqrt(vrx * vrx + vry * vry + vrz * vrz);
	vr[1] = acos(vrz);
	vr[2] = atan2(vry, vrz);

	return 0;
}
