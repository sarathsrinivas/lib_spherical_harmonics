#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include "lebedev.h"
#include "ylm.h"
#include "tests.h"
#include "../lib_arr/lib_arr.h"

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
