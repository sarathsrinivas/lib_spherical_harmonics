#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib_spherical.h"

#define PI (3.1415926535897)

static int sph_to_cart(const double *v_sph, double *v_cart)
{

	v_cart[0] = v_sph[0] * cos(v_sph[2]) * sin(v_sph[1]);
	v_cart[1] = v_sph[0] * sin(v_sph[2]) * sin(v_sph[1]);
	v_cart[2] = v_sph[0] * cos(v_sph[1]);

	return 0;
}

static double get_distance(const double *v1, const double *v2)
{
	double dist;

	dist = (v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1])
	       + (v1[2] - v2[2]) * (v1[2] - v2[2]);

	return sqrt(dist);
}

static double get_dot_prod(const double *v1, const double *v2)
{
	double dotp;

	dotp = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

	return dotp;
}

int rodrigues_rot(const double *v, double th, const double *nc, double *vr)
{
	double cos_th, sin_th, vrc[3], nxv[3], vc[3], ndotv;

	vc[0] = cos(v[2]) * sin(v[1]);
	vc[1] = sin(v[2]) * sin(v[1]);
	vc[2] = cos(v[1]);

	cos_th = cos(th);
	sin_th = sin(th);

	nxv[0] = nc[1] * vc[2] - nc[2] * vc[1];
	nxv[1] = -nc[0] * vc[2] + nc[2] * vc[0];
	nxv[2] = nc[0] * vc[1] - nc[1] * vc[0];

	ndotv = nc[0] * vc[0] + nc[1] * vc[1] + nc[2] * vc[2];

	vrc[0] = vc[0] * cos_th + nxv[0] * sin_th + nc[0] * ndotv * (1 - cos_th);
	vrc[1] = vc[1] * cos_th + nxv[1] * sin_th + nc[1] * ndotv * (1 - cos_th);
	vrc[2] = vc[2] * cos_th + nxv[2] * sin_th + nc[2] * ndotv * (1 - cos_th);

	vr[0] = v[0] * sqrt(vrc[0] * vrc[0] + vrc[1] * vrc[1] + vrc[2] * vrc[2]);
	vr[1] = acos(vrc[2]);
	vr[2] = atan2(vrc[1], vrc[0]);

	return 0;
}

double test_rodrigues_rot(unsigned long nt)
{
	double p[3], k[3], kp[3], pc[3], kc[3], kpc[3], nc[3], pk[3], kk[3], kpk[3], pkc[3], kkc[3], kpkc[3],
	    a_k_kp, a_p_kp, a_k_p, a_kk_kpk, a_kk_pk, a_pk_kpk, d_k_kp, d_p_kp, d_k_p, d_kk_kpk, d_kk_pk,
	    d_pk_kpk, err;

	unsigned long i;

	err = 0;
	for (i = 0; i < nt; i++) {
		p[0] = 2.0;
		p[1] = 0;
		p[2] = 0;

		k[0] = 2.0;
		k[1] = 0.6 * PI;
		k[2] = 2 * PI * rand() / (1.0 + RAND_MAX);

		kp[0] = 2.0;
		kp[1] = PI * rand() / (1.0 + RAND_MAX);
		kp[2] = 2 * PI * rand() / (1.0 + RAND_MAX);

		sph_to_cart(k, kc);
		sph_to_cart(p, pc);
		sph_to_cart(kp, kpc);

		nc[0] = -sin(k[2]);
		nc[1] = cos(k[2]);
		nc[2] = 0;

		a_k_kp = get_dot_prod(kc, kpc);
		a_p_kp = get_dot_prod(pc, kpc);
		a_k_p = get_dot_prod(kc, pc);

		d_k_kp = get_distance(kc, kpc);
		d_p_kp = get_distance(pc, kpc);
		d_k_p = get_distance(kc, pc);

		rodrigues_rot(p, -k[1], nc, pk);
		rodrigues_rot(k, -k[1], nc, kk);
		rodrigues_rot(kp, -k[1], nc, kpk);

		sph_to_cart(kk, kkc);
		sph_to_cart(pk, pkc);
		sph_to_cart(kpk, kpkc);

		a_kk_kpk = get_dot_prod(kkc, kpkc);
		a_pk_kpk = get_dot_prod(pkc, kpkc);
		a_kk_pk = get_dot_prod(kkc, pkc);

		d_kk_kpk = get_distance(kkc, kpkc);
		d_pk_kpk = get_distance(pkc, kpkc);
		d_kk_pk = get_distance(kkc, pkc);

		err += fabs(a_k_kp - a_kk_kpk);
		err += fabs(a_k_p - a_kk_pk);
		err += fabs(a_p_kp - a_pk_kpk);

		err += fabs(d_k_kp - d_kk_kpk);
		err += fabs(d_k_p - d_kk_pk);
		err += fabs(d_p_kp - d_pk_kpk);

		printf("%+.15E %+.15E %+.15E %+.15E %+.15E\n", k[1] / PI, k[2] / PI, kp[1] / PI, kp[2] / PI,
		       err);
	}

	return err;
}
