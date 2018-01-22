#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib_spherical.h"

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
