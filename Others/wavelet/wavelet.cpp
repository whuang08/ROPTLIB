#include "Others/wavelet/wavelet.h"

void haarFWT_1d(int n, realdpcomplex *v)
{
	realdp r2 = static_cast<realdp> (sqrt(2.0));
	realdp *tmp = new realdp[2 * n];
	for (integer i = 0; i < 2 * n; i++)
		tmp[i] = 0;

	integer j = 1;
	while (j * 2 <= n)
	{
		j = j * 2;
	}

	while (1 < j)
	{
		j = j / 2;
		for (integer i = 0; i < j; i++)
		{
			tmp[2 * i] = (v[2 * i].r + v[2 * i + 1].r) / r2;
			tmp[2 * i + 1] = (v[2 * i].i + v[2 * i + 1].i) / r2;
			tmp[2 * (i + j)] = (v[2 * i].r - v[2 * i + 1].r) / r2;
			tmp[2 * (i + j) + 1] = (v[2 * i].i - v[2 * i + 1].i) / r2;
		}
		for (integer i = 0; i < 2 * j; i++)
		{
			v[i].r = tmp[2 * i];
			v[i].i = tmp[2 * i + 1];
		}
	}

	delete[] tmp;

	return;
};

void haarFWT_1d_inverse(int n, realdpcomplex *v)
{
	realdp r2 = static_cast<realdp> (sqrt(2.0));
	realdp *tmp = new realdp[2 * n];
	for (integer i = 0; i < 2 * n; i++)
		tmp[i] = 0;

	integer j = 1;
	while (j * 2 <= n)
	{
		for (integer i = 0; i < j; i++)
		{
			tmp[2 * (2 * i)] = (v[i].r + v[i + j].r) / r2;
			tmp[2 * (2 * i) + 1] = (v[i].i + v[i + j].i) / r2;
			tmp[2 * (2 * i + 1)] = (v[i].r - v[i + j].r) / r2;
			tmp[2 * (2 * i + 1) + 1] = (v[i].i - v[i + j].i) / r2;
		}
		for (integer i = 0; i < j * 2; i++)
		{
			v[i].r = tmp[2 * i];
			v[i].i = tmp[2 * i + 1];
		}
		j = j * 2;
	}

	delete [] tmp;

	return;
};

void haarFWT_2d(int n1, int n2, realdpcomplex *vv)
{
	realdp r2 = static_cast<realdp> (sqrt(2.0));

	realdpcomplex *tmp = new realdpcomplex[n1 * n2];

	for (integer i = 0; i < n1 * n2; i++)
	{
		tmp[i].r = vv[i].r;
		tmp[i].i = vv[i].i;
	}

	integer k = 1;
	while (k * 2 <= n1)
	{
		k = k * 2;
	}
	while (1 < k)
	{
		k = k / 2;

		for (integer j = 0; j < n2; j++)
		{
			for (integer i = 0; i < k; i++)
			{
				tmp[i + j * n1].r = (vv[2 * i + j * n1].r + vv[2 * i + 1 + j * n1].r) / r2;
				tmp[i + j * n1].i = (vv[2 * i + j * n1].i + vv[2 * i + 1 + j * n1].i) / r2;
				tmp[k + i + j * n1].r = (vv[2 * i + j*n1].r - vv[2 * i + 1 + j * n1].r) / r2;
				tmp[k + i + j * n1].i = (vv[2 * i + j*n1].i - vv[2 * i + 1 + j * n1].i) / r2;
			}
		}
		for (integer j = 0; j < n2; j++)
		{
			for (integer i = 0; i < 2 * k; i++)
			{
				vv[i + j * n1].r = tmp[i + j * n1].r;
				vv[i + j * n1].i = tmp[i + j * n1].i;
			}
		}
	}
	k = 1;
	while (k * 2 <= n2)
	{
		k = k * 2;
	}
	while (1 < k)
	{
		k = k / 2;

		for (integer j = 0; j < k; j++)
		{
			for (integer i = 0; i < n1; i++)
			{
				tmp[i + j * n1].r = (vv[i + 2 * j * n1].r + vv[i + (2 * j + 1) * n1].r) / r2;
				tmp[i + j * n1].i = (vv[i + 2 * j * n1].i + vv[i + (2 * j + 1) * n1].i) / r2;
				tmp[i + (k + j) * n1].r = (vv[i + 2 * j*n1].r - vv[i + (2 * j + 1) * n1].r) / r2;
				tmp[i + (k + j) * n1].i = (vv[i + 2 * j*n1].i - vv[i + (2 * j + 1) * n1].i) / r2;
			}
		}

		for (integer j = 0; j < 2 * k; j++)
		{
			for (integer i = 0; i < n1; i++)
			{
				vv[i + j*n1].r = tmp[i + j*n1].r;
				vv[i + j*n1].i = tmp[i + j*n1].i;
			}
		}
	}
	delete[] tmp;
};

void haarFWT_2d_inverse(int n1, int n2, realdpcomplex *vv)
{
	realdp r2 = static_cast<realdp> (sqrt(2.0));

	realdpcomplex *tmp = new realdpcomplex[n1 * n2];

	for (integer j = 0; j < n2; j++)
	{
		for (integer i = 0; i < n1; i++)
		{
			tmp[i + j * n1] = vv[i + j * n1];
		}
	}
	integer k = 1;

	while (k * 2 <= n2)
	{
		for (integer j = 0; j < k; j++)
		{
			for (integer i = 0; i < n1; i++)
			{
				tmp[i + (2 * j) * n1].r = (vv[i + j * n1].r + vv[i + (k + j)*n1].r) / r2;
				tmp[i + (2 * j) * n1].i = (vv[i + j * n1].i + vv[i + (k + j)*n1].i) / r2;
				tmp[i + (2 * j + 1) * n1].r = (vv[i + j * n1].r - vv[i + (k + j)*n1].r) / r2;
				tmp[i + (2 * j + 1) * n1].i = (vv[i + j * n1].i - vv[i + (k + j)*n1].i) / r2;
			}
		}

		for (integer j = 0; j < 2 * k; j++)
		{
			for (integer i = 0; i < n1; i++)
			{
				vv[i + j * n1].r = tmp[i + j * n1].r;
				vv[i + j * n1].i = tmp[i + j * n1].i;
			}
		}
		k = k * 2;
	}

	k = 1;
	while (k * 2 <= n1)
	{
		for (integer j = 0; j < n2; j++)
		{
			for (integer i = 0; i < k; i++)
			{
				tmp[2 * i + j * n1].r = (vv[i + j * n1].r + vv[k + i + j * n1].r) / r2;
				tmp[2 * i + j * n1].i = (vv[i + j * n1].i + vv[k + i + j * n1].i) / r2;
				tmp[2 * i + 1 + j * n1].r = (vv[i + j * n1].r - vv[k + i + j * n1].r) / r2;
				tmp[2 * i + 1 + j * n1].i = (vv[i + j * n1].i - vv[k + i + j * n1].i) / r2;
			}
		}

		for (integer j = 0; j < n2; j++)
		{
			for (integer i = 0; i < 2 * k; i++)
			{
				vv[i + j * n1].r = tmp[i + j * n1].r;
				vv[i + j * n1].i = tmp[i + j * n1].i;
			}
		}
		k = k * 2;
	}
	delete[] tmp;
};
