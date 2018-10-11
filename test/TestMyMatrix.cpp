
#include "test/TestMyMatrix.h"

using namespace ROPTLIB;

void testEigenSymmetricM(void)
{
	double *M = new double[16 + 4 + 16];
	double *eigvalues = M + 16;
	double *eigvectors = eigvalues + 4;
	for (integer i = 0; i < 4; i++)
	{
		for (integer j = i; j < 4; j++)
		{
			M[i + j * 4] = genrandnormal();
			M[j + i * 4] = M[i + j * 4];
		}
	}
	Matrix A(M, 4, 4), E(eigvalues, 4, 1), V(eigvectors, 4, 4);
	Matrix::EigenSymmetricM(GLOBAL::U, A, E, V);
	delete[] M;
};

void testExpSymmetricM(void)
{
	double *M = new double[16 + 16];
	double *ExpM = M + 16;
	for (integer i = 0; i < 4; i++)
	{
		for (integer j = i; j < 4; j++)
		{
			M[i + j * 4] = genrandnormal();
			M[j + i * 4] = M[i + j * 4];
		}
	}
	Matrix A(M, 4, 4), B(ExpM, 4, 4);
	Matrix::ExpSymmetricM(GLOBAL::U, A, B);
	delete[] M;
};

void testLogSymmetricM(void)
{
	integer N = 4;
	double *M = new double[16 + 16 + 16];
	double *MMt = M + 16;
	double *LogM = MMt + 16;
	for (integer i = 0; i < 4; i++)
	{
		for (integer j = i; j < 4; j++)
		{
			M[i + j * 4] = genrandnormal();
			M[j + i * 4] = M[i + j * 4];
		}
	}
	dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, M, &N, M, &N, &GLOBAL::DZERO, MMt, &N);
	Matrix A(M, 4, 4), B(LogM, 4, 4), C(MMt, 4, 4);
	Matrix::LogSymmetricM(GLOBAL::U, C, B);
	delete[] M;
};

void testDSYL(void)
{
	integer N = 4;
	double *A = new double[16 + 16 + 16];
	double *B = A + 16;
	double *C = B + 16;
	for (integer i = 0; i < 4; i++)
	{
		for (integer j = 0; j < 4; j++)
		{
			A[i + j * 4] = genrandnormal();
			B[i + j * 4] = genrandnormal();
			C[i + j * 4] = genrandnormal();
		}
	}

	Matrix MA(A, 4, 4), MB(B, 4, 4), MC(C, 4, 4);
	std::cout << "A:" << MA << std::endl;
	std::cout << "B:" << MB << std::endl;
	std::cout << "C:" << MC << std::endl;
	Matrix::DSYL(MA, MA, MC);
	std::cout << "X:" << MC << std::endl;
	delete[] A;
};
