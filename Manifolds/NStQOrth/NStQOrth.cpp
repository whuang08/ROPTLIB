
#include "Manifolds/NStQOrth/NStQOrth.h"

/*Define the namespace*/
namespace ROPTLIB{

	NStQOrth::NStQOrth(integer r, integer c)
	{
		n = r;
		p = c;
		IsIntrApproach = true;
		HasHHR = false;
        HasLockCon = false;
		UpdBetaAlone = false;
		name.assign("Noncompact Stiefel manifold quotient orthogonal group");
		IntrinsicDim = r * c - c * (c - 1) / 2;
		ExtrinsicDim = r * c;
		metric = NSOHGZ;
		VecTran = NSOPARALLELIZATION;
		EMPTYEXTR = new NSOVector(r, c);
		EMPTYINTR = new NSOVector(IntrinsicDim);
	};

	// Choose the default parameters
	void NStQOrth::ChooseNSOParamsSet1(void)
	{
		metric = NSOHGZ;
		VecTran = NSOPARALLELIZATION;
	};

	void NStQOrth::ChooseNSOParamsSet2(void)
	{
		metric = NSOEUC;
		VecTran = NSOPARALLELIZATION;
	};

	void NStQOrth::ChooseNSOParamsSet3(void)
	{
		metric = NSOQEUC;
		VecTran = NSOPARALLELIZATION;
	};

	void NStQOrth::ChooseNSOParamsSet4(void)
	{
		metric = NSOHGZ;
		VecTran = NSOPROJECTION;
	};

	void NStQOrth::ChooseNSOParamsSet5(void)
	{
		metric = NSOEUC;
		VecTran = NSOPROJECTION;
	};

	void NStQOrth::ChooseNSOParamsSet6(void)
	{
		metric = NSOQEUC;
		VecTran = NSOPROJECTION;
	};

	NStQOrth::~NStQOrth(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	double NStQOrth::Metric(Variable *x, Vector *etax, Vector *xix) const
	{
		//Vector *exetax = EMPTYEXTR->ConstructEmpty();
		//Vector *exxix = EMPTYEXTR->ConstructEmpty();
		//ObtainExtr(x, etax, exetax);
		//ObtainExtr(x, xix, exxix);
		//double result = Manifold::Metric(x, exetax, exxix);
		//delete exetax;
		//delete exxix;
		if (IsIntrApproach)
			return Manifold::Metric(x, etax, xix);
		printf("Warning: Metric for extrinsic representation has not been done!\n");
		return 0;
	};

	void NStQOrth::CheckParams(void) const
	{
		std::string NSOMetricnames[NSOMETRICLENGTH] = { "NSOEUC", "NSOQEUC", "NSOHGZ" };
		std::string NSOVTnames[NSOVECTORTRANSPORTLENGTH] = { "NSOPARALLELIZATION", "NSOPROJECTION" };
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		if (p == 1)
		{
			printf("n           :%15d,\t", n);
			printf("metric      :%15s\n", NSOMetricnames[metric].c_str());
			printf("VecTran     :%15s\n", NSOVTnames[VecTran].c_str());
		}
		else
		{
			printf("n           :%15d,\t", n);
			printf("p           :%15d\n", p);
			printf("metric      :%15s\t", NSOMetricnames[metric].c_str());
			printf("VecTran     :%15s\n", NSOVTnames[VecTran].c_str());
		}
	};

	void NStQOrth::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (metric == NSOEUC)
			EucGradToGradEUC(x, egf, gf, prob);
		else
			if (metric == NSOHGZ)
				EucGradToGradHGZ(x, egf, gf, prob);
			else
			if (metric == NSOQEUC)
				EucGradToGradQEUC(x, egf, gf, prob);
			else
				printf("warning: NStQOrth::EucGradToGrad has not been done for this metric!\n");
	};

	void NStQOrth::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		if (metric == NSOEUC)
			EucHvToHvEUC(x, etax, exix, xix, prob);
		else
			if (metric == NSOHGZ)
				EucHvToHvHGZ(x, etax, exix, xix, prob);
			else
				if (metric == NSOQEUC)
					EucHvToHvQEUC(x, etax, exix, xix, prob);
			else
				printf("warning: NStQOrth::EucHvToHv has not been done for this metric!\n");
	};

	void NStQOrth::Projection(Variable *x, Vector *v, Vector *result) const
	{
		if (IsIntrApproach)
			IntrProjection(x, v, result);
		else
			ExtrProjection(x, v, result);
	};

	void NStQOrth::Retraction(Variable *x, Vector *etax, Variable *result, double stepsize) const
	{
		if (IsIntrApproach)
		{
			Vector *exetax = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(x, etax, exetax);
			Manifold::Retraction(x, exetax, result, 1);
			delete exetax;
		}
		else
		{
			Manifold::Retraction(x, etax, result, 1);
		}
	};

	void NStQOrth::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		printf("warning:NStQOrth::coTangentVector has not been done!\n");
		xiy->CopyTo(result);
	};

	void NStQOrth::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsIntrApproach)
		{
			Vector *exxix = EMPTYEXTR->ConstructEmpty();
			Vector *exresult = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(x, xix, exxix);
			ExtrProjection(y, exxix, exresult);
			ObtainIntr(y, exresult, result);
			delete exxix;
			delete exresult;
		}
		else
		{
			ExtrProjection(y, xix, result);
		}
	};

	double NStQOrth::Beta(Variable *x, Vector *etax) const
	{
		return 1;
	};


	void NStQOrth::ComputeHHR(Variable *x)
	{
		const double *xM = x->ObtainReadData();
		SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
		double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
		SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1]);
		double *tau = HHRTau->ObtainWriteEntireData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1;
		// ptrHHR <- xM, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		integer Length = N * P;
		dcopy_(&Length, const_cast<double *> (xM), &inc, ptrHHR, &inc);
		integer *jpvt = new integer[P];
		integer info;
		integer lwork = -1;
		double lworkopt;
		// compute the size of space required in the dgeqp3
#ifndef MATLAB_MEX_FILE
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
#else
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		for (integer i = 0; i < P; i++)
			jpvt[i] = i + 1;
		// QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
		// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
#ifndef MATLAB_MEX_FILE
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
#else
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
#endif

		x->AddToTempData("HHR", HouseHolderResult);
		x->AddToTempData("HHRTau", HHRTau);
		if (info < 0)
			printf("Error in qr decomposition in NStQOrth::ComputeHHR!\n");
		for (integer i = 0; i < P; i++)
		{
			if (jpvt[i] != (i + 1))
				printf("Error in qf retraction in NStQOrth::ComputeHHR!\n");
		}
		delete[] jpvt;
		delete[] work;
	};

	void NStQOrth::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		if (metric == NSOEUC)
			ObtainIntrEUC(x, etax, result);
		else
			if (metric == NSOHGZ)
				ObtainIntrHGZ(x, etax, result);
			else
				if (metric == NSOQEUC)
					ObtainIntrQEUC(x, etax, result);
			else
				printf("warning: NStQOrth::ObtainIntr has not been done for this metric!\n");
	};

	void NStQOrth::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		if (metric == NSOEUC)
			ObtainExtrEUC(x, intretax, result);
		else
			if (metric == NSOHGZ)
				ObtainExtrHGZ(x, intretax, result);
			else
				if (metric == NSOQEUC)
					ObtainExtrQEUC(x, intretax, result);
			else
				printf("warning: NStQOrth::ObtainExtr has not been done for this metric!\n");
	};

	void NStQOrth::IntrProjection(Variable *x, Vector *v, Vector *result) const
	{
		v->CopyTo(result);
	};

	void NStQOrth::ExtrProjection(Variable *x, Vector *v, Vector *result) const
	{
		if (metric == NSOEUC || metric == NSOHGZ)
			ExtrProjectionEUCorHGZ(x, v, result);
		else
			if (metric == NSOQEUC)
				ExtrProjectionQEUC(x, v, result);
			else
				printf("Warning: NStQOrth::ExtrProjection has not been done for this Riemannian metric!\n");
	};

	void NStQOrth::ExtrProjectionQEUC(Variable *x, Vector *v, Vector *result) const
	{
		const double *xM = x->ObtainReadData();
		const double *vptr = v->ObtainReadData();

		double *XTX = new double[p * p * 2];
		double *XTV = XTX + p * p;

		// XTX <- XM^T XM
		Matrix MxM(xM, n, p), MXTX(XTX, p, p), MV(vptr, n, p), MXTV(XTV, p, p);
		Matrix::DGEMM(GLOBAL::DONE, MxM, true, MxM, false, GLOBAL::DZERO, MXTX);
		// XTV <- XM^T V
		Matrix::DGEMM(GLOBAL::DONE, MxM, true, MV, false, GLOBAL::DZERO, MXTV);

		// XTV <- XTV - XTV^T
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				XTV[i + j * p] -= XTV[j + i * p];
				XTV[j + i * p] = -XTV[i + j * p];
			}
		}
		v->CopyTo(result);
		double *resultTV = result->ObtainWritePartialData();
		Matrix Mresult(resultTV, n, p);
		Matrix::DSYL(MXTX, MXTX, MXTV);
		Matrix::DGEMM(GLOBAL::DNONE, MxM, false, MXTV, false,
			GLOBAL::DONE, Mresult);

		delete[] XTX;
	};

	void NStQOrth::ExtrProjectionEUCorHGZ(Variable *x, Vector *v, Vector *result) const
	{
		const double *xM = x->ObtainReadData();
		const double *V = v->ObtainReadData();

		double *XTX = new double[2 * p * p];
		double *XTV = XTX + p * p;
		Matrix MxM(xM, n, p), MXTX(XTX, p, p), MV(V, n, p), MXTV(XTV, p, p);
		// XHX <- XM^T XM
		Matrix::DGEMM(GLOBAL::DONE, MxM, true, MxM, false, GLOBAL::DZERO, MXTX);
		// XHV <- XM^H V
		Matrix::DGEMM(GLOBAL::DONE, MxM, true, MV, false, GLOBAL::DZERO, MXTV);

		integer N = n, P = p, info;
#ifndef MATLAB_MEX_FILE
		// solve for (XTX)^{-1} XTV
		dpotrf_(GLOBAL::L, &P, XTX, &P, &info);
		dpotrs_(GLOBAL::L, &P, &P, XTX, &P, XTV, &P, &info);
#else
		// solve for (XTX)^{-1} XTV
		dpotrf_(GLOBAL::L, &P, (double *)XTX, &P, &info);
		dpotrs_(GLOBAL::L, &P, &P, (double *)XTX, &P, (double *)XTV, &P, &info);
#endif

		if (info != 0)
		{
			printf("warning: dpotrs failed in NStQOrth::ExtrProjection with info:%d!\n", info);
	}

		// XTV <- (XTV - XTV^T)/2
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				XTV[i + j * p] -= XTV[j + i * p];
				XTV[j + i * p] = -XTV[i + j * p];
			}
		}
		for (integer i = 0; i < p * p; i++)
		{
			XTV[i] /= 2;
		}
		v->CopyTo(result);
		double *resultTV = result->ObtainWritePartialData();
		Matrix MresultTV(resultTV, n, p);
		Matrix::DGEMM(GLOBAL::DNONE, MxM, false, MXTV, false, GLOBAL::DONE, MresultTV);

		delete[] XTX;
	};

	void NStQOrth::EucGradToGradEUC(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}

		egf->CopyTo(gf);
		double *gfv = gf->ObtainWriteEntireData();

		//Y / (x'*x) use compute HHR in NstQorth
		if (!x->TempDataExist("HHR"))
			ComputeHHR(x);
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const double *ptrHHR = HHR->ObtainReadData();

		double *YT = new double[n * p];

		/*YT = Y^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				YT[j + i * p] = gfv[i + j * n];
			}
		}
		integer info, N = n, P = p;
#ifndef MATLAB_MEX_FILE
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
#else
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
#endif
		/*YT = Y^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				gfv[i + j * n] = YT[j + i * p];
			}
		}
		delete[] YT;
		double half = 0.5;
		integer length = n * p;
		dscal_(&length, &half, gfv, &GLOBAL::IONE);

		double *tempspace = new double[p * p];

		const double *Y = x->ObtainReadData();

		/*tempspace2 = Y^T * gfv */
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (Y), &N, gfv, &N, &GLOBAL::DZERO, tempspace, &P);

		/*compute (Y^T Y)^{-1} tempspace2. It is also R^{-1} R^{-T} tempspace2 */
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, tempspace, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, tempspace, &P, &info);

		double nhalf = -0.5;
		/*final result: gfv = (I - 1/2 * Y * (Y^T Y)^{-1}Y^T) gfv */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &nhalf, const_cast<double *> (Y), &N, tempspace, &P, &GLOBAL::DONE, gfv, &N);

		delete[] tempspace;
		//ExtrProjection(x, gf, gf); Not necessary since gf is already in the horizontal space.
	};

	void NStQOrth::EucHvToHvEUC(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		const double *Y = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();

		const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
		Vector *egfx = Sharedegf->GetSharedElement();
		const double *egfTV = egfx->ObtainReadData();

		if (!x->TempDataExist("HHR"))
			ComputeHHR(x);
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const double *ptrHHR = HHR->ObtainReadData();

		double *tmp = new double[n * p * 2 + p * p * 4];
		double *YYinvYTEH = tmp + n * p;
		double *EGYYinv = YYinvYTEH + p * p;
		double *YTeta = EGYYinv + n * p;
		double *EtaTEGYYinv = YTeta + p * p;
		double *YYinvYTEGYYinv = EtaTEGYYinv + p * p;
		/*tmp = Egrad^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				tmp[j + i * p] = egfTV[i + j * n];
			}
		}
		integer info, N = n, P = p;
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);
		/*EGYYinv = Egrad * (Y^T * Y)^{-1}*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				EGYYinv[i + j * n] = tmp[j + i * p];
			}
		}

		/*YYinvYTEH = (Y^T * Y)^{-1} * Y^T * EHess */
		const double *exixTV = exix->ObtainReadData();
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (Y), &N, const_cast<double *> (exixTV), &N, &GLOBAL::DZERO, YYinvYTEH, &P);
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, YYinvYTEH, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, YYinvYTEH, &P, &info);

		/*YYinvYTEGYYinv = (Y^T * Y)^{-1} * Y^T * Egrad * (Y^T * Y)^{-1}*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (Y), &N, EGYYinv, &N, &GLOBAL::DZERO, YYinvYTEGYYinv, &P);
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, YYinvYTEGYYinv, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, YYinvYTEGYYinv, &P, &info);

		/*EtaTEGYYinv = eta^T * Egrad * (Y^T * Y)^{-1}*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (etaxTV), &N, EGYYinv, &N, &GLOBAL::DZERO, EtaTEGYYinv, &P);

		/*YTeta = Y^T * eta */
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (Y), &N, const_cast<double *> (etaxTV), &N, &GLOBAL::DZERO, YTeta, &P);

		exix->CopyTo(xix);
		double *xixTV = xix->ObtainWritePartialData();
		double half = 0.5, nhalf = -0.5;

		/*xix = EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Y^T * EucHess */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &nhalf, const_cast<double *> (Y), &N, YYinvYTEH, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Y^T * EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Egrad^T * eta */
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &P, &P, &nhalf, const_cast<double *>(Y), &N, EtaTEGYYinv, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Y^T * EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Egrad^T * eta - Egrad * (Y^T * Y)^{-1} * Y * eta */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DNONE, const_cast<double *>(EGYYinv), &N, YTeta, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Y^T * EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Egrad^T * eta - Egrad * (Y^T * Y)^{-1} * Y * eta 
		         + Y * (Y^T * Y)^{-1} * Y^T * Egrad * (Y^T * Y)^{-1} * Y^T * eta */
		dgemm_(GLOBAL::N, GLOBAL::N, &P, &P, &P, &GLOBAL::DONE, YYinvYTEGYYinv, &P, YTeta, &P, &GLOBAL::DZERO, tmp, &P);
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DONE, const_cast<double *>(Y), &N, tmp, &P, &GLOBAL::DONE, xixTV, &N);

		integer length = n * p;
		dscal_(&length, &half, xixTV, &GLOBAL::IONE);

		/*tmp = xix^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				tmp[j + i * p] = xixTV[i + j * n];
			}
		}
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);

		/*xix = xix * (Y^T * Y)^{-1}*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				xixTV[i + j * n] = tmp[j + i * p];
			}
		}

		delete[] tmp;

		ExtrProjection(x, xix, xix);
	};

	void NStQOrth::ObtainIntrEUC(Variable *x, Vector *etax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		double *resultTV = result->ObtainWriteEntireData();
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		integer info;
		integer lwork = -1;
		double lworkopt;
		double *tempspace = new double[n * p];
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, &lworkopt, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// tempspace <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (etaxTV), &inc, tempspace, &inc);
		// tempspace <- Q^T * tempspace, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, work, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, work, &lwork, &info);
#endif
		double sign = 0;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
#ifndef MATLAB_MEX_FILE
			dscal_(&P, &sign, tempspace + i, &N);
#else
			dscal_(&P, &sign, tempspace + i, &N);
#endif
		}

		double *L = new double[p * p];
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				L[j + i * p] = 0;
			}
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			for (integer j = i; j < p; j++)
			{
				L[j + i * p] = ptrHHR[i + n * j] * sign;
			}
		}

		double *tempspaceL = new double[n * p];

		Matrix MtL(tempspaceL, n, p), ML(L, p, p), Mtempspace(tempspace, n, p);
		Matrix::DGEMM(GLOBAL::DONE, Mtempspace, false, ML, false, GLOBAL::DZERO, MtL);

		delete[] L;
		delete[] tempspace;

		/*Matrix MxM(xM, n, p), MetaxTV(etaxTV, n, p), Mtempspace((double *)tempspace, p, p, n);
		Matrix::CGEMM(GLOBAL::ZONE, MxM, true, MetaxTV, false, GLOBAL::ZZERO, Mtempspace);*/
		double r2 = sqrt(2.0);
		double factor = 1;//-- sqrt(Manifold::Metric(x, x, x));
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			resultTV[idx] = 2 * tempspaceL[i + i * n] / factor;
			idx++;
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[idx] = 2 * r2 * tempspaceL[j + i * n] / factor;
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[idx] = r2 * tempspaceL[j + i * n];
				idx++;
			}
		}
		delete[] work;
		delete[] tempspaceL;
	};

	void NStQOrth::ObtainExtrEUC(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		integer info;
		integer idx = 0;
		//		doublecomplex *S = new doublecomplex[p * p];
		double r2 = sqrt(2.0);
		double factor = 1;//-- sqrt(Manifold::Metric(x, x, x));
		for (integer i = 0; i < p; i++)
		{
			resultTV[i + i * n] = intretaxTV[idx] * factor / 2;
			idx++;
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx] / r2 / 2 * factor;
				resultTV[i + j * n] = intretaxTV[idx] / r2 / 2 * factor;
				idx++;
			}
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx] / r2;
				idx++;
			}
		}

		double sign = 0;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			// result(i, :) <- sign * result(i, :), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
#ifndef MATLAB_MEX_FILE
			dscal_(&P, &sign, resultTV + i, &N);
#else
			dscal_(&P, &sign, resultTV + i, &N);
#endif
		}
		integer lwork = -1;
		double lworkopt;
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, &lworkopt, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt);

		double *work = new double[lwork];
		// resultTV <- Q * resultTV, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, work, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, work, &lwork, &info);
#endif
		delete[] work;

		double *L = new double[p * p + n * p];
		double *r_T = L + p * p;
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				L[j + i * p] = 0;
			}
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			for (integer j = i; j < p; j++)
			{
				L[j + i * p] = ptrHHR[i + n * j] * sign;
			}
		}

		/*r_T <-  resultTV transpose*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				r_T[j + i * p] = resultTV[i + j * n];
			}
		}
#ifndef MATLAB_MEX_FILE
		/*solve linear system L^H M = r_T, the solution M is stored in r_T*/
		dtrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &P, &N, L, &P, r_T, &P, &info);
#else
		dtrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &P, &N, L, &P, r_T, &P, &info);
#endif
		/*resultTV <-  r_T transpose*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				resultTV[i + j * n] = r_T[j + i * p];
			}
		}
		delete[] L;
	};

	void NStQOrth::EucGradToGradQEUC(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		ExtrProjection(x, egf, gf);
	};

	void NStQOrth::EucHvToHvQEUC(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		ExtrProjection(x, exix, xix);
	};

	void NStQOrth::ObtainIntrQEUC(Variable *x, Vector *etax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		double *resultTV = result->ObtainWriteEntireData();
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		integer info;
		integer lwork = -1;
		double lworkopt;
		double *tempspace = new double[n * p];
		// compute the size of space required in the dormqr
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// tempspace <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (etaxTV), &inc, tempspace, &inc);
		// tempspace <- Q^T * tempspace, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, work, &lwork, &info);
		delete[] work;

		double sign = 0;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			dscal_(&P, &sign, tempspace + i, &N);
		}

		double *L = new double[p * p];
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				L[j + i * p] = 0;
			}
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			for (integer j = i; j < p; j++)
			{
				L[j + i * p] = ptrHHR[i + n * j] * sign;
			}
		}

		/* Ltempspace = L* tempspace(1:p, 1:p) */

		double *Ltempspace = new double[p * p];
		dgemm_(GLOBAL::N, GLOBAL::N, &P, &P, &P, &GLOBAL::DONE, L, &P, tempspace, &N, &GLOBAL::DZERO, Ltempspace, &P);
		for (integer i = 0; i < P; i++)
		{
			for (integer j = 0; j < P; j++)
				tempspace[j + i * n] = Ltempspace[j + i * p];
		}
		delete[] Ltempspace;

		double *resultptr = result->ObtainWriteEntireData();
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				resultptr[idx] = tempspace[j + i * n];
				idx++;
			}
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultptr[idx] = tempspace[j + i * n];
				idx++;
			}
		}
		delete[] tempspace;
		/* Above: result is then R^{-T} S + K. Next, we first find the coefficients of the etax under a non-orthonormal basis.
		Then orthonormalize the basis and update the coefficient. To the end, we find the non-orthogonal basis and store it in B.
		Then use qr decomposition to find the matrix R. Then use the matrix R to find the coefficient vector under the orthonormalized
		basis. total complexity is O(p^6) + O(n p)*/


		/* cholesky decomposition */
		double *B = new double[p * p * p * (p + 1) / 2];
		idx = 0;
		for (integer i = 0; i < p * p * p * (p + 1) / 2; i++)
			B[i] = 0;

		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				B[j + i * p + idx * p * p] = 1;
				B[i + j * p + idx * p * p] = 1;
				dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &P, &P, L, &P, B + idx * p * p, &P, &info);
				idx++;
			}
		}
		delete[] L;

		integer *jpvt = new integer[p * (p + 1) / 2];
		double *tau = new double[p * (p + 1) / 2];
		lwork = -1;
		integer dim = p * (p + 1) / 2;
		integer Psquare = p * p;
		lworkopt;
		// compute the size of space required in the dgeqp3
		dgeqp3_(&Psquare, &dim, B, &Psquare, jpvt, tau, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		work = new double[lwork];
		for (integer i = 0; i < p * (p + 1) / 2; i++)
			jpvt[i] = i + 1;
		// QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
		// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
		dgeqp3_(&Psquare, &dim, B, &Psquare, jpvt, tau, work, &lwork, &info);
		delete tau;
		delete[] work;
		delete[]jpvt;


		SharedSpace *SharedBL = new SharedSpace(2, dim, dim);
		double *BL = SharedBL->ObtainWriteEntireData();
		for (integer i = 0; i < dim; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				BL[j + i * dim] = 0;
			}
			sign = (B[i + p * p * i] >= 0) ? 1 : -1;
			for (integer j = i; j < dim; j++)
			{
				BL[j + i * dim] = B[i + p * p * j] * sign;
			}
		}

		delete[] B;

		tempspace = new double[dim];
		dgemm_(GLOBAL::T, GLOBAL::N, &dim, &GLOBAL::IONE, &dim, &GLOBAL::DONE, BL, &dim, resultptr, &dim, &GLOBAL::DZERO, tempspace, &dim);
		dcopy_(&dim, tempspace, &GLOBAL::IONE, resultptr, &GLOBAL::IONE);
		delete[] tempspace;
		x->AddToTempData("BL", SharedBL);
	};

	void NStQOrth::ObtainExtrQEUC(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		const SharedSpace *BL = x->ObtainReadTempData("BL");
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();
		const double *ptrBL = BL->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		integer dim = p * (p + 1) / 2, info;
		double *tempspace = new double[dim];
		dcopy_(&dim, const_cast<double *> (intretaxTV), &GLOBAL::IONE, tempspace, &GLOBAL::IONE);
		dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &dim, &GLOBAL::IONE, const_cast<double *>(ptrBL), &dim, tempspace, &dim, &info);

		integer N = n, P = p, Length = N * P;
		integer idx = 0;

		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				resultTV[j + i * n] = tempspace[idx];
				resultTV[i + j * n] = resultTV[j + i * n];
				idx++;
			}
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx];
				idx++;
			}
		}
		delete[] tempspace;
		double *L = new double[p * p];
		double sign = 0;
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				L[j + i * p] = 0;
			}
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			for (integer j = i; j < p; j++)
			{
				L[j + i * p] = ptrHHR[i + n * j] * sign;
			}
		}

		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &P, &P, L, &P, resultTV, &N, &info);

		delete[] L;

		/*=========*/
		sign = 0;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			// result(i, :) <- sign * result(i, :), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
#ifndef MATLAB_MEX_FILE
			dscal_(&P, &sign, resultTV + i, &N);
#else
			dscal_(&P, &sign, resultTV + i, &N);
#endif
		}
		integer lwork = -1;
		double lworkopt;
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, &lworkopt, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt);

		double *work = new double[lwork];
		// resultTV <- Q * resultTV, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, work, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, work, &lwork, &info);
#endif
		delete[] work;
	};

	void NStQOrth::EucGradToGradHGZ(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		egf->CopyTo(gf);
		double *gfv = gf->ObtainWriteEntireData();

		//Y / (x'*x) use compute HHR in NstQorth
		if (!x->TempDataExist("HHR"))
			ComputeHHR(x);
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const double *ptrHHR = HHR->ObtainReadData();

		double *YT = new double[n * p];

		/*YT = Y^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				YT[j + i * p] = gfv[i + j * n];
			}
		}
		integer info, N = n, P = p;
#ifndef MATLAB_MEX_FILE
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
#else
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, YT, &P, &info);
#endif
		/*YT = Y^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				gfv[i + j * n] = YT[j + i * p];
			}
		}
		delete[] YT;

		//ExtrProjection(x, gf, gf); Not necessary since gf is already in the horizontal space.
	};

	void NStQOrth::EucHvToHvHGZ(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		const double *Y = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();

		const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
		Vector *egfx = Sharedegf->GetSharedElement();
		const double *egfTV = egfx->ObtainReadData();

		if (!x->TempDataExist("HHR"))
			ComputeHHR(x);
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const double *ptrHHR = HHR->ObtainReadData();

		double *tmp = new double[n * p * 2 + p * p * 3];
		double *EGYYinv = tmp + n * p;
		double *YTEGYYinv = EGYYinv + n * p;
		double *EtaTEGYYinv = YTEGYYinv + p * p;
		double *EtaTY = EtaTEGYYinv + p * p;
		/*tmp = Egrad^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				tmp[j + i * p] = egfTV[i + j * n];
			}
		}
		integer info, N = n, P = p;
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);

		/*EGYYinv = Egrad * (Y^T * Y)^{-1}*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				EGYYinv[i + j * n] = tmp[j + i * p];
			}
		}

		/*YTEGYYinv = Y^T * Egrad * (Y^T * Y)^{-1}*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (Y), &N, EGYYinv, &N, &GLOBAL::DZERO, YTEGYYinv, &P);

		/*EtaTEGYYinv = eta^T * Egrad * (Y^T * Y)^{-1}*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (etaxTV), &N, EGYYinv, &N, &GLOBAL::DZERO, EtaTEGYYinv, &P);

		/*EtaTY = eta^T * Y */
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (etaxTV), &N, const_cast<double *> (Y), &N, &GLOBAL::DZERO, EtaTY, &P);

		exix->CopyTo(xix);
		double *xixTV = xix->ObtainWritePartialData();
		double half = 0.5, nhalf = -0.5;

		/*xix = EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta */
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &P, &P, &nhalf, EGYYinv, &N, EtaTY, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta - 0.5 * Y * (Y^T * Y)^{-1} Egrad^T * eta */
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &P, &P, &nhalf, const_cast<double *>(Y), &N, EtaTEGYYinv, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta - 0.5 * Y * (Y^T * Y)^{-1} Egrad^T * eta + 0.5 * eta * Y^T * Egrad * (Y^T * Y)^{-1} */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &half, const_cast<double *>(etaxTV), &N, YTEGYYinv, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta - 0.5 * Y * (Y^T * Y)^{-1} Egrad^T * eta + 0.5 * eta * Y^T * Egrad * (Y^T * Y)^{-1} 
		         - 0.5 * Y * eta^T * Egrad * (Y^T * Y)^{-1} */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &nhalf, const_cast<double *>(Y), &N, EtaTEGYYinv, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta - 0.5 * Y * (Y^T * Y)^{-1} Egrad^T * eta + 0.5 * eta * Y^T * Egrad * (Y^T * Y)^{-1}
				 - 0.5 * Y * eta^T * Egrad * (Y^T * Y)^{-1} + 0.5 * eta *(Y^T * Y)^{-1} Egrad^T Y */
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &P, &P, &half, const_cast<double *>(etaxTV), &N, YTEGYYinv, &P, &GLOBAL::DONE, xixTV, &N);

		/*xix =  EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta - 0.5 * Y * (Y^T * Y)^{-1} Egrad^T * eta + 0.5 * eta * Y^T * Egrad * (Y^T * Y)^{-1}
				 - 0.5 * Y * eta^T * Egrad * (Y^T * Y)^{-1} + 0.5 * eta *(Y^T * Y)^{-1} Egrad^T Y - 0.5 * Egrad * (Y^T * Y)^{-1} * eta^T * Y */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &nhalf, EGYYinv, &N, EtaTY, &P, &GLOBAL::DONE, xixTV, &N);
		
		/*tmp = xix^T*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				tmp[j + i * p] = xixTV[i + j * n];
			}
		}
		dtrtrs_(GLOBAL::U, GLOBAL::T, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);
		dtrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &N, const_cast<double *> (ptrHHR), &N, tmp, &P, &info);

		/*xix = xix * (Y^T * Y)^{-1}*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				xixTV[i + j * n] = tmp[j + i * p];
			}
		}

		delete[] tmp;
		
		ExtrProjection(x, xix, xix);
	};

	void NStQOrth::ObtainIntrHGZ(Variable *x, Vector *etax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		double *resultTV = result->ObtainWriteEntireData();
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		integer info;
		integer lwork = -1;
		double lworkopt;
		double *tempspace = new double[n * p];
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, &lworkopt, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// tempspace <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (etaxTV), &inc, tempspace, &inc);
		// tempspace <- Q^T * tempspace, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, work, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::T, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), tempspace, &N, work, &lwork, &info);
#endif
		double sign = 0;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
#ifndef MATLAB_MEX_FILE
			dscal_(&P, &sign, tempspace + i, &N);
#else
			dscal_(&P, &sign, tempspace + i, &N);
#endif
		}

		double *L = new double[p * p];
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				L[j + i * p] = 0;
			}
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			for (integer j = i; j < p; j++)
			{
				L[j + i * p] = ptrHHR[i + n * j] * sign;
			}
		}

		double *tempspaceL = new double[n * p];

		Matrix MtL(tempspaceL, n, p), ML(L, p, p), Mtempspace(tempspace, n, p);
		Matrix::DGEMM(GLOBAL::DONE, Mtempspace, false, ML, false, GLOBAL::DZERO, MtL);

		delete[] L;
		delete[] tempspace;

		/*Matrix MxM(xM, n, p), MetaxTV(etaxTV, n, p), Mtempspace((double *)tempspace, p, p, n);
		Matrix::CGEMM(GLOBAL::ZONE, MxM, true, MetaxTV, false, GLOBAL::ZZERO, Mtempspace);*/
		double r2 = sqrt(2.0);
		double factor = 1;//-- sqrt(Manifold::Metric(x, x, x));
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			resultTV[idx] = tempspaceL[i + i * n] / factor;
			idx++;
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[idx] = r2 * tempspaceL[j + i * n] / factor;
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[idx] = tempspaceL[j + i * n];
				idx++;
			}
		}
		delete[] work;
		delete[] tempspaceL;
	};

	void NStQOrth::ObtainExtrHGZ(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		integer info;
		integer idx = 0;
		//		doublecomplex *S = new doublecomplex[p * p];
		double r2 = sqrt(2.0);
		double factor = 1;//-- sqrt(Manifold::Metric(x, x, x));
		for (integer i = 0; i < p; i++)
		{
			resultTV[i + i * n] = intretaxTV[idx] * factor;
			idx++;
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx] / r2 * factor;
				resultTV[i + j * n] = intretaxTV[idx] / r2 * factor;
				idx++;
			}
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx];
				idx++;
			}
		}

		double sign = 0;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			// result(i, :) <- sign * result(i, :), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
#ifndef MATLAB_MEX_FILE
			dscal_(&P, &sign, resultTV + i, &N);
#else
			dscal_(&P, &sign, resultTV + i, &N);
#endif
		}
		integer lwork = -1;
		double lworkopt;
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, &lworkopt, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt);

		double *work = new double[lwork];
		// resultTV <- Q * resultTV, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, work, &lwork, &info);
#else
		dormqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, work, &lwork, &info);
#endif
		delete[] work;

		double *L = new double[p * p + n * p];
		double *r_T = L + p * p;
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				L[j + i * p] = 0;
			}
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			for (integer j = i; j < p; j++)
			{
				L[j + i * p] = ptrHHR[i + n * j] * sign;
			}
		}

		/*r_T <-  resultTV transpose*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				r_T[j + i * p] = resultTV[i + j * n];
			}
		}
#ifndef MATLAB_MEX_FILE
		/*solve linear system L^H M = r_T, the solution M is stored in r_T*/
		dtrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &P, &N, L, &P, r_T, &P, &info);
#else
		dtrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &P, &N, L, &P, r_T, &P, &info);
#endif
		/*resultTV <-  r_T transpose*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				resultTV[i + j * n] = r_T[j + i * p];
			}
		}
		delete[] L;
	};

	void NStQOrth::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		if (VecTran == NSOPARALLELIZATION && !HasHHR)
		{
			return Manifold::VectorTransport(x, etax, y, xix, result);
		}

		if (VecTran == NSOPROJECTION && !HasHHR)
		{
			if (IsIntrApproach)
			{
				VectorTransportProj(x, etax, y, xix, result);
			}
			else
				printf("Warning: Vector transport by projection has not been done for extrinsic approach!\n");
			return;
		}

		if (HasHHR)
			return LCVectorTransport(x, etax, y, xix, result);

		printf("Error: VectorTransport has not been done!\n");
	};

	void NStQOrth::VectorTransportProj(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		Vector *exxix = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, xix, exxix);
		ObtainIntr(y, exxix, result);
		delete exxix;
	};

}; /*end of ROPTLIB namespace*/
