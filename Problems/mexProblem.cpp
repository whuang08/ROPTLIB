
#include "Problems/mexProblem.h"

#ifdef MATLAB_MEX_FILE

/*Define the namespace*/
namespace ROPTLIB{

	mexProblem::mexProblem(const mxArray *inf, const mxArray *ingf, const mxArray *inHess, const mxArray *inPreCon)
	{
		mxf = inf;
		mxgf = ingf;
		mxHess = inHess;
		mxPreCon = inPreCon;
        NumGradHess = false;
	};

	mexProblem::~mexProblem(void)
	{
	};

	realdp mexProblem::f(const Variable &x) const
	{
		mxArray *Xmx;
		ObtainMxArrayFromElement(Xmx, &x);
		mxArray *lhs[2], *rhs[2];
		rhs[0] = const_cast<mxArray *> (mxf);
		rhs[1] = const_cast<mxArray *> (Xmx);
		mexCallMATLAB(2, lhs, 2, rhs, "feval");
        AddToElementFromMxArray(const_cast<Variable *> (&x), lhs[1]);
		realdp result = mxGetScalar(lhs[0]);
		mxDestroyArray(Xmx);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
		return result;
	};

	Vector &mexProblem::EucGrad(const Variable &x, Vector *result) const
	{
        if(! mxIsClass(mxgf, "function_handle"))
            return Problem::EucGrad(x, result);
        
		mxArray *Xmx;
		ObtainMxArrayFromElement(Xmx, &x);
		mxArray *lhs[2], *rhs[2];
		rhs[0] = const_cast<mxArray *> (mxgf);
		rhs[1] = const_cast<mxArray *> (Xmx);
		mexCallMATLAB(2, lhs, 2, rhs, "feval");
        AddToElementFromMxArray(const_cast<Variable *> (&x), lhs[1]);
		ObtainElementFromMxArray(result, lhs[0]);
		mxDestroyArray(Xmx);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
        return *result;
	};

	Vector &mexProblem::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        if(! mxIsClass(mxHess, "function_handle"))
            return Problem::EucHessianEta(x, etax, result);
        
		mxArray *Xmx, *Etaxmx;
		ObtainMxArrayFromElement(Xmx, &x);
		ObtainMxArrayFromElement(Etaxmx, &etax);
		mxArray *lhs[2], *rhs[3];
		rhs[0] = const_cast<mxArray *> (mxHess);
		rhs[1] = const_cast<mxArray *> (Xmx);
		rhs[2] = const_cast<mxArray *> (Etaxmx);
		mexCallMATLAB(2, lhs, 3, rhs, "feval");
		AddToElementFromMxArray(const_cast<Variable *> (&x), lhs[1]);
		ObtainElementFromMxArray(result, lhs[0]);
		mxDestroyArray(Xmx);
		mxDestroyArray(Etaxmx);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
        return *result;
	};
	
	Vector &mexProblem::PreConditioner(const Variable &x, const Vector &eta, Vector *result) const
	{
        if(! mxIsClass(mxPreCon, "function_handle"))
            return Problem::PreConditioner(x, eta, result);
        
        Vector Exeta(Domain->GetEMPTYEXTR()), Exresult(Domain->GetEMPTYEXTR());
        if(Domain->GetIsIntrinsic())
            Domain->ObtainExtr(x, eta, &Exeta);
        else
            Exeta = eta;
        
		mxArray *Xmx, *eta1mx;
		mexProblem::ObtainMxArrayFromElement(Xmx, &x);
        mexProblem::ObtainMxArrayFromElement(eta1mx, &Exeta);
        
		mxArray *lhs[2], *rhs[3];
		rhs[0] = const_cast<mxArray *> (mxPreCon);
		rhs[1] = const_cast<mxArray *> (Xmx);
		rhs[2] = const_cast<mxArray *> (eta1mx);
		mexCallMATLAB(2, lhs, 3, rhs, "feval");
		AddToElementFromMxArray(const_cast<Variable *> (&x), lhs[1]);
        
        
        ObtainElementFromMxArray(&Exresult, lhs[0]);
        if(Domain->GetIsIntrinsic())
            Domain->ObtainIntr(x, Exresult, result);
        else
            *result = Exresult;
        
		mxDestroyArray(Xmx);
		mxDestroyArray(eta1mx);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
        return *result;
	};

	void mexProblem::ObtainMxArrayFromElement(mxArray *&Xmx, const Element *X)
	{
		integer sizeoftempdata = X->GetSizeofFields();
		integer nfields = sizeoftempdata + 1;
		std::string *fnames = new std::string[nfields];
		X->ObtainTempNames(fnames);
		fnames[sizeoftempdata].assign("main");
		char **names = new char *[nfields];
		for (integer i = 0; i < nfields; i++)
		{
			names[i] = new char[30];
			if (fnames[i].size() >= 30)
			{
				mexErrMsgTxt("The lengths of field names should be less than 30!");
			}
			strcpy(names[i], fnames[i].c_str());
		}

		mxArray *tmp;
		const char *name;
        integer lengthtmp;
		const realdp *Sharedtmpptr;
        integer nsize;
        bool iscomplex;
        const integer *size;
        mwSize mwndim;
        mwSize *mwdims;
		double *tmpptr;
		Xmx = mxCreateStructMatrix(1, 1, nfields, const_cast<const char **> (names));
		for (integer i = 0; i < nfields; i++)
		{
			name = mxGetFieldNameByNumber(Xmx, i);
			if (strcmp(name, "main") != 0) /* if not main */
			{
                Sharedtmpptr = X->Field(name).ObtainReadData();
                nsize = X->Field(name).Getls();
                size = X->Field(name).Getsize();
                iscomplex = X->Field(name).Getiscomplex();
                lengthtmp = X->Field(name).Getlength();
			}
			else
			{
				Sharedtmpptr = X->ObtainReadData();
                nsize = X->Getls();
                size = X->Getsize();
                iscomplex = X->Getiscomplex();
				lengthtmp = X->Getlength();
			}
            
            mwndim = nsize;
            mwdims = new mwSize[mwndim];
            for(integer j = 0; j < mwndim; j++)
                mwdims[j] = size[j];
            if(iscomplex)
                mwdims[0] /= 2;
            
            if(iscomplex)
                tmp = mxCreateNumericArray(mwndim, mwdims, mxDOUBLE_CLASS, mxCOMPLEX);
            else
                tmp = mxCreateNumericArray(mwndim, mwdims, mxDOUBLE_CLASS, mxREAL);
            
            delete[] mwdims;
            
            if(iscomplex)
                tmpptr = (realdp *) mxGetComplexDoubles(tmp);
            else
                tmpptr = mxGetDoubles(tmp);
            
			for (integer i = 0; i < lengthtmp; i++)
				tmpptr[i] = Sharedtmpptr[i];
			mxSetFieldByNumber(Xmx, 0, i, tmp);
		}

		for (integer i = 0; i < nfields; i++)
			delete[] names[i];
		delete[] names;
		delete[] fnames;
	};

	void mexProblem::ObtainElementFromMxArray(Element *X, const mxArray *Xmx)
	{
		/* copy the main data from mxArray to X */
		mxArray *XmxA = GetFieldbyName(Xmx, 0, "main");
		if (XmxA != nullptr)
		{
            realdp *Xmxptr = nullptr;
            if(mxIsComplex(XmxA))
                Xmxptr = (realdp *) mxGetComplexDoubles(XmxA);
            else
                Xmxptr = mxGetDoubles(XmxA);
//			double *Xmxptr = mxGetPr(XmxA);
			integer lengthX = X->Getlength();
			realdp *Xptr = X->ObtainWriteEntireData();
			for (integer i = 0; i < lengthX; i++)
				Xptr[i] = Xmxptr[i];
		}

		/* copy temp data from mxArray to X */
		integer nfields = mxGetNumberOfFields(Xmx);
		const char **fnames;
		mxArray *tmp;

		fnames = static_cast<const char **> (mxCalloc(nfields, sizeof(*fnames)));
		for (integer i = 0; i < nfields; i++)
		{
//            std::cout << "i:" << i << std::endl;//---
			fnames[i] = mxGetFieldNameByNumber(Xmx, i);
			if (strcmp(fnames[i], "main") != 0)
			{
				tmp = GetFieldbyName(Xmx, 0, fnames[i]);
                realdp *tmpptr = nullptr;
                if(mxIsComplex(tmp))
                    tmpptr = (realdp *) mxGetComplexDoubles(tmp);
                else
                    tmpptr = mxGetDoubles(tmp);
                
                mwSize mwndim = mxGetNumberOfDimensions(tmp);
                if(mwndim > 3)
                    mexErrMsgTxt("ROPTLIB does not support a 4-th or higher order tensor!");
                const mwSize *mwdims = mxGetDimensions(tmp);
                integer row = mxGetM(tmp);
                integer col = mxGetN(tmp);
                integer num = (mwndim < 3) ? 1 : mwdims[2];
                Vector Sharedtmp(row, col, num, mxIsComplex(tmp));
                realdp *Sharedtmpptr = Sharedtmp.ObtainWriteEntireData();
                for (integer i = 0; i < Sharedtmp.Getlength(); i++)
                    Sharedtmpptr[i] = tmpptr[i];
				X->AddToFields(fnames[i], Sharedtmp);
			}
		}
	};

	void mexProblem::AddToElementFromMxArray(Element *X, const mxArray *Xmx)
	{
		/* copy temp data--which does not exist in X--from mxArray to X */
		integer nfields = mxGetNumberOfFields(Xmx);
		const char **fnames;
		mxArray *tmp;

		fnames = static_cast<const char **> (mxCalloc(nfields, sizeof(*fnames)));
		for (integer i = 0; i < nfields; i++)
		{
			fnames[i] = mxGetFieldNameByNumber(Xmx, i);
			if (strcmp(fnames[i], "main") != 0 && ! X->FieldsExist(fnames[i]))
			{
                tmp = GetFieldbyName(Xmx, 0, fnames[i]);
                realdp *tmpptr = nullptr;
                if(mxIsComplex(tmp))
                    tmpptr = (realdp *) mxGetComplexDoubles(tmp);
                else
                    tmpptr = mxGetDoubles(tmp);
                
                mwSize mwndim = mxGetNumberOfDimensions(tmp);
                const mwSize *mwdims = mxGetDimensions(tmp);
                if(mwndim > 3)
                    mexErrMsgTxt("ROPTLIB does not support a 4-th or higher order tensor!");
                integer row = mxGetM(tmp);
                integer col = mxGetN(tmp);
                integer num = (mwndim < 3) ? 1 : mwdims[2];
                Vector Sharedtmp(row, col, num, mxIsComplex(tmp));
                realdp *Sharedtmpptr = Sharedtmp.ObtainWriteEntireData();
                for (integer i = 0; i < Sharedtmp.Getlength(); i++)
                    Sharedtmpptr[i] = tmpptr[i];
                X->AddToFields(fnames[i], Sharedtmp);
			}
		}
	};

	mxArray *mexProblem::GetFieldbyName(const mxArray *S, integer idxstruct, const char *name)
	{
		integer nfields = mxGetNumberOfFields(S);
		for (integer i = 0; i < nfields; i++)
		{
			if (strcmp(mxGetFieldNameByNumber(S, i), name) == 0)
			{
				return mxGetFieldByNumber(S, idxstruct, i);
			}
		}
		return nullptr;
	};
}; /*end of ROPTLIB namespace*/

#endif
