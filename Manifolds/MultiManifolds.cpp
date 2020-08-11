#include "Manifolds/MultiManifolds.h"

/*Define the namespace*/
namespace ROPTLIB{

	MultiManifolds::MultiManifolds(integer numberofmanifolds, ...)
	{
		numoftypes = numberofmanifolds;
		powsinterval = new integer[numoftypes + 1];
		manifolds = new Manifold *[numoftypes];
		va_list argptr;
		va_start(argptr, numberofmanifolds);
		powsinterval[0] = 0;
		for (integer i = 0; i < numoftypes; i++)
		{
			manifolds[i] = va_arg(argptr, Manifold *);
			powsinterval[i + 1] = powsinterval[i] + va_arg(argptr, integer);
		}
		va_end(argptr);

		HasHHR = false;
		numoftotalmani = 0;
		ExtrinsicDim = 0;
		IntrinsicDim = 0;
		for (integer i = 0; i < numoftypes; i++)
		{
			ExtrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetExtrDim();
			IntrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetIntrDim();
			numoftotalmani += (powsinterval[i + 1] - powsinterval[i]);
		}
		name.assign("Multi-manifolds");

        /*If IsIntrApproach is true, then individual manifolds use their own IsIntrApproach parameters to
         specify the manifolds. Otherwise, extrinsic representations are used.*/
		IsIntrApproach = true;

		Element *elements = new Element [numoftypes];
		for (integer i = 0; i < numoftypes; i++)
		{
			if (manifolds[i]->GetIsIntrinsic())
			{
                elements[i] = manifolds[i]->GetEMPTYINTR();
			}
			else
			{
                elements[i] = manifolds[i]->GetEMPTYEXTR();
			}
		}
        
        EMPTYINTR = Element(numoftypes, elements, powsinterval);
        
		for (integer i = 0; i < numoftypes; i++)
		{
            elements[i] = manifolds[i]->GetEMPTYEXTR();
		}
        EMPTYEXTR = Element(numoftypes, elements, powsinterval);
		delete[] elements;
	};

    MultiManifolds::MultiManifolds(Manifold **inmanifolds, integer innumoftypes, integer *inpowsinterval)
    {
        numoftypes = innumoftypes;
        numoftotalmani = inpowsinterval[innumoftypes];
        powsinterval = new integer[numoftypes + 1];
        manifolds = new Manifold *[numoftypes];
        powsinterval[0] = 0;
        for (integer i = 0; i < numoftypes; i++)
        {
            manifolds[i] = inmanifolds[i];
            powsinterval[i + 1] = inpowsinterval[i + 1];
        }

        HasHHR = false;
        ExtrinsicDim = 0;
        IntrinsicDim = 0;
        for (integer i = 0; i < numoftypes; i++)
        {
            ExtrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetExtrDim();
            IntrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetIntrDim();
        }
        name.assign("Multi-manifolds");
        /*If IsIntrApproach is true, then individual manifolds use their own IsIntrApproach parameters to
         specify the manifolds. Otherwise, extrinsic representations are used.*/
        IsIntrApproach = true;
        
        Element *elements = new Element [numoftypes];
        for (integer i = 0; i < numoftypes; i++)
        {
            if (manifolds[i]->GetIsIntrinsic())
            {
                elements[i] = manifolds[i]->GetEMPTYINTR();
            }
            else
            {
                elements[i] = manifolds[i]->GetEMPTYEXTR();
            }
        }
        
        EMPTYINTR = Element(numoftypes, elements, powsinterval);
        
        for (integer i = 0; i < numoftypes; i++)
        {
            elements[i] = manifolds[i]->GetEMPTYEXTR();
        }
        EMPTYEXTR = Element(numoftypes, elements, powsinterval);
        delete[] elements;
    };

	MultiManifolds::~MultiManifolds(void)
	{
		delete[] manifolds;
		delete[] powsinterval;
	};

    Variable MultiManifolds::RandominManifold(void) const
    {
        Variable result(EMPTYEXTR);
        result.NewMemoryOnWrite();
        for(integer i = 0; i < numoftypes; i++)
        {
            for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                result.GetElement(j) = manifolds[i]->RandominManifold();
            }
        }
        return result;
    };

	realdp MultiManifolds::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
	{
        realdp result = 0;
        for(integer i = 0; i < numoftypes; i++)
        {
            for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                result += manifolds[i]->Metric(x.GetElement(j), etax.GetElement(j), xix.GetElement(j));
            }
        }
        return result;
	};

	Vector &MultiManifolds::LinearOPEEta(const Variable &x, const LinearOPE &Hx, const Vector &etax, Vector *result) const
	{
        return Manifold::LinearOPEEta(x, Hx, etax, result);
	};

    Vector &MultiManifolds::ScalarTimesVector(const Variable &x, const realdp &scalar, const Vector &etax, Vector *result) const
    {
        return Manifold::ScalarTimesVector(x, scalar, etax, result);
    };

    Vector &MultiManifolds::ScalarVectorAddVector(const Variable &x, const realdp &scalar, const Vector &etax, const Vector &xix, Vector *result) const
    {
        return Manifold::ScalarVectorAddVector(x, scalar, etax, xix, result);
    };

	Vector &MultiManifolds::VectorLinearCombination(const Variable &x, realdp scalar1, const Vector &etax, realdp scalar2, const Vector &xix, Vector *result) const
	{
		return Manifold::VectorLinearCombination(x, scalar1, etax, scalar2, xix, result);
	};

	Vector &MultiManifolds::Projection(const Variable &x, const Vector &etax, Vector *result) const
	{
        /*etax and result may be the same. Store etax first then new memory for result. Otherwise, the memory in etax is also changed.*/
        Vector inetax = etax;
        result->NewMemoryOnWrite();
        
        for(integer i = 0; i < numoftypes; i++)
        {
            for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->Projection(x.GetElement(j), inetax.GetElement(j), &(result->GetElement(j)));
            }
        }
        
        return *result;
	};

    Vector &MultiManifolds::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        /*etax and result may be the same. Store etax first then new memory for result. Otherwise, the memory in etax is also changed.*/
        Vector inetax = etax;
        result->NewMemoryOnWrite();
        
        for(integer i = 0; i < numoftypes; i++)
        {
            for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->ExtrProjection(x.GetElement(j), inetax.GetElement(j), &(result->GetElement(j)));
            }
        }
        
        return *result;
    };

	Variable &MultiManifolds::Retraction(const Variable &x, const Vector &etax, Variable *result) const
	{
        result->NewMemoryOnWrite();
        
        if(IsIntrApproach)
        {
            for(integer i = 0; i < numoftypes; i++)
            {
                for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
                {
                    manifolds[i]->Retraction(x.GetElement(j), etax.GetElement(j), &(result->GetElement(j)));
                }
            }
            return *result;
        }
        
        /*if using extrinsic representation, then*/
        for(integer i = 0; i < numoftypes; i++)
        {
            for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->SetIsIntrApproach(false);
                manifolds[i]->Retraction(x.GetElement(j), etax.GetElement(j), &(result->GetElement(j)));
                manifolds[i]->SetIsIntrApproach(true);
            }
        }
        return *result;
	};

    Vector &MultiManifolds::InvRetraction(const Variable &x, const Variable &y, Vector *result) const
    {
        result->NewMemoryOnWrite();
        if(IsIntrApproach)
        {
            for(integer i = 0; i < numoftypes; i++)
            {
                for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
                {
                    manifolds[i]->InvRetraction(x.GetElement(j), y.GetElement(j), &(result->GetElement(j)));
                }
            }
            return *result;
        }
        
        for(integer i = 0; i < numoftypes; i++)
        {
            for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->SetIsIntrApproach(false);
                manifolds[i]->InvRetraction(x.GetElement(j), y.GetElement(j), &(result->GetElement(j)));
                manifolds[i]->SetIsIntrApproach(true);
            }
        }
        
        return *result;
    };

	Vector &MultiManifolds::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
        /*xiy and result may be the same. Store xiy first then allocate memory for result. Otherwise, the memory in xiy may be also changed.*/
        Vector inxiy = xiy;
        result->NewMemoryOnWrite();
        
        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->coTangentVector(x.GetElement(j), etax.GetElement(j), y.GetElement(j), inxiy.GetElement(j), &(result->GetElement(j)));
            }
        }
        return *result;
	};

	realdp MultiManifolds::Beta(const Variable &x, Vector &etax) const
	{
		if (!HasHHR)
			return 1;

		if (etax.FieldsExist("beta"))
		{
            const realdp *betav = etax.Field("beta").ObtainReadData();
			return betav[0];
		}

		const realdp *betav;
		realdp numerator = 0, denominator = 0, tmp;
		for (integer i = 0; i < numoftypes; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				if (x.GetElement(j).FieldsExist("beta"))
				{
					betav = x.GetElement(j).Field("beta").ObtainReadData();
					numerator += betav[1];
					denominator += betav[2];
				}
				else
				{
                    tmp = manifolds[i]->Metric(x.GetElement(j), etax.GetElement(j), etax.GetElement(j));
					numerator += tmp;
					denominator += tmp;
				}
			}
		}
		return sqrt(numerator / denominator);
	};

	Vector &MultiManifolds::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
        /*xix and result may be the same. Store xix first then allocate memory for result. Otherwise, the memory in xix may be also changed.*/
        Vector inxix = xix;
        result->NewMemoryOnWrite();
        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->DiffRetraction(x.GetElement(j), etax.GetElement(j), y.GetElement(j), inxix.GetElement(j), &(result->GetElement(j)), IsEtaXiSameDir);
            }
        }
        
        if (IsEtaXiSameDir && HasHHR)
        {
            Vector beta(3);
            realdp *betaptr = beta.ObtainWriteEntireData();
            realdp EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
            betaptr[0] = std::sqrt(Metric(x, etax, etax) / Metric(x, *result, *result)) / EtatoXi;
            betaptr[1] = Metric(x, etax, etax);
            betaptr[2] = Metric(x, *result, *result) * EtatoXi * EtatoXi;
            etax.AddToFields("beta", beta);
            
            if (HasHHR)
            {
                etax.AddToFields("betaTReta", (*result) * (betaptr[0] * EtatoXi));
            }
        }
        
        return *result;
	};

	realdp MultiManifolds::Dist(const Variable &x1, const Variable &x2) const
	{
		realdp result = 0, tmp = 0;
		for (integer i = 0; i < numoftypes; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				tmp = manifolds[i]->Dist(x1.GetElement(j), x2.GetElement(j));
				result += tmp * tmp;
			}
		}
		return std::sqrt(result);
	};

	Vector &MultiManifolds::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
		if (HasHHR)
			return LCVectorTransport(x, etax, y, xix, result);
        
        /*xix and result may be the same. Store xix first then allocate memory for result. Otherwise, the memory in xix may be also changed.*/
        Vector inxix = xix;
        result->NewMemoryOnWrite();
        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->VectorTransport(x.GetElement(j), etax.GetElement(j), y.GetElement(j), inxix.GetElement(j), &(result->GetElement(j)));
            }
        }
        return *result;
	};

	Vector &MultiManifolds::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		if (HasHHR)
			return LCInverseVectorTransport(x, etax, y, xiy, result);
        
        /*xiy and result may be the same. Store xiy first then allocate memory for result. Otherwise, the memory in xiy may be also changed.*/
        Vector inxiy = xiy;
        result->NewMemoryOnWrite();
        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->InverseVectorTransport(x.GetElement(j), etax.GetElement(j), y.GetElement(j), inxiy.GetElement(j), &(result->GetElement(j)));
            }
        }
        return *result;
	};

	LinearOPE &MultiManifolds::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
	{
		if (HasHHR)
			return LCTranHInvTran(x, etax, y, Hx, result);

        result->CopyOnWrite();
        
		integer start, end = 0;
		for (integer i = 0; i < numoftypes; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				start = end;
				end = start + etax.GetElement(j).Getlength();
				manifolds[i]->HInvTran(x.GetElement(j), etax.GetElement(j), y.GetElement(j), *result, start, end, result);
				manifolds[i]->TranH(x.GetElement(j), etax.GetElement(j), y.GetElement(j), *result, start, end, result);
			}
		}
        
        return *result;
	};

	LinearOPE &MultiManifolds::HaddScaledRank1OPE(const Variable &x, const LinearOPE &Hx, realdp scalar, const Vector &etax, const Vector &xix, LinearOPE *result) const
	{
        Vector xixflat(xix);
        xixflat.NewMemoryOnWrite();
        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                manifolds[i]->ObtainEtaxFlat(x.GetElement(j), xix.GetElement(j), &(xixflat.GetElement(j)));
            }
        }
        return Manifold::HaddScaledRank1OPE(x, Hx, scalar, etax, xixflat, result);
	};

	Vector &MultiManifolds::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
        result->NewMemoryOnWrite();
        
        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                if (manifolds[i]->GetIsIntrinsic())
                {
                    manifolds[i]->ObtainIntr(x.GetElement(j), etax.GetElement(j), &(result->GetElement(j)));
                }
                else
                {
                    result->GetElement(j) = etax.GetElement(j);
                }
            }
        }
        return *result;
	};

	Vector &MultiManifolds::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
        result->NewMemoryOnWrite();
        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                if (manifolds[i]->GetIsIntrinsic())
                {
                    manifolds[i]->ObtainExtr(x.GetElement(j), intretax.GetElement(j), &(result->GetElement(j)));
                }
                else
                {
                    result->GetElement(j) = intretax.GetElement(j);
                }
            }
        }
        return *result;
	};

	Vector &MultiManifolds::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
	{
        if(prob->GetUseHess())
        {
            /*The copy on write is necessary. The reason is that the egf may be from a component in a product of elements.
            Therefore, if CopyOnWrite is not used, then the attached data in x and the product of elements share the same
            memory. This may cause an issue: if the product of elements are released before the attached data in x, then
            release the attached data in x would attempt to delete memory that has been released. This is an error!*/
            Vector Sharedegf(egf);
            Sharedegf.CopyOnWrite();
            x.AddToFields("EGrad", Sharedegf);
        }
        
        /*egf and result may be the same. Store egf first then allocate memory for result. Otherwise, the memory in egf may be also changed.*/
        Vector inegf = egf;
        result->NewMemoryOnWrite();

        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                egf.CopyFieldsTo(egf.GetElement(j));
                manifolds[i]->EucGradToGrad(x.GetElement(j), inegf.GetElement(j), prob, &(result->GetElement(j)));
            }
        }
        return *result;
	};

	Vector &MultiManifolds::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        /*exix and result may be the same. Store exix first then allocate memory for result. Otherwise, the memory in exix may be also changed.*/
        Vector inexix = exix;
        result->NewMemoryOnWrite();

        for (integer i = 0; i < numoftypes; i++)
        {
            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                inexix.CopyFieldsTo(inexix.GetElement(j));
                manifolds[i]->EucHvToHv(x.GetElement(j), etax.GetElement(j), inexix.GetElement(j), prob, &(result->GetElement(j)));
            }
        }
        return *result;
	};

	void MultiManifolds::CheckParams(void) const
	{
		if (numoftotalmani == 1)
		{
			manifolds[0]->CheckParams();
		}
		else
		{
			Manifold::CheckParams();
			for (integer i = 0; i < numoftypes; i++)
			{
				printf("%d-th manifold parameters (the number is %d):\n", i, powsinterval[i + 1] - powsinterval[i]);
				manifolds[i]->CheckParams();
			}
		}
	};
}; /*end of ROPTLIB namespace*/
