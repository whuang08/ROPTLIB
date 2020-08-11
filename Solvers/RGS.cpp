
#include "Solvers/RGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	RGS::RGS(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RGS::SetDefaultParams(void)
	{
		SolversNSMSubLS::SetDefaultParams();
		SolverName.assign("RGS");
	};

	RGS::~RGS(void)
	{
	};

    void RGS::UpdateData(void)
    {
    };

	void RGS::GetSearchDir(void)
	{
        gfs[0] = gf1;
        for(integer i = 1; i < Lengthgfs; i++)
        {
            gfs[i] = gf1;
            gfs[i].RandGaussian();
            Mani->Projection(x1, gfs[i], &gfs[i]);
//            std::cout << "i:" << i << std::endl;//---
//            gfs[i].GetTranspose().Print("gfi transpose:");//----
        }
        
		/*normalized all the tangent vectors such that their norms equal Eps*/
		realdp tmp = 0;
		for (integer i = 1; i < Lengthgfs; i++)
		{
			tmp = sqrt(Mani->Metric(x1, gfs[i], gfs[i]));
			Mani->ScalarTimesVector(x1, Eps / tmp, gfs[i], &gfs[i]);
		}
        
		/*Apply retraction and obtain points around x1*/
		/*Compute gradients at Xs[i] and transport them to the tangent space at x1*/
        Vector gfXsi(gf1);
		for (integer i = 1; i < Lengthgfs; i++)
		{
			Mani->Retraction(x1, gfs[i], &Xs[i]); nR++;
			Prob->f(Xs[i]); nf++;
            Prob->Grad(Xs[i], &gfXsi);
            Mani->InverseVectorTransport(x1, gfs[i], Xs[i], gfXsi, &gfs[i]); ng++; nVp++;
		}
        
		ndir1 = sqrt(MinPNormConHull(Mani, x1, gfs, Lengthgfs, nullptr, nullptr, minPv));
		subprobtimes++;
        
//        Mani->ScalarTimesVector(x1, -1.0, gf1, &eta1); //---
//        minPv = eta1;//---
//        ndir1 = 1;//---
//        return;//----
        
		/*eta1 is viewed as the search direction*/
		Mani->ScalarTimesVector(x1, static_cast<realdp> (-1) / ndir1, minPv, &eta1);
	};
}; /*end of ROPTLIB namespace*/
