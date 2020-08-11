/*
This file defines the class for the problem 
min_{X1 \in St(p, n), X2 \in St(p, n), X3 \in St(q, m)}     tr(X1^T B1 X1 D1) + tr(X2^T B2 X2 D2) + tr(X3^T B3 X3 D3), 
where B1, B2 are n by n symmetric matrices, B3 is m by m symmetric matrix, D1 and D2 are p by p diagonal matrices, and
D3 is q by q diagonal matrix.

Problem --> ProdStieSumBrockett

---- WH
*/

#ifndef PRODSTIESUMBROCKETT_H
#define PRODSTIESUMBROCKETT_H

#include "Manifolds/Stiefel.h"
#include "Manifolds/MultiManifolds.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ProdStieSumBrockett : public Problem{
	public:
		ProdStieSumBrockett(Vector inB1, Vector inD1, Vector inB2, Vector inD2, Vector inB3, Vector inD3);
		virtual ~ProdStieSumBrockett();
        virtual realdp f(const Variable &x) const;

        virtual Vector &EucGrad(const Variable &x, Vector *result) const;
        virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;
        
        Vector B1;
        Vector D1;
        Vector B2;
        Vector D2;
        Vector B3;
        Vector D3;
		integer n;
		integer p;
		integer m;
		integer q;
	};
}; /*end of ROPTLIB namespace*/

#endif /* end of PRODSTIESUMBROCKETT_H */
