/*
This file defines the class for the unit sphere S^{n-1} = \{x \in R^{n} | x^T x = 1\}
It defines the common properties and features of the manifold.

Manifold --> Stiefel --> Sphere

---- WH
*/

#ifndef SPHERE_H
#define SPHERE_H

#include "Manifolds/Stiefel.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Sphere : public Stiefel{
	public:
		/*Construct the unit sphere S^{n-1}*/
		Sphere(integer n);

		/*Delete the sphere*/
		virtual ~Sphere();

		/* choose qf retraction, parallelization and intrinsic approach and no householder reflections */
		void ChooseParamsSet1();

		/* choose exponential map, parallel translation and extrinsic approach and no householder reflections
		Even though the Householder reflections are not used, the locking condition is satisfied.*/
		void ChooseParamsSet2();

		/* choose qf, parallel translation and extrinsic approach and no householder reflections
		The locking conidition is not satisfied*/
		void ChooseParamsSet3();

		/* choose qf, parallel translation and extrinsic approach and no householder reflections
		Beta \neq 1 is used and the locking conidition is satisfied*/
		void ChooseParamsSet4();

		/*Beside the exponential mapping of the sphere, the retractions defined in Stiefel.h also can be used.*/
		virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

		/*Beside the cotangent vector of exponential mapping of the sphere, the retractions defined in Stiefel.h also can be used.*/
		virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Beside the vector transport by differeitiated the exponential mapping, the differentiated retraction defined in Stiefel.h
		also can be used.*/
		virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Beside the parallel translation of the unit sphere, the vector transport defined in Stiefel.h also can be used.*/
		virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/*Beside the inverse parallel translation of the unit sphere, the inverse vector transport defined in Stiefel.h also can be used.*/
		virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*If parallel translation is used, then call member function "ExpHInvTran",
		otherwise, call function in Stiefel::HInvTran*/
        virtual LinearOPE &HInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*If parallel translation is used, then call member function "ExpTranH",
		otherwise, call function in Stiefel::TranH*/
        virtual LinearOPE &TranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*If parallel translation is used, then call member function "ExpTranHInvTran",
		otherwise, call function in Stiefel::TranHInvTran*/
        virtual LinearOPE &TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);
	private:

		/* exponential mapping using extrinsic approach*/
		virtual Variable &ExpRetraction(const Variable &x, const Vector &etax, Variable *result) const;

		/* the cotangent vector using exponential mapping and extrinsic approach*/
		virtual Vector &ExpcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/* the vector transport by differentiated exponential mapping using extrinsic representation*/
		virtual Vector &ExpDiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/* parallel translation using extrinsic representation*/
		virtual Vector &ExpVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/* inverse parallel translation using extrinsic representation */
		virtual Vector &ExpInverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Compute result = H(:, start : end) * \mathcal{T}^{-1}, where H(:, start : end) denotes the matrix formed by columns from "start" to "end".
		\mathcal{T}^{-1} is the inverse parallel translation*/
		virtual LinearOPE &ExpHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H(start : end, :), where H(start : end, :) denotes the matrix formed by cows from "start" to "end".
		\mathcal{T} is the parallel translation*/
		virtual LinearOPE &ExpTranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}. \mathcal{T} is the parallel translation*/
		virtual LinearOPE &ExpTranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of SPHERE_H */
