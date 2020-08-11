#include "test/TestElement.h"

using namespace ROPTLIB;

void testElement(void)
{ /*This is used to test whether memory leakage happens in certain cases.*/

	/*Set the random seed*/
	unsigned tt = (unsigned)time(NULL);
    tt = 1; /*The following test is only for random seed zero*/
	std::cout << "seed:" << tt << std::endl;
	genrandseed(tt);
    
    {
        Vector M(4, 4, "complex"); M.RandGaussian();
        M.GetHaarFWT();
        M.GetInvHaarFWT();
        
//        M.Print("M1:");
//        M.GetHaarFWT().Print("M2:");
//        M.GetInvHaarFWT().Print("M3:");
//        M.GetHaarFWT().GetInvHaarFWT().Print("M4:");
    }
    
#ifdef ROPTLIB_WITH_FFTW
    {
        Vector M(4, 4, "complex"); M.RandGaussian();
        M.GetFFT2D(FFTW_FORWARD);
        M.GetFFT2D(FFTW_BACKWARD);
        
//        realdpcomplex coef = {16, 0};
//        M.Print("M1:");
//        M.GetFFT2D(FFTW_FORWARD).Print("M2:");
//        M.GetFFT2D(FFTW_BACKWARD).Print("M3:");
//        M.GetFFT2D(FFTW_FORWARD).GetFFT2D(FFTW_BACKWARD).Print("M4:");
//        (M * coef - M.GetFFT2D(FFTW_FORWARD).GetFFT2D(FFTW_BACKWARD)).Print("diff");
    }
#endif
    
    {
        integer *ir = new integer[10];
        integer *jc = ir + 5;
        realdp *vals = new realdp[5];
        ir[0] = 0; ir[1] = 0; ir[2] = 1; ir[3] = 2; ir[4] = 3;
        jc[0] = 0; jc[1] = 2; jc[2] = 1; jc[3] = 2; jc[4] = 1;
        vals[0] = 1; vals[1] = 2; vals[2] = 1; vals[3] = 1; vals[4] = -1;
        /* SM = [1, 0, 2; 0 1 0; 0 0 1; 0 -1 0] */
        SparseMatrix SM(4, 3, ir, jc, vals, 5);
        delete [] ir;
        delete[] vals;
        
        Vector A(3, 2);
        A.RandGaussian();
//        A.Print("A:");
//        (SM * A).Print("SM * A:");
        
        Vector B(2, 4);
        B.RandGaussian();
//        B.Print("B:");
//        (B * SM).Print("B * SM:");
    }
    
    {
        integer *ir = new integer[10];
        integer *jc = ir + 5;
        realdpcomplex *vals = new realdpcomplex[5];
        ir[0] = 0; ir[1] = 0; ir[2] = 1; ir[3] = 2; ir[4] = 3;
        jc[0] = 0; jc[1] = 2; jc[2] = 1; jc[3] = 2; jc[4] = 1;
        vals[0].r = 1; vals[0].i = -1; vals[1].r = 2; vals[1].i = 1; vals[2].r = 1; vals[2].i = 2; vals[3].r = 1; vals[3].i = 0; vals[4].r = -1; vals[4].i = 5;
        /* SM = [1-i, 0, 2+i; 0, 1 + 2 * i, 0; 0, 0, 1; 0, -1 + 5 * i, 0] */
        SparseMatrix SM(4, 3, ir, jc, vals, 5);
        delete [] ir;
        delete[] vals;
        
        Vector A(3, 2, "complex");
        A.RandGaussian();
//        A.Print("A:");
//        (SM * A).Print("SM * A:");
        
        Vector B(2, 4, "complex");
        B.RandGaussian();
//        B.Print("B:");
//        (B * SM).Print("B * SM:");
    }
    
    {
        Vector X(2, 3), Y;
        X = Y;
    }

    {
        Vector X(2, 3), Y;
        Y = X;
    }

    {
        Vector X(3, 2), Y(2);
        Vector Z(2, &X, 2, &Y, 1), W;
        W = Z;
    }

    {
        Vector X(3, 2), Y(2);
        Vector Z(2, &X, 2, &Y, 1), W;
        Z = W;
    }

    {
        Vector X(3, 2), Y(2);
        Vector Z(2, &X, 2, &Y, 1), W(2, &X, 1, &Y, 2);
        W = Z;
    }

    {
        Vector X(3, 2), Y(2);
        Vector Z(2, &X, 2, &Y, 1), W(2, &X, 1, &Y, 2);
        Z = W;
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX1(2, &X, 2, &Y, 1), ProdX2(2, &ProdX1, 1, &Y, 2), Z;
        ProdX2.RandGaussian();
        Z = ProdX2;
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX1(2, &X, 2, &Y, 1), ProdX2(2, &ProdX1, 1, &Y, 2), Z;
        ProdX2.RandGaussian();
        ProdX2 = Z;
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX1(2, &X, 2, &Y, 1), ProdX2(2, &ProdX1, 1, &Y, 2), ProdX3(2, &ProdX1, 1, &Y, 2);
        ProdX2.RandGaussian();
        ProdX3 = ProdX2;
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX1(2, &X, 2, &Y, 1), ProdX2(2, &ProdX1, 1, &Y, 2), ProdX3(2, &ProdX1, 1, &Y, 2);
        ProdX2.RandGaussian();
        ProdX2 = ProdX3;
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX(2, &X, 2, &Y, 1);
        ProdX.RandGaussian();
        X.SetToZeros();
        ProdX.GetElement(0) = X;
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX(2, &X, 2, &Y, 1);
        ProdX.RandGaussian();
        //Note that this variable Z must be declared after the construction of ProdX.
        //In this case, the variable Z will be deleted before the variable ProdX.
        //Otherwise, the Space used in ProdX has been deleted but still referred by Z.
        //It follows that deleting Z would cause memory problems.
        Vector Z(3, 2);
        Z = ProdX.GetElement(0);
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX(2, &X, 2, &Y, 1);
        ProdX.NewMemoryOnWrite();
		Vector Z(3, 2); Z.SetToZeros();
        Vector W(2); W.RandUnform();
        ProdX.GetElement(0) = Z;
        ProdX.GetElement(1) = Z + 1;
        ProdX.GetElement(2) = W;
    }

    {
        Vector X(3, 2), Y(2);
        Vector ProdX(2, &X, 2, &Y, 1);
        ProdX.RandGaussian();
        Vector Z(3, 2);
        //New memory for Z is necessary, otherwise, there would be a memory error!
        //If new memory is not created for Z, then the space in ProdX.GetElement(0) need be shared twice.
        //Then assigning W to ProdX.GetElement(0) will not make the memory in ProdX consecutive. This undesired
        //feature makes a memroy error!
        Z.NewMemoryOnWrite();
        Z = ProdX.GetElement(0);
        Vector W(3, 2); W.RandUnform();
        ProdX.GetElement(0) = W;
    }
    return;
}
