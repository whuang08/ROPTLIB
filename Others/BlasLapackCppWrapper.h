/*
Some blas and lapack functions are not easy to use. This class defines wrapper functions
which can be used to solve some problems or decomposition, e.g., SVD decomposition, Sylevster equation.
More functions will be added in this class.

-----WH
*/
#pragma   once  
#ifndef BLASLAPACKCPPWRAPPER_H
#define BLASLAPACKCPPWRAPPER_H

#include "Others/SparseBLAS/blas_sparse.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include "matrix.h"
#define logical integer
#define integer ptrdiff_t
#define real float
#define doublereal double
typedef integer (*L_fp)();
#endif

#ifndef MATLAB_MEX_FILE
/* blas and lapack related
BLAS */
#include <caxpy.h>
#include <ccopy.h>
#include <cdotc.h>
#include <cdotu.h>
#include <cgbmv.h>
#include <cgemm.h>
#include <cgemv.h>
#include <cgerc.h>
#include <cgeru.h>
#include <chbmv.h>
#include <chemm.h>
#include <chemv.h>
#include <cher.h>
#include <cher2.h>
#include <cher2k.h>
#include <cherk.h>
#include <chpmv.h>
#include <chpr.h>
#include <chpr2.h>
#include <crotg.h>
#include <cscal.h>
#include <csscal.h>
#include <cswap.h>
#include <csymm.h>
#include <csyr2k.h>
#include <csyrk.h>
#include <ctbmv.h>
#include <ctbsv.h>
#include <ctpmv.h>
#include <ctpsv.h>
#include <ctrmm.h>
#include <ctrmv.h>
#include <ctrsm.h>
#include <ctrsv.h>
#include <dasum.h>
#include <daxpy.h>
#include <dcabs1.h>
#include <dcopy.h>
#include <ddot.h>
#include <dgbmv.h>
#include <dgemm.h>
#include <dgemv.h>
#include <dger.h>
#include <dnrm2.h>
#include <drot.h>
#include <drotg.h>
#include <dsbmv.h>
#include <dscal.h>
#include <dspmv.h>
#include <dspr.h>
#include <dspr2.h>
#include <dswap.h>
#include <dsymm.h>
#include <dsymv.h>
#include <dsyr.h>
#include <dsyr2.h>
#include <dsyr2k.h>
#include <dsyrk.h>
#include <dtbmv.h>
#include <dtbsv.h>
#include <dtpmv.h>
#include <dtpsv.h>
#include <dtrmm.h>
#include <dtrmv.h>
#include <dtrsm.h>
#include <dtrsv.h>
#include <dzasum.h>
#include <dznrm2.h>
#include <f2c.h>
#include <icamax.h>
#include <idamax.h>
#include <isamax.h>
#include <izamax.h>
#include <lsame.h>
#include <sasum.h>
#include <saxpy.h>
#include <scasum.h>
#include <scnrm2.h>
#include <scopy.h>
#include <sdot.h>
#include <sgbmv.h>
#include <sgemm.h>
#include <sgemv.h>
#include <sger.h>
#include <snrm2.h>
#include <srot.h>
#include <srotg.h>
#include <ssbmv.h>
#include <sscal.h>
#include <sspmv.h>
#include <sspr.h>
#include <sspr2.h>
#include <sswap.h>
#include <ssymm.h>
#include <ssymv.h>
#include <ssyr.h>
#include <ssyr2.h>
#include <ssyr2k.h>
#include <ssyrk.h>
#include <stbmv.h>
#include <stbsv.h>
#include <stpmv.h>
#include <stpsv.h>
#include <strmm.h>
#include <strmv.h>
#include <strsm.h>
#include <strsv.h>
#include <xerbla.h>
#include <zaxpy.h>
#include <zcopy.h>
#include <zdotc.h>
#include <zdotu.h>
#include <zdscal.h>
#include <zgbmv.h>
#include <zgemm.h>
#include <zgemv.h>
#include <zgerc.h>
#include <zgeru.h>
#include <zhbmv.h>
#include <zhemm.h>
#include <zhemv.h>
#include <zher.h>
#include <zher2.h>
#include <zher2k.h>
#include <zherk.h>
#include <zhpmv.h>
#include <zhpr.h>
#include <zhpr2.h>
#include <zrotg.h>
#include <zscal.h>
#include <zswap.h>
#include <zsymm.h>
#include <zsyr2k.h>
#include <zsyrk.h>
#include <ztbmv.h>
#include <ztbsv.h>
#include <ztpmv.h>
#include <ztpsv.h>
#include <ztrmm.h>
#include <ztrmv.h>
#include <ztrsm.h>
#include <ztrsv.h>




/* LAPACK */
#include <cbdsqr.h>
#include <cgbbrd.h>
#include <cgbcon.h>
#include <cgbequ.h>
#include <cgbrfs.h>
#include <cgbsv.h>
#include <cgbsvx.h>
#include <cgbtf2.h>
#include <cgbtrf.h>
#include <cgbtrs.h>
#include <cgebak.h>
#include <cgebal.h>
#include <cgebd2.h>
#include <cgebrd.h>
#include <cgecon.h>
#include <cgeequ.h>
#include <cgees.h>
#include <cgeesx.h>
#include <cgeev.h>
#include <cgeevx.h>
#include <cgegs.h>
#include <cgegv.h>
#include <cgehd2.h>
#include <cgehrd.h>
#include <cgelq2.h>
#include <cgelqf.h>
#include <cgels.h>
#include <cgelsd.h>
#include <cgelss.h>
#include <cgelsx.h>
#include <cgelsy.h>
#include <cgeql2.h>
#include <cgeqlf.h>
#include <cgeqp3.h>
#include <cgeqpf.h>
#include <cgeqr2.h>
#include <cgeqrf.h>
#include <cgerfs.h>
#include <cgerq2.h>
#include <cgerqf.h>
#include <cgesc2.h>
#include <cgesdd.h>
#include <cgesv.h>
#include <cgesvd.h>
#include <cgesvx.h>
#include <cgetc2.h>
#include <cgetf2.h>
#include <cgetrf.h>
#include <cgetri.h>
#include <cgetrs.h>
#include <cggbak.h>
#include <cggbal.h>
#include <cgges.h>
#include <cggesx.h>
#include <cggev.h>
#include <cggevx.h>
#include <cggglm.h>
#include <cgghrd.h>
#include <cgglse.h>
#include <cggqrf.h>
#include <cggrqf.h>
#include <cggsvd.h>
#include <cggsvp.h>
#include <cgtcon.h>
#include <cgtrfs.h>
#include <cgtsv.h>
#include <cgtsvx.h>
#include <cgttrf.h>
#include <cgttrs.h>
#include <cgtts2.h>
#include <chbev.h>
#include <chbevd.h>
#include <chbevx.h>
#include <chbgst.h>
#include <chbgv.h>
#include <chbgvd.h>
#include <chbgvx.h>
#include <chbtrd.h>
#include <checon.h>
#include <cheev.h>
#include <cheevd.h>
#include <cheevr.h>
#include <cheevx.h>
#include <chegs2.h>
#include <chegst.h>
#include <chegv.h>
#include <chegvd.h>
#include <chegvx.h>
#include <cherfs.h>
#include <chesv.h>
#include <chesvx.h>
#include <chetd2.h>
#include <chetf2.h>
#include <chetrd.h>
#include <chetrf.h>
#include <chetri.h>
#include <chetrs.h>
#include <chgeqz.h>
#include <chpcon.h>
#include <chpev.h>
#include <chpevd.h>
#include <chpevx.h>
#include <chpgst.h>
#include <chpgv.h>
#include <chpgvd.h>
#include <chpgvx.h>
#include <chprfs.h>
#include <chpsv.h>
#include <chpsvx.h>
#include <chptrd.h>
#include <chptrf.h>
#include <chptri.h>
#include <chptrs.h>
#include <chsein.h>
#include <chseqr.h>
#include <clabrd.h>
#include <clacgv.h>
#include <clacon.h>
#include <clacp2.h>
#include <clacpy.h>
#include <clacrm.h>
#include <clacrt.h>
#include <cladiv.h>
#include <claed0.h>
#include <claed7.h>
#include <claed8.h>
#include <claein.h>
#include <claesy.h>
#include <claev2.h>
#include <clags2.h>
#include <clagtm.h>
#include <clahef.h>
#include <clahqr.h>
#include <clahrd.h>
#include <claic1.h>
#include <clals0.h>
#include <clalsa.h>
#include <clalsd.h>
#include <clangb.h>
#include <clange.h>
#include <clangt.h>
#include <clanhb.h>
#include <clanhe.h>
#include <clanhp.h>
#include <clanhs.h>
#include <clanht.h>
#include <clansb.h>
#include <clansp.h>
#include <clansy.h>
#include <clantb.h>
#include <clantp.h>
#include <clantr.h>
#include <clapll.h>
#include <clapmt.h>
#include <claqgb.h>
#include <claqge.h>
#include <claqhb.h>
#include <claqhe.h>
#include <claqhp.h>
#include <claqp2.h>
#include <claqps.h>
#include <claqsb.h>
#include <claqsp.h>
#include <claqsy.h>
#include <clar1v.h>
#include <clar2v.h>
#include <clarcm.h>
#include <clarf.h>
#include <clarfb.h>
#include <clarfg.h>
#include <clarft.h>
#include <clarfx.h>
#include <clargv.h>
#include <clarnv.h>
#include <clarrv.h>
#include <clartg.h>
#include <clartv.h>
#include <clarz.h>
#include <clarzb.h>
#include <clarzt.h>
#include <clascl.h>
#include <claset.h>
#include <clasr.h>
#include <classq.h>
#include <claswp.h>
#include <clasyf.h>
#include <clatbs.h>
#include <clatdf.h>
#include <clatps.h>
#include <clatrd.h>
#include <clatrs.h>
#include <clatrz.h>
#include <clatzm.h>
#include <clauu2.h>
#include <clauum.h>
#include <cpbcon.h>
#include <cpbequ.h>
#include <cpbrfs.h>
#include <cpbstf.h>
#include <cpbsv.h>
#include <cpbsvx.h>
#include <cpbtf2.h>
#include <cpbtrf.h>
#include <cpbtrs.h>
#include <cpocon.h>
#include <cpoequ.h>
#include <cporfs.h>
#include <cposv.h>
#include <cposvx.h>
#include <cpotf2.h>
#include <cpotrf.h>
#include <cpotri.h>
#include <cpotrs.h>
#include <cppcon.h>
#include <cppequ.h>
#include <cpprfs.h>
#include <cppsv.h>
#include <cppsvx.h>
#include <cpptrf.h>
#include <cpptri.h>
#include <cpptrs.h>
#include <cptcon.h>
#include <cpteqr.h>
#include <cptrfs.h>
#include <cptsv.h>
#include <cptsvx.h>
#include <cpttrf.h>
#include <cpttrs.h>
#include <cptts2.h>
#include <crot.h>
#include <cspcon.h>
#include <cspmv.h>
#include <cspr.h>
#include <csprfs.h>
#include <cspsv.h>
#include <cspsvx.h>
#include <csptrf.h>
#include <csptri.h>
#include <csptrs.h>
#include <csrot.h>
#include <csrscl.h>
#include <cstedc.h>
#include <cstegr.h>
#include <cstein.h>
#include <csteqr.h>
#include <csycon.h>
#include <csymv.h>
#include <csyr.h>
#include <csyrfs.h>
#include <csysv.h>
#include <csysvx.h>
#include <csytf2.h>
#include <csytrf.h>
#include <csytri.h>
#include <csytrs.h>
#include <ctbcon.h>
#include <ctbrfs.h>
#include <ctbtrs.h>
#include <ctgevc.h>
#include <ctgex2.h>
#include <ctgexc.h>
#include <ctgsen.h>
#include <ctgsja.h>
#include <ctgsna.h>
#include <ctgsy2.h>
#include <ctgsyl.h>
#include <ctpcon.h>
#include <ctprfs.h>
#include <ctptri.h>
#include <ctptrs.h>
#include <ctrcon.h>
#include <ctrevc.h>
#include <ctrexc.h>
#include <ctrrfs.h>
#include <ctrsen.h>
#include <ctrsna.h>
#include <ctrsyl.h>
#include <ctrti2.h>
#include <ctrtri.h>
#include <ctrtrs.h>
#include <ctzrqf.h>
#include <ctzrzf.h>
#include <cung2l.h>
#include <cung2r.h>
#include <cungbr.h>
#include <cunghr.h>
#include <cungl2.h>
#include <cunglq.h>
#include <cungql.h>
#include <cungqr.h>
#include <cungr2.h>
#include <cungrq.h>
#include <cungtr.h>
#include <cunm2l.h>
#include <cunm2r.h>
#include <cunmbr.h>
#include <cunmhr.h>
#include <cunml2.h>
#include <cunmlq.h>
#include <cunmql.h>
#include <cunmqr.h>
#include <cunmr2.h>
#include <cunmr3.h>
#include <cunmrq.h>
#include <cunmrz.h>
#include <cunmtr.h>
#include <cupgtr.h>
#include <cupmtr.h>
#include <dbdsdc.h>
#include <dbdsqr.h>
#include <ddisna.h>
#include <dgbbrd.h>
#include <dgbcon.h>
#include <dgbequ.h>
#include <dgbrfs.h>
#include <dgbsv.h>
#include <dgbsvx.h>
#include <dgbtf2.h>
#include <dgbtrf.h>
#include <dgbtrs.h>
#include <dgebak.h>
#include <dgebal.h>
#include <dgebd2.h>
#include <dgebrd.h>
#include <dgecon.h>
#include <dgeequ.h>
#include <dgees.h>
#include <dgeesx.h>
#include <dgeev.h>
#include <dgeevx.h>
#include <dgegs.h>
#include <dgegv.h>
#include <dgehd2.h>
#include <dgehrd.h>
#include <dgelq2.h>
#include <dgelqf.h>
#include <dgels.h>
#include <dgelsd.h>
#include <dgelss.h>
#include <dgelsx.h>
#include <dgelsy.h>
#include <dgeql2.h>
#include <dgeqlf.h>
#include <dgeqp3.h>
#include <dgeqpf.h>
#include <dgeqr2.h>
#include <dgeqrf.h>
#include <dgerfs.h>
#include <dgerq2.h>
#include <dgerqf.h>
#include <dgesc2.h>
#include <dgesdd.h>
#include <dgesv.h>
#include <dgesvd.h>
#include <dgesvx.h>
#include <dgetc2.h>
#include <dgetf2.h>
#include <dgetrf.h>
#include <dgetri.h>
#include <dgetrs.h>
#include <dggbak.h>
#include <dggbal.h>
#include <dgges.h>
#include <dggesx.h>
#include <dggev.h>
#include <dggevx.h>
#include <dggglm.h>
#include <dgghrd.h>
#include <dgglse.h>
#include <dggqrf.h>
#include <dggrqf.h>
#include <dggsvd.h>
#include <dggsvp.h>
#include <dgtcon.h>
#include <dgtrfs.h>
#include <dgtsv.h>
#include <dgtsvx.h>
#include <dgttrf.h>
#include <dgttrs.h>
#include <dgtts2.h>
#include <dhgeqz.h>
#include <dhsein.h>
#include <dhseqr.h>
#include <dlabad.h>
#include <dlabrd.h>
#include <dlacon.h>
#include <dlacpy.h>
#include <dladiv.h>
#include <dlae2.h>
#include <dlaebz.h>
#include <dlaed0.h>
#include <dlaed1.h>
#include <dlaed2.h>
#include <dlaed3.h>
#include <dlaed4.h>
#include <dlaed5.h>
#include <dlaed6.h>
#include <dlaed7.h>
#include <dlaed8.h>
#include <dlaed9.h>
#include <dlaeda.h>
#include <dlaein.h>
#include <dlaev2.h>
#include <dlaexc.h>
#include <dlag2.h>
#include <dlags2.h>
#include <dlagtf.h>
#include <dlagtm.h>
#include <dlagts.h>
#include <dlagv2.h>
#include <dlahqr.h>
#include <dlahrd.h>
#include <dlaic1.h>
#include <dlaln2.h>
#include <dlals0.h>
#include <dlalsa.h>
#include <dlalsd.h>
#include <dlamch.h>
#include <dlamrg.h>
#include <dlangb.h>
#include <dlange.h>
#include <dlangt.h>
#include <dlanhs.h>
#include <dlansb.h>
#include <dlansp.h>
#include <dlanst.h>
#include <dlansy.h>
#include <dlantb.h>
#include <dlantp.h>
#include <dlantr.h>
#include <dlanv2.h>
#include <dlapll.h>
#include <dlapmt.h>
#include <dlapy2.h>
#include <dlapy3.h>
#include <dlaqgb.h>
#include <dlaqge.h>
#include <dlaqp2.h>
#include <dlaqps.h>
#include <dlaqsb.h>
#include <dlaqsp.h>
#include <dlaqsy.h>
#include <dlaqtr.h>
#include <dlar1v.h>
#include <dlar2v.h>
#include <dlarf.h>
#include <dlarfb.h>
#include <dlarfg.h>
#include <dlarft.h>
#include <dlarfx.h>
#include <dlargv.h>
#include <dlarnv.h>
#include <dlarrb.h>
#include <dlarre.h>
#include <dlarrf.h>
#include <dlarrv.h>
#include <dlartg.h>
#include <dlartv.h>
#include <dlaruv.h>
#include <dlarz.h>
#include <dlarzb.h>
#include <dlarzt.h>
#include <dlas2.h>
#include <dlascl.h>
#include <dlasd0.h>
#include <dlasd1.h>
#include <dlasd2.h>
#include <dlasd3.h>
#include <dlasd4.h>
#include <dlasd5.h>
#include <dlasd6.h>
#include <dlasd7.h>
#include <dlasd8.h>
#include <dlasd9.h>
#include <dlasda.h>
#include <dlasdq.h>
#include <dlasdt.h>
#include <dlaset.h>
#include <dlasq1.h>
#include <dlasq2.h>
#include <dlasq3.h>
#include <dlasq4.h>
#include <dlasq5.h>
#include <dlasq6.h>
#include <dlasr.h>
#include <dlasrt.h>
#include <dlassq.h>
#include <dlasv2.h>
#include <dlaswp.h>
#include <dlasy2.h>
#include <dlasyf.h>
#include <dlatbs.h>
#include <dlatdf.h>
#include <dlatps.h>
#include <dlatrd.h>
#include <dlatrs.h>
#include <dlatrz.h>
#include <dlatzm.h>
#include <dlauu2.h>
#include <dlauum.h>
#include <dopgtr.h>
#include <dopmtr.h>
#include <dorg2l.h>
#include <dorg2r.h>
#include <dorgbr.h>
#include <dorghr.h>
#include <dorgl2.h>
#include <dorglq.h>
#include <dorgql.h>
#include <dorgqr.h>
#include <dorgr2.h>
#include <dorgrq.h>
#include <dorgtr.h>
#include <dorm2l.h>
#include <dorm2r.h>
#include <dormbr.h>
#include <dormhr.h>
#include <dorml2.h>
#include <dormlq.h>
#include <dormql.h>
#include <dormqr.h>
#include <dormr2.h>
#include <dormr3.h>
#include <dormrq.h>
#include <dormrz.h>
#include <dormtr.h>
#include <dpbcon.h>
#include <dpbequ.h>
#include <dpbrfs.h>
#include <dpbstf.h>
#include <dpbsv.h>
#include <dpbsvx.h>
#include <dpbtf2.h>
#include <dpbtrf.h>
#include <dpbtrs.h>
#include <dpocon.h>
#include <dpoequ.h>
#include <dporfs.h>
#include <dposv.h>
#include <dposvx.h>
#include <dpotf2.h>
#include <dpotrf.h>
#include <dpotri.h>
#include <dpotrs.h>
#include <dppcon.h>
#include <dppequ.h>
#include <dpprfs.h>
#include <dppsv.h>
#include <dppsvx.h>
#include <dpptrf.h>
#include <dpptri.h>
#include <dpptrs.h>
#include <dptcon.h>
#include <dpteqr.h>
#include <dptrfs.h>
#include <dptsv.h>
#include <dptsvx.h>
#include <dpttrf.h>
#include <dpttrs.h>
#include <dptts2.h>
#include <drscl.h>
#include <dsbev.h>
#include <dsbevd.h>
#include <dsbevx.h>
#include <dsbgst.h>
#include <dsbgv.h>
#include <dsbgvd.h>
#include <dsbgvx.h>
#include <dsbtrd.h>
#include <dsecnd.h>
#include <dspcon.h>
#include <dspev.h>
#include <dspevd.h>
#include <dspevx.h>
#include <dspgst.h>
#include <dspgv.h>
#include <dspgvd.h>
#include <dspgvx.h>
#include <dsprfs.h>
#include <dspsv.h>
#include <dspsvx.h>
#include <dsptrd.h>
#include <dsptrf.h>
#include <dsptri.h>
#include <dsptrs.h>
#include <dstebz.h>
#include <dstedc.h>
#include <dstegr.h>
#include <dstein.h>
#include <dsteqr.h>
#include <dsterf.h>
#include <dstev.h>
#include <dstevd.h>
#include <dstevr.h>
#include <dstevx.h>
#include <dsycon.h>
#include <dsyev.h>
#include <dsyevd.h>
#include <dsyevr.h>
#include <dsyevx.h>
#include <dsygs2.h>
#include <dsygst.h>
#include <dsygv.h>
#include <dsygvd.h>
#include <dsygvx.h>
#include <dsyrfs.h>
#include <dsysv.h>
#include <dsysvx.h>
#include <dsytd2.h>
#include <dsytf2.h>
#include <dsytrd.h>
#include <dsytrf.h>
#include <dsytri.h>
#include <dsytrs.h>
#include <dtbcon.h>
#include <dtbrfs.h>
#include <dtbtrs.h>
#include <dtgevc.h>
#include <dtgex2.h>
#include <dtgexc.h>
#include <dtgsen.h>
#include <dtgsja.h>
#include <dtgsna.h>
#include <dtgsy2.h>
#include <dtgsyl.h>
#include <dtpcon.h>
#include <dtprfs.h>
#include <dtptri.h>
#include <dtptrs.h>
#include <dtrcon.h>
#include <dtrevc.h>
#include <dtrexc.h>
#include <dtrrfs.h>
#include <dtrsen.h>
#include <dtrsna.h>
#include <dtrsyl.h>
#include <dtrti2.h>
#include <dtrtri.h>
#include <dtrtrs.h>
#include <dtzrqf.h>
#include <dtzrzf.h>
#include <dzsum1.h>
#include <f2c.h>
#include <icmax1.h>
#include <ieeeck.h>
#include <ilaenv.h>
#include <izmax1.h>
#include <lsame.h>
#include <lsamen.h>
#include <sbdsdc.h>
#include <sbdsqr.h>
#include <scsum1.h>
#include <sdisna.h>
#include <second.h>
#include <sgbbrd.h>
#include <sgbcon.h>
#include <sgbequ.h>
#include <sgbrfs.h>
#include <sgbsv.h>
#include <sgbsvx.h>
#include <sgbtf2.h>
#include <sgbtrf.h>
#include <sgbtrs.h>
#include <sgebak.h>
#include <sgebal.h>
#include <sgebd2.h>
#include <sgebrd.h>
#include <sgecon.h>
#include <sgeequ.h>
#include <sgees.h>
#include <sgeesx.h>
#include <sgeev.h>
#include <sgeevx.h>
#include <sgegs.h>
#include <sgegv.h>
#include <sgehd2.h>
#include <sgehrd.h>
#include <sgelq2.h>
#include <sgelqf.h>
#include <sgels.h>
#include <sgelsd.h>
#include <sgelss.h>
#include <sgelsx.h>
#include <sgelsy.h>
#include <sgeql2.h>
#include <sgeqlf.h>
#include <sgeqp3.h>
#include <sgeqpf.h>
#include <sgeqr2.h>
#include <sgeqrf.h>
#include <sgerfs.h>
#include <sgerq2.h>
#include <sgerqf.h>
#include <sgesc2.h>
#include <sgesdd.h>
#include <sgesv.h>
#include <sgesvd.h>
#include <sgesvx.h>
#include <sgetc2.h>
#include <sgetf2.h>
#include <sgetrf.h>
#include <sgetri.h>
#include <sgetrs.h>
#include <sggbak.h>
#include <sggbal.h>
#include <sgges.h>
#include <sggesx.h>
#include <sggev.h>
#include <sggevx.h>
#include <sggglm.h>
#include <sgghrd.h>
#include <sgglse.h>
#include <sggqrf.h>
#include <sggrqf.h>
#include <sggsvd.h>
#include <sggsvp.h>
#include <sgtcon.h>
#include <sgtrfs.h>
#include <sgtsv.h>
#include <sgtsvx.h>
#include <sgttrf.h>
#include <sgttrs.h>
#include <sgtts2.h>
#include <shgeqz.h>
#include <shsein.h>
#include <shseqr.h>
#include <slabad.h>
#include <slabrd.h>
#include <slacon.h>
#include <slacpy.h>
#include <sladiv.h>
#include <slae2.h>
#include <slaebz.h>
#include <slaed0.h>
#include <slaed1.h>
#include <slaed2.h>
#include <slaed3.h>
#include <slaed4.h>
#include <slaed5.h>
#include <slaed6.h>
#include <slaed7.h>
#include <slaed8.h>
#include <slaed9.h>
#include <slaeda.h>
#include <slaein.h>
#include <slaev2.h>
#include <slaexc.h>
#include <slag2.h>
#include <slags2.h>
#include <slagtf.h>
#include <slagtm.h>
#include <slagts.h>
#include <slagv2.h>
#include <slahqr.h>
#include <slahrd.h>
#include <slaic1.h>
#include <slaln2.h>
#include <slals0.h>
#include <slalsa.h>
#include <slalsd.h>
#include <slamch.h>
#include <slamrg.h>
#include <slangb.h>
#include <slange.h>
#include <slangt.h>
#include <slanhs.h>
#include <slansb.h>
#include <slansp.h>
#include <slanst.h>
#include <slansy.h>
#include <slantb.h>
#include <slantp.h>
#include <slantr.h>
#include <slanv2.h>
#include <slapll.h>
#include <slapmt.h>
#include <slapy2.h>
#include <slapy3.h>
#include <slaqgb.h>
#include <slaqge.h>
#include <slaqp2.h>
#include <slaqps.h>
#include <slaqsb.h>
#include <slaqsp.h>
#include <slaqsy.h>
#include <slaqtr.h>
#include <slar1v.h>
#include <slar2v.h>
#include <slarf.h>
#include <slarfb.h>
#include <slarfg.h>
#include <slarft.h>
#include <slarfx.h>
#include <slargv.h>
#include <slarnv.h>
#include <slarrb.h>
#include <slarre.h>
#include <slarrf.h>
#include <slarrv.h>
#include <slartg.h>
#include <slartv.h>
#include <slaruv.h>
#include <slarz.h>
#include <slarzb.h>
#include <slarzt.h>
#include <slas2.h>
#include <slascl.h>
#include <slasd0.h>
#include <slasd1.h>
#include <slasd2.h>
#include <slasd3.h>
#include <slasd4.h>
#include <slasd5.h>
#include <slasd6.h>
#include <slasd7.h>
#include <slasd8.h>
#include <slasd9.h>
#include <slasda.h>
#include <slasdq.h>
#include <slasdt.h>
#include <slaset.h>
#include <slasq1.h>
#include <slasq2.h>
#include <slasq3.h>
#include <slasq4.h>
#include <slasq5.h>
#include <slasq6.h>
#include <slasr.h>
#include <slasrt.h>
#include <slassq.h>
#include <slasv2.h>
#include <slaswp.h>
#include <slasy2.h>
#include <slasyf.h>
#include <slatbs.h>
#include <slatdf.h>
#include <slatps.h>
#include <slatrd.h>
#include <slatrs.h>
#include <slatrz.h>
#include <slatzm.h>
#include <slauu2.h>
#include <slauum.h>
#include <sopgtr.h>
#include <sopmtr.h>
#include <sorg2l.h>
#include <sorg2r.h>
#include <sorgbr.h>
#include <sorghr.h>
#include <sorgl2.h>
#include <sorglq.h>
#include <sorgql.h>
#include <sorgqr.h>
#include <sorgr2.h>
#include <sorgrq.h>
#include <sorgtr.h>
#include <sorm2l.h>
#include <sorm2r.h>
#include <sormbr.h>
#include <sormhr.h>
#include <sorml2.h>
#include <sormlq.h>
#include <sormql.h>
#include <sormqr.h>
#include <sormr2.h>
#include <sormr3.h>
#include <sormrq.h>
#include <sormrz.h>
#include <sormtr.h>
#include <spbcon.h>
#include <spbequ.h>
#include <spbrfs.h>
#include <spbstf.h>
#include <spbsv.h>
#include <spbsvx.h>
#include <spbtf2.h>
#include <spbtrf.h>
#include <spbtrs.h>
#include <spocon.h>
#include <spoequ.h>
#include <sporfs.h>
#include <sposv.h>
#include <sposvx.h>
#include <spotf2.h>
#include <spotrf.h>
#include <spotri.h>
#include <spotrs.h>
#include <sppcon.h>
#include <sppequ.h>
#include <spprfs.h>
#include <sppsv.h>
#include <sppsvx.h>
#include <spptrf.h>
#include <spptri.h>
#include <spptrs.h>
#include <sptcon.h>
#include <spteqr.h>
#include <sptrfs.h>
#include <sptsv.h>
#include <sptsvx.h>
#include <spttrf.h>
#include <spttrs.h>
#include <sptts2.h>
#include <srscl.h>
#include <ssbev.h>
#include <ssbevd.h>
#include <ssbevx.h>
#include <ssbgst.h>
#include <ssbgv.h>
#include <ssbgvd.h>
#include <ssbgvx.h>
#include <ssbtrd.h>
#include <sspcon.h>
#include <sspev.h>
#include <sspevd.h>
#include <sspevx.h>
#include <sspgst.h>
#include <sspgv.h>
#include <sspgvd.h>
#include <sspgvx.h>
#include <ssprfs.h>
#include <sspsv.h>
#include <sspsvx.h>
#include <ssptrd.h>
#include <ssptrf.h>
#include <ssptri.h>
#include <ssptrs.h>
#include <sstebz.h>
#include <sstedc.h>
#include <sstegr.h>
#include <sstein.h>
#include <ssteqr.h>
#include <ssterf.h>
#include <sstev.h>
#include <sstevd.h>
#include <sstevr.h>
#include <sstevx.h>
#include <ssycon.h>
#include <ssyev.h>
#include <ssyevd.h>
#include <ssyevr.h>
#include <ssyevx.h>
#include <ssygs2.h>
#include <ssygst.h>
#include <ssygv.h>
#include <ssygvd.h>
#include <ssygvx.h>
#include <ssyrfs.h>
#include <ssysv.h>
#include <ssysvx.h>
#include <ssytd2.h>
#include <ssytf2.h>
#include <ssytrd.h>
#include <ssytrf.h>
#include <ssytri.h>
#include <ssytrs.h>
#include <stbcon.h>
#include <stbrfs.h>
#include <stbtrs.h>
#include <stgevc.h>
#include <stgex2.h>
#include <stgexc.h>
#include <stgsen.h>
#include <stgsja.h>
#include <stgsna.h>
#include <stgsy2.h>
#include <stgsyl.h>
#include <stpcon.h>
#include <stprfs.h>
#include <stptri.h>
#include <stptrs.h>
#include <strcon.h>
#include <strevc.h>
#include <strexc.h>
#include <strrfs.h>
#include <strsen.h>
#include <strsna.h>
#include <strsyl.h>
#include <strti2.h>
#include <strtri.h>
#include <strtrs.h>
#include <stzrqf.h>
#include <stzrzf.h>
#include <xerbla.h>
#include <zbdsqr.h>
#include <zdrot.h>
#include <zdrscl.h>
#include <zgbbrd.h>
#include <zgbcon.h>
#include <zgbequ.h>
#include <zgbrfs.h>
#include <zgbsv.h>
#include <zgbsvx.h>
#include <zgbtf2.h>
#include <zgbtrf.h>
#include <zgbtrs.h>
#include <zgebak.h>
#include <zgebal.h>
#include <zgebd2.h>
#include <zgebrd.h>
#include <zgecon.h>
#include <zgeequ.h>
#include <zgees.h>
#include <zgeesx.h>
#include <zgeev.h>
#include <zgeevx.h>
#include <zgegs.h>
#include <zgegv.h>
#include <zgehd2.h>
#include <zgehrd.h>
#include <zgelq2.h>
#include <zgelqf.h>
#include <zgels.h>
#include <zgelsd.h>
#include <zgelss.h>
#include <zgelsx.h>
#include <zgelsy.h>
#include <zgeql2.h>
#include <zgeqlf.h>
#include <zgeqp3.h>
#include <zgeqpf.h>
#include <zgeqr2.h>
#include <zgeqrf.h>
#include <zgerfs.h>
#include <zgerq2.h>
#include <zgerqf.h>
#include <zgesc2.h>
#include <zgesdd.h>
#include <zgesv.h>
#include <zgesvd.h>
#include <zgesvx.h>
#include <zgetc2.h>
#include <zgetf2.h>
#include <zgetrf.h>
#include <zgetri.h>
#include <zgetrs.h>
#include <zggbak.h>
#include <zggbal.h>
#include <zgges.h>
#include <zggesx.h>
#include <zggev.h>
#include <zggevx.h>
#include <zggglm.h>
#include <zgghrd.h>
#include <zgglse.h>
#include <zggqrf.h>
#include <zggrqf.h>
#include <zggsvd.h>
#include <zggsvp.h>
#include <zgtcon.h>
#include <zgtrfs.h>
#include <zgtsv.h>
#include <zgtsvx.h>
#include <zgttrf.h>
#include <zgttrs.h>
#include <zgtts2.h>
#include <zhbev.h>
#include <zhbevd.h>
#include <zhbevx.h>
#include <zhbgst.h>
#include <zhbgv.h>
#include <zhbgvd.h>
#include <zhbgvx.h>
#include <zhbtrd.h>
#include <zhecon.h>
#include <zheev.h>
#include <zheevd.h>
#include <zheevr.h>
#include <zheevx.h>
#include <zhegs2.h>
#include <zhegst.h>
#include <zhegv.h>
#include <zhegvd.h>
#include <zhegvx.h>
#include <zherfs.h>
#include <zhesv.h>
#include <zhesvx.h>
#include <zhetd2.h>
#include <zhetf2.h>
#include <zhetrd.h>
#include <zhetrf.h>
#include <zhetri.h>
#include <zhetrs.h>
#include <zhgeqz.h>
#include <zhpcon.h>
#include <zhpev.h>
#include <zhpevd.h>
#include <zhpevx.h>
#include <zhpgst.h>
#include <zhpgv.h>
#include <zhpgvd.h>
#include <zhpgvx.h>
#include <zhprfs.h>
#include <zhpsv.h>
#include <zhpsvx.h>
#include <zhptrd.h>
#include <zhptrf.h>
#include <zhptri.h>
#include <zhptrs.h>
#include <zhsein.h>
#include <zhseqr.h>
#include <zlabrd.h>
#include <zlacgv.h>
#include <zlacon.h>
#include <zlacp2.h>
#include <zlacpy.h>
#include <zlacrm.h>
#include <zlacrt.h>
#include <zladiv.h>
#include <zlaed0.h>
#include <zlaed7.h>
#include <zlaed8.h>
#include <zlaein.h>
#include <zlaesy.h>
#include <zlaev2.h>
#include <zlags2.h>
#include <zlagtm.h>
#include <zlahef.h>
#include <zlahqr.h>
#include <zlahrd.h>
#include <zlaic1.h>
#include <zlals0.h>
#include <zlalsa.h>
#include <zlalsd.h>
#include <zlangb.h>
#include <zlange.h>
#include <zlangt.h>
#include <zlanhb.h>
#include <zlanhe.h>
#include <zlanhp.h>
#include <zlanhs.h>
#include <zlanht.h>
#include <zlansb.h>
#include <zlansp.h>
#include <zlansy.h>
#include <zlantb.h>
#include <zlantp.h>
#include <zlantr.h>
#include <zlapll.h>
#include <zlapmt.h>
#include <zlaqgb.h>
#include <zlaqge.h>
#include <zlaqhb.h>
#include <zlaqhe.h>
#include <zlaqhp.h>
#include <zlaqp2.h>
#include <zlaqps.h>
#include <zlaqsb.h>
#include <zlaqsp.h>
#include <zlaqsy.h>
#include <zlar1v.h>
#include <zlar2v.h>
#include <zlarcm.h>
#include <zlarf.h>
#include <zlarfb.h>
#include <zlarfg.h>
#include <zlarft.h>
#include <zlarfx.h>
#include <zlargv.h>
#include <zlarnv.h>
#include <zlarrv.h>
#include <zlartg.h>
#include <zlartv.h>
#include <zlarz.h>
#include <zlarzb.h>
#include <zlarzt.h>
#include <zlascl.h>
#include <zlaset.h>
#include <zlasr.h>
#include <zlassq.h>
#include <zlaswp.h>
#include <zlasyf.h>
#include <zlatbs.h>
#include <zlatdf.h>
#include <zlatps.h>
#include <zlatrd.h>
#include <zlatrs.h>
#include <zlatrz.h>
#include <zlatzm.h>
#include <zlauu2.h>
#include <zlauum.h>
#include <zpbcon.h>
#include <zpbequ.h>
#include <zpbrfs.h>
#include <zpbstf.h>
#include <zpbsv.h>
#include <zpbsvx.h>
#include <zpbtf2.h>
#include <zpbtrf.h>
#include <zpbtrs.h>
#include <zpocon.h>
#include <zpoequ.h>
#include <zporfs.h>
#include <zposv.h>
#include <zposvx.h>
#include <zpotf2.h>
#include <zpotrf.h>
#include <zpotri.h>
#include <zpotrs.h>
#include <zppcon.h>
#include <zppequ.h>
#include <zpprfs.h>
#include <zppsv.h>
#include <zppsvx.h>
#include <zpptrf.h>
#include <zpptri.h>
#include <zpptrs.h>
#include <zptcon.h>
#include <zpteqr.h>
#include <zptrfs.h>
#include <zptsv.h>
#include <zptsvx.h>
#include <zpttrf.h>
#include <zpttrs.h>
#include <zptts2.h>
#include <zrot.h>
#include <zspcon.h>
#include <zspmv.h>
#include <zspr.h>
#include <zsprfs.h>
#include <zspsv.h>
#include <zspsvx.h>
#include <zsptrf.h>
#include <zsptri.h>
#include <zsptrs.h>
#include <zstedc.h>
#include <zstegr.h>
#include <zstein.h>
#include <zsteqr.h>
#include <zsycon.h>
#include <zsymv.h>
#include <zsyr.h>
#include <zsyrfs.h>
#include <zsysv.h>
#include <zsysvx.h>
#include <zsytf2.h>
#include <zsytrf.h>
#include <zsytri.h>
#include <zsytrs.h>
#include <ztbcon.h>
#include <ztbrfs.h>
#include <ztbtrs.h>
#include <ztgevc.h>
#include <ztgex2.h>
#include <ztgexc.h>
#include <ztgsen.h>
#include <ztgsja.h>
#include <ztgsna.h>
#include <ztgsy2.h>
#include <ztgsyl.h>
#include <ztpcon.h>
#include <ztprfs.h>
#include <ztptri.h>
#include <ztptrs.h>
#include <ztrcon.h>
#include <ztrevc.h>
#include <ztrexc.h>
#include <ztrrfs.h>
#include <ztrsen.h>
#include <ztrsna.h>
#include <ztrsyl.h>
#include <ztrti2.h>
#include <ztrtri.h>
#include <ztrtrs.h>
#include <ztzrqf.h>
#include <ztzrzf.h>
#include <zung2l.h>
#include <zung2r.h>
#include <zungbr.h>
#include <zunghr.h>
#include <zungl2.h>
#include <zunglq.h>
#include <zungql.h>
#include <zungqr.h>
#include <zungr2.h>
#include <zungrq.h>
#include <zungtr.h>
#include <zunm2l.h>
#include <zunm2r.h>
#include <zunmbr.h>
#include <zunmhr.h>
#include <zunml2.h>
#include <zunmlq.h>
#include <zunmql.h>
#include <zunmqr.h>
#include <zunmr2.h>
#include <zunmr3.h>
#include <zunmrq.h>
#include <zunmrz.h>
#include <zunmtr.h>
#include <zupgtr.h>
#include <zupmtr.h>
#endif /* end of ifndef MATLAB_MEX_FILE */
#ifdef MATLAB_MEX_FILE

#include "mex.h"
#include "blas.h"
#include "lapack.h"

#define caxpy_ caxpy
#define ccopy_ ccopy
#define cdotc_ cdotc
#define cdotu_ cdotu
#define cgbmv_ cgbmv
#define cgemm_ cgemm
#define cgemv_ cgemv
#define cgerc_ cgerc
#define cgeru_ cgeru
#define chbmv_ chbmv
#define chemm_ chemm
#define chemv_ chemv
#define cher_ cher
#define cher2_ cher2
#define cher2k_ cher2k
#define cherk_ cherk
#define chpmv_ chpmv
#define chpr_ chpr
#define chpr2_ chpr2
#define crotg_ crotg
#define cscal_ cscal
#define csscal_ csscal
#define cswap_ cswap
#define csymm_ csymm
#define csyr2k_ csyr2k
#define csyrk_ csyrk
#define ctbmv_ ctbmv
#define ctbsv_ ctbsv
#define ctpmv_ ctpmv
#define ctpsv_ ctpsv
#define ctrmm_ ctrmm
#define ctrmv_ ctrmv
#define ctrsm_ ctrsm
#define ctrsv_ ctrsv
#define dasum_ dasum
#define daxpy_ daxpy
#define dcabs1_ dcabs1
#define dcopy_ dcopy
#define ddot_ ddot
#define dgbmv_ dgbmv
#define dgemm_ dgemm
#define dgemv_ dgemv
#define dger_ dger
#define dnrm2_ dnrm2
#define drot_ drot
#define drotg_ drotg
#define dsbmv_ dsbmv
#define dscal_ dscal
#define dspmv_ dspmv
#define dspr_ dspr
#define dspr2_ dspr2
#define dswap_ dswap
#define dsymm_ dsymm
#define dsymv_ dsymv
#define dsyr_ dsyr
#define dsyr2_ dsyr2
#define dsyr2k_ dsyr2k
#define dsyrk_ dsyrk
#define dtbmv_ dtbmv
#define dtbsv_ dtbsv
#define dtpmv_ dtpmv
#define dtpsv_ dtpsv
#define dtrmm_ dtrmm
#define dtrmv_ dtrmv
#define dtrsm_ dtrsm
#define dtrsv_ dtrsv
#define dzasum_ dzasum
#define dznrm2_ dznrm2
#define f2c_ f2c
#define icamax_ icamax
#define idamax_ idamax
#define isamax_ isamax
#define izamax_ izamax
#define lsame_ lsame
#define sasum_ sasum
#define saxpy_ saxpy
#define scasum_ scasum
#define scnrm2_ scnrm2
#define scopy_ scopy
#define sdot_ sdot
#define sgbmv_ sgbmv
#define sgemm_ sgemm
#define sgemv_ sgemv
#define sger_ sger
#define snrm2_ snrm2
#define srot_ srot
#define srotg_ srotg
#define ssbmv_ ssbmv
#define sscal_ sscal
#define sspmv_ sspmv
#define sspr_ sspr
#define sspr2_ sspr2
#define sswap_ sswap
#define ssymm_ ssymm
#define ssymv_ ssymv
#define ssyr_ ssyr
#define ssyr2_ ssyr2
#define ssyr2k_ ssyr2k
#define ssyrk_ ssyrk
#define stbmv_ stbmv
#define stbsv_ stbsv
#define stpmv_ stpmv
#define stpsv_ stpsv
#define strmm_ strmm
#define strmv_ strmv
#define strsm_ strsm
#define strsv_ strsv
#define xerbla_ xerbla
#define zaxpy_ zaxpy
#define zcopy_ zcopy
#define zdotc_ zdotc
#define zdotu_ zdotu
#define zdscal_ zdscal
#define zgbmv_ zgbmv
#define zgemm_ zgemm
#define zgemv_ zgemv
#define zgerc_ zgerc
#define zgeru_ zgeru
#define zhbmv_ zhbmv
#define zhemm_ zhemm
#define zhemv_ zhemv
#define zher_ zher
#define zher2_ zher2
#define zher2k_ zher2k
#define zherk_ zherk
#define zhpmv_ zhpmv
#define zhpr_ zhpr
#define zhpr2_ zhpr2
#define zrotg_ zrotg
#define zscal_ zscal
#define zswap_ zswap
#define zsymm_ zsymm
#define zsyr2k_ zsyr2k
#define zsyrk_ zsyrk
#define ztbmv_ ztbmv
#define ztbsv_ ztbsv
#define ztpmv_ ztpmv
#define ztpsv_ ztpsv
#define ztrmm_ ztrmm
#define ztrmv_ ztrmv
#define ztrsm_ ztrsm
#define ztrsv_ ztrsv


/* LAPACK */
#define cbdsqr_ cbdsqr
#define cgbbrd_ cgbbrd
#define cgbcon_ cgbcon
#define cgbequ_ cgbequ
#define cgbrfs_ cgbrfs
#define cgbsv_ cgbsv
#define cgbsvx_ cgbsvx
#define cgbtf2_ cgbtf2
#define cgbtrf_ cgbtrf
#define cgbtrs_ cgbtrs
#define cgebak_ cgebak
#define cgebal_ cgebal
#define cgebd2_ cgebd2
#define cgebrd_ cgebrd
#define cgecon_ cgecon
#define cgeequ_ cgeequ
#define cgees_ cgees
#define cgeesx_ cgeesx
#define cgeev_ cgeev
#define cgeevx_ cgeevx
#define cgegs_ cgegs
#define cgegv_ cgegv
#define cgehd2_ cgehd2
#define cgehrd_ cgehrd
#define cgelq2_ cgelq2
#define cgelqf_ cgelqf
#define cgels_ cgels
#define cgelsd_ cgelsd
#define cgelss_ cgelss
#define cgelsx_ cgelsx
#define cgelsy_ cgelsy
#define cgeql2_ cgeql2
#define cgeqlf_ cgeqlf
#define cgeqp3_ cgeqp3
#define cgeqpf_ cgeqpf
#define cgeqr2_ cgeqr2
#define cgeqrf_ cgeqrf
#define cgerfs_ cgerfs
#define cgerq2_ cgerq2
#define cgerqf_ cgerqf
#define cgesc2_ cgesc2
#define cgesdd_ cgesdd
#define cgesv_ cgesv
#define cgesvd_ cgesvd
#define cgesvx_ cgesvx
#define cgetc2_ cgetc2
#define cgetf2_ cgetf2
#define cgetrf_ cgetrf
#define cgetri_ cgetri
#define cgetrs_ cgetrs
#define cggbak_ cggbak
#define cggbal_ cggbal
#define cgges_ cgges
#define cggesx_ cggesx
#define cggev_ cggev
#define cggevx_ cggevx
#define cggglm_ cggglm
#define cgghrd_ cgghrd
#define cgglse_ cgglse
#define cggqrf_ cggqrf
#define cggrqf_ cggrqf
#define cggsvd_ cggsvd
#define cggsvp_ cggsvp
#define cgtcon_ cgtcon
#define cgtrfs_ cgtrfs
#define cgtsv_ cgtsv
#define cgtsvx_ cgtsvx
#define cgttrf_ cgttrf
#define cgttrs_ cgttrs
#define cgtts2_ cgtts2
#define chbev_ chbev
#define chbevd_ chbevd
#define chbevx_ chbevx
#define chbgst_ chbgst
#define chbgv_ chbgv
#define chbgvd_ chbgvd
#define chbgvx_ chbgvx
#define chbtrd_ chbtrd
#define checon_ checon
#define cheev_ cheev
#define cheevd_ cheevd
#define cheevr_ cheevr
#define cheevx_ cheevx
#define chegs2_ chegs2
#define chegst_ chegst
#define chegv_ chegv
#define chegvd_ chegvd
#define chegvx_ chegvx
#define cherfs_ cherfs
#define chesv_ chesv
#define chesvx_ chesvx
#define chetd2_ chetd2
#define chetf2_ chetf2
#define chetrd_ chetrd
#define chetrf_ chetrf
#define chetri_ chetri
#define chetrs_ chetrs
#define chgeqz_ chgeqz
#define chpcon_ chpcon
#define chpev_ chpev
#define chpevd_ chpevd
#define chpevx_ chpevx
#define chpgst_ chpgst
#define chpgv_ chpgv
#define chpgvd_ chpgvd
#define chpgvx_ chpgvx
#define chprfs_ chprfs
#define chpsv_ chpsv
#define chpsvx_ chpsvx
#define chptrd_ chptrd
#define chptrf_ chptrf
#define chptri_ chptri
#define chptrs_ chptrs
#define chsein_ chsein
#define chseqr_ chseqr
#define clabrd_ clabrd
#define clacgv_ clacgv
#define clacon_ clacon
#define clacp2_ clacp2
#define clacpy_ clacpy
#define clacrm_ clacrm
#define clacrt_ clacrt
#define cladiv_ cladiv
#define claed0_ claed0
#define claed7_ claed7
#define claed8_ claed8
#define claein_ claein
#define claesy_ claesy
#define claev2_ claev2
#define clags2_ clags2
#define clagtm_ clagtm
#define clahef_ clahef
#define clahqr_ clahqr
#define clahrd_ clahrd
#define claic1_ claic1
#define clals0_ clals0
#define clalsa_ clalsa
#define clalsd_ clalsd
#define clangb_ clangb
#define clange_ clange
#define clangt_ clangt
#define clanhb_ clanhb
#define clanhe_ clanhe
#define clanhp_ clanhp
#define clanhs_ clanhs
#define clanht_ clanht
#define clansb_ clansb
#define clansp_ clansp
#define clansy_ clansy
#define clantb_ clantb
#define clantp_ clantp
#define clantr_ clantr
#define clapll_ clapll
#define clapmt_ clapmt
#define claqgb_ claqgb
#define claqge_ claqge
#define claqhb_ claqhb
#define claqhe_ claqhe
#define claqhp_ claqhp
#define claqp2_ claqp2
#define claqps_ claqps
#define claqsb_ claqsb
#define claqsp_ claqsp
#define claqsy_ claqsy
#define clar1v_ clar1v
#define clar2v_ clar2v
#define clarcm_ clarcm
#define clarf_ clarf
#define clarfb_ clarfb
#define clarfg_ clarfg
#define clarft_ clarft
#define clarfx_ clarfx
#define clargv_ clargv
#define clarnv_ clarnv
#define clarrv_ clarrv
#define clartg_ clartg
#define clartv_ clartv
#define clarz_ clarz
#define clarzb_ clarzb
#define clarzt_ clarzt
#define clascl_ clascl
#define claset_ claset
#define clasr_ clasr
#define classq_ classq
#define claswp_ claswp
#define clasyf_ clasyf
#define clatbs_ clatbs
#define clatdf_ clatdf
#define clatps_ clatps
#define clatrd_ clatrd
#define clatrs_ clatrs
#define clatrz_ clatrz
#define clatzm_ clatzm
#define clauu2_ clauu2
#define clauum_ clauum
#define cpbcon_ cpbcon
#define cpbequ_ cpbequ
#define cpbrfs_ cpbrfs
#define cpbstf_ cpbstf
#define cpbsv_ cpbsv
#define cpbsvx_ cpbsvx
#define cpbtf2_ cpbtf2
#define cpbtrf_ cpbtrf
#define cpbtrs_ cpbtrs
#define cpocon_ cpocon
#define cpoequ_ cpoequ
#define cporfs_ cporfs
#define cposv_ cposv
#define cposvx_ cposvx
#define cpotf2_ cpotf2
#define cpotrf_ cpotrf
#define cpotri_ cpotri
#define cpotrs_ cpotrs
#define cppcon_ cppcon
#define cppequ_ cppequ
#define cpprfs_ cpprfs
#define cppsv_ cppsv
#define cppsvx_ cppsvx
#define cpptrf_ cpptrf
#define cpptri_ cpptri
#define cpptrs_ cpptrs
#define cptcon_ cptcon
#define cpteqr_ cpteqr
#define cptrfs_ cptrfs
#define cptsv_ cptsv
#define cptsvx_ cptsvx
#define cpttrf_ cpttrf
#define cpttrs_ cpttrs
#define cptts2_ cptts2
#define crot_ crot
#define cspcon_ cspcon
#define cspmv_ cspmv
#define cspr_ cspr
#define csprfs_ csprfs
#define cspsv_ cspsv
#define cspsvx_ cspsvx
#define csptrf_ csptrf
#define csptri_ csptri
#define csptrs_ csptrs
#define csrot_ csrot
#define csrscl_ csrscl
#define cstedc_ cstedc
#define cstegr_ cstegr
#define cstein_ cstein
#define csteqr_ csteqr
#define csycon_ csycon
#define csymv_ csymv
#define csyr_ csyr
#define csyrfs_ csyrfs
#define csysv_ csysv
#define csysvx_ csysvx
#define csytf2_ csytf2
#define csytrf_ csytrf
#define csytri_ csytri
#define csytrs_ csytrs
#define ctbcon_ ctbcon
#define ctbrfs_ ctbrfs
#define ctbtrs_ ctbtrs
#define ctgevc_ ctgevc
#define ctgex2_ ctgex2
#define ctgexc_ ctgexc
#define ctgsen_ ctgsen
#define ctgsja_ ctgsja
#define ctgsna_ ctgsna
#define ctgsy2_ ctgsy2
#define ctgsyl_ ctgsyl
#define ctpcon_ ctpcon
#define ctprfs_ ctprfs
#define ctptri_ ctptri
#define ctptrs_ ctptrs
#define ctrcon_ ctrcon
#define ctrevc_ ctrevc
#define ctrexc_ ctrexc
#define ctrrfs_ ctrrfs
#define ctrsen_ ctrsen
#define ctrsna_ ctrsna
#define ctrsyl_ ctrsyl
#define ctrti2_ ctrti2
#define ctrtri_ ctrtri
#define ctrtrs_ ctrtrs
#define ctzrqf_ ctzrqf
#define ctzrzf_ ctzrzf
#define cung2l_ cung2l
#define cung2r_ cung2r
#define cungbr_ cungbr
#define cunghr_ cunghr
#define cungl2_ cungl2
#define cunglq_ cunglq
#define cungql_ cungql
#define cungqr_ cungqr
#define cungr2_ cungr2
#define cungrq_ cungrq
#define cungtr_ cungtr
#define cunm2l_ cunm2l
#define cunm2r_ cunm2r
#define cunmbr_ cunmbr
#define cunmhr_ cunmhr
#define cunml2_ cunml2
#define cunmlq_ cunmlq
#define cunmql_ cunmql
#define cunmqr_ cunmqr
#define cunmr2_ cunmr2
#define cunmr3_ cunmr3
#define cunmrq_ cunmrq
#define cunmrz_ cunmrz
#define cunmtr_ cunmtr
#define cupgtr_ cupgtr
#define cupmtr_ cupmtr
#define dbdsdc_ dbdsdc
#define dbdsqr_ dbdsqr
#define ddisna_ ddisna
#define dgbbrd_ dgbbrd
#define dgbcon_ dgbcon
#define dgbequ_ dgbequ
#define dgbrfs_ dgbrfs
#define dgbsv_ dgbsv
#define dgbsvx_ dgbsvx
#define dgbtf2_ dgbtf2
#define dgbtrf_ dgbtrf
#define dgbtrs_ dgbtrs
#define dgebak_ dgebak
#define dgebal_ dgebal
#define dgebd2_ dgebd2
#define dgebrd_ dgebrd
#define dgecon_ dgecon
#define dgeequ_ dgeequ
#define dgees_ dgees
#define dgeesx_ dgeesx
#define dgeev_ dgeev
#define dgeevx_ dgeevx
#define dgegs_ dgegs
#define dgegv_ dgegv
#define dgehd2_ dgehd2
#define dgehrd_ dgehrd
#define dgelq2_ dgelq2
#define dgelqf_ dgelqf
#define dgels_ dgels
#define dgelsd_ dgelsd
#define dgelss_ dgelss
#define dgelsx_ dgelsx
#define dgelsy_ dgelsy
#define dgeql2_ dgeql2
#define dgeqlf_ dgeqlf
#define dgeqp3_ dgeqp3
#define dgeqpf_ dgeqpf
#define dgeqr2_ dgeqr2
#define dgeqrf_ dgeqrf
#define dgerfs_ dgerfs
#define dgerq2_ dgerq2
#define dgerqf_ dgerqf
#define dgesc2_ dgesc2
#define dgesdd_ dgesdd
#define dgesv_ dgesv
#define dgesvd_ dgesvd
#define dgesvx_ dgesvx
#define dgetc2_ dgetc2
#define dgetf2_ dgetf2
#define dgetrf_ dgetrf
#define dgetri_ dgetri
#define dgetrs_ dgetrs
#define dggbak_ dggbak
#define dggbal_ dggbal
#define dgges_ dgges
#define dggesx_ dggesx
#define dggev_ dggev
#define dggevx_ dggevx
#define dggglm_ dggglm
#define dgghrd_ dgghrd
#define dgglse_ dgglse
#define dggqrf_ dggqrf
#define dggrqf_ dggrqf
#define dggsvd_ dggsvd
#define dggsvp_ dggsvp
#define dgtcon_ dgtcon
#define dgtrfs_ dgtrfs
#define dgtsv_ dgtsv
#define dgtsvx_ dgtsvx
#define dgttrf_ dgttrf
#define dgttrs_ dgttrs
#define dgtts2_ dgtts2
#define dhgeqz_ dhgeqz
#define dhsein_ dhsein
#define dhseqr_ dhseqr
#define dlabad_ dlabad
#define dlabrd_ dlabrd
#define dlacon_ dlacon
#define dlacpy_ dlacpy
#define dladiv_ dladiv
#define dlae2_ dlae2
#define dlaebz_ dlaebz
#define dlaed0_ dlaed0
#define dlaed1_ dlaed1
#define dlaed2_ dlaed2
#define dlaed3_ dlaed3
#define dlaed4_ dlaed4
#define dlaed5_ dlaed5
#define dlaed6_ dlaed6
#define dlaed7_ dlaed7
#define dlaed8_ dlaed8
#define dlaed9_ dlaed9
#define dlaeda_ dlaeda
#define dlaein_ dlaein
#define dlaev2_ dlaev2
#define dlaexc_ dlaexc
#define dlag2_ dlag2
#define dlags2_ dlags2
#define dlagtf_ dlagtf
#define dlagtm_ dlagtm
#define dlagts_ dlagts
#define dlagv2_ dlagv2
#define dlahqr_ dlahqr
#define dlahrd_ dlahrd
#define dlaic1_ dlaic1
#define dlaln2_ dlaln2
#define dlals0_ dlals0
#define dlalsa_ dlalsa
#define dlalsd_ dlalsd
#define dlamch_ dlamch
#define dlamrg_ dlamrg
#define dlangb_ dlangb
#define dlange_ dlange
#define dlangt_ dlangt
#define dlanhs_ dlanhs
#define dlansb_ dlansb
#define dlansp_ dlansp
#define dlanst_ dlanst
#define dlansy_ dlansy
#define dlantb_ dlantb
#define dlantp_ dlantp
#define dlantr_ dlantr
#define dlanv2_ dlanv2
#define dlapll_ dlapll
#define dlapmt_ dlapmt
#define dlapy2_ dlapy2
#define dlapy3_ dlapy3
#define dlaqgb_ dlaqgb
#define dlaqge_ dlaqge
#define dlaqp2_ dlaqp2
#define dlaqps_ dlaqps
#define dlaqsb_ dlaqsb
#define dlaqsp_ dlaqsp
#define dlaqsy_ dlaqsy
#define dlaqtr_ dlaqtr
#define dlar1v_ dlar1v
#define dlar2v_ dlar2v
#define dlarf_ dlarf
#define dlarfb_ dlarfb
#define dlarfg_ dlarfg
#define dlarft_ dlarft
#define dlarfx_ dlarfx
#define dlargv_ dlargv
#define dlarnv_ dlarnv
#define dlarrb_ dlarrb
#define dlarre_ dlarre
#define dlarrf_ dlarrf
#define dlarrv_ dlarrv
#define dlartg_ dlartg
#define dlartv_ dlartv
#define dlaruv_ dlaruv
#define dlarz_ dlarz
#define dlarzb_ dlarzb
#define dlarzt_ dlarzt
#define dlas2_ dlas2
#define dlascl_ dlascl
#define dlasd0_ dlasd0
#define dlasd1_ dlasd1
#define dlasd2_ dlasd2
#define dlasd3_ dlasd3
#define dlasd4_ dlasd4
#define dlasd5_ dlasd5
#define dlasd6_ dlasd6
#define dlasd7_ dlasd7
#define dlasd8_ dlasd8
#define dlasd9_ dlasd9
#define dlasda_ dlasda
#define dlasdq_ dlasdq
#define dlasdt_ dlasdt
#define dlaset_ dlaset
#define dlasq1_ dlasq1
#define dlasq2_ dlasq2
#define dlasq3_ dlasq3
#define dlasq4_ dlasq4
#define dlasq5_ dlasq5
#define dlasq6_ dlasq6
#define dlasr_ dlasr
#define dlasrt_ dlasrt
#define dlassq_ dlassq
#define dlasv2_ dlasv2
#define dlaswp_ dlaswp
#define dlasy2_ dlasy2
#define dlasyf_ dlasyf
#define dlatbs_ dlatbs
#define dlatdf_ dlatdf
#define dlatps_ dlatps
#define dlatrd_ dlatrd
#define dlatrs_ dlatrs
#define dlatrz_ dlatrz
#define dlatzm_ dlatzm
#define dlauu2_ dlauu2
#define dlauum_ dlauum
#define dopgtr_ dopgtr
#define dopmtr_ dopmtr
#define dorg2l_ dorg2l
#define dorg2r_ dorg2r
#define dorgbr_ dorgbr
#define dorghr_ dorghr
#define dorgl2_ dorgl2
#define dorglq_ dorglq
#define dorgql_ dorgql
#define dorgqr_ dorgqr
#define dorgr2_ dorgr2
#define dorgrq_ dorgrq
#define dorgtr_ dorgtr
#define dorm2l_ dorm2l
#define dorm2r_ dorm2r
#define dormbr_ dormbr
#define dormhr_ dormhr
#define dorml2_ dorml2
#define dormlq_ dormlq
#define dormql_ dormql
#define dormqr_ dormqr
#define dormr2_ dormr2
#define dormr3_ dormr3
#define dormrq_ dormrq
#define dormrz_ dormrz
#define dormtr_ dormtr
#define dpbcon_ dpbcon
#define dpbequ_ dpbequ
#define dpbrfs_ dpbrfs
#define dpbstf_ dpbstf
#define dpbsv_ dpbsv
#define dpbsvx_ dpbsvx
#define dpbtf2_ dpbtf2
#define dpbtrf_ dpbtrf
#define dpbtrs_ dpbtrs
#define dpocon_ dpocon
#define dpoequ_ dpoequ
#define dporfs_ dporfs
#define dposv_ dposv
#define dposvx_ dposvx
#define dpotf2_ dpotf2
#define dpotrf_ dpotrf
#define dpotri_ dpotri
#define dpotrs_ dpotrs
#define dppcon_ dppcon
#define dppequ_ dppequ
#define dpprfs_ dpprfs
#define dppsv_ dppsv
#define dppsvx_ dppsvx
#define dpptrf_ dpptrf
#define dpptri_ dpptri
#define dpptrs_ dpptrs
#define dptcon_ dptcon
#define dpteqr_ dpteqr
#define dptrfs_ dptrfs
#define dptsv_ dptsv
#define dptsvx_ dptsvx
#define dpttrf_ dpttrf
#define dpttrs_ dpttrs
#define dptts2_ dptts2
#define drscl_ drscl
#define dsbev_ dsbev
#define dsbevd_ dsbevd
#define dsbevx_ dsbevx
#define dsbgst_ dsbgst
#define dsbgv_ dsbgv
#define dsbgvd_ dsbgvd
#define dsbgvx_ dsbgvx
#define dsbtrd_ dsbtrd
#define dsecnd_ dsecnd
#define dspcon_ dspcon
#define dspev_ dspev
#define dspevd_ dspevd
#define dspevx_ dspevx
#define dspgst_ dspgst
#define dspgv_ dspgv
#define dspgvd_ dspgvd
#define dspgvx_ dspgvx
#define dsprfs_ dsprfs
#define dspsv_ dspsv
#define dspsvx_ dspsvx
#define dsptrd_ dsptrd
#define dsptrf_ dsptrf
#define dsptri_ dsptri
#define dsptrs_ dsptrs
#define dstebz_ dstebz
#define dstedc_ dstedc
#define dstegr_ dstegr
#define dstein_ dstein
#define dsteqr_ dsteqr
#define dsterf_ dsterf
#define dstev_ dstev
#define dstevd_ dstevd
#define dstevr_ dstevr
#define dstevx_ dstevx
#define dsycon_ dsycon
#define dsyev_ dsyev
#define dsyevd_ dsyevd
#define dsyevr_ dsyevr
#define dsyevx_ dsyevx
#define dsygs2_ dsygs2
#define dsygst_ dsygst
#define dsygv_ dsygv
#define dsygvd_ dsygvd
#define dsygvx_ dsygvx
#define dsyrfs_ dsyrfs
#define dsysv_ dsysv
#define dsysvx_ dsysvx
#define dsytd2_ dsytd2
#define dsytf2_ dsytf2
#define dsytrd_ dsytrd
#define dsytrf_ dsytrf
#define dsytri_ dsytri
#define dsytrs_ dsytrs
#define dtbcon_ dtbcon
#define dtbrfs_ dtbrfs
#define dtbtrs_ dtbtrs
#define dtgevc_ dtgevc
#define dtgex2_ dtgex2
#define dtgexc_ dtgexc
#define dtgsen_ dtgsen
#define dtgsja_ dtgsja
#define dtgsna_ dtgsna
#define dtgsy2_ dtgsy2
#define dtgsyl_ dtgsyl
#define dtpcon_ dtpcon
#define dtprfs_ dtprfs
#define dtptri_ dtptri
#define dtptrs_ dtptrs
#define dtrcon_ dtrcon
#define dtrevc_ dtrevc
#define dtrexc_ dtrexc
#define dtrrfs_ dtrrfs
#define dtrsen_ dtrsen
#define dtrsna_ dtrsna
#define dtrsyl_ dtrsyl
#define dtrti2_ dtrti2
#define dtrtri_ dtrtri
#define dtrtrs_ dtrtrs
#define dtzrqf_ dtzrqf
#define dtzrzf_ dtzrzf
#define dzsum1_ dzsum1
#define f2c_ f2c
#define icmax1_ icmax1
#define ieeeck_ ieeeck
#define ilaenv_ ilaenv
#define izmax1_ izmax1
#define lsame_ lsame
#define lsamen_ lsamen
#define sbdsdc_ sbdsdc
#define sbdsqr_ sbdsqr
#define scsum1_ scsum1
#define sdisna_ sdisna
#define second_ second
#define sgbbrd_ sgbbrd
#define sgbcon_ sgbcon
#define sgbequ_ sgbequ
#define sgbrfs_ sgbrfs
#define sgbsv_ sgbsv
#define sgbsvx_ sgbsvx
#define sgbtf2_ sgbtf2
#define sgbtrf_ sgbtrf
#define sgbtrs_ sgbtrs
#define sgebak_ sgebak
#define sgebal_ sgebal
#define sgebd2_ sgebd2
#define sgebrd_ sgebrd
#define sgecon_ sgecon
#define sgeequ_ sgeequ
#define sgees_ sgees
#define sgeesx_ sgeesx
#define sgeev_ sgeev
#define sgeevx_ sgeevx
#define sgegs_ sgegs
#define sgegv_ sgegv
#define sgehd2_ sgehd2
#define sgehrd_ sgehrd
#define sgelq2_ sgelq2
#define sgelqf_ sgelqf
#define sgels_ sgels
#define sgelsd_ sgelsd
#define sgelss_ sgelss
#define sgelsx_ sgelsx
#define sgelsy_ sgelsy
#define sgeql2_ sgeql2
#define sgeqlf_ sgeqlf
#define sgeqp3_ sgeqp3
#define sgeqpf_ sgeqpf
#define sgeqr2_ sgeqr2
#define sgeqrf_ sgeqrf
#define sgerfs_ sgerfs
#define sgerq2_ sgerq2
#define sgerqf_ sgerqf
#define sgesc2_ sgesc2
#define sgesdd_ sgesdd
#define sgesv_ sgesv
#define sgesvd_ sgesvd
#define sgesvx_ sgesvx
#define sgetc2_ sgetc2
#define sgetf2_ sgetf2
#define sgetrf_ sgetrf
#define sgetri_ sgetri
#define sgetrs_ sgetrs
#define sggbak_ sggbak
#define sggbal_ sggbal
#define sgges_ sgges
#define sggesx_ sggesx
#define sggev_ sggev
#define sggevx_ sggevx
#define sggglm_ sggglm
#define sgghrd_ sgghrd
#define sgglse_ sgglse
#define sggqrf_ sggqrf
#define sggrqf_ sggrqf
#define sggsvd_ sggsvd
#define sggsvp_ sggsvp
#define sgtcon_ sgtcon
#define sgtrfs_ sgtrfs
#define sgtsv_ sgtsv
#define sgtsvx_ sgtsvx
#define sgttrf_ sgttrf
#define sgttrs_ sgttrs
#define sgtts2_ sgtts2
#define shgeqz_ shgeqz
#define shsein_ shsein
#define shseqr_ shseqr
#define slabad_ slabad
#define slabrd_ slabrd
#define slacon_ slacon
#define slacpy_ slacpy
#define sladiv_ sladiv
#define slae2_ slae2
#define slaebz_ slaebz
#define slaed0_ slaed0
#define slaed1_ slaed1
#define slaed2_ slaed2
#define slaed3_ slaed3
#define slaed4_ slaed4
#define slaed5_ slaed5
#define slaed6_ slaed6
#define slaed7_ slaed7
#define slaed8_ slaed8
#define slaed9_ slaed9
#define slaeda_ slaeda
#define slaein_ slaein
#define slaev2_ slaev2
#define slaexc_ slaexc
#define slag2_ slag2
#define slags2_ slags2
#define slagtf_ slagtf
#define slagtm_ slagtm
#define slagts_ slagts
#define slagv2_ slagv2
#define slahqr_ slahqr
#define slahrd_ slahrd
#define slaic1_ slaic1
#define slaln2_ slaln2
#define slals0_ slals0
#define slalsa_ slalsa
#define slalsd_ slalsd
#define slamch_ slamch
#define slamrg_ slamrg
#define slangb_ slangb
#define slange_ slange
#define slangt_ slangt
#define slanhs_ slanhs
#define slansb_ slansb
#define slansp_ slansp
#define slanst_ slanst
#define slansy_ slansy
#define slantb_ slantb
#define slantp_ slantp
#define slantr_ slantr
#define slanv2_ slanv2
#define slapll_ slapll
#define slapmt_ slapmt
#define slapy2_ slapy2
#define slapy3_ slapy3
#define slaqgb_ slaqgb
#define slaqge_ slaqge
#define slaqp2_ slaqp2
#define slaqps_ slaqps
#define slaqsb_ slaqsb
#define slaqsp_ slaqsp
#define slaqsy_ slaqsy
#define slaqtr_ slaqtr
#define slar1v_ slar1v
#define slar2v_ slar2v
#define slarf_ slarf
#define slarfb_ slarfb
#define slarfg_ slarfg
#define slarft_ slarft
#define slarfx_ slarfx
#define slargv_ slargv
#define slarnv_ slarnv
#define slarrb_ slarrb
#define slarre_ slarre
#define slarrf_ slarrf
#define slarrv_ slarrv
#define slartg_ slartg
#define slartv_ slartv
#define slaruv_ slaruv
#define slarz_ slarz
#define slarzb_ slarzb
#define slarzt_ slarzt
#define slas2_ slas2
#define slascl_ slascl
#define slasd0_ slasd0
#define slasd1_ slasd1
#define slasd2_ slasd2
#define slasd3_ slasd3
#define slasd4_ slasd4
#define slasd5_ slasd5
#define slasd6_ slasd6
#define slasd7_ slasd7
#define slasd8_ slasd8
#define slasd9_ slasd9
#define slasda_ slasda
#define slasdq_ slasdq
#define slasdt_ slasdt
#define slaset_ slaset
#define slasq1_ slasq1
#define slasq2_ slasq2
#define slasq3_ slasq3
#define slasq4_ slasq4
#define slasq5_ slasq5
#define slasq6_ slasq6
#define slasr_ slasr
#define slasrt_ slasrt
#define slassq_ slassq
#define slasv2_ slasv2
#define slaswp_ slaswp
#define slasy2_ slasy2
#define slasyf_ slasyf
#define slatbs_ slatbs
#define slatdf_ slatdf
#define slatps_ slatps
#define slatrd_ slatrd
#define slatrs_ slatrs
#define slatrz_ slatrz
#define slatzm_ slatzm
#define slauu2_ slauu2
#define slauum_ slauum
#define sopgtr_ sopgtr
#define sopmtr_ sopmtr
#define sorg2l_ sorg2l
#define sorg2r_ sorg2r
#define sorgbr_ sorgbr
#define sorghr_ sorghr
#define sorgl2_ sorgl2
#define sorglq_ sorglq
#define sorgql_ sorgql
#define sorgqr_ sorgqr
#define sorgr2_ sorgr2
#define sorgrq_ sorgrq
#define sorgtr_ sorgtr
#define sorm2l_ sorm2l
#define sorm2r_ sorm2r
#define sormbr_ sormbr
#define sormhr_ sormhr
#define sorml2_ sorml2
#define sormlq_ sormlq
#define sormql_ sormql
#define sormqr_ sormqr
#define sormr2_ sormr2
#define sormr3_ sormr3
#define sormrq_ sormrq
#define sormrz_ sormrz
#define sormtr_ sormtr
#define spbcon_ spbcon
#define spbequ_ spbequ
#define spbrfs_ spbrfs
#define spbstf_ spbstf
#define spbsv_ spbsv
#define spbsvx_ spbsvx
#define spbtf2_ spbtf2
#define spbtrf_ spbtrf
#define spbtrs_ spbtrs
#define spocon_ spocon
#define spoequ_ spoequ
#define sporfs_ sporfs
#define sposv_ sposv
#define sposvx_ sposvx
#define spotf2_ spotf2
#define spotrf_ spotrf
#define spotri_ spotri
#define spotrs_ spotrs
#define sppcon_ sppcon
#define sppequ_ sppequ
#define spprfs_ spprfs
#define sppsv_ sppsv
#define sppsvx_ sppsvx
#define spptrf_ spptrf
#define spptri_ spptri
#define spptrs_ spptrs
#define sptcon_ sptcon
#define spteqr_ spteqr
#define sptrfs_ sptrfs
#define sptsv_ sptsv
#define sptsvx_ sptsvx
#define spttrf_ spttrf
#define spttrs_ spttrs
#define sptts2_ sptts2
#define srscl_ srscl
#define ssbev_ ssbev
#define ssbevd_ ssbevd
#define ssbevx_ ssbevx
#define ssbgst_ ssbgst
#define ssbgv_ ssbgv
#define ssbgvd_ ssbgvd
#define ssbgvx_ ssbgvx
#define ssbtrd_ ssbtrd
#define sspcon_ sspcon
#define sspev_ sspev
#define sspevd_ sspevd
#define sspevx_ sspevx
#define sspgst_ sspgst
#define sspgv_ sspgv
#define sspgvd_ sspgvd
#define sspgvx_ sspgvx
#define ssprfs_ ssprfs
#define sspsv_ sspsv
#define sspsvx_ sspsvx
#define ssptrd_ ssptrd
#define ssptrf_ ssptrf
#define ssptri_ ssptri
#define ssptrs_ ssptrs
#define sstebz_ sstebz
#define sstedc_ sstedc
#define sstegr_ sstegr
#define sstein_ sstein
#define ssteqr_ ssteqr
#define ssterf_ ssterf
#define sstev_ sstev
#define sstevd_ sstevd
#define sstevr_ sstevr
#define sstevx_ sstevx
#define ssycon_ ssycon
#define ssyev_ ssyev
#define ssyevd_ ssyevd
#define ssyevr_ ssyevr
#define ssyevx_ ssyevx
#define ssygs2_ ssygs2
#define ssygst_ ssygst
#define ssygv_ ssygv
#define ssygvd_ ssygvd
#define ssygvx_ ssygvx
#define ssyrfs_ ssyrfs
#define ssysv_ ssysv
#define ssysvx_ ssysvx
#define ssytd2_ ssytd2
#define ssytf2_ ssytf2
#define ssytrd_ ssytrd
#define ssytrf_ ssytrf
#define ssytri_ ssytri
#define ssytrs_ ssytrs
#define stbcon_ stbcon
#define stbrfs_ stbrfs
#define stbtrs_ stbtrs
#define stgevc_ stgevc
#define stgex2_ stgex2
#define stgexc_ stgexc
#define stgsen_ stgsen
#define stgsja_ stgsja
#define stgsna_ stgsna
#define stgsy2_ stgsy2
#define stgsyl_ stgsyl
#define stpcon_ stpcon
#define stprfs_ stprfs
#define stptri_ stptri
#define stptrs_ stptrs
#define strcon_ strcon
#define strevc_ strevc
#define strexc_ strexc
#define strrfs_ strrfs
#define strsen_ strsen
#define strsna_ strsna
#define strsyl_ strsyl
#define strti2_ strti2
#define strtri_ strtri
#define strtrs_ strtrs
#define stzrqf_ stzrqf
#define stzrzf_ stzrzf
#define xerbla_ xerbla
#define zbdsqr_ zbdsqr
#define zdrot_ zdrot
#define zdrscl_ zdrscl
#define zgbbrd_ zgbbrd
#define zgbcon_ zgbcon
#define zgbequ_ zgbequ
#define zgbrfs_ zgbrfs
#define zgbsv_ zgbsv
#define zgbsvx_ zgbsvx
#define zgbtf2_ zgbtf2
#define zgbtrf_ zgbtrf
#define zgbtrs_ zgbtrs
#define zgebak_ zgebak
#define zgebal_ zgebal
#define zgebd2_ zgebd2
#define zgebrd_ zgebrd
#define zgecon_ zgecon
#define zgeequ_ zgeequ
#define zgees_ zgees
#define zgeesx_ zgeesx
#define zgeev_ zgeev
#define zgeevx_ zgeevx
#define zgegs_ zgegs
#define zgegv_ zgegv
#define zgehd2_ zgehd2
#define zgehrd_ zgehrd
#define zgelq2_ zgelq2
#define zgelqf_ zgelqf
#define zgels_ zgels
#define zgelsd_ zgelsd
#define zgelss_ zgelss
#define zgelsx_ zgelsx
#define zgelsy_ zgelsy
#define zgeql2_ zgeql2
#define zgeqlf_ zgeqlf
#define zgeqp3_ zgeqp3
#define zgeqpf_ zgeqpf
#define zgeqr2_ zgeqr2
#define zgeqrf_ zgeqrf
#define zgerfs_ zgerfs
#define zgerq2_ zgerq2
#define zgerqf_ zgerqf
#define zgesc2_ zgesc2
#define zgesdd_ zgesdd
#define zgesv_ zgesv
#define zgesvd_ zgesvd
#define zgesvx_ zgesvx
#define zgetc2_ zgetc2
#define zgetf2_ zgetf2
#define zgetrf_ zgetrf
#define zgetri_ zgetri
#define zgetrs_ zgetrs
#define zggbak_ zggbak
#define zggbal_ zggbal
#define zgges_ zgges
#define zggesx_ zggesx
#define zggev_ zggev
#define zggevx_ zggevx
#define zggglm_ zggglm
#define zgghrd_ zgghrd
#define zgglse_ zgglse
#define zggqrf_ zggqrf
#define zggrqf_ zggrqf
#define zggsvd_ zggsvd
#define zggsvp_ zggsvp
#define zgtcon_ zgtcon
#define zgtrfs_ zgtrfs
#define zgtsv_ zgtsv
#define zgtsvx_ zgtsvx
#define zgttrf_ zgttrf
#define zgttrs_ zgttrs
#define zgtts2_ zgtts2
#define zhbev_ zhbev
#define zhbevd_ zhbevd
#define zhbevx_ zhbevx
#define zhbgst_ zhbgst
#define zhbgv_ zhbgv
#define zhbgvd_ zhbgvd
#define zhbgvx_ zhbgvx
#define zhbtrd_ zhbtrd
#define zhecon_ zhecon
#define zheev_ zheev
#define zheevd_ zheevd
#define zheevr_ zheevr
#define zheevx_ zheevx
#define zhegs2_ zhegs2
#define zhegst_ zhegst
#define zhegv_ zhegv
#define zhegvd_ zhegvd
#define zhegvx_ zhegvx
#define zherfs_ zherfs
#define zhesv_ zhesv
#define zhesvx_ zhesvx
#define zhetd2_ zhetd2
#define zhetf2_ zhetf2
#define zhetrd_ zhetrd
#define zhetrf_ zhetrf
#define zhetri_ zhetri
#define zhetrs_ zhetrs
#define zhgeqz_ zhgeqz
#define zhpcon_ zhpcon
#define zhpev_ zhpev
#define zhpevd_ zhpevd
#define zhpevx_ zhpevx
#define zhpgst_ zhpgst
#define zhpgv_ zhpgv
#define zhpgvd_ zhpgvd
#define zhpgvx_ zhpgvx
#define zhprfs_ zhprfs
#define zhpsv_ zhpsv
#define zhpsvx_ zhpsvx
#define zhptrd_ zhptrd
#define zhptrf_ zhptrf
#define zhptri_ zhptri
#define zhptrs_ zhptrs
#define zhsein_ zhsein
#define zhseqr_ zhseqr
#define zlabrd_ zlabrd
#define zlacgv_ zlacgv
#define zlacon_ zlacon
#define zlacp2_ zlacp2
#define zlacpy_ zlacpy
#define zlacrm_ zlacrm
#define zlacrt_ zlacrt
#define zladiv_ zladiv
#define zlaed0_ zlaed0
#define zlaed7_ zlaed7
#define zlaed8_ zlaed8
#define zlaein_ zlaein
#define zlaesy_ zlaesy
#define zlaev2_ zlaev2
#define zlags2_ zlags2
#define zlagtm_ zlagtm
#define zlahef_ zlahef
#define zlahqr_ zlahqr
#define zlahrd_ zlahrd
#define zlaic1_ zlaic1
#define zlals0_ zlals0
#define zlalsa_ zlalsa
#define zlalsd_ zlalsd
#define zlangb_ zlangb
#define zlange_ zlange
#define zlangt_ zlangt
#define zlanhb_ zlanhb
#define zlanhe_ zlanhe
#define zlanhp_ zlanhp
#define zlanhs_ zlanhs
#define zlanht_ zlanht
#define zlansb_ zlansb
#define zlansp_ zlansp
#define zlansy_ zlansy
#define zlantb_ zlantb
#define zlantp_ zlantp
#define zlantr_ zlantr
#define zlapll_ zlapll
#define zlapmt_ zlapmt
#define zlaqgb_ zlaqgb
#define zlaqge_ zlaqge
#define zlaqhb_ zlaqhb
#define zlaqhe_ zlaqhe
#define zlaqhp_ zlaqhp
#define zlaqp2_ zlaqp2
#define zlaqps_ zlaqps
#define zlaqsb_ zlaqsb
#define zlaqsp_ zlaqsp
#define zlaqsy_ zlaqsy
#define zlar1v_ zlar1v
#define zlar2v_ zlar2v
#define zlarcm_ zlarcm
#define zlarf_ zlarf
#define zlarfb_ zlarfb
#define zlarfg_ zlarfg
#define zlarft_ zlarft
#define zlarfx_ zlarfx
#define zlargv_ zlargv
#define zlarnv_ zlarnv
#define zlarrv_ zlarrv
#define zlartg_ zlartg
#define zlartv_ zlartv
#define zlarz_ zlarz
#define zlarzb_ zlarzb
#define zlarzt_ zlarzt
#define zlascl_ zlascl
#define zlaset_ zlaset
#define zlasr_ zlasr
#define zlassq_ zlassq
#define zlaswp_ zlaswp
#define zlasyf_ zlasyf
#define zlatbs_ zlatbs
#define zlatdf_ zlatdf
#define zlatps_ zlatps
#define zlatrd_ zlatrd
#define zlatrs_ zlatrs
#define zlatrz_ zlatrz
#define zlatzm_ zlatzm
#define zlauu2_ zlauu2
#define zlauum_ zlauum
#define zpbcon_ zpbcon
#define zpbequ_ zpbequ
#define zpbrfs_ zpbrfs
#define zpbstf_ zpbstf
#define zpbsv_ zpbsv
#define zpbsvx_ zpbsvx
#define zpbtf2_ zpbtf2
#define zpbtrf_ zpbtrf
#define zpbtrs_ zpbtrs
#define zpocon_ zpocon
#define zpoequ_ zpoequ
#define zporfs_ zporfs
#define zposv_ zposv
#define zposvx_ zposvx
#define zpotf2_ zpotf2
#define zpotrf_ zpotrf
#define zpotri_ zpotri
#define zpotrs_ zpotrs
#define zppcon_ zppcon
#define zppequ_ zppequ
#define zpprfs_ zpprfs
#define zppsv_ zppsv
#define zppsvx_ zppsvx
#define zpptrf_ zpptrf
#define zpptri_ zpptri
#define zpptrs_ zpptrs
#define zptcon_ zptcon
#define zpteqr_ zpteqr
#define zptrfs_ zptrfs
#define zptsv_ zptsv
#define zptsvx_ zptsvx
#define zpttrf_ zpttrf
#define zpttrs_ zpttrs
#define zptts2_ zptts2
#define zrot_ zrot
#define zspcon_ zspcon
#define zspmv_ zspmv
#define zspr_ zspr
#define zsprfs_ zsprfs
#define zspsv_ zspsv
#define zspsvx_ zspsvx
#define zsptrf_ zsptrf
#define zsptri_ zsptri
#define zsptrs_ zsptrs
#define zstedc_ zstedc
#define zstegr_ zstegr
#define zstein_ zstein
#define zsteqr_ zsteqr
#define zsycon_ zsycon
#define zsymv_ zsymv
#define zsyr_ zsyr
#define zsyrfs_ zsyrfs
#define zsysv_ zsysv
#define zsysvx_ zsysvx
#define zsytf2_ zsytf2
#define zsytrf_ zsytrf
#define zsytri_ zsytri
#define zsytrs_ zsytrs
#define ztbcon_ ztbcon
#define ztbrfs_ ztbrfs
#define ztbtrs_ ztbtrs
#define ztgevc_ ztgevc
#define ztgex2_ ztgex2
#define ztgexc_ ztgexc
#define ztgsen_ ztgsen
#define ztgsja_ ztgsja
#define ztgsna_ ztgsna
#define ztgsy2_ ztgsy2
#define ztgsyl_ ztgsyl
#define ztpcon_ ztpcon
#define ztprfs_ ztprfs
#define ztptri_ ztptri
#define ztptrs_ ztptrs
#define ztrcon_ ztrcon
#define ztrevc_ ztrevc
#define ztrexc_ ztrexc
#define ztrrfs_ ztrrfs
#define ztrsen_ ztrsen
#define ztrsna_ ztrsna
#define ztrsyl_ ztrsyl
#define ztrti2_ ztrti2
#define ztrtri_ ztrtri
#define ztrtrs_ ztrtrs
#define ztzrqf_ ztzrqf
#define ztzrzf_ ztzrzf
#define zung2l_ zung2l
#define zung2r_ zung2r
#define zungbr_ zungbr
#define zunghr_ zunghr
#define zungl2_ zungl2
#define zunglq_ zunglq
#define zungql_ zungql
#define zungqr_ zungqr
#define zungr2_ zungr2
#define zungrq_ zungrq
#define zungtr_ zungtr
#define zunm2l_ zunm2l
#define zunm2r_ zunm2r
#define zunmbr_ zunmbr
#define zunmhr_ zunmhr
#define zunml2_ zunml2
#define zunmlq_ zunmlq
#define zunmql_ zunmql
#define zunmqr_ zunmqr
#define zunmr2_ zunmr2
#define zunmr3_ zunmr3
#define zunmrq_ zunmrq
#define zunmrz_ zunmrz
#define zunmtr_ zunmtr
#define zupgtr_ zupgtr
#define zupmtr_ zupmtr
#endif /* end of ifdef MATLAB_MEX_FILE */

/*Define the namespace*/
namespace ROPTLIB{
//#include <dasum.h>
//#include <daxpy.h>
	void axpy_(integer *n, complex *ca, complex *cx, integer *incx, complex *cy, integer *incy);
	void axpy_(integer *n, doublereal *da, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
	void axpy_(integer *n, real *sa, real *sx, integer *incx, real *sy, integer *incy);
	void axpy_(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);

//#include <dcabs1.h>
//#include <dcopy.h>
	void copy_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy);
	void copy_(integer *n, doublereal *cx, integer *incx, doublereal *cy, integer *incy);
	void copy_(integer *n, real *cx, integer *incx, real *cy, integer *incy);
	void copy_(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy);

//#include <ddot.h>
	doublereal dot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
	real dot_(integer *n, real *dx, integer *incx, real *dy, integer *incy);
	void dotu_(complex * ret_val, integer *n, complex *cx, integer *incx, complex *cy, integer *incy);
	void dotu_(doublecomplex * ret_val, integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy);
	void dotc_(complex * ret_val, integer *n, complex *cx, integer *incx, complex *cy, integer *incy);
	void dotc_(doublecomplex * ret_val, integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy);

//#include <dgbmv.h>
//#include <dgemm.h>
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, complex *alpha, complex *a, integer *lda, complex *b, integer *ldb, complex *beta, complex *c__, integer *ldc);
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc);
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta, real *c__, integer *ldc);
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *ldc);

//#include <dgemv.h>
	void gemv_(char *trans, integer *m, integer *n, complex *alpha, complex *a, integer *lda, complex *x, integer *incx, complex *beta, complex *y, integer *incy);
	void gemv_(char *trans, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);
	void gemv_(char *trans, integer *m, integer *n, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy);
	void gemv_(char *trans, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);


//#include <dger.h>
	void ger_(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *a, integer *lda);
	void ger_(integer *m, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy, real *a, integer *lda);

	void gerc_(integer *m, integer *n, complex *alpha, complex *x, integer *incx, complex *y, integer *incy, complex *a, integer *lda);
	void gerc_(integer *m, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda);

	void ger_(integer *m, integer *n, complex *alpha, complex *x, integer *incx, complex *y, integer *incy, complex *a, integer *lda);
	void ger_(integer *m, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda);

//#include <dnrm2.h>
	doublereal nrm2_(integer *n, doublereal *x, integer *incx);
	real nrm2_(integer *n, real *x, integer *incx);

//#include <drot.h>
//#include <drotg.h>
//#include <dsbmv.h>
//#include <dscal.h>
	void scal_(integer *n, complex *ca, complex *cx, integer *incx);
	void scal_(integer *n, doublereal *ca, doublereal *cx, integer *incx);
	void scal_(integer *n, real *ca, real *cx, integer *incx);
	void scal_(integer *n, doublecomplex *ca, doublecomplex *cx, integer *incx);

//#include <dspmv.h>
//#include <dspr.h>
//#include <dspr2.h>
//#include <dswap.h>
//#include <dsymm.h>
//#include <dsymv.h>
	void symv_(char *uplo, integer *n, complex *alpha, complex *a, integer *lda, complex *x, integer *incx, complex *beta, complex *y, integer *incy);
	void symv_(char *uplo, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);
	void symv_(char *uplo, integer *n, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy);
	void symv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

//#include <dsyr.h>
//#include <dsyr2.h>
//#include <dsyr2k.h>
//#include <dsyrk.h>
//#include <dtbmv.h>
//#include <dtbsv.h>
//#include <dtpmv.h>
//#include <dtpsv.h>
//#include <dtrmm.h>
//#include <dtrmv.h>
//#include <dtrsm.h>
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb);
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, integer *ldb);
    void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);
    void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, complex *alpha, complex *a, integer *lda, complex *b, integer *ldb);
//#include <dtrsv.h>
//#include <dzasum.h>
//#include <dznrm2.h>



//#include <dbdsdc.h>
//#include <dbdsqr.h>
//#include <ddisna.h>
//#include <dgbbrd.h>
//#include <dgbcon.h>
//#include <dgbequ.h>
//#include <dgbrfs.h>
//#include <dgbsv.h>
//#include <dgbsvx.h>
//#include <dgbtf2.h>
//#include <dgbtrf.h>
//#include <dgbtrs.h>
//#include <dgebak.h>
//#include <dgebal.h>
//#include <dgebd2.h>
//#include <dgebrd.h>
//#include <dgecon.h>
//#include <dgeequ.h>
//#include <dgees.h>
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, complex *a, integer *lda, integer *sdim, complex *w, complex *vs, integer *ldvs, complex *work, integer *lwork, real *rwork, logical *bwork, integer *info);
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, doublereal *a, integer *lda, integer *sdim, doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, integer *lwork, logical *bwork, integer *info);
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, integer *ldvs, real *work, integer *lwork, logical *bwork, integer *info);
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

//#include <dgeesx.h>
//#include <dgeev.h>
//#include <dgeevx.h>
//#include <dgegs.h>
//#include <dgegv.h>
//#include <dgehd2.h>
//#include <dgehrd.h>
//#include <dgelq2.h>
//#include <dgelqf.h>
//#include <dgels.h>
//#include <dgelsd.h>
//#include <dgelss.h>
//#include <dgelsx.h>
//#include <dgelsy.h>
//#include <dgeql2.h>
//#include <dgeqlf.h>
//#include <dgeqp3.h>
	void geqp3_(integer *m, integer *n, complex *a, integer *lda, integer *jpvt, complex *tau, complex *work, integer *lwork, real *rwork, integer *info);
	void geqp3_(integer *m, integer *n, doublereal *a, integer *lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork, integer *info);
	void geqp3_(integer *m, integer *n, real *a, integer *lda, integer *jpvt, real *tau, real *work, integer *lwork, integer *info);
	void geqp3_(integer *m, integer *n, doublecomplex *a, integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);
//#include <dgeqpf.h>
//#include <dgeqr2.h>
//#include <dgeqrf.h>
//#include <dgerfs.h>
//#include <dgerq2.h>
//#include <dgerqf.h>
//#include <dgesc2.h>
//#include <dgesdd.h>
	void gesdd_(char *jobz, integer *m, integer *n, complex *a, integer *lda, real *s, complex *u, integer *ldu, complex *vt, integer *ldvt, complex *work, integer *lwork, real *rwork, integer *iwork, integer *info);
	void gesdd_(char *jobz, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *iwork, integer *info);
	void gesdd_(char *jobz, integer *m, integer *n, real *a, integer *lda, real *s, real *u, integer *ldu, real *vt, integer *ldvt, real *work, integer *lwork, integer *iwork, integer *info);
	void gesdd_(char *jobz, integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, integer *lwork, doublereal *rwork, integer *iwork, integer *info);
//#include <dgesv.h>
	void gesv_(integer *n, integer *nrhs, complex *a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer *info);
	void gesv_(integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
	void gesv_(integer *n, integer *nrhs, real *a, integer *lda, integer *ipiv, real *b, integer *ldb, integer *info);
	void gesv_(integer *n, integer *nrhs, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *info);

//#include <dgesvd.h>
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, complex *a, integer *lda, real *s, complex *u, integer *ldu, complex *vt, integer *ldvt, complex *work, integer *lwork, real *rwork, integer *info);
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *info);
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, real *a, integer *lda, real *s, real *u, integer *ldu, real *vt, integer *ldvt, real *work, integer *lwork, integer *info);
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

//#include <dgesvx.h>
//#include <dgetc2.h>
//#include <dgetf2.h>
//#include <dgetrf.h>
	void getrf_(integer *m, integer *n, complex *a, integer *lda, integer *ipiv, integer *info);
	void getrf_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);
	void getrf_(integer *m, integer *n, real *a, integer *lda, integer *ipiv, integer *info);
	void getrf_(integer *m, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info);

//#include <dgetri.h>
	void getri_(integer *n, complex *a, integer *lda, integer *ipiv, complex *work, integer *lwork, integer *info);
	void getri_(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);
	void getri_(integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *lwork, integer *info);
	void getri_(integer *n, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, integer *info);

//#include <dgetrs.h>
	void getrs_(char *trans, integer *n, integer *nrhs, complex *a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer *info);
	void getrs_(char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
	void getrs_(char *trans, integer *n, integer *nrhs, real *a, integer *lda, integer *ipiv, real *b, integer *ldb, integer *info);
	void getrs_(char *trans, integer *n, integer *nrhs, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *info);
//#include <dggbak.h>
//#include <dggbal.h>
//#include <dgges.h>
//#include <dggesx.h>
//#include <dggev.h>
//#include <dggevx.h>
//#include <dggglm.h>
//#include <dgghrd.h>
//#include <dgglse.h>
//#include <dggqrf.h>
//#include <dggrqf.h>
//#include <dggsvd.h>
//#include <dggsvp.h>
//#include <dgtcon.h>
//#include <dgtrfs.h>
//#include <dgtsv.h>
//#include <dgtsvx.h>
//#include <dgttrf.h>
//#include <dgttrs.h>
//#include <dgtts2.h>
//#include <dhgeqz.h>
//#include <dhsein.h>
//#include <dhseqr.h>
//#include <dlabad.h>
//#include <dlabrd.h>
//#include <dlacon.h>
//#include <dlacpy.h>
//#include <dladiv.h>
//#include <dlae2.h>
//#include <dlaebz.h>
//#include <dlaed0.h>
//#include <dlaed1.h>
//#include <dlaed2.h>
//#include <dlaed3.h>
//#include <dlaed4.h>
//#include <dlaed5.h>
//#include <dlaed6.h>
//#include <dlaed7.h>
//#include <dlaed8.h>
//#include <dlaed9.h>
//#include <dlaeda.h>
//#include <dlaein.h>
//#include <dlaev2.h>
//#include <dlaexc.h>
//#include <dlag2.h>
//#include <dlags2.h>
//#include <dlagtf.h>
//#include <dlagtm.h>
//#include <dlagts.h>
//#include <dlagv2.h>
//#include <dlahqr.h>
//#include <dlahrd.h>
//#include <dlaic1.h>
//#include <dlaln2.h>
//#include <dlals0.h>
//#include <dlalsa.h>
//#include <dlalsd.h>
//#include <dlamch.h>
//#include <dlamrg.h>
//#include <dlangb.h>
//#include <dlange.h>
//#include <dlangt.h>
//#include <dlanhs.h>
//#include <dlansb.h>
//#include <dlansp.h>
//#include <dlanst.h>
//#include <dlansy.h>
//#include <dlantb.h>
//#include <dlantp.h>
//#include <dlantr.h>
//#include <dlanv2.h>
//#include <dlapll.h>
//#include <dlapmt.h>
	void lapmt_(logical *forwrd, integer *m, integer *n, complex *x, integer *ldx, integer *k);
	void lapmt_(logical *forwrd, integer *m, integer *n, doublereal *x, integer *ldx, integer *k);
	void lapmt_(logical *forwrd, integer *m, integer *n, real *x, integer *ldx, integer *k);
	void lapmt_(logical *forwrd, integer *m, integer *n, doublecomplex *x, integer *ldx, integer *k);

//#include <dlapy2.h>
//#include <dlapy3.h>
//#include <dlaqgb.h>
//#include <dlaqge.h>
//#include <dlaqp2.h>
//#include <dlaqps.h>
//#include <dlaqsb.h>
//#include <dlaqsp.h>
//#include <dlaqsy.h>
//#include <dlaqtr.h>
//#include <dlar1v.h>
//#include <dlar2v.h>
//#include <dlarf.h>
//#include <dlarfb.h>
//#include <dlarfg.h>
//#include <dlarft.h>
//#include <dlarfx.h>
	void larfx_(char *side, integer *m, integer *n, complex *v, complex *tau, complex *c__, integer *ldc, complex *work);
	void larfx_(char *side, integer *m, integer *n, doublereal *v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);
	void larfx_(char *side, integer *m, integer *n, real *v, real *tau, real *c__, integer *ldc, real *work);
	void larfx_(char *side, integer *m, integer *n, doublecomplex *v, doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *work);
//#include <dlargv.h>
//#include <dlarnv.h>
//#include <dlarrb.h>
//#include <dlarre.h>
//#include <dlarrf.h>
//#include <dlarrv.h>
//#include <dlartg.h>
//#include <dlartv.h>
//#include <dlaruv.h>
//#include <dlarz.h>
//#include <dlarzb.h>
//#include <dlarzt.h>
//#include <dlas2.h>
//#include <dlascl.h>
//#include <dlasd0.h>
//#include <dlasd1.h>
//#include <dlasd2.h>
//#include <dlasd3.h>
//#include <dlasd4.h>
//#include <dlasd5.h>
//#include <dlasd6.h>
//#include <dlasd7.h>
//#include <dlasd8.h>
//#include <dlasd9.h>
//#include <dlasda.h>
//#include <dlasdq.h>
//#include <dlasdt.h>
//#include <dlaset.h>
//#include <dlasq1.h>
//#include <dlasq2.h>
//#include <dlasq3.h>
//#include <dlasq4.h>
//#include <dlasq5.h>
//#include <dlasq6.h>
//#include <dlasr.h>
//#include <dlasrt.h>
//#include <dlassq.h>
//#include <dlasv2.h>
//#include <dlaswp.h>
//#include <dlasy2.h>
//#include <dlasyf.h>
//#include <dlatbs.h>
//#include <dlatdf.h>
//#include <dlatps.h>
//#include <dlatrd.h>
//#include <dlatrs.h>
//#include <dlatrz.h>
//#include <dlatzm.h>
//#include <dlauu2.h>
//#include <dlauum.h>
//#include <dopgtr.h>
//#include <dopmtr.h>
//#include <dorg2l.h>
//#include <dorg2r.h>
//#include <dorgbr.h>
//#include <dorghr.h>
//#include <dorgl2.h>
//#include <dorglq.h>
//#include <dorgql.h>
//#include <dorgqr.h>
	void orgqr_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
	void orgqr_(integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info);

//#include <dorgr2.h>
//#include <dorgrq.h>
//#include <dorgtr.h>
//#include <dorm2l.h>
//#include <dorm2r.h>
//#include <dormbr.h>
//#include <dormhr.h>
//#include <dorml2.h>
//#include <dormlq.h>
//#include <dormql.h>
//#include <dormqr.h>
	void ormqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
	void ormqr_(char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, real *work, integer *lwork, integer *info);

//#include <dormr2.h>
//#include <dormr3.h>
//#include <dormrq.h>
//#include <dormrz.h>
//#include <dormtr.h>
//#include <dpbcon.h>
//#include <dpbequ.h>
//#include <dpbrfs.h>
//#include <dpbstf.h>
//#include <dpbsv.h>
//#include <dpbsvx.h>
//#include <dpbtf2.h>
//#include <dpbtrf.h>
//#include <dpbtrs.h>
//#include <dpocon.h>
//#include <dpoequ.h>
//#include <dporfs.h>
//#include <dposv.h>
//#include <dposvx.h>
//#include <dpotf2.h>
//#include <dpotrf.h>
	void potrf_(char *uplo, integer *n, complex *a, integer *lda, integer *info);
	void potrf_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info);
	void potrf_(char *uplo, integer *n, real *a, integer *lda, integer *info);
	void potrf_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info);

//#include <dpotri.h>
//#include <dpotrs.h>
	void potrs_(char *uplo, integer *n, integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, integer *info);
	void potrs_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *info);
	void potrs_(char *uplo, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer *info);
	void potrs_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *info);

//#include <dppcon.h>
//#include <dppequ.h>
//#include <dpprfs.h>
//#include <dppsv.h>
//#include <dppsvx.h>
//#include <dpptrf.h>
//#include <dpptri.h>
//#include <dpptrs.h>
//#include <dptcon.h>
//#include <dpteqr.h>
//#include <dptrfs.h>
//#include <dptsv.h>
//#include <dptsvx.h>
//#include <dpttrf.h>
//#include <dpttrs.h>
//#include <dptts2.h>
//#include <drscl.h>
//#include <dsbev.h>
//#include <dsbevd.h>
//#include <dsbevx.h>
//#include <dsbgst.h>
//#include <dsbgv.h>
//#include <dsbgvd.h>
//#include <dsbgvx.h>
//#include <dsbtrd.h>
//#include <dsecnd.h>
//#include <dspcon.h>
//#include <dspev.h>
//#include <dspevd.h>
//#include <dspevx.h>
//#include <dspgst.h>
//#include <dspgv.h>
//#include <dspgvd.h>
//#include <dspgvx.h>
//#include <dsprfs.h>
//#include <dspsv.h>
//#include <dspsvx.h>
//#include <dsptrd.h>
//#include <dsptrf.h>
//#include <dsptri.h>
//#include <dsptrs.h>
//#include <dstebz.h>
//#include <dstedc.h>
//#include <dstegr.h>
//#include <dstein.h>
//#include <dsteqr.h>
//#include <dsterf.h>
//#include <dstev.h>
//#include <dstevd.h>
//#include <dstevr.h>
//#include <dstevx.h>
//#include <dsycon.h>
//#include <dsyev.h>
	void syev_(char *jobz, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *w, doublereal *work, integer *lwork, integer *info);
	void syev_(char *jobz, char *uplo, integer *n, real *a, integer *lda, real *w, real *work, integer *lwork, integer *info);

//#include <dsyevd.h>
	void syevd_(char *jobz, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *w, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
	void syevd_(char *jobz, char *uplo, integer *n, real *a, integer *lda, real *w, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
//#include <dsyevr.h>
	void syevr_(char *jobz, char *range, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
	void syevr_(char *jobz, char *range, char *uplo, integer *n, real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, real *z__, integer *ldz, integer *isuppz, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
//#include <dsyevx.h>
	void syevx_(char *jobz, char *range, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *ifail, integer *info);
	void syevx_(char *jobz, char *range, char *uplo, integer *n, real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, real *z__, integer *ldz, real *work, integer *lwork, integer *iwork, integer *ifail, integer *info);
//#include <dsygs2.h>
//#include <dsygst.h>
//#include <dsygv.h>
//#include <dsygvd.h>
//#include <dsygvx.h>
//#include <dsyrfs.h>
//#include <dsysv.h>
//#include <dsysvx.h>
//#include <dsytd2.h>
//#include <dsytf2.h>
//#include <dsytrd.h>
//#include <dsytrf.h>
//#include <dsytri.h>
//#include <dsytrs.h>
//#include <dtbcon.h>
//#include <dtbrfs.h>
//#include <dtbtrs.h>
//#include <dtgevc.h>
//#include <dtgex2.h>
//#include <dtgexc.h>
//#include <dtgsen.h>
//#include <dtgsja.h>
//#include <dtgsna.h>
//#include <dtgsy2.h>
//#include <dtgsyl.h>
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *c__, integer *ldc, complex *d__, integer *ldd, complex *e, integer *lde, complex *f, integer *ldf, real *scale, real *dif, complex *work, integer *lwork, integer *iwork, integer *info);
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *scale, doublereal *dif, doublereal *work, integer *lwork, integer *iwork, integer *info);
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, real *a, integer *lda, real *b, integer *ldb, real *c__, integer *ldc, real *d__, integer *ldd, real *e, integer *lde, real *f, integer *ldf, real *scale, real *dif, real *work, integer *lwork, integer *iwork, integer *info);
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, doublereal *scale, doublereal *dif, doublecomplex *work, integer *lwork, integer *iwork, integer *info);
//#include <dtpcon.h>
//#include <dtprfs.h>
//#include <dtptri.h>
//#include <dtptrs.h>
//#include <dtrcon.h>
//#include <dtrevc.h>
//#include <dtrexc.h>
//#include <dtrrfs.h>
//#include <dtrsen.h>
//#include <dtrsna.h>
//#include <dtrsyl.h>
//#include <dtrti2.h>
//#include <dtrtri.h>
//#include <dtrtrs.h>
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, integer *info);
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *info);
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer *info);
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *info);
//#include <dtzrqf.h>
//#include <dtzrzf.h>
//#include <dzsum1.h>

	//zgegs_
	void gegs_(char *jobvsl, char *jobvsr, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, complex *work, integer *lwork, real *rwork, integer *info);
	void gegs_(char *jobvsl, char *jobvsr, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

	//zunmqr_
	void unmqr_(char *side, char *trans, integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info);
	void unmqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, integer *info);

	//zpotri_
	void potri_(char *uplo, integer *n, complex *a, integer *lda, integer *info);
	void potri_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info);

    //zheev_
    void heev_(char *jobz, char *uplo, integer *n, complex *A, integer *LDA, real *W, complex *work, integer *lwork, real *rwork, integer *info);
    void heev_(char *jobz, char *uplo, integer *n, doublecomplex *A, integer *LDA, doublereal *W, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

    //zungqr_
    void ungqr_(integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *work, integer *lwork, integer *info);
    void ungqr_(integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info);

	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, int nrhs, float alpha, blas_sparse_matrix A, const float *b, int ldb, float *c, int ldc);
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa,	int nrhs, double alpha, blas_sparse_matrix A, const double *b, int ldb, double *c, int ldc);
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa,	int nrhs, const complex *alpha, blas_sparse_matrix A, const complex *b,	int ldb, complex *c, int ldc);
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa,	int nrhs, const doublecomplex *alpha, blas_sparse_matrix A, const doublecomplex *b,	int ldb, doublecomplex *c, int ldc);

//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, float alpha, blas_sparse_matrix A, const float *b, integer ldb, float *c, integer ldc);
//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, double alpha, blas_sparse_matrix A, const double *b, integer ldb, double *c, integer ldc);
//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, const complex *alpha, blas_sparse_matrix A, const complex *b, integer ldb, complex *c, integer ldc);
//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, const doublecomplex *alpha, blas_sparse_matrix A, const doublecomplex *b, integer ldb, doublecomplex *c, integer ldc);


	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const float *val, const integer *indx, const integer *jndx);
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const double *val, const integer *indx, const integer *jndx);
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const complex *val, const integer *indx, const integer *jndx);
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const doublecomplex *val, const integer *indx, const integer *jndx);


//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const float *val, const int *indx, const int *jndx);
//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const double *val, const int *indx, const int *jndx);
//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const complex *val, const int *indx, const int *jndx);
//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const doublecomplex *val, const int *indx, const int *jndx);

}; /*end of ROPTLIB namespace*/

#endif
