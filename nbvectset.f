      subroutine nbvectset(x,y,z)
c
c     setup for vectorized routines
c     currently for nonperiodic only
*     assumes vector.and.nbmeth.eq.'SPLINE'
c
      implicit none
      include '../includes/maxparm.dc'
      include '../includes/eopt.dc'
      include '../includes/nbond.dc'
      include '../includes/gparm.dc'
!      include '../includes/periodc.xtl'  !natmcel
      include '../includes/machine.dc'
      include '../includes/ener.dc'  !lsetchrg,lupchrg
      include '../includes/nbchkcm.dc'  !ismall
      include '../includes/edrnb.dc'    !c2onnb
      include '../includes/misc.dc'    !coulen
!      include '../includes/atodat.dc'
!      include '../includes/param.dc'
!      include '../includes/pertth.dc'
!      include '../includes/shell.dc'
!      include '../includes/nonbond.xtl'  !inbx
c
************************************************************************
*arguments
*       type          argument         in/out     description
*-----------------------------------------------------------------------
        real          x(*)
        real          y(*)
        real          z(*)
c
c90     double precision da(1)
c90      real ra(1)
c90      integer ia(1)
c90      equivalence (ia,ra,da)
c
      integer junk
c
c     ----- memory offsets -----
c
c90      integer locia
c90      integer locra
c     90      integer locda
      c
c     vectorized section (nonperiodic), use nblst4
c
      mxsw=max(mxsw,1)
      mxrn=max(mxrn,1)
      mxicrn=max(mxicrn,1)
c
c           NBLST4
c
c90 -- the following arrays are local so they could be dynamically allocated here and deallocated at the end of the routine
c90 for now they are just allocated with fixed dimension
     
      integer jsave(mxatpr)
      integer ijsave(mxatpr)
      integer jnbrn(mxrn)
      integer nnbrn(mxrn)
      integer icptrn(mxrn)
      integer icnbrn(mxicrn)
      integer inbsw(mxsw)
      integer jnbsw(mxsw)
      integer icnbsw(mxsw)
      integer icnb(mxjnpa)
      integer irnlo(mxatpr+1)
      integer iswlo(mxatpr+1)
      integer nicrn
      integer cnbain(mxjnpa)
      integer cnbbin(mxjnpa)
      integer cnbasw(mxsw)
      integer cnbbsw(mxsw)
      integer cnbarn(mxicrn)
      integer cnbbrn(mxicrn)
c
      logical memerr
c
      integer loc
      integer memsys
c
c     ----- location of memory anchors -----
c
c90      locia=loc(ia)
c90      locra=loc(ra)
c90      locda=loc(da)
*
 1    continue
c
c     ensure that memory is not already allocated -----
c
c90      if (.not.allnb.or.lnblim.or.
c90     $        memsys('pointer inbsw', 0,inbsw, locia).ne.0.or.
c90     $        memsys('pointer jnbsw', 0,jnbsw, locia).ne.0.or.
c90     $        memsys('pointer icnbsw',0,icnbsw,locia).ne.0.or.
c90     $        memsys('pointer icnb',  0,icnb,  locia).ne.0.or.
c90     $        memsys('pointer iswlo', 0,iswlo, locia).ne.0.or.
c90     $        memsys('pointer icnbrn',0,icnbrn,locia).ne.0.or.
c90     $        memsys('pointer jnbrn', 0,jnbrn, locia).ne.0.or.
c90     $        memsys('pointer nnbrn', 0,nnbrn, locia).ne.0.or.
c90     $        memsys('pointer icptrn',0,icptrn,locia).ne.0.or.
c90     $        memsys('pointer irnlo', 0,irnlo, locia).ne.0.or.
c90     $        memsys('pointer cnbain',0,cnbain,locra).ne.0.or.
c90     $        memsys('pointer cnbbin',0,cnbbin,locra).ne.0.or.
c90     $        memsys('pointer cnbasw',0,cnbasw,locra).ne.0.or.
c90     $        memsys('pointer cnbbsw',0,cnbbsw,locra).ne.0.or.
c90     $        memsys('pointer cnbarn',0,cnbarn,locra).ne.0.or.
c90     $        memsys('pointer cnbbrn',0,cnbbrn,locra).ne.0) then
c90        call nbfree
c
c       Allocate memory -----
c
c90        if (memsys('create integer nbjsave',nmovatm,jsave,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer nbijsave',nmovatm,ijsave,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer jnbrn',mxrn,jnbrn,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer nnbrn',mxrn,nnbrn,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer icptrn',mxrn,icptrn,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer icnbrn',mxicrn,icnbrn,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer inbsw',mxsw,inbsw,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer jnbsw',mxsw,jnbsw,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer icnbsw',mxsw,icnbsw,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer icnb',MXJNPA,icnb,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer irnlo',nmovatm+1,irnlo,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create integer iswlo',nmovatm+1,iswlo,locia)
c90     $           .ne.0) go to 9000
c90        if (memsys('create real cnbain',MXJNPA,cnbain,locra)
c90     $           .ne.0) go to 9000
c90        if (memsys('create real cnbbin',MXJNPA,cnbbin,locra)
c90     $           .ne.0) go to 9000
c90        if (memsys('create real cnbasw',mxsw,cnbasw,locra)
c90     $           .ne.0) go to 9000
c90        if (memsys('create real cnbbsw',mxsw,cnbbsw,locra)
c90     $           .ne.0) go to 9000
c90        if (memsys('create real cnbarn',mxicrn,cnbarn,locra)
c90     $           .ne.0) go to 9000
c90        if (memsys('create real cnbbrn',mxicrn,cnbbrn,locra)
c90     $           .ne.0) go to 9000
c
        ctnbin=ctonnb-(cutnb-ctofnb)
c
        call nblst4(x,y,z,nmovatm,jsave,ijsave,
     $           jnbrn,nnbrn,icptrn,icnbrn,
     $           inbsw,jnbsw,icnbsw,
     $           icnb,irnlo,iswlo,memerr,nicrn,
     $           cnbain,cnbbin,cnbasw,cnbbsw,
     $           cnbarn,cnbbrn)
c
c90        junk=memsys('free nbjsave',0,0,0)
c90        junk=memsys('free nbijsave',0,0,0)
c
c90        if (memerr) then
c90              call nbfree
c90               mxrn=nnnbrn*1.1
c90               mxicrn=nicrn*1.1
c90               mxsw=nnnbsw*1.1
c90               go to 1
c90        end if
c        finished vectorized non periodic case
c      endif
c
*     ON SUCCESSFUL CONCLUSION NBLST4 SETS NBFIN=.TRUE.
*     TO INDICATE THAT NBLIST HAS BEEN CONSTRUCTED
*     OTHERWISE NBFIN=.FALSE.
*
      return
c
 9000 continue
      call error('NBLIST: Can not create data arrays for NBLST4 '//
     $     '(VECTOR)')
c90      call nbfree
      call die
      return
c
      end
