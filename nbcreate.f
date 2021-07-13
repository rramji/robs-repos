      SUBROUTINE nbcreate
czmc   this routine allocate dynamic allocation for hugh arrays
c     in COMMON blocks, NOTICE that arrays are NOT LOCAL!
c     this should be called whenever energy expression is setup
      IMPLICIT NONE
      include '../includes/maxparm.dc'
      include '../includes/gparm.dc'
      include '../includes/nbond.dc'
      include '../includes/pertth.dc'
      include '../includes/vibrat.dc'
      include '../includes/phon.xtl'
 
c
      integer memsys,memsys3
      integer loc
      integer locra
      integer junk
      integer locia
      integer ntem
czmc
      vdwadd = 0
      elcadd = 0
      jnbadd = 0
      jnbptadd = 0
      jnbceladd = 0
!      call gsinf1(' Initial run-time memory allocation for NB listing')
!!      ntem = MXJNPA
czmc  initial assignment
      NBALOC = MXJNPA
cmem   kludge to skip memory allocation for nbonds for now  
c     mem
      if( ntem .lt. 100000000) then
         return
      endif

czmc 2 tc,Wag,BDO
c     nbond limit is not n(n-1)/2 for ewald and defcel case
c     while JNB used to be in common where data-overrun will
c     result in crash, JNB in memsys will result in UNcertain
c     results!! the operation of twice in memsys may be used
c     to adjust the limit.  DEC,1992
      if (ntem .gt. 0 ) then
c     1)
         locra = loc(vdwpr)
cmem         junk = memsys('free vdwpr.nbond',0,0,0)
cmem         if (memsys('create real vdwpr.nbond',ntem,vdwadd,locra)
cmem     1        .ne. 0) then
cmem            call error(
cmem     $           'Warning In Nbcreate: VDWPR array is not created')
cmem            vdwadd = 0
         end if
c     2)
         locra = loc(elcpr)
         junk = memsys('free elcpr.nbond',0,0,0)
         if (memsys('create real elcpr.nbond',ntem,elcadd,locra)
     1        .ne. 0) then
            call error(
     $           'Warning In Nbcreate: ELCPR array is not created')
            elcadd = 0
         end if
c     3)
         locia = loc(jnb)
         junk = memsys('free jnb.nbond',0,0,0)
         if (memsys('create integer jnb.nbond',ntem,jnbadd,locia)
     1        .ne. 0) then
            call error
     $           ('Warning In Nbcreate: ELCPR array is not created')
            jnbadd = 0
         end if
c     4)
         if (PRTTHRM) then
            locia = loc(jnbprt)
            junk = memsys('free jnbprt.nbond',0,0,0)
            if (memsys('create integer jnbprt.nbond',ntem,jnbptadd
     $           ,locia).ne. 0) then
               call error
     $         ('Warning In Nbcreate: JNBPRT array is not created')
               jnbptadd = 0
               PRTTHRM = .false.
            end if
         end if
c     5)
         locia = loc(jnbcel)
         junk = memsys('free jnbcel.nbond',0,0,0)
         if (memsys('create integer jnbcel.nbond',ntem,jnbceladd,locia)
     1        .ne. 0) then
            call error(
     $           'Warning In Nbcreate: JNBCEL array is not created')
            jnbceladd = 0
         end if
      end if
cmem   end of kludge to skip memory allocation for nbonds for now  
cmem      

c     6)  mxvibmod may be changed in run-time
czmc 2 zmc : this is not space saving measure
c     the best way is to allocate DDM whenever NSECD = .false.
c     and free them NSECD = .true.
      ddmadd = 0
      locia = loc(ddm(1))
      junk = memsys('free ddm.vibrat',0,0,0)
      if (memsys('create double precision ddm.vibrat',
     1     MXVIBSYM,ddmadd,locia) .ne. 0) then
         call error('Warning In Nbcreate: DDM array is not created')
         call error(
     $'Warning In Nbcreate: Second Deriv will not work')
         ddmadd = 0
      end if
      if (vdwadd .le.0) call error(
     $'Warning In Nbcreate: VDWPR is not created in energy expression')
      if (elcadd .le.0) call error(
     $     'Warning In Nbcreate: ELCPR is not created in EEX')
czmc  failure to create VDWPR and ELCPR can only screw up analyze
      if (VDWADD .le. 0 .or. ELCADD .le.0 ) call error(
     $     'Warning In Nbcreate: Energy-Force in Analyze module'//
     $' may not work')
c
      if (jnbadd .le.0 .or. jnbceladd .le.0) then
czmc  can not do nbond listing, fatal for dynamics
         call error('Warning In Nbcreate: non-bond listing is not
     $        created in energy expression')
         leeset = .false.
      end if
      if (ddmadd .le. 0)
     1     call error('Warning In Nbcreate: Second Deriv may not work')
czmc  phonon allocation
      save_PHREAL = 0
      save_PHIMAG = 0
      if (phonon) then
         junk = memsys3('free phreal.array',0,0,0)
         if (memsys3('create double precision phreal.array',
     $        mxvibmod*mxvibmod,pphreal) .ne. 0) then
            write(*,*)
     $           'warning in nbcreate: phreal array is not created'
            phonon = .false.
         else
            save_phreal = pphreal
         end if
         junk = memsys3('free phimag.array',0,0,0)
         if (memsys3('create double precision phimag.array',
     $        mxvibmod*mxvibmod,pphimag) .ne. 0) then
            write(*,*)
     $           'warning in nbcreate: phimag array is not created'
            phonon = .false.
         else
            save_phimag = pphimag
         end if
      end if
      return
      end
 
 
 
 
 
 
 
 
