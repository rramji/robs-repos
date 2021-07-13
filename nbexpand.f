      SUBROUTINE nbexpand(ntem)
czmc   this routine allocate dynamic allocation for hugh arrays
c     in COMMON blocks, NOTICE that arrays are NOT LOCAL!
c     this should be called whenever energy expression is setup
      IMPLICIT NONE
      include '../includes/maxparm.dc'
      include '../includes/gparm.dc'
      include '../includes/nbond.dc'
      include '../includes/pertth.dc'
      include '../includes/vibrat.dc'
 
c
      integer memsys
      integer loc
      integer locra
      integer junk
      integer locia
      integer ntem
c      call gsinf1('Allocating extra run-time memory for NB listing')
       write(ibout,*)'Allocating extra run-time memory for NB listing'
czmc
c     ntem = MXJNPA/10 in nbcreate
cmem   kludge by returning to skip nbond array expansion -- no longer dynamically allocated
c      if( ntem .tr. 0 ) then 
       
       if (ntem .lt. 10000000) then
          return
       endif
cmem end of kludge
       
czmc  energy/force never works for period
         if ( (.not. period) .or. fixmolec) then
c     1)
            locra = loc(vdwpr)
            if (memsys('reallocate real vdwpr.nbond',ntem,vdwadd,locra)
     $           .ne. 0) then
               call error(
     $          'Warning In Nbexpand: VDWPR array is not reallocated')
               vdwadd = 0
            end if
c     2)
            locra = loc(elcpr)
            if (memsys('reallocate real elcpr.nbond',ntem,elcadd,locra)
     $           .ne. 0) then
               call error(
     $        'Warning In Nbexpand: ELCPR array is not reallocated')
               elcadd = 0
            end if
         end if
c     3)
         locia = loc(jnb)
         if (memsys('reallocate integer jnb.nbond',ntem,jnbadd,locia)
     1        .ne. 0) then
            call error(
     $        'Warning In Nbexpand: ELCPR array is not reallocated')
            jnbadd = 0
         end if
c     5)
         if (period) then
            locia = loc(jnbcel)
            if (memsys('reallocate integer jnbcel.nbond',
     $           ntem,jnbceladd,locia).ne. 0) then
               call error(
     $        'Warning In Nbexpand: JNBCEL array is not reallocated')
               jnbceladd = 0
            end if
         end if
      end if
      if (jnbadd .le.0 .or. jnbceladd .le.0) then
czmc  can not do nbond listing, fatal for dynamics
         call error('Warning In Nbexpand: non-bond listing is not
     $        reallocated in energy expression')
         leeset = .false.
      end if
      RETURN
      END
 
 
 
 
 
 
 
 
