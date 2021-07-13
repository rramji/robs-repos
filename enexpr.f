C***********************************************************************
C***********************************************************************
      SUBROUTINE ENEXPR
C***********************************************************************
*     CONTROL FOR CALCULATING THE VALENCE ENERGY EXPRESSION
*
*@v250 NEED TO CHECK WHETHER LP DONE CORRECTLY IF LPHI ETC IS OFF
*MUST add subr for doing prt on constraints
      IMPLICIT NONE
      include '../includes/maxparm.dc'
      include '../includes/ener.dc'
      include '../includes/contrl.dc'
      include '../includes/param.dc'  !ffindx
      include '../includes/eeint.dc'
      include '../includes/atodat.dc'
      include '../includes/eparm.dc'
      include '../includes/number.dc'
      include '../includes/fildat.dc'
      include '../includes/gparm.dc'
      include '../includes/misc.dc'
      include '../includes/eopt.dc'
      include '../includes/resdat.dc'
      include '../includes/timer.dc'
      include '../includes/constraint.dc'
      include '../includes/nbond.dc'
      include '../includes/hbond.dc'
      include '../includes/element.dc'
      include '../includes/periodc.xtl'
      include '../includes/nbewald.xtl' !lautoew
      include '../includes/coord.dc'    !x,y,z needed for setewald call
      include '../includes/shell.dc'
      include '../includes/pertth.dc'
c
c wag 24/nov/91
c CM option
c accumulate nb exclusions from calls to ex123cm etc
c need SCRAT2 for time between calls to eexang and exclnb
c
      LOGICAL LSAVE
*
      INTEGER II,I,IH,ICARB,IUNTRA,maxknb,knb,j,i1
      INTEGER limatom,limmovfix
      LOGICAL NEEDLP
c
      integer memsys 
      integer loc
      integer locra
      integer junk
      integer ia(1)
      integer locia
      integer molpt
      integer molnatom
      integer molatoms
      integer centmass
      real ra(1)
c nk 19 aug 91
      integer npole
      double precision scale
*
      locia=loc(ia)
      locra=loc(ra)
*
*     RESET COUNTERS FOR MINIMIZATION (CNTE) AND DYNAMICS (ITSTEP,CURTIM
*
      CNTE=0
      ITSTEP=0
      CURTIM=0.0
      goodcm=.false.
*
      NRES   = 1
      NMOVATM = 0
      NBOND  = 0
      NTHETA = 0
      NIMPHI = 0
      NANGANG = 0
      NBXA    = 0
      NPHI   = 0
      NPITW = 0
      NNB    = 0
      NDON   = 0
      NACC   = 0
C
      CHKMISS=.FALSE.
      LPATOM=.FALSE.
C
C     FRSTENR IS SET TO .F. AFTER COMPLETION OF THE FIRST ENERGY
C       CALCULATION (WHETHER CALLED BY MIN OR DYN)
      FRSTENR=.TRUE.
      FRSTDYN=.TRUE.
      FRSTMIN=.TRUE.
C
C     NBFIN=F INDICATES THAT THE NB LIST MUST BE CONSTRUCTED
C       NBFIN IS RESET TO .F. AFTER COMPLETION OF NBLIST
      NBFIN=.FALSE.
C
C     HBFIN=F INDICATES THAT THE HB LIST MUST BE CONSTRUCTED
C       HBFIN IS RESET TO .F. AFTER COMPLETION OF HBLIST
      HBFIN=.FALSE.
C
c     LSETCHRG IS RESET TO .T. AFTER CHARGES ARE CALCULATED
      LSETCHRG=.FALSE.
C
C     LNBLIM=T INDICATE that the nonbond lists must be redone
      if (period) then
        if (NBMETHX.EQ.'NOLIST') go to 2532
      else
        if (NBMETH.EQ.'NOLIST') go to 2532
      endif
      lnblim=.true.
cmem -- skip for G90      call nbfree
 2532 continue
*
      IF (TIMER.GT.0) THEN
C       THE INTEGER VARIABLE TIMER IS SET TO CALL THE TIMING ROUTINES.
        WRITE (IBOUT,2000) 'TIME GOING INTO ENEXPR:'
 2000   FORMAT(1X,A)
        CALL TIMERE
        CALL TIMERB
      ENDIF
c
c     limatom is total of ALL atoms (mov + fix + notineex)
c     limmovfix is total of all mov+fix atoms
      if (period) then
        limatom=natmcel
        limmovfix=movfixcel
      else
        limatom=natom
        limmovfix=ntotatm
      endif
c
c     find all "molecules" for periodic cases -----
c
      if (period .or. FIXMOLEC) then
c90         junk=memsys('free molecules.pointer',0,0,0)
c90         junk=memsys('free molecules.natoms',0,0,0)
c90         junk=memsys('free molecules.atoms',0,0,0)
*
c90         if (memsys('create integer molecules.pointer',limatom,molpt,
c90     $        locia).ne.0.or.
c90     $        memsys('create integer molecules.natoms',limatom,
c90     $        molnatom,locia).ne.0.or.
c90     $        memsys('create integer molecules.atoms',limatom,
c90     $        molatoms,locia).ne.0) then
c90           call error2('Ran out of memory finding discrete '//
c90     $           'molecules setting up the energy expression')
c90           return
c90         end if
c
         call fmolecules(1,limatom,nmolecules,ia(molpt),
     $        ia(molnatom),ia(molatoms))
 
czmc  pointer is allocated with size of limatom, may lock up too much
c      memory. Since nmolecules is obtained after calling fmolecules,
c      redo again to minimize
c      memory allocation. fmolecules does not cost much CPU time.
c      I  will use reallocate later. zmc
c      nmolecules = mmovmol + n of fixed molecules
c90         junk=memsys('free molecules.atoms',0,0,0)
c90         junk=memsys('free molecules.natoms',0,0,0)
c90         junk=memsys('free molecules.pointer',0,0,0)
*
c90         if (memsys('create integer molecules.pointer',nmolecules,molpt,
c90     $        locia).ne.0.or.
c90     $        memsys('create integer molecules.natoms',nmolecules,
c90     $        molnatom,locia).ne.0.or.
c90     $        memsys('create integer molecules.atoms',limatom,
c90     $        molatoms,locia).ne.0) then
c90           call error2('Ran out of memory finding discrete '//
c90     $           'molecules setting up the energy expression')
c90           return
c90         end if
c90         call fmolecules(1,limatom,nmolecules,ia(molpt),
c90     $        ia(molnatom),ia(molatoms))
c90      endif
C
c**** FIND MOVABLE AND FIXED ATOMS
C
C     EEXATM OR EEXATMX sets NATOM=IBIG
C       and sets up the MUNAI AND MOVATM lists
C         (this contains a list of all atoms allowed to move in the
C         energy calculations)
*
      IF (PERIOD) THEN
        CALL  EEXATMX
      ELSE
        CALL EEXATM
      ENDIF
*
***** check for shell model, lone pairs etc
c
c     limmovfix is total of all mov+fix atoms
      if (period) then
        limmovfix=movfixcel
      else
        limmovfix=ntotatm
      endif
c
c     zero zcg and icore arrays for shell model
      do 2343 ii=1,limmovfix
        i=movatm(ii)
*       note that zcg and icore for non mov+fix atoms are not redefined
        zcg(i)=cg(i)
        icore(i)=0
2343  continue
c
      lshell=.false.
      NEEDLP=.FALSE.
      do 2344 ii=1,limmovfix
        i=movatm(ii)
        if (fflabel(iac(i))(1:2).eq.'Lw') then
          if (period) then
            j=iconx(munai(i),1)
          else
            j=icons(munai(i),1)
          end if
          if (mtavom(j).eq.0) then
            write(ibout,2345)i,j
2345        format(2x,'error found in enexpr:shell atom ',i6,' is ',
     $        'bonded to atom ',i6,' that is not mov or fix')
            return
          endif
c
c         j is the core atom corresponding to shell i
c         zcg is the total charge for shell and core
c         icore(j) is the shell for this core (=0 if no shell)
          zcg(j)=cg(i)+cg(j)
          icore(j)=i
          lshell=.true.
        end if
c
        IF (ffneedlp(iac(i))) CALL CHKLP(I,NEEDLP)
!@v250  check whether works correctly for cases of LP atoms and for
!         ____L atoms
*@v250  add check for other types of special atoms
2344  continue
c
      if (lshell.and.prtthrm) then
       write(ibout,*)' shell model for pert therm not programmed yet.'
      end if
c
c     finished shell model setup
*     ALSO CHECKED TO SEE THAT LONE PAIRS ARE ALL PRESENT
C
***** CHECK FOR FAILURE ( MISSING PARAMETERS OR NMOVATM=0)
*
      IF (NEEDLP) THEN
        WRITE(IBOUT,33)
   33   FORMAT(/2X,'*** MUST ADD LP FOR AT LEAST ONE ATOM:')
        CHKMISS=.TRUE.
      ENDIF
      IF (CHKMISS)  THEN
        CALL ERROR2('**SOME ATOMS ARE NOT DEFINED,'
     1    //' DO SO NOW (USING MOD_PARMS)')
        CALL ERROR2( ' THEN TOGGLE UPDATE EE TO REGENERATE THE ENERGY'
     1    //' EXPRESSION')
        IF (BTCH) CALL DIE
      ELSE
        IF (NMOVATM.GT.0) THEN
          GO TO 30
        ELSE
           call error2('WARNING - No movable atoms, will skip '//
     $          'optimization')
        ENDIF
      ENDIF
      RETURN
*
   30 CONTINUE
!@v250 get rid over nmove,nmovat
*     NMOVE (NUMBER.DC) USED IN MINMIB,TRAJCT,BENER,BIOBATCH
      NMOVE=NMOVATM
*     NMOVAT (ATODAT.DC) USED IN MANY BIO1.2 ROUTINES
*     NMOVATM (EOPT.DC) USED IN XTL1.1 ROUTINES
      NMOVAT=NMOVATM
*
*     CHECK FOR SHRUNK CH BOND
*
      IF (SHRNKCH) THEN
        NUMBSCH=0
        DO 96 II=1,LIMATOM
          I=MOVATM(II)
          IF (fflabel(IAC(I)).EQ.ATCSCH) THEN
            IH=MUNAI(I)
            ICARB=ICONS(IH,1)
            IF (ffatno(IAC(ICARB)).EQ.CATNO) THEN
              NUMBSCH=NUMBSCH+1
              LSTSCH(NUMBSCH)=I
            ENDIF
          ENDIF
   96   CONTINUE
        IF (NUMBSCH.GT.MAXSCH) THEN
          WRITE(IBOUT,97)NUMBSCH,MAXSCH
   97     FORMAT(2X,'NUMBER OF SHRUNK CH BONDS (',I5,') EXCEEDS THE ',
     1      'MAXIMUM NUMBER (',I5,'), REDO CALCULATION')
          RETURN
        ENDIF
      ENDIF
c
c     initialize for CM option (eg for center of Cp ring) (tey option)
c
      if (ndepend.ne.0) then
!@v250  change initbzermas so that limtot input is deleted
!       also change so that only do for movable atoms
        call initbzermas(ndepend,lendepend,iscmatm,ldepend,
     .            maxdepend,mdepend,btrn,limmovfix)
*       quit if exceeded dimensions
        if (ndepend.eq.0) return
      end if
*
*  DO FIXED MOLECULE OPTION
*
      IF (FIXMOLEC) THEN
        CALL EEXFIXMOL(ia(molpt),ia(molnatom),ia(molatoms))
*       SKIP VALENCE FF STUFF
        GO TO 355
      ENDIF
*
C  Find all bonds
*
      IF (LBOND) CALL EEXBND
*
C  Find all angle interactions (two bonds sharing a common atom)
*
      IF (LTHETA) CALL EEXANG
*
C  Find all torsions
c
c     cases where two bonds are bonded to opposite ends of a third bond
      IF ((LPHI.OR.LBNDTOR.OR.LANGTOR) .OR. LPITWIST) THEN
        CALL EEXTOR
      ENDIF
*
C  Find all inversions
c
c     cases with three bonds to a common atom
      IF (LIMPHI) THEN
        CALL EEXINV
      ENDIF
*
C  Find all cases of 3 center angle-angle interactions
c
      IF (LANGANG) THEN
        CALL EEXAAI
      ENDIF
*
C  Find all bnd-x-angs
c
c     cases with three bonds to a common atom
      IF (LBNDXANG) THEN
        CALL EEXBXA
      ENDIF
*
C  Find all one center ang-x-angs
c
c     cases with no common bond
      IF (LANGXANG) THEN
        CALL EEXAXA
      ENDIF
*
C  Eliminate bond,angle,torsion,and inversion terms not allowed for
*    PERIODIC BOUNDARY CONDITIONS
*
      IF (PERIOD .AND. (.NOT.MOLXTL) ) CALL ENEXPRX
*
*  ADD PARAMETERS FOR GENFF (DREIDING B) FF
*
!@v250 replace with logical to use genff to unknown variables
      IF (FFINDX.EQ.5) THEN
        IF (.NOT.LALLTOR) THEN
*         GENFF (DREIDING B) REQUIRES EEFULL=.T.
          WRITE(IBOUT,2354)
 2354     format(2x,'warning:genff generally expects to include all',
     $      ' POSSIBLE TORSION '/
     $      4X,'PARAMETERS, HOWEVER LALLTOR=.FALSE., ',
     $      'PLEASE CHECK WHETHER'/
     $      2X'THIS IS INTENDED')
        ENDIF
        IF (.NOT.LALLINV) THEN
*         GENFF (DREIDING B) REQUIRES EEFULL=.T.
          WRITE(IBOUT,2355)
 2355     format(2x,'warning:genff generally expects to include all',
     $      ' POSSIBLE INVERSION'/
     $      4X,'PARAMETERS, HOWEVER LALLINV=.FALSE., ',
     $      'PLEASE CHECK WHETHER'/
     $      2X'THIS IS INTENDED')
        ENDIF
        CALL GENFF
      ENDIF
*
C  Construct list of interactions (bonds, angles) to be excluded from
C    nonbond calculations
*
      IF (LNBOND) THEN
        CALL CLEXCLNB
      ENDIF
*
C  Find h-bond donor cases
 355  IF (LHBOND) THEN
        CALL EEXHBD
      ENDIF
*
C  Find h-bond acceptor cases
      IF (LHBOND) THEN
        CALL EEXHBA
      ENDIF
*
C  Find energy expression for user type force field
*
      IF (.NOT.FIXMOLEC) CALL UGETEE
*
C  Do SHAKE EEX
*
      IF (SHAKE) THEN
        CALL EEXSHK
        IF (NBNDSHK.LE.0 .AND.  NANGSHK.LE.0) THEN
          CALL ERROR2(' *** NO BONDS OR ANGLES FOUND FOR RIGID '//
     1      'CONSTRAINT SET, WILL IGNORE RIGID')
          SHAKE= .FALSE.
        ENDIF
      ELSE
        NBNDSHK=0
        NANGSHK=0
      ENDIF
*
*     SETUP NONBOND QUANTITIES
*
      CALL NBSPLSET
      IF (PERIOD) THEN
        IF (NBMETHX.EQ.'EWALD ') THEN
*         THIS DOES EXTRA WORK, ONLY NEED ETA,RCUT,HCUT
          LSAVE=LAUTOEW
          LAUTOEW=.TRUE.
          CALL SETEWALD
          LAUTOEW=LSAVE
        elseif ( nbmethx.eq.'CELLMM') then
c         nk 1 july 91
          call cmmsetup
          scale = 0.8d0
          npole = 35
          call redsetup(scale,npole)
          call mnmaglst
        endif
      ELSE
        IF (NBMETH.EQ.'CELLMM') THEN
          CALL CMMSETUP
        ENDIF
      ENDIF
*
!@v250 reintroduce the res,resid,ibase options to automatically group
!      into sets with zero charge
C  Create the SEGID AND IBASE parameters
*     SEGID AND IBASE ARE PASSED IN .EEX FILE, RES,RESID ARE NOT
      RES(1)='GRD '
      RESID(1)='GRD '
!2/90      SEGID(1)='GRD '
      IBASE(1)=0
      IBASE(2)=IBIG
* (WAG) IBASE CAN BE USED IN NBLIST TO DO A PRESORT OVER RESIDUES BEFORE
*    CHECKING PAIRWISE DISTANCES OF EACH ATOM PAIR
*  THERE SHOULD BE NRES + 1 VALUES OF IBASE
*
      IF (.NOT.FIXMOLEC) THEN
*
*       CHECK THAT IORDER ARE CONSISTENT FOR LPOLYENE CALC
*
        IF (LPOLYENE) CALL CHKIORDER
*
C       Write all of the above to the output file
*
        CALL EEXOUT
      ENDIF
*
*  SET PROGRAM SWITCHES
*  ____________________
*
      leeset=.true.
czmc  create hugh array in memsys
      call nbcreate
c
      IF( LWRTTJ )  THEN
         iuntra=29
         if (.not.lftrj) then
            CLOSE(IUNTRA)
            WRITE(IGOUT,9854)
 9854       FORMAT(/2X,'The current trajectory file has been closed.')
         endif
      ENDIF
*
      IF (TIMER.GT.0) THEN
C       THE INTEGER VARIABLE TIMER IS SET TO CALL THE TIMING ROUTINES.
        WRITE (IBOUT,2000) 'TIMES IN ENEXPR:'
        CALL TIMERE
        CALL TIMERB
      ENDIF
*
      RETURN
      END
