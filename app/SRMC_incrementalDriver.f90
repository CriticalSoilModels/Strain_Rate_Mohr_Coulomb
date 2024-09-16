! Copyright (C)  2007 ... 2023  Andrzej Niemunis
!
! incrementalDriver is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! incrementalDriver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.


! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,USA.

! August 2007: utility routines added
! Februar 2008 :   *Repetition repaired again
! March 2008 :    *ObeyRestrictions added
! April 2008 : command line options: test= , param= , ini= , out= , verbose=  added
! November 2008: warn if  divergent equil. iter.  $\| u\_dstress \|$ too large
! February 2009: *ObeyRestrictions works with *Repetitions
! August 2010 increases output precision   +  keyword(3) error interception
! December 2015: fragments of niemunis\_tools\_lt  unsymmetric\_module incorporated into a single file incrementalDriver.f
! January 2016 if nstatev = 0 no statev( ) values will be read in but we set nstatv= 1 and statev(1)= 0.0d0
! October 2016 exit from step on inequality condition, import step loading data from a file, write every n-th state only
! November 2016 alignment of stress
! Jan 2017 exit on inequality  corrected  twice
! June 2017 undo changes in stress and state from ZERO call of  umat  (just for jacobian, with zero dstran and zero dtime )
! Dec  2017 parser disregards comments beyond \#
! Sept 2019   c\_dstran(:) = 0   in line  590  otherwise c\_dstran may be used  before being initialized.
! Sept 2019   random walk
! 2020    exit on inequality condition extended to mixed components
! Jan 2023  the initial temperature and the increase of temperature per step can be prescribed
!           (tested with thermal expansion or thermal stress in elasticity  under *ObeyRestrictions only )

!  Main program that\_calls\_umat ( performs calculation writing  to output.txt).
PROGRAM that_calls_umat   ! written by  A.Niemunis  2007 - 2023
   use incrementalDriver_funcs, only: splitaLine, ReadStepCommons, PARSER, get_increment,&
      USOLVER, EXITNOW
   use MOD_NAMC_ESM           , only: UMAT => UMAT_NAMC

   implicit none

   character*80  cmname,rebarn
   integer ndi,nshr,ntens,nstatv,nprops,ncrds
   integer noel,npt,layer,kspt,lrebar,kinc,i
   real(8),parameter,dimension(3,3):: delta = reshape([1,0,0,0,1,0,0,0,1],[3,3])

   parameter(ntens=6,ndi=3,nshr=3,ncrds=3) ! same ntens as in SOLVER
   parameter( noel=1 ,npt=1,layer=1,kspt=1,lrebar=1)
   parameter( rebarn ='xxx')
   real*8 dtime,temp,dtemp,sse,spd,scd,rpl,drpldt,pnewdt,celent
   real*8 stress(ntens),&
      ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),&
      stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
      coords(ncrds),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
   character(len=1) :: aChar                                       ! AN 2016
   character(len=40):: keywords(10), outputfilename,&
      parametersfilename,&
      initialconditionsfilename, testfilename, outputfilename1,&
      exitCond, ImportFileName,mString, keyword2, &                 ! AN 2016
      aShortLine, leftLine, rightLine
   character(len=260) ::  inputline(6), aLine, heading
   character(len=520) :: hugeLine

   character(len=10):: timeHead(2), stranHead(6), stressHead(6)      ! AN 2016
   character(len=15), allocatable :: statevHead(:)                   ! AN 2016

   logical :: verbose
   ! logical :: EXITNOW, existCond,okSplit                             ! AN 2016 ! WaveHello: ExitNow (bool) conflicts with the function
   logical :: existCond,okSplit
   real(8), dimension(6,6)  :: cMt , cMe
   real(8), dimension(6)  :: mb, mbinc


   integer :: mImport, columnsInFile(7),every,ievery                 ! AN 2016
   real(8) ::  importFactor(7)
   real(8),dimension(20) :: oldState, newState,dState
   real(8), allocatable :: props(:), statev(:), r_statev(:)

   real(8),dimension(3,3):: Qb33,eps33,T33

   integer:: ifstress(ntens), maxiter, ninc,kiter, ikeyword, &
      iRepetition, nRepetitions, kStep,iStep,nSteps,ntens_in

   real(8):: r_stress(ntens),a_dstress(ntens),u_dstress(ntens),&
      stress_Rosc(ntens),r_stress_Rosc(ntens),           &
      ddstress(ntens), c_dstran(ntens) ,                 &
      deltaLoadCirc(6),phase0(6),deltaLoad(9),           &
      dstran_Cart(6), ddsdde_bar(6,6), deltaTime,        &
      deltaTemp
   real(8),parameter :: sq3=1.7320508075688772935d0,&
      sq6=2.4494897427831780982d0,&
      sq2=1.4142135623730950488d0,&
      Pi =3.1415926535897932385d0
   real(8),parameter :: i3=0.3333333333333333333d0,&
      i2=0.5d0,                  &
      isq2=1/sq2,                &
      isq3=1.0d0/sq3,            &
      isq6=1.0d0/sq6

   !  M for isomorphic Roscoe variables P,Q,Z,....
   real(8), parameter,dimension(1:6,1:6)::MRoscI=reshape( &
      (/-isq3,-2.0d0*isq6,0.0d0,  0.0d0, 0.0d0, 0.0d0,     &
      -isq3, isq6,      -isq2,  0.0d0, 0.0d0, 0.0d0,  &
      -isq3, isq6,       isq2,  0.0d0, 0.0d0, 0.0d0,  &
      0.0d0, 0.0d0,     0.0d0, 1.0d0, 0.0d0, 0.0d0,  &
      0.0d0, 0.0d0,     0.0d0, 0.0d0, 1.0d0, 0.0d0,  &
      0.0d0, 0.0d0,     0.0d0, 0.0d0, 0.0d0, 1.0d0/),&
      (/6,6/))

   real(8), parameter,dimension(1:6,1:6)::MRoscImT=MRoscI            !  latest $\cM^{-T}$ (is orthogonal)

   !  M for isomorphic Rendulic $ sigma_{11}= -T_{11}$,  $sigma_{22}  = -(T_{22} + T_{33}) / \sqrt(2) $,  $ Z= \dots$
   real(8), parameter,dimension(1:6,1:6)::MRendul=reshape( &
      (/ -1.0d0, 0.0d0,  0.0d0,  0.0d0,0.0d0,0.0d0,         &
      0.0d0, -isq2,  -isq2,  0.0d0,0.0d0,0.0d0,       &
      0.0d0, -isq2,   isq2,  0.0d0,0.0d0,0.0d0,       &
      0.0d0,  0.0d0, 0.0d0,  1.0d0,0.0d0,0.0d0,       &
      0.0d0,  0.0d0, 0.0d0,  0.0d0,1.0d0,0.0d0,       &
      0.0d0,  0.0d0, 0.0d0,  0.0d0,0.0d0,1.0d0/),     &
      (/6,6/))

   real(8), parameter,dimension(1:6,1:6)::MRendulmT=MRendul          !  latest  $\cM^{-T}$   (is orthogonal)

   !  M for Roscoe variables $p,q,z,....$
   real(8), parameter,dimension(1:6,1:6)::MRosc=reshape( &
      (/-i3,-1.0d0, 0.0d0,    0.0d0,0.0d0,0.0d0,          &
      -i3, i2, -1.0d0,      0.0d0,0.0d0,0.0d0,       &
      -i3,i2, 1.0d0,       0.0d0,0.0d0,0.0d0,       &
      0.0d0, 0.0d0,0.0d0,  1.0d0,0.0d0,0.0d0,       &
      0.0d0, 0.0d0, 0.0d0, 0.0d0,1.0d0,0.0d0,       &
      0.0d0, 0.0d0, 0.0d0, 0.0d0,0.0d0,1.0d0/),     &
      (/6,6/))

   !  latest  $\cM^{-T}$   (is not orthogonal)
   real(8), parameter,dimension(1:6,1:6)::MRoscmT=reshape( &
      (/-1.0d0, -2.0d0*i3, 0.0d0,  0.0d0, 0.0d0,0.0d0,     &
      -1.0d0,   i3,      -i2,    0.0d0,0.0d0,0.0d0,       &
      -1.0d0,   i3,       i2,    0.0d0,0.0d0,0.0d0,       &
      0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0,           &
      0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,           &
      0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/),         &
      (/6,6/))

   !  M for Cartesian coords $T_{11}, T_{22}, T_{33}, T_{12},.....$
   real(8), parameter,dimension(1:6,1:6)::MCart=reshape( &
      (/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,    &
      0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,   &
      0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,   &
      0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0,   &
      0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,   &
      0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/), &
      (/6,6/))

   real(8), parameter,dimension(1:6,1:6)::MCartmT=MCart              !  latest  $\cM^{-T}$  (is orthogonal)


   real(8),dimension(1:6,1:6)::M,MmT                                 !  currrent $\cM$ and $\cM^{-T}$  for a given iStep
   real(8) :: aux1,aux2

   type descriptionOfStep
      integer:: ninc,maxiter, ifstress(ntens),columnsInFile(7),mImport ! AN 2016
      real(8) :: deltaLoadCirc(ntens),phase0(ntens),deltaLoad(9),    &
         dfgrd0(3,3), dfgrd1(3,3),deltaTime, importFactor(7),&
         deltaTemp                                            ! AN 2023 temperat
      character(40) :: keyword2, keyword3, exitCond,ImportFileName    ! AN 2016
      real(8),dimension(1:6,1:6) :: cMt, cMe
      real(8),dimension(1:6) :: mbinc
      logical::existCond                                              ! AN 2016
   end type  descriptionOfStep

   type StressAlignment
      logical:: active
      character(len=40) :: ImportFileName
      integer:: kblank,nrec,kReversal, ncol
      integer,dimension(100) :: Reversal
      integer,dimension(6):: isig
      real(8),dimension(6) :: sigFac
   end type StressAlignment

   type(StressAlignment) :: align

   type(descriptionOfStep) :: ofStep(30)            !  stores descriptions of up to 30 steps which are repeated


   ! [1]  read the command-line parameters to set the file names  *************
   continue

   parametersfilename = 'parameters.inp'
   initialconditionsfilename  = 'initialconditions.inp'
   testfilename = 'test.inp'
   outputfilename = '--'
   verbose = .true.
   call  get_command_line_arguments() ! command line  can override the above file names


!  [2] read the material parameters ************************************************
   open(1,err=901,file=parametersfilename,status='old')
   read(1,'(a)') cmname
   i = index(cmname, '#')
   if(i == 0) then
      cmname = trim(cmname)
   else
      cmname = cmname(:i-1)
      cmname = trim(cmname)
   endif
   read(1,*) nprops
   allocate( props(nprops) )
   do i=1,nprops
      read(1,*) props(i)
   enddo
   close(1)

![3]  read initial conditions and  initialize everything **************************************
   open(1,err=902,file=initialconditionsfilename,status='old')
   read(1,*) ntens_in
   stress(:) = 0.0d0
   time(:) = 0.0d0
   stran(:)=0.0d0
   dtime = 0.0d0
   do i=1,ntens_in
      read(1,*) stress(i)
   enddo

   temp  = 0.0d0                    !  default initial temperature               AN 2023
   read(1,'(a)') aLine              ! AN 2023  aLine may be *temperature= 30.0 or   nstatv
   i = index(aLine,'=')             ! AN 2023
   if(i==0) then                   ! AN 2023 '=' is absent so read nstatv
      read(aLine,*) nstatv          ! AN 2023
   else                              ! AN 2023
      aLine = trim(aLine(i+1:))
      read(aLine,*) temp             ! AN 2023
      read(1,*) nstatv                ! AN 2023
   endif                             ! AN 2023

   if(nstatv >= 1) then
      allocate(statev(nstatv), r_statev(nstatv), statevHead(nstatv)) !  AN 2016
      statev(:) = 0.0d0
      do i=1,nstatv
         read(1,*,end=500) statev(i)   !
      enddo
   else
      allocate( statev(1) , r_statev(1), statevHead(1)  )             !  AN 2016 formal placeholder not really used
      statev(:) = 0.0d0
      nstatv = 1
   endif
500 continue
   close(1)

!      ntens\_in = 6
!      nstatv = 300
!      allocate( statev(nstatv) , r\_statev(nstatv) )
!      statev = 0.0d0
!      call ParaelasticInitialCondition(statev)        ! Loads two states into the stack

![4] read a piece from the loading path ***********************************************************





   open(1,err=903,file=testfilename,status='old')
![4.1] read the outputfilename from test.inp, create/open this file and write the tablehead, heading(if any)  and the first line = initial conditions
   read(1,'(a)') aLine
   i = index(aLine,'#')
   if(i==0) then
      outputfilename1=trim(aLine)
      heading = '#'
   else
      outputfilename1=trim(aLine(:i-1))
      heading = trim( aLine(i+1:) )
   endif
   if(outputfilename == '--') outputfilename = outputfilename1
   open(2,err=904,file=outputfilename)

   do i=1,2
      write(timeHead(i),'(a,i1,a)')  'time(',i, ')'
   enddo
   do i=1,6
      write( stranHead(i), '(a,i1,a)' )   'stran(',i, ')'
      write(stressHead(i), '(a,i1,a)' )  'stress(',i, ')'
   enddo
   do i=1,nstatv
      write(statevHead(i), '(a,i0,a)' )  '  statev(',i, ')'
   enddo
   write(2,'(a14,500a20)') timeHead,stranHead,stressHead,statevHead


   if(heading(1:1) /= '#') write(2,*) trim(heading)
   write(2,'(500(g17.10,3h    ))') time+(/dtime,dtime/), stran, stress, statev

![4.2]  loop over keywords(1) unless keyword(1) = *Repetition  it is copied to keyword(2) which is the true type of loading
   kStep = 0  ! kStep = counter over all steps whereas  iStep = counter over steps within a *Repetition
   do 200 ikeyword=1,10000
      read(1,'(a)',end=999) keywords(1)
      keywords(1) = trim( keywords(1) )
      if(keywords(1) == '*Repetition') then
         read(1,*) nSteps, nRepetitions
      else
         nRepetitions=1
         nSteps=1
         keywords(2) = keywords(1)
      endif

      do 130  iRepetition  = 1,nRepetitions
         do 120 iStep = 1,nSteps
            kStep = kStep + 1

            ! ! write to the screen before the first increment
            ! write(*,'(12H ikeyword = ,i3, 8H kstep = ,i3,7H kinc = ,i5, 9H kiter = ,i2, 8H TEMP = ,f9.4, 8H TIME = ,f9.4)'), ikeyword , kstep, kinc, kiter ,TEMP,TIME(1)                    ! AN 2023 temperat
            ! write to the screen before the first increment
            ! WaveHello: Replacing this format
            write(*,'(A,I3,A,I3,A,I5,A,I2,A,F9.4,A,F9.4)') &
               ' ikeyword = ', ikeyword, &
               ' kstep = ', kstep, &
               ' kinc = ', kinc, &
               ' kiter = ', kiter, &
               ' TEMP = ', TEMP, &
               ' TIME = ', TIME(1)
            if(iRepetition > 1) then  ! while repeating  recall the loading parameters of the repeated step read in during the first iRepetition
               ninc             = ofStep(istep)%ninc
               maxiter          = ofStep(istep)%maxiter
               ifstress         = ofStep(istep)%ifstress
               deltaLoadCirc    = ofStep(istep)%deltaLoadCirc
               phase0           = ofStep(istep)%phase0
               deltaLoad        = ofStep(istep)%deltaLoad
               dfgrd0           = ofStep(istep)%dfgrd0
               dfgrd1           = ofStep(istep)%dfgrd1
               deltaTime        = ofStep(istep)%deltaTime
               keywords(2)      = ofStep(istep)%keyword2
               keywords(3)      = ofStep(istep)%keyword3
               cMe              = ofStep(istep)%cMe
               cMt              = ofStep(istep)%cMt
               mbinc            = ofStep(istep)%mbinc
               exitCond         = ofStep(istep)%exitCond                       ! AN 2016
               existCond        = ofStep(istep)%existCond                      ! AN 2016
               ImportFileName   = ofStep(istep)%ImportFileName                 ! AN 2016
               mImport          = ofStep(istep)%mImport                        ! AN 2016
               columnsInFile    = ofStep(istep)%columnsInFile                  ! AN 2016   7 integers with numbers of columns  (or value = 0)
               importFactor     = ofStep(istep)%importFactor                   ! AN 2016   7 real factors to be multiplied with columns
               goto 10  ! jump over reading, because reading of steps is performed only on the first loop, when iRepetition==1
            endif

            if(keywords(1) == '*Repetition') read(1,'(a)') keywords(2)          ! = *LinearLoad  or *CirculatingLoad or *ObeyRestrictions...
            ! otherwise keywords(2) = keywords(1)

            call splitaLine(keywords(2),'?', keywords(2), exitCond, existCond)  ! AN 2016 look for exit condition in keywords(2)

            keywords(2)  = trim(keywords(2))

            ifstress(:)=0                                                      ! default strain control
            deltaLoadCirc(:)=0.0d0                                             ! default zero step increment
            phase0(:)=0.0d0                                                    ! default no phase shift
            deltaLoad(:) = 0.0d0
            deltaTemp = 0.0d0    !  AN 2023  temperature increase per step
            dfgrd0 = delta
            dfgrd1 = delta


            if(keywords(2) == '*DeformationGradient') then
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016 ! AN 2023 temperat
               keywords(3) = '*Cartesian'
               do i=1,9
                  read(1,*)  deltaLoad(i)                                      !  dload means total change in the whole step here
               enddo
               goto 10
            endif
            if (keywords(2) == '*CirculatingLoad') then
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)        ! AN 2016  read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

               read(1,*) keywords(3)                                           !  = Cartesian or Roscoe  or RoscoeIsomorph or Rendulic
               keywords(3)  = trim(keywords(3))
               do i=1,6
                  read(1,*) ifstress(i),deltaLoadCirc(i),phase0(i),deltaLoad(i)   !  dload means amplitude here
               enddo
               goto 10
            endif
            if(keywords(2) == '*LinearLoad') then
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               read(1,'(a)') keywords(3)
               do i=1,6
                  read(1,*) ifstress(i), deltaLoad(i)                           !  dload means total change in the whole step here
               enddo
               goto 10
            endif

            keyword2 = keywords(2)
            if(keyword2(1:11) == '*ImportFile') then
               keywords(2) = '*ImportFile'; keyword2 = keyword2(12:)
               call splitaLine( keyword2,'|',ImportFileName, mString,okSplit)
               if(.not.okSplit)    stop 'missing | in line *ImportFile'
               read(mString,*) mImport
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016  read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               read(1,'(a)') keywords(3)

               columnsInFile(:) = 0; importFactor(:) = 1

               do i=1,6
                  read(1,'(a)') aShortLine
                  call splitaLine(aShortLine,'*',leftLine,rightLine,okSplit )
                  read(leftLine,*)  ifstress(i), columnsInFile(i)
                  if(okSplit)  read(rightLine,*)  ImportFactor(i)
                  !         read(1,*) ifstress(i), columnsInFile(i) , ImportFactor(i)                           !  dload means total change in the whole step here
               enddo
               if(deltaTime <= 0) then
                  read(1,'(a)') aShortLine
                  call splitaLine(aShortLine,'*',leftLine,rightLine,okSplit)
                  read(leftLine,*) columnsInFile(7)
                  if(okSplit)  read(rightLine,*)  ImportFactor(7)
                  !         read(1,*) columnsInFile(7),  ImportFac(7)
               endif !deltaTime

               !**********************************************************
               call readAlignment(align, ImportFileName )
               !**********************************************************

               goto 10
            endif

            if(keywords(2) == '*OedometricE1') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Cartesian'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               read(1,*)   deltaLoad(1)
               goto 10
            endif
            if(keywords(2) == '*OedometricS1') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Cartesian'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               ifstress(1) = 1
               read(1,*)   deltaLoad(1)
               goto 10
            endif
            if(keywords(2) == '*TriaxialE1') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Cartesian'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               read(1,*) deltaLoad(1)
               ifstress(2:3) = 1
               goto 10
            endif
            if(keywords(2) == '*TriaxialS1') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Cartesian'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               read(1,*)   deltaLoad(1)
               ifstress(1:3) = 1
               goto 10
            endif
            if(keywords(2) == '*TriaxialUEq') then
               keywords(2) = '*LinearLoad'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               keywords(3) ='*Roscoe'
               read(1,*)   deltaLoad(2)                                      ! = deviatoric strain
               goto 10
            endif
            if(keywords(2) == '*TriaxialUq') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Roscoe'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016    read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               read(1,*)  deltaLoad(2)                                       ! = deviatoric stress
               ifstress(2) = 1
               goto 10
            endif
            if(keywords(2) == '*PureRelaxation') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Cartesian'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016  read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               goto 10
            endif
            if(keywords(2) == '*PureCreep') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Cartesian'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)     ! AN 2016 read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               ifstress(:) =  1
               goto 10
            endif
            if(keywords(2) == '*UndrainedCreep') then
               keywords(2) = '*LinearLoad'
               keywords(3) ='*Roscoe'
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)      ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               ifstress(2:6) =  1
               goto 10
            endif
            if(keywords(2) == '*ObeyRestrictions') then  ! ======================= *ObeyRestrictions ==================================
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)      ! AN 2016    read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               do i=1,6
                  read(1,'(a)')  inputline(i)  !=== a line of form  ''-sd1 + sd2 + 3.0*sd3 = -10  ! a comment '' is expected
                  if(index(inputline(i),'=')== 0) stop 'restr without "=" '
               enddo
               call parser(inputline, cMt,cMe,mb )
               mbinc = mb/ninc
               keywords(3) ='*Cartesian'
               ifstress(1:6) = 1            !=== because  we solve (cMt.ddsdde + cMe).dstran = mbinc for dstran
               goto 10
            endif
            if(keywords(2) == '*PerturbationsS') then
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)       ! AN 2016     read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               deltaTemp = 0   ! AN 2023 temperat
               read(1,*) keywords(3)  ! = *Rendulic  or *RoscoeIsomorph
               keywords(3) = trim( keywords(3) )
               if(keywords(3) .ne. '*Rendulic' .and. &
                  keywords(3) .ne. '*RoscoeIsomorph') &
                  write(*,*) 'warning: non-Isomorphic perturburbation'
               read(1,*)  deltaLoad(1)
               ifstress(1:6) =  1
               goto 10
            endif
            if(keywords(2) == '*PerturbationsE') then
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)    ! AN 2016      read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
               deltaTemp = 0  ! AN 2023 temperat
               read(1,*) keywords(3)
               keywords(3) = trim( keywords(3) )
               if(keywords(3) .ne. '*Rendulic' .and.  &
                  keywords(3) .ne. '*RoscoeIsomorph') &
                  write(*,*) 'warning: Anisomorphic perturburbation'
               read(1,*) deltaLoad(1)
               goto 10
            endif

            if(keywords(2) == '*RandomWalk') then                         ! AN 2019
               call ReadStepCommons(1,ninc,maxiter,deltaTime,deltaTemp,every)   ! AN 2023 temperat
               deltaTemp = 0  ! AN 2023 temperat
               read(1,*) keywords(3)
               keywords(3) = trim( keywords(3) )
               do i=1,6
                  read(1,*) ifstress(i),deltaLoad(i)    !  dload means max abs value of to be multiplied by random in (-1,1)
               enddo
               goto 10
            endif


            if(keywords(2) == '*End') stop '*End encountered in test.inp'
            write(*,*) 'error: unknown keywords(2)=',keywords(2)
            stop 'stopped by unknown keyword(2) in test.inp'

10          keywords(3) = trim(keywords(3))

            if(keywords(1) == '*Repetition' .and. iRepetition == 1) then      !  remember the description of step for the next repetition
               ofStep(istep)%ninc          =    ninc
               ofStep(istep)%maxiter       =    maxiter
               ofStep(istep)%ifstress      =    ifstress
               ofStep(istep)%deltaLoadCirc =    deltaLoadCirc
               ofStep(istep)%phase0        =    phase0
               ofStep(istep)%deltaLoad     =    deltaLoad
               ofStep(istep)%dfgrd0        =    dfgrd0
               ofStep(istep)%dfgrd1        =    dfgrd1
               ofStep(istep)%deltaTime     =    deltaTime
               ofStep(istep)%keyword2      =    keywords(2)
               ofStep(istep)%keyword3      =    keywords(3)
               ofStep(istep)%cMe           =    cMe
               ofStep(istep)%cMt           =    cMt
               ofStep(istep)%mbinc         =    mbinc
               ofStep(istep)%deltaTemp     =    deltaTemp               ! AN 2023 temperat
               ofStep(istep)%exitCond      =    exitCond                ! AN 2016
               ofStep(istep)%existCond     =    existCond                ! AN 2016
               ofStep(istep)%ImportFileName =   ImportFileName            ! AN 2016
               ofStep(istep)%mImport        =   mImport                    ! AN 2016
               ofStep(istep)%columnsInFile  =   columnsInFile              ! AN 2016    7 integers with numbers of columns  (or value = 0)
               ofStep(istep)%importFactor   =   importFactor               ! AN 2016
            endif

            if(any(ifstress==1)) maxiter = max(maxiter,5)                    ! at least 5 iterations
            if(all(ifstress==0) .and. keywords(2) .ne. '*ObeyRestrictions') maxiter = 1  ! no iterations are necessary

!     start the current step with zero-load call of umat() just to get the stiffness
            dstran(:)=0
            dtime=0
            dtemp=0
            kinc=0
            r_statev(:)=statev(:);  r_stress(:)=stress(:)     ! AN 21.06.2017 remember the initial state and stress
            !=== first call umat with dstrain=0 dtime=0 just for stiffness (=jacobian ddsdde)
            call  UMAT(stress,statev,ddsdde,sse,spd,scd,                       &
               rpl,ddsddt,drplde,drpldt,                               &
               stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
               ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,  &
               celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,0,kinc)   !=== some constitutive models require kStep=0 other do not

            statev(:)=r_statev(:);  stress(:)=r_stress(:)   !  AN 21.06.2017 recover stress and state  although the ZERO call of umat should not modify them

            select case( keywords(3) )
             case('*Cartesian' ) ;      M =  MCart   ;    MmT = MCartmT
             case('*Roscoe')     ;      M = MRosc  ; MmT = MRoscmT
             case('*RoscoeIsomorph');   M = MRoscI ;  MmT = MRoscImT
             case('*Rendulic')      ;   M = MRendul;   MmT = MRendulmT
             case default ;   write(*,*) 'Unknown keyword = ', keywords(3)
               stop  ' stopped by unknown keywords(3) in test.inp'
            end select


            if(keywords(2) == '*ImportFile' ) then                           ! AN 2016
               open(3,file=ImportFileName, status ='old', err=905)            ! AN 2016
               do                                                             ! AN 2016
                  read(3,'(a)',err=906) hugeLine;                             ! AN 2016
                  hugeLine= adjustL(hugeLine) ; aChar = hugeLine(1:1)         ! AN 2016
                  if(index('1234567890+-.',aChar) > 0) exit                  ! preceding non-numeric lines in ImportFile will be ignored
               enddo
               read(hugeLine,*,err=907) oldState(1:mImport)                  ! AN 2016
            endif

            ievery=1
            do 100 kinc=1,ninc

               if(keywords(2) == '*ImportFile') then                             ! AN 2016
                  deltaTemp = 0                                                  ! AN 2023 temperat
                  dTemp = 0                                                      ! AN 2023 temperat
                  read(3,*,iostat=i)  newState(1:mImport)                         ! AN 2016
                  if(i > 0) then                                                  ! AN 2016
                     write(*,*) 'error Import file',ImportFileName, 'line=', kinc+1  ! AN 2016
                     stop                                                            ! AN 2016
                  endif                                                           ! AN 2016
                  if(i < 0) then                                                  ! AN 2016
                     close(3)                                                    ! AN 2016
                     write(*,*) 'finished reading file', ImportFileName          ! AN 2016
                     exit                                                        ! AN 2016
                  endif                                                           ! AN 2016
                  dState = newState(:) -  oldState(:)                             ! AN 2016
                  do i=1,6                                                        ! AN 2016
                     if  (columnsInFile(i) == 0) cycle                               ! AN 2016
                     if (ifstress(i)==1  )  ddstress(i)= dState(columnsInFile(i))* ImportFactor(i)      ! AN 2016
                     if (ifstress(i)==0)    dstran(i)  = dState(columnsInFile(i))* ImportFactor(i)      ! AN 2016
                  enddo                                                            ! AN 2016
                  if(columnsInFile(7)/= 0) deltaTime= dState(columnsInFile(7)) * ImportFactor(i)   ! AN 2016
                  dtime = deltaTime                                             ! AN 2016
                  oldState(:) = newState(:)                                    ! AN 2016
               endif                                                             ! AN 2016

               if(keywords(2) /= '*ImportFile') then                            ! AN 2016

                  call get_increment(keywords, time, deltaTime, ifstress, ninc,  &    ! get inc. in terms of Rosc. variables
                     deltaLoadCirc,phase0,deltaLoad,deltaTemp,&
                     dtime, ddstress,  dstran, dTemp,  Qb33,  &    ! AB 2023 deltaTemp and dTemp added
                     dfgrd0, dfgrd1,drot )   ! to be called in each increment
               endif


               a_dstress(:)= 0.0d0    ! approximated Roscoe's dstress
               r_statev(:)=statev(:)  ! remember the initial state and stress till the iteration is completed
               r_stress(:)= stress    ! remembered Cartesian stress

               do 95 kiter=1, maxiter  !--------Equilibrium Iteration--------
                  c_dstran(:) = 0
                  if(keywords(2)== '*ObeyRestrictions'  ) then  ! ======================= ObeyRestrictions =========
                     ddsdde_bar = matmul(cMt,ddsdde) + cMe
                     u_dstress = - matmul(cMt,a_dstress)-matmul(cMe,dstran)+ mbinc
                     call  USOLVER(ddsdde_bar,c_dstran,u_dstress,ifstress,ntens)
                     dstran = dstran + c_dstran
                     print *, "Inc driver dstran", dstran
                     call  UMAT(stress,statev,ddsdde,sse,spd,scd,                       &
                        rpl,ddsddt,drplde,drpldt,                               &
                        stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
                        ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,  &
                        celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kStep,kinc)

                     if(kiter.lt.maxiter) then                                      ! continue iteration
                        statev(:)=r_statev(:)                                       ! 1) undo the update of state (done by umat)
                        a_dstress  = stress  - r_stress                             ! 2) compute the new approximation of dstress
                        stress(:)=r_stress(:)                                       ! 3) undo the update of (stress done by umat)
                     else
                        stran(:)=stran(:)+dstran(:)                                 !  accept  the updated state and stress (Cartesian)
                     endif
                  endif  ! ==== obey-restrictions

                  if(keywords(2) /= '*ObeyRestrictions'  ) then   ! ======================= disObeyRestrictions ==========
                     u_dstress = 0.0d0
                     where (ifstress == 1)  u_dstress =ddstress -a_dstress           ! undesired Roscoe stress
                     ddsdde_bar = matmul(matmul(M,ddsdde),transpose(M))              ! Roscoe-Roscoe stiffness

                     call  USOLVER(ddsdde_bar,c_dstran,u_dstress,ifstress,ntens)     ! get Rosc. correction  c\_dstran() caused by undesired Rosc. dstress
                     where (ifstress == 1) dstran = dstran + c_dstran                ! corrected Rosc. dstran where stress-controlled
                     dstran_Cart = matmul( transpose(M),dstran )                     ! transsform Rosc. to Cartesian dstran
                     call  UMAT(stress,statev,ddsdde,sse,spd,scd,                       &
                        rpl,ddsddt,drplde,drpldt,                               &
                        stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
                        ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,  &
                        celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kStep,kinc)

                     if (kiter.lt.maxiter) then                                      ! continue iteration
                        statev(:)=r_statev(:)                                        ! 1) forget the changes of state done in umat
                        stress_Rosc = matmul(M,stress)                               !    output from umat transform to Roscoe ?
                        r_stress_Rosc = matmul(M,r_stress)
                        where (ifstress ==1) a_dstress = stress_Rosc - r_stress_Rosc ! 2) compute the new approximation of stress
                        stress(:)=r_stress(:)                                        ! 3) forget the changes of stress done in umat
                     else
                        stran(:)=stran(:)+dstran_Cart(:)                              !  accept  the updated state and stress (Cartesian)
                     endif
                  endif  ! ==== disObey-restrictions

! 94                continue ! WaveHello: Label not used

                  if((kiter==maxiter) .and. mod(kinc,10)==0 .and. verbose ) then    ! write to screen after each increment
                     write(*, '(A,I3,A,I3,A,I5,A,I2,A,F9.4,A,F9.4)') &
                        ' ikeyword = ', ikeyword, &
                        ' kstep = ', kStep, &
                        ' kinc = ', kinc, &
                        ' kiter = ', kiter, &
                        ' TEMP = ', TEMP, &
                        ' TIME = ', TIME(1)
                  endif

95             continue !--------------------end of Equilibrium Iteration

               aux1 =  dot_product(a_dstress,a_dstress)
               aux2 =  dot_product(u_dstress,u_dstress)
               if((aux1>1.d-10 .and. aux2/aux1 > 1.0d-2) .or. (aux1<1.d-10 .and. aux2 > 1.0d-12)  ) then
                  write(*,*) 'I cannot apply the prescribed stress components,'// &
                     '||u_dstress|| too large.'                 ! check  the Rosc.stress error < toler
               endif

               if(keywords(2) =='*DeformationGradient' ) then                    !  rigid rotation of stress
                  T33 = map2T(stress,6)
                  T33 = matmul( matmul(Qb33,T33),transpose(Qb33))
                  stress=map2stress(T33,6)                                          ! rigid rotation of strain
                  eps33 = map2D(stran,6)
                  eps33 = matmul( matmul(Qb33,eps33),transpose(Qb33))
                  stran=map2stran(eps33,6)
               endif

               where(abs(time) < 1.0d-99)   time = 0.0d0   ! prevents fortran error write 1.3E-391
               where(abs(stran) < 1.0d-99)  stran = 0.0d0
               where(abs(stress) < 1.0d-99) stress = 0.0d0
               where(abs(statev) < 1.0d-99) statev = 0.0d0
               if(ievery==1) then
                  write(2,'(500(g17.10,3h    ))') time+(/dtime,dtime/), stran, stress, statev
               endif
               if(keywords(2) =='*PerturbationsS' .or. keywords(2) =='*PerturbationsE' ) then ! having plotted everything undo the increment
                  stran(:)=stran(:) - dstran_Cart(:)
                  statev(:)=r_statev(:)
                  stress(:)=r_stress(:)
               endif

               time(1)=time(1)+dtime     !  step time at the beginning of the next increment
               time(2)=time(2)+dtime     !  total time at the beginning of the next increment
               Temp = Temp + dTemp       !  AN 2023 total temperature at the beginning of the increment

               if( existCond ) then                                             ! AN 2016 only if a condition exists
                  if(  EXITNOW(exitCond, stress,stran,statev,nstatv)  ) exit   ! AN 2016 depending on  exitCond go to next step
               endif                                                            ! AN 2016
               ievery = ievery+1; if(ievery > every) ievery = 1

!*****************************************************
               if(keywords(2) == '*ImportFile' ) then
                  call  tryAlignStress(align, kinc, newState, mImport,stress,ntens)
               endif
               !***********************************************

100         continue  ! next kinc
120      continue  ! next iStep
130   continue  ! next iRepetition
200 continue  ! next keyword




! 998 stop 'End or record encountered in test.inp' ! WaveHello: Label not used
999 close(1)
   close(2)
   stop 'I have reached end of the file test.inp'
901 stop 'I cannot open the file parameters.inp'
902 stop 'I cannot open the file initialconditions.inp'
903 stop 'I cannot open the file test.inp'
904 stop 'I cannot open the outputfile'
905 stop 'I cannot open the ImportFile'
906 stop 'Error reading ImportFile in the first non-numeric records '
907 stop 'Error reading ImportFile in the first numeric record '

contains !========================================================
   !  contained in program\_that\_calls\_umat that reads the command line
   subroutine get_command_line_arguments()
      implicit none            ! ===file names in the command line override defaults
      integer :: iarg,narg, is, iargc
      integer, parameter :: argLength=40
      character(argLength) :: anArgument, argType, argValue
      narg = iargc()
      do iarg = 1,narg
         call getarg(iarg,anArgument)
         is = index(anArgument,'=')
         if(is == 0) stop 'error: a command line argument without "=" '
         argType = anArgument(:is-1)
         argValue =  anArgument(is+1:)
         select case (argType)
          case ('param')
            parametersfilename = argValue
          case ('ini')
            initialconditionsfilename = argValue
          case ('test')
            testfilename = argValue
          case ('out')
            outputfilename = argValue
          case ('verbose')
            if (argValue == 'true') verbose = .true.
            if (argValue == 'false') verbose = .false.
         end select
      enddo
   end  subroutine get_command_line_arguments

   !   contained in program\_that\_calls\_umat writes a 6x6 matrix for debugging with Mma
   subroutine write66(a)
      implicit none
      real(8),dimension(6,6) :: a,aT
      aT = Transpose(a)
      open(12,file='nic.m',access='append')
      write(12,'(6ha66={ ,( 2h{  ,5(f15.4,2h,  ),f15.4, 3h}, ))' ) aT
      close(12)
   end subroutine write66

   !   contained in program\_that\_calls\_umat writes a 6x1 matrix  for debugging with Mma
   subroutine write6(a)
      implicit none
      real(8), dimension(6) :: a
      open(12,file='nic.m',access='append')
      write(12,'( 5hx6={ , 5(f15.4,2h,  ),f15.4, 3h}  )' ) a
      close(12)
   end subroutine write6


   !   contained in  program\_that\_calls\_umat converts D(3,3)  to stran(6)
   function map2stran(a,ntens)
      implicit none             !===converts D(3,3)  to stran(6) with $\gamma_{12} = 2 \epsilon_{12}$ etc.
      real(8), intent(in), dimension(1:3,1:3) :: a
      integer, intent(in) :: ntens
      real(8),  dimension(1:ntens) :: map2stran
      real(8), dimension(1:6) :: b
      b =[a(1,1),a(2,2),a(3,3),2*a(1,2),2*a(1,3),2*a(2,3)]
      map2stran(1:ntens)=b(1:ntens)
   end function map2stran

   !   contained in  program\_that\_calls\_umat converts strain rate from vector dstran(1:ntens) to  D(3,3)
   function map2D(a,ntens)
      implicit none
      real(8),  dimension(1:3,1:3) :: map2D
      integer, intent(in) :: ntens
      real(8), intent(in), dimension(:) :: a
      real(8),dimension(1:6) :: b = 0
      b(1:ntens) = a(1:ntens)
      map2D = reshape( [b(1), b(4)/2, b(5)/2, b(4)/2,b(2),b(6)/2, b(5)/2,b(6)/2, b(3)],[3,3] )
   end function map2D

   !   contained in  program\_that\_calls\_umat converts tensor T(3,3)  to matrix stress(ntens)
   function map2stress(a,ntens)
      implicit none
      real(8), intent(in), dimension(1:3,1:3) :: a
      integer, intent(in) :: ntens
      real(8),  dimension(1:ntens) :: map2stress
      real(8), dimension(1:6) :: b
      b = [a(1,1),a(2,2),a(3,3),a(1,2),a(1,3),a(2,3)]
      map2stress = b(1:ntens)
   end function map2stress

   !   contained in  program\_that\_calls\_umat converts matrix stress(1:ntens)  to tensor T(3,3)
   function map2T(a,ntens)
      implicit none
      real(8),  dimension(1:3,1:3) :: map2T
      integer, intent(in) :: ntens
      real(8), intent(in), dimension(:) :: a
      real(8), dimension(1:6) :: b= 0
      b(1:ntens) = a(1:ntens)
      map2T = reshape( [b(1),b(4),b(5), b(4),b(2),b(6),  b(5),b(6),b(3) ],[3,3] )
   end function map2T

   !   contained in  program\_that\_calls\_umat  reads a file with instructions for stress alignment
   subroutine  readAlignment(align, ImportFileName )
      implicit none
      character(len=40) ImportFileName, trunc, extension
      character(len=80) ReversalFileName
      logical ::  okSplit
      type(StressAlignment) :: align
      call splitaLine( ImportFileName ,'.',trunc, extension, okSplit )
      if(.not. okSplit) stop 'error  readAlignment FileName without . '
      reversalFileName = Trim(trunc) // 'rev'
      open(22, file=reversalFileName,status ='old', err=555 )
      align%active=.True.
      align%reversal(:) = 0
      read(22,*,err=556)  align%kblank,align%nrec,align%kReversal, align%ncol
      read(22,*,err=557)  align%reversal(1:align%kReversal)
      read(22,*,err=558)  align%isig(1:6)
      read(22,*,err=559)  align%sigFac(1:6)
      return
555   align%active=.False.
      return
556   stop 'error   readAlignment  cannot read kblank... '
557   stop 'error   readAlignment  cannot read reversal() '
558   stop 'error   readAlignment  cannot read sigCol() '
559   stop 'error   readAlignment  cannot read factor() '
   end subroutine  readAlignment

!   contained in  program\_that\_calls\_umat tries to align stress to values from aState(1:mImport)
   subroutine  tryAlignStress(align, kinc, aState, mImport,stress,ntens)
      implicit none
      integer:: mImport,kinc,ntens,ie
      real(8) :: aState(mImport)
      real(8) :: stress(ntens)
      type(StressAlignment) :: align

      if(.not. align%active) return
      if(.not. any(align%Reversal == kinc)) return

      ! only  stress components for which isig(ie) /= 0 will be aligned
      forall(ie=1:ntens, align%isig(ie) /= 0) stress(ie)= aState( align%isig(ie))*align%sigFac(ie)
      return
   end subroutine  tryAlignStress

end program that_calls_umat

