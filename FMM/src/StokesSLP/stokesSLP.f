      subroutine stokesSLP(s1,s2,xs,ys,ns,xt,yt,nt,riequal,riprec,
     $      u1,u2)
c     Evaluates the single-layer potential for Stokes' equation
c     (s1,s2) is the source strength of length ns
c     xs and ys are the source points and xt and yt are the target 
c     locations of length ns
c     iequal is 1 if sources .eq. targest
c     iequal is 0 if sources .ne. targets
c     iprec controls the accuracy of the FMM
c         -2 => tolerance =.5d0
c         -1 => tolerance =.5d-1
c          0 => tolerance =.5d-2
c          1 => tolerance =.5d-3
c          2 => tolerance =.5d-6
c          3 => tolerance =.5d-9
c          4 => tolerance =.5d-12
c          5 => tolerance =.5d-15
c     (u1,u2) are the two components of the velocity field
      implicit real*8 (a-h,o-z)

      integer error
      real *8 s1(ns),s2(ns)
c     strength of Stokes SLP
      real *8 xs(ns),ys(ns)
c     location of the source points
      real *8 xt(nt),yt(nt)
c     location of the target points
      real *8 u1(nt),u2(nt)
c     x and y components of the velocity field

      real *8, allocatable :: sources(:,:)
c     location of sources
      real *8, allocatable :: targets(:,:)
c     location of targets
      real *8, allocatable :: charges(:)
c     charge strength of single-layer potential term
      real *8, allocatable :: dipstr(:),dipvec(:,:)
c     charge strength and direction of the 
c     double-layer potential term
      real *8, allocatable :: pot(:)
c     room for two components that have to be summed
c     to form the velocity field

      iequal = int(riequal)
      iprec = int(riprec)

      if (iequal .eq. 0) then
        allocate(targets(2,nt),stat=error)
      endif

      allocate(sources(2,ns),stat=error)
      allocate(charges(ns),stat=error)
      allocate(dipstr(ns),stat=error)
      allocate(dipvec(2,ns),stat=error)
      allocate(pot(nt),stat=error)
c     allocate memory for temporary variables

      if (iequal .eq. 1) then
        ifpot = 1 ! need the potential since sources .eq. targets
      else
        ifpot = 0 ! don't need the potential since sources ~= targets
        ifpottarg = 1 ! need the potential at the targets
      endif
      ifgrad = 0 ! don't need the gradient
      ifhess = 0 ! don't need the Hessian
      ifgradtarg = 0 ! don't need the gradient at the targets
      ifhesstarg = 0 ! don't need the Hessian at the targets
c     set flags for what output values we require

      do i=1,ns
        sources(1,i) = xs(i)
        sources(2,i) = ys(i)
      enddo
c     set source locations
      if (iequal .ne. 1) then
        do i=1,nt
          targets(1,i) = xt(i)
          targets(2,i) = yt(i)
        enddo
      endif
c     set target locations


      ifcharge = 1 ! need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        charges(i) = -s1(i)
        dipstr(i) = xs(i)
        dipvec(1,i) = s1(i)
        dipvec(2,i) = s2(i)
      enddo

      if (iequal .eq. 1) then
c        call system_clock(iclock1)
        call rfmm2dpart(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,u1,ifgrad,trash,ifhess,trash)
c        call system_clock(iclock2,irate)
c        write(6,*) '****************'
c        write(6,*) 'FMM CALL 1'
c        write(6,1000) real(iclock2 - iclock1)/real(irate)
      else
        call rfmm2dparttarg(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,trash,ifgrad,thras,ifhess,trash,
     $     nt,targets,ifpottarg,u1,ifgradtarg,trash,
     $     ifhesstarg,trash)
      endif
c     compute the first part of the first component in the Stokes SLP

      do i=1,ns
        charges(i) = -s2(i)
        dipstr(i) = ys(i)
      enddo

      if (iequal .eq. 1) then
c        call system_clock(iclock1)
        call rfmm2dpart(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,u2,ifgrad,trash,ifhess,trash)
c        call system_clock(iclock2,irate)
c        write(6,*) '****************'
c        write(6,*) 'FMM CALL 2'
c        write(6,1000) real(iclock2 - iclock1)/real(irate)
      else
        call rfmm2dparttarg(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,trash,ifgrad,trash,ifhess,trash,
     $     nt,targets,ifpottarg,u2,ifgradtarg,trash,
     $     ifhesstarg,trash)
      endif
c     compute the first part of the second component in the Stokes SLP


c     START OF COMPUTING THE SECOND PART OF THE STOKES SLP
      ifcharge = 0 ! don't need a charge component

      do i=1,ns
        dipstr(i) = -1.d0
      enddo

      if (iequal .eq. 1) then
c        call system_clock(iclock1)
        call rfmm2dpart(ierr,iprec,ns,sources,
     $     ifcharge,trash,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,trash,ifhess,trash)
c        call system_clock(iclock2,irate)
      else
        call rfmm2dparttarg(ierr,iprec,ns,sources,
     $     ifcharge,trash,ifdipole,dipstr,dipvec,
     $     ifpot,trash,ifgrad,trash,ifhess,trash,
     $     nt,targets,ifpottarg,pot,ifgradtarg,trash,
     $     ifhesstarg,trash)
      endif
c     compute the second part of the first and second  component in the
c     Stokes SLP
c      write(6,*) '****************'
c      write(6,*) 'FMM CALL 3'
c      write(6,1000) real(iclock2 - iclock1)/real(irate)
c      write(6,*) '****************'


c      call system_clock(iclock1)
      do i=1,nt
        u1(i) = u1(i) + xt(i)*pot(i)
        u1(i) = 7.957747154594767e-2*u1(i)
        u2(i) = u2(i) + yt(i)*pot(i)
        u2(i) = 7.957747154594767e-2*u2(i)
      enddo
c     constant is 1/(4*pi) which is the constant multiplying the SLP
c      call system_clock(iclock2,irate)
c      write(6,*) '****************'
c      write(6,*) 'ADDING TERMS AND SCALING'
c      write(6,1000) real(iclock2 - iclock1)/real(irate)
c      write(6,*) '****************'
          

c 1000 format(ES12.2)

      end
 
