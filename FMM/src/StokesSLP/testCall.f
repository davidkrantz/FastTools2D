      program testCall
      

      implicit real*8 (a-h,o-z)

      parameter (ns = 2**10)
      parameter (nt = ns)

      dimension xs(ns),ys(ns)
      dimension q1(ns),q2(ns)
      dimension xt(nt),yt(nt)
      dimension u1(nt),v1(nt)
      dimension u1Exact(nt),v1Exact(nt)
      dimension u2(nt),v2(nt)
      dimension u2Exact(nt),v2Exact(nt)

      twopi = 8.d0*datan(1.d0)

      write(*,*) "Starting test program..."

      dx = twopi/dble(ns)
      write(*,*) 'ns = ', ns
      write(*,*) 'dx = ', dx

      do i=1,ns
        theta = dble(i-1)*dx
        xs(i) = dcos(theta)
        ys(i) = dsin(theta)
        q1(i) = 1.d0/dble(ns)
        q2(i) = 1.d0/dble(ns)
      enddo
      dx = twopi/dble(nt)
      do i=1,nt
        theta = dble(i-1)*dx
        xt(i) = dcos(theta) + 2.d0
        yt(i) = dsin(theta) + 1.d0
      enddo
	
      write(*,*) 'Calling FMM...'
      call system_clock(iclock1)
      call stokesSLP(q1,q2,xs,ys,ns,xs,ys,ns,1.d0,6.d0,u1,u2)
      call system_clock(iclock2,irate)
      write(*,*) '***********'
      write(*,*) 'FMM CPU TIME IS'
      write(*,1000) real(iclock2 - iclock1)/real(irate)
      write(*,*) '***********'

      write(*,*) 'Calling direct sum...'
      call system_clock(iclock1)
      call stokesSLPdirect(q1,q2,xs,ys,ns,xs,ys,ns,1,u1Exact,u2Exact)
      call system_clock(iclock2,irate)
      write(*,*) '***********'
      write(*,*) 'DIRECT CPU TIME IS'
      write(*,1000) real(iclock2 - iclock1)/real(irate)
      write(*,*) '***********'

      r_err = 0.d0
      c_err = 0.d0
      do i = 1,ns
        r_err = r_err + (u1(i) - u1Exact(i))**2
        c_err = c_err + (u2(i) - u2Exact(i))**2
      enddo
      write(*,*) '***********'
      write(*,*) 'COMPONENT 1 ERROR'
      write(*,1000) r_err
      write(*,*) 'COMPONENT 2 ERROR'
      write(*,1000) c_err
      write(*,*) '***********'





c      call cpu_time(t0)
c      call stokesSLP(q1,q2,xs,ys,ns,xt,yt,nt,0,u1,u2)
c      call cpu_time(t1)
c      print*,t1-t0
cc      call cpu_time(t0)
cc      call stokesSLPdirect(q1,q2,xs,ys,ns,xt,yt,nt,0,u1Exact,u2Exact)
cc      call cpu_time(t1)
c      print*,t1-t0
c
c      r_err = 0.d0
c      c_err = 0.d0
c      do i = 1,ns
c        r_err = r_err + (u1(i) - u1Exact(i))**2
c        c_err = c_err + (u2(i) - u2Exact(i))**2
c      enddo
c      print*,r_err
c      print*,c_err

1000  format(ES12.4)
      end





