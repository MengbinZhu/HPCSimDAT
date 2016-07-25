      subroutine temper(temp,psi)
!c
!c Computes temperatures from the layers thickness. Returns a 
!C three levels field with T at the 200-500 level, T at the 
!c 500-800 level and T at 1000mb as a fraction of that at the 
!c above level.
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real psi(nvaria2,nlevels),costa(nvaria2)
      real temp(nvaria2,nlevels),gt(n2long,nlat)
      real fco(nlat), fcof(nlat,n2long),sfcof(nvaria2)
      real phi(nvaria2,nlevels),prof(32)
!c
      data prof/ 25.0, 25.0, 25.0, 24.2, 23.4, 22.7, 21.9, 21.1, 20.4, &
                19.6, 18.8, 18.1, 17.3, 16.5, 15.8, 15.0,             &
                25.0, 25.0, 25.0, 24.2, 23.4, 22.7, 21.9, 21.1, 20.4, &
                19.6, 18.8, 18.1, 17.3, 16.5, 15.8, 15.0/             
!c
!c coefficient -1/R 
!c 
!c      coe=-f0/(RAIR)
      coe=-1./rair
!c
!c compute phi at the three levels
!c
!c  linear balance
      do nl=1,nlevels
         call psi2phi(phi(1,nl),psi(1,nl))         
      enddo
!c  add prescribed spatial mean (from US standard)
      phi(243,3)=(1949.*g)
      phi(243,2)=(5574.*g)
      phi(243,1)=(11784.*g)         
!c      print*,phi
!c
!c level 500-800
!c
      coe1=alog(8.)-alog(5.)
!c      Print*,coe1,coe,coe1*coe,cost
      call scopy(nvaria2,phi(1,3),1,temp(1,2),1)
      call saxpy(nvaria2,-1.,phi(1,2),1,temp(1,2),1)
!c      print*,phi(243,2),phi(243,3),temp(243,2),coe/coe1
      call sscal(nvaria2,coe/coe1,temp(1,2),1)
!c
!c level 200-500
!c
      coe1=alog(5.)-alog(2.)
      call scopy(nvaria2,phi(1,2),1,temp(1,1),1)
      call saxpy(nvaria2,-1.,phi(1,1),1,temp(1,1),1)
      call sscal(nvaria2,coe/coe1,temp(1,1),1) 
!c
!c reduction to 1000mb temperature, as propotional to T at 650 mb
!c     
      call scopy(nvaria2,temp(1,2),1,temp(1,3),1)
      call spectogrid(gt,temp(1,3))
!c      rk=0.09
      do ila=1,nlat
         do ilo=1,nlong
!c             print*,gt(ilo,ila),gt(ilo,ila)*(1.+rk),(1.+rk),rk
             gt(ilo,ila)=gt(ilo,ila)*(1.+rktemp)
         enddo
      enddo

!ccc      temp(243,3)=temp(243,2)+24.     
!ccc      call sscal(nvaria2,(1.+k),temp(1,3),1)

      call gridtospec(temp(1,3),gt)
!c
!c      do nl=1,nlevels
!cc         temp(243,nl)=temp(243,nl)-273.15
!c         print*,temp(243,nl)
!c      enddo
!c
      return
!c
      end
!c
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!c
!c
      Subroutine psi2phi(phi,psi)
!c
!c computes the geopotential from streanfunction by solving the linear 
!c balance equation. NO INTEGRATION CONSTANT ADDED (do it yourself).
!c
      include 'sphectra.h'
      include 'paramod.h'
!c      
      real phi(nvaria2),psi(nvaria2),lpsi(nvaria2)
      real fco(nlat),fcof(n2long,nlat),sfcof(nvaria2)
      real phi1(nvaria2),phi2(nvaria2),wrk(n2long,nlat)
      real wrk1(n2long,nlat)
!c
      data fco/                  &
        0.0000070,0.0000211,0.0000349,0.0000484,0.0000615,    &
        0.0000739,0.0000857,0.0000967,0.0001068,0.0001159,0.0001239, &
         0.0001307,0.0001363,0.0001407,0.0001437,0.0001454,   &
        -0.0000070,-0.0000211,-0.0000349,-0.0000484,-0.0000615,  &
        -0.0000739,-0.0000857,-0.0000967,-0.0001068,-0.0001159,-0.0001239,&
           -0.0001307,-0.0001363,-0.0001407,-0.0001437,-0.0001454/
!c
!c
!c  definition of f
!c
      do ila=1,nlat
         do ilo=1,nlong
            fcof(ilo,ila)=fco(ila)
         enddo
      enddo
      call gridtospec(sfcof,fcof)
!c
!c term f lap(psi)
!c
      call laplace(lpsi,psi,2)
      call spectogrid(wrk,lpsi)
      call shprod(n2long*nlat,fcof,1,wrk,1,wrk1,1)
      call gridtospec(phi1,wrk1)
!c      print*,(phi1(i),i=1,10) 
!c
!c term grad(f)*grad(psi)
!c
      call progras(phi2,sfcof,psi)
!c      print*,(phi2(i),i=1,10)
!c      stop
!c
!c sum and inverse laplacian
!c
      call saxpy(nvaria2,1.,phi1,1,phi2,1)
      call laplace(phi,phi2,-2)
!c      print*,(phi2(i),i=1,10)
!c      print*,(phi(i),i=1,10)
!c      
      return
      end
!c
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine wind(u,v,psi)
!c
!c  Computes the wind from streamfunction.
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real psi(nvaria2,nlevels),u(nvaria2,nlevels)
      real v(nvaria2,nlevels),chi(nvaria2)
!c
      call sset(nvaria2,0.,chi,1)
!c
!c  compute wind at the three levels
!c
      do nl=1,nlevels
         call calvent(u(1,nl),v(1,nl),psi(1,nl),chi,20)
      enddo
!c
      return
      end
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c
      real function ssum(n,sx,incx)
      integer n,incx
      real sx(1)
      integer ix,iy,i

      ssum = 0.0e0
      sum=0.
      do ii=incx,n
         sum=sum+sx(ii)
      enddo
      ssum=sum
      return
      end
