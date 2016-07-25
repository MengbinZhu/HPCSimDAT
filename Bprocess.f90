
!--transfer the array, before output
!--input:wrk(n2long,nlat)
!--output:wrk(n2long,nlat)
    subroutine output_trans(wrk)
	include 'sphectra.h'

	real wrk(nlong,nlat),wrk_out(nlong,nlat)
	wrk_out=0.
	wrk_out(1:nlong,1:nlat/2)=wrk(1:nlong,nlat:nlat/2+1:-1)
	wrk_out(1:nlong,nlat/2+1:nlat)=wrk(1:nlong,1:nlat/2)
	wrk=wrk_out
	return
    end subroutine


	!--transfer the array, before output
!--input:wrk(n2long,nlat)
!--output:wrk(n2long,nlat)
        subroutine trans(wrk)
        include 'sphectra.h'

        real wrk(nlong,nlat),wrk_out(nlong,nlat)
        wrk_out=0.
        wrk_out(1:nlong,1:nlat/2)=wrk(1:nlong,nlat/2+1:nlat)
        wrk_out(1:nlong,nlat/2+1:nlat)=wrk(1:nlong,nlat/2:1:-1)
        wrk=wrk_out
        return
        end subroutine
!------------------------right

!---transfer gpsi to ght
        subroutine gpsi2ght(wrk)
        include 'sphectra.h'
        include 'paramod.h'

        real wrk(nlong,nlat),wrk_out(nlong,nlat)
        do i=1,nlat/2
             wrk_out(:,i)=-wrk(:,i)*f0/g
        end do
        do i=nlat/2+1,nlat
             wrk_out(:,i)=wrk(:,i)*f0/g
        end do
        wrk=wrk_out
        return
        end subroutine

!integrate for nt_day given an initial state
!input: 
!--nt_hour:total hours of integration
!--nt_inter:interval of output
!--nt_times=nt_hour/nt_inter
!output:
!--out_gpsi(nlong,nlat,nlevels,0:nt_times)
        subroutine integrate(nt_hour,nt_inter,nt_times,out_gpsi)
        include 'sphectra.h'
        include 'paramod.h'       
 
        integer::nt_hour,nt_inter,nt_times
        real::out_gpsi(nlong,nlat,nlevels,0:nt_times)
        real::gpsi_wrk(n2long,nlat),spsi_wrk(nvaria2,nlevels)
        real::sq_wrk(nvaria2,nlevels),sq1_wrk(nvaria2,nlevels)
        integer::nl,nst,iflag
!------------------
        spsi_wrk=0.
        sq_wrk=0.
        do nl=1,nlevels
        gpsi_wrk=0.
        gpsi_wrk(1:nlong,:)=out_gpsi(1:nlong,:,nl,0)
        call gridtospec(spsi_wrk(:,nl),gpsi_wrk)
        end do    !input each level of gpsi to spsi
!--------------
        call psi2q(sq_wrk,spsi_wrk)  !spsi --> sq
!--------------
        iflag=0
        do nst=1,nt_hour
              call onestep(sq1_wrk,sq_wrk)
!
              if(mod(nst,nt_inter).eq.0.and.nst.gt.nsteq) then
              iflag=iflag+1
              call q2psi(spsi_wrk,sq1_wrk)
! 
              do nl=1,nlevels
                   call spectogrid(gpsi_wrk,spsi_wrk(:,nl))
                   out_gpsi(1:nlong,:,nl,iflag)=gpsi_wrk(1:nlong,:)
              end do
              endif
              call scopy(nvaria2*nlevels,sq1_wrk,1,sq_wrk,1)
        end do
        return
        end subroutine
        !------------------right----

!transfer 3-d array to 1-d array
        subroutine three2one(nx,ny,nl,x,nn,y)
        implicit none
        integer::nx,ny,nl,nn
        real::x(nx,ny,nl)
        real::y(nn)
        integer::flag
        integer::ii,jj,zz
        flag=0
        do zz=1,nl
             do jj=1,ny
                 do ii=1,nx
                   flag=flag+1
                   y(flag)=x(ii,jj,zz)
                 end do
             end do
        end do
        return
        end subroutine
!-------------------------------------------

!transfer 3-d array to 1-d array
        subroutine one2three(nn,y,nx,ny,nl,x)
        implicit none
        integer::nx,ny,nl,nn
        real::x(nx,ny,nl)
        real::y(nn)
        integer::flag
        integer::ii,jj,zz
        flag=0
        do zz=1,nl
             do jj=1,ny
                 do ii=1,nx
                   flag=flag+1
                   x(ii,jj,zz)=y(flag)
                 end do
             end do
        end do
        return
        end subroutine        


        !-----------------------------------------------
        !!!!!field rescale
        subroutine field_rescale_2d(nx,ny,nl,x,stan,y)
        implicit none
        integer::nx,ny,nl
        real :: x(nx,ny,nl),y(nx,ny,nl)
        real :: ave_error,stan
        call weight_rmse(.false.,nx,ny,nl,x,x,ave_error)
        y=x*stan/ave_error
        return
        end subroutine
        !-----------------------------------------right

!-----------------------------------------------
	!--regionally rescale-----
    !-input:
    !--smth_r:smoothing radius
    !--nx: ngrid on x direction
    !--ny: ngrid on y direction
    !--nl: number of vertical levels
    !--x(nx,ny,nl): orginal error field
    !--mask: used for regional rescaling
    !-output:
    !--y(nx,ny,nl)
    subroutine region_rescale(smth_r,nx,ny,nl,x,mask,y)
    implicit none
    integer::i,j,k
    integer::nx,ny,nl,nxy,smth_r,npoint
    real::x(nx,ny,nl),y(nx,ny,nl)
    real::xwrk(1-nx:2*nx,1-ny:2*ny)
    real::mask(nx,ny,nl)
    real::ewrk2(-smth_r:smth_r,-smth_r:smth_r)
    real, allocatable::ewrk1(:)
    npoint = (2*smth_r+1)*(2*smth_r+1)
    allocate (ewrk1(npoint))
    nxy = nx*ny

    do k = 1,nl
        xwrk=0.
        call output_trans(x(:,:,k))    !correct the order of latitude
        call extend_field(nx,ny,x(:,:,k),xwrk)  
        do j = 1,ny,smth_r*2+1
            !----------------
            do i = 1,nx,smth_r*2+1
            
            ewrk2=0.
            ewrk1=0.
            ewrk2(:,:)=xwrk(i-smth_r:i+smth_r,j-smth_r:j+smth_r)
            call two2one(2*smth_r+1,2*smth_r+1,ewrk2,npoint,ewrk1)
            call region_standard(npoint,ewrk1,mask(i,j,k)) 
            call one2two(npoint,ewrk1,2*smth_r+1,2*smth_r+1,ewrk2)
            xwrk(i-smth_r:i+smth_r,j-smth_r:j+smth_r) = ewrk2(:,:)

            end do
        end do
        y(:,:,k)=xwrk(1:nx,1:ny)
        call trans(y(:,:,k))
    end do
    deallocate(ewrk1)

    return
    end subroutine
!--------------------right-----


subroutine extend_field(mx,my,x,y)
implicit none
integer::mx,my,mm
real :: x(0:mx-1,0:my-1)
real :: y(-mx:2*mx-1,-my:2*my-1)
    y(0:mx-1,0:my-1)=x(0:mx-1,0:my-1)      !c
    y(-mx:-1,0:my-1)=x(0:mx-1,0:my-1)       !west
    y(mx:2*mx-1,0:my-1)=x(0:mx-1,0:my-1)    !east
    
    y(0:mx-1,-my:-1)=y(0+mx/2:mx-1+mx/2,my-1:0:-1)       !south
    y(0:mx-1,my:2*my-1)=y(0+mx/2:mx-1+mx/2,my-1:0:-1)    !north
    
    y(-mx:-1,-my:-1)=y(-mx+mx/2:-1+mx/2,my-1:0:-1)       !south west
    y(-mx:-1,my:2*my-1)=y(-mx+mx/2:-1+mx/2,my-1:0:-1)    !north west
    y(mx:2*mx-1,-my:-1)=y(mx-mx/2:2*mx-1-mx/2,my-1:0:-1)    !south east
    y(mx:2*mx-1,my:2*my-1)=y(mx-mx/2:2*mx-1-mx/2,my-1:0:-1) !north east
return
end subroutine
!---------------------------------------

        !-----------------------------------------------------
        !!!transfer two-dimension field to one-dimension series
        subroutine two2one(nx,ny,x,n,y)
        implicit none

        integer::nx,ny,n
        real ::x(nx,ny),y(n)
        integer::i,j,flag
        flag=0
        do j=1,ny
                do i=1,nx
                flag=flag+1
                y(flag)=x(i,j)
                end do
        end do
        return
        end subroutine
        !--------------------------------------------------right
 
        !----------------------------------------------------
        !!!transfer one-dimension series to two-dimen field
        subroutine one2two(n,x,nx,ny,y)
        implicit none
        
        integer :: n,nx,ny
        real :: x(n),y(nx,ny)
        integer::flag,i,j
        flag=0
        do j=1,ny
                do i=1,nx
                flag=flag+1
                y(i,j)=x(flag)
                end do
        end do
        return
        end subroutine
        !-----------------------------------------------right

        !!!!make standerlised process on x(n), the final magnitude is dl.The output is y(n)
        subroutine region_standard(n,x,stan)
        implicit none
        
        integer::n,i
        real :: x(n),y(n)
        real :: stan,dd
        dd=0.0
        do i=1,n
                dd=dd+x(i)*x(i)
        end do  
        dd=dd/real(n)
        dd=sqrt(dd)
        
            if (dd>=1.0*stan) then  !if larger than stan, than normalized
            x = x*stan*1.0/dd
                       !else unchanged     
            end if
        return
        end subroutine
        !-------------------------------------------------right

        !---calculate the rms error (weighted)---
        !ang0rad::angle or radian (radian then true)
        !lat: the latitude
        subroutine weight_rmse(ang0rad,nx,ny,nl,error1,error2,ave_error)
        include 'sphectra.h'
        include 'paramod.h'      
 
        logical::ang0rad
        integer::nx,ny,nl
        real::wgt,total_wgt
        real::lat(ny)
        real::error1(nx,ny,nl),error2(nx,ny,nl)
        real::error1_temp(nx,ny,nl),error2_temp(nx,ny,nl)
        real::ave_error,summ
        integer::ix,iy,il
        summ=0.
        total_wgt=0.
        ave_error=0.
!------------------------
        error1_temp=error1*f0/g
        error2_temp=error2*f0/g
        
        if (ang0rad==.false.) then
              lat=lat_angle/180.0*PI
        end if

        do il=1,nl
            do iy=1,ny
                    wgt=cos(lat(iy))
                do ix=1,nx
                    summ=summ+error1_temp(ix,iy,il)*error2_temp(ix,iy,il)*wgt
                    total_wgt=total_wgt+wgt
                end do
            end do
        end do
        ave_error=sqrt(summ/total_wgt)
        ave_error=ave_error*g/f0
        return
        end subroutine
!---------------------------------right!

!---calculate the cross product (weighted)---
        !ang0rad::angle or radian (radian then true)
        !lat: the latitude
        subroutine weight_cross_product(ang0rad,nx,ny,nl,error1,error2,ave_error)
        include 'sphectra.h'
        include 'paramod.h'

        logical::ang0rad
        integer::nx,ny,nl
        real::wgt,total_wgt
        real::lat(ny)
        real::error1(nx,ny,nl),error2(nx,ny,nl)
        real::error1_temp(nx,ny,nl),error2_temp(nx,ny,nl)
        real::ave_error,summ
        integer::ix,iy,il
        summ=0.
        total_wgt=0.
        ave_error=0.
!------------------------
        error1_temp=error1*f0/g
        error2_temp=error2*f0/g

        if (ang0rad==.false.) then
              lat=lat_angle/180.0*PI
        end if

        do il=1,nl
            do iy=1,ny
                    wgt=cos(lat(iy))
                do ix=1,nx
                    summ=summ+error1_temp(ix,iy,il)*error2_temp(ix,iy,il)*wgt
                    total_wgt=total_wgt+wgt
                end do
            end do
        end do
        ave_error=summ/total_wgt
        ave_error=ave_error*g*g/(f0*f0)
        return
        end subroutine
        !-----------------right------------------

!-----orthoginalize the perts using GSR
!input: 
!---error(dimen,nl,nper)
!---or2or3: if .true., orthogonalizing each level; 
!          .false. orthogonalizing all levels
!output:
!---errorout(dimen,nl,nout)
!---nout<=nper
subroutine gsr_gb(or2or3,nx,ny,nl,nper,nout,error,stan,errorout)  !if the global scale is used
!change the gsr(or2or3,nx,ny,nl,nper,nout,stan,error,errorout)  stan(nlevels)
implicit none
logical::or2or3
integer::nx,ny,nl,nper,nout
integer::i,j,il
real::error(nx,ny,nl,nper)
real::errorout(nx,ny,nl,nout)
real::stan(nl)
real::yout(nx,ny,nl)
!-----------------
 errorout(:,:,:,:)=error(:,:,:,1:nout)
!-----------------
if (or2or3) then   !orthogonalization each level
    do il=1,nl
        do i=2,nout
            do j=1,i-1
			yout=0.
            call project(nx,ny,1,errorout(:,:,il,j),&
			    error(:,:,il,i),yout(:,:,1))
            errorout(:,:,il,i)=errorout(:,:,il,i)-yout(:,:,1)
            end do
        end do
    end do
else            !orthogonalization three levels
    do i=2,nout
        do j=1,i-1
		yout=0.
        call project(nx,ny,nl,errorout(:,:,:,j),&
		    error(:,:,:,i),yout(:,:,:))
		errorout(:,:,:,i)=errorout(:,:,:,i)-yout(:,:,:)
        end do
    end do
end if
!-----rescale------
do i=1,nout
    do il=1,nl
    call field_rescale_2d(nx,ny,1,errorout(:,:,il,i),&
    	stan(il),errorout(:,:,il,i))
    end do
end do

return
end subroutine
!-------------------------------------------------

!------------------------------
 !independent unit vectors
subroutine onebase(m,n,bs)
implicit none
integer::m,n
real::bs(m,n),x(m)
real,parameter :: norm=1.0
integer::i
!------------
do i = 1,n
    call error_R(norm,m,x)
    bs(:,i) = x(:)
end do
call gsr_stan(m,n,bs,norm)
return
end subroutine
!----------------------------------

        !------------------------------------------------
        !!!!form random error of random distribution
        !!!stan is the rescaling size
        subroutine error_R(stan,n,x)
        implicit none
        real::a,stan
        integer::i,n
        real::xtemp(n),x(n)
        
        do i=1,n
                call random_number(a)
                a=a-0.5
                xtemp(i)=a
        end do
        x = xtemp

        !call standard(n,xtemp,stan,x) 
        return
        end subroutine
        !---------------------------------------------right

        !!!input: 
!!!--dimen : dimension
!!!--error(dimen,dimen): the error vectors to be GSR
!!!output:
!!!--errorout(dimen,dimen): output result
subroutine gsr_stan(m,n,error,dl)
implicit none
integer::i,j,k,m,n
real::dl
real::yout(m)
real::error(m,n),errorout(m,n)
do i=1,m
        do j=1,n
                errorout(i,j)=error(i,j)
        end do
end do
do i=2,n
        do j=1,i-1
                call project_1d(m,errorout(:,j),error(:,i),yout)
                do k=1,m
                        errorout(k,i)=errorout(k,i)-yout(k)
                end do
        end do
end do

do i=1,n
        call standard(m,errorout(:,i),dl,yout)
        error(:,i)=yout(:)
end do
return
end subroutine
!-------------------------------------------------

!-----------------------------------------------------------
!!!calculate projection of y on x
!----output:
!-----yout(dimen): projection of y on x
subroutine project_1d(m,x,y,yout)
implicit none
integer::i,m
real::x(m),y(m),yout(m)
real::d2,d1
d1=0.
d2=0.
do i=1,m
        d2=d2+x(i)*y(i)
        d1=d1+x(i)*x(i)
end do
do i=1,m
        yout(i)=d2*x(i)/d1
end do
return
end subroutine
!---------------------------------------------------

!!!------------------------------------------------------------------------
        !!!!make standerlised process on x(n), the final magnitude is dl.The output is y(n)
        subroutine standard(n,x,stan,y)
        implicit none

        integer::n,i
        real::x(n),y(n)
        real::stan,dd
        dd=0.0
        do i=1,n
                dd=dd+x(i)*x(i)
        end do
        dd=dd/real(n)
        dd=sqrt(dd)
        do i=1,n
                y(i)=x(i)*stan/dd
        end do
        return
        end subroutine
        !-------------------------------------------------right

!-----------------------------------------------------------
!!!calculate projection of y on x
!----output:
!-----yout(dimen): projection of y on x
subroutine project(nx,ny,nl,x,y,yout)
implicit none
integer::nx,ny,nl
real::x(nx,ny,nl),y(nx,ny,nl),yout(nx,ny,nl)
real::d2,d1
d1=0.
d2=0.
call weight_cross_product(.false.,nx,ny,nl,x,x,d1)
call weight_cross_product(.false.,nx,ny,nl,x,y,d2)
yout=d2*x/d1
return
end subroutine
!---------------------------------------------------


!!=========calculate the ensemble covariace matrix of x
!input:
!---m:the number of ensmeble samples
!---n:the number of variables
!---x(m,n):original matrix
!output:
!---cov_M(n,n):covariance matrix
subroutine cov_matrix(nn,ne,x,cov_x)
implicit none
    integer::nn,ne,i,j
	real::x(nn,ne),y(nn,ne)
	real::cov_x(nn,nn)
	real::cov
	real::a1(ne),a2(ne)
!-----------------	
    do i=1,nn
        do j=i,nn
        cov=0.
		a1(:)=x(i,:)
		a2(:)=x(j,:)
        call covariance(ne,a1,a2,cov)
        cov_x(i,j)=cov
        cov_x(j,i)=cov
        end do
        y(i,:)=a1(:)
    end do
    x=y
    return
    end subroutine
	!!========================right===================


!calculate the average,variance,standard deviation!
    subroutine meanvar(m,x,ax,vx,mx)
    implicit none
	integer::m,i
	real::x(m)
	real::ax,vx,mx
	real::summ
	summ=0.0
    do i=1,m
        	summ=summ+x(i)
    end do
	ax=summ/real(m)
	summ=0.0
    do i=1,m
		summ=(x(i)-ax)**2+summ
    end do
	vx=summ/real(m)
	mx=sqrt(vx)
    return
    end subroutine
!--------------------------right
	
!!====================calculate the covariance
!input:
!---m:the number of series
!---a1,a2:two series
!output:
!---variance:covariance
subroutine covariance(m,a1,a2,cov)
implicit none
	integer::m,i
	real::a1(m),a2(m)
	real::ave1,ave2
	real::cov,summ
	summ=0.0
	ave1=sum(a1)/real(m)
    ave2=sum(a2)/real(m)
	
    a1(:)=a1(:)-ave1
    a2(:)=a2(:)-ave2
        
    a1(:)=a1(:)*(1.0+0.28)  !inflation factor
    a2(:)=a2(:)*(1.0+0.28)
	
    do i=1,m
		summ=summ+a1(i)*a2(i)
    end do
	cov=summ/real(m-1)
    
    a1(:)=a1(:)+ave1
    a2(:)=a2(:)+ave2
    return
    end subroutine
	!!==========================right===============




!-----orthoginalize the perts using GSR
!input: 
!---error(dimen,nl,nper)
!---or2or3: if .true., orthogonalizing each level; 
!          .false. orthogonalizing all levels
!output:
!---errorout(dimen,nl,nout)
!---nout<=nper
subroutine gsr_region(or2or3,nx,ny,nl,nper,nout,error,errorout)  !if the global scale is used
!change the gsr(or2or3,nx,ny,nl,nper,nout,stan,error,errorout)  stan(nlevels)
implicit none
logical::or2or3
integer::nx,ny,nl,nper,nout
integer::i,j,il
real::error(nx,ny,nl,nper)
real::errorout(nx,ny,nl,nout)
real::stan(nl)
real::yout(nx,ny,nl)
!-----------------
 errorout(:,:,:,:)=error(:,:,:,1:nout)
!-----------------
if (or2or3) then   !orthogonalization each level
    do il=1,nl
        do i=2,nout
            do j=1,i-1
			yout=0.
            call project(nx,ny,1,errorout(:,:,il,j),&
			    error(:,:,il,i),yout(:,:,1))
            errorout(:,:,il,i)=errorout(:,:,il,i)-yout(:,:,1)
            end do
        end do
    end do
else            !orthogonalization three levels
    do i=2,nout
        do j=1,i-1
		yout=0.
        call project(nx,ny,nl,errorout(:,:,:,j),&
		    error(:,:,:,i),yout(:,:,:))
		errorout(:,:,:,i)=errorout(:,:,:,i)-yout(:,:,:)
        end do
    end do
end if
!-----rescale------
!do i=1,nout
!    do il=1,nl
!    call field_rescale_2d(nx,ny,1,errorout(:,:,il,i),&
!    	stan(il),errorout(:,:,il,i))
!    end do
!end do

return
end subroutine
!-------------------------------------------------

!!!-------------------------------------------------------
!!This subroutine rank x(n) according to their size.
!--input:
!----x(n):original array
!--output:
!----y(n):final array from largest to smallest.
!----array:final sequence number.
subroutine order(x,n,y,array)
implicit none   
integer::n,ne,i,j    
integer::array(n)
real::x(n),y(n),t
y(:)=x(:)    
do i=1,n-1    
    do j=1,n-i
    if (y(j)<y(j+1)) then
        t=y(j+1)
        y(j+1)=y(j)    
        y(j)=t    
!!---------------------
        ne=array(j+1)
        array(j+1)=array(j)
        array(j)=ne
    end if
    end do
end do
return
end subroutine
!---------------------------------------------------------------
