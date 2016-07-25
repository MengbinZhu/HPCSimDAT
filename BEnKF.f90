
!assimilate from the initial state id0 for two months!
!1-ob-2ob-ob-3ob.....
subroutine barotropic_ENKF(M_mdl2ob,start_gpsi)
use my_inverse
include 'sphectra.h'
include 'paramod.h'

    integer::i,j,ii,jj,k,flag,iflag,inumber
    real::gpsi_std(nlong,nlat,nlevels)
    real::covini_constant(nlevels,nlevels)
    real::ocov_constant(nlevels,nlevels)
    real::random_mtrx(nlatlonglv,nens+1)
    !real::random_mtrx(allnstn,nens+1)
    real::backfc_field_3d(nlong,nlat,nlevels,nens)
    real::backfc_field_1d(nlatlonglv,nens)
    real::BCOV(nlatlonglv,nlatlonglv)
    real::OCOV(allnstn,allnstn)
    !real(8)::BOCOV_8kind(allnstn,allnstn)
    real::BOCOV(allnstn,allnstn)
    real::WCOV(nlatlonglv,allnstn)
    real::true_mdl_3d(nlong,nlat,nlevels)
    real::true_mdl_1d(nlatlonglv),true_stn(allnstn)
    real::obser_stn(allnstn),obser_per(allnstn)
    real::xana_1d(nlatlonglv)
    real::xana_field(nlong,nlat,nlevels,nens)
    real::work_field(nlong,nlat,nlevels,0:1)
    real::cov
    !inumber=(iflag-1)*nassim_times
    !write(*,*)inumber
!------------------------------------------------
    data covini_constant /390.,80.,20.,80.,120.,20.,20.,20.,50./
    data ocov_constant   /260.,80.,20.,80.,80.,20.,20.,20.,30./
    !ocov_constant=ocov_constant*
    !==========generate first guess states===========
    !---generate the initial forecast-error covariance--------
    !random pers for first guess        
    true_mdl_3d = start_gpsi
    random_mtrx=0.
    call genVGaussRandom(nlevels,nlatlonglv,0.,covini_constant,nens,random_mtrx(:,1:nens))
    random_mtrx=random_mtrx*g/f0   !unit: h--->fi

    !first guess=true state+ random pers
    do i=1,nens
        call one2three(nlatlonglv,random_mtrx(:,i),&
            nlong,nlat,nlevels,backfc_field_3d(:,:,:,i))
        backfc_field_3d(:,:,:,i)=true_mdl_3d(:,:,:)+backfc_field_3d(:,:,:,i)
    end do

    !================================================
    !===============generate observation variance=========
    OCOV=0.
    do i=1,nlevels
        do j=i,nlevels
            do k=1,nstn
            ii=k+(i-1)*nstn
            jj=k+(j-1)*nstn
            OCOV(ii,jj)=ocov_constant(i,j)
            OCOV(jj,ii)=ocov_constant(i,j)
            end do
        end do
    end do
    !OCOV=OCOV*g*g/(f0*f0)
    !----------------------
    !===============assimilation=====================
    do i=0,nassim_times
        write(*,*) "ASSIM TIME: ",i

        if (i==0) then
            xana_field = backfc_field_3d

        else

        call cov_matrix(nlatlonglv,nens,backfc_field_1d,BCOV)  !get Pb
        !call schur(BCOV)
        BCOV = BCOV*f0*f0/(g*g)
        !write(*,*)BCOV(1,1:10)
        !pause
        !-------H*Pb*HT+OB-------
        BOCOV = MATMUL(MATMUL(M_mdl2ob,BCOV),TRANSPOSE(M_mdl2ob))+OCOV
        !----------BOCOV(-1)-----
        !BOCOV_8kind = BOCOV
        call inverse(BOCOV,flag)
        !BOCOV = real(BOCOV_8kind)
        !-------Pb*HT*(BOCOV(-1))---
        WCOV = 0.
        WCOV = MATMUL(BCOV,TRANSPOSE(M_mdl2ob))
        WCOV = MATMUL(WCOV,BOCOV)
        !-------integrate the true state forward
        work_field=0.
        work_field(:,:,:,0)=true_mdl_3d(:,:,:)
        call integrate(bcycle,bcycle,1,work_field)
        true_mdl_3d(:,:,:)=work_field(:,:,:,1)

        call three2one(nlong,nlat,nlevels,&
            true_mdl_3d,nlatlonglv,true_mdl_1d)

        true_stn(:)= MATMUL(M_mdl2ob,true_mdl_1d)
        !--------------------------
        random_mtrx=0.
        call genVGaussRandom(nlevels,allnstn,0.,&
            ocov_constant,nens+1,random_mtrx(1:allnstn,:))
        
        random_mtrx=random_mtrx*g/f0
        obser_stn(:)=true_stn(:)+random_mtrx(1:allnstn,nens+1)

        do j=1,nens
            obser_per(:)=obser_stn(:)+&
                random_mtrx(1:allnstn,j)  !observation perturbation
   
            xana_1d(:)=backfc_field_1d(:,j)&
                +MATMUL(WCOV,(obser_per(:)-&
                MATMUL(M_mdl2ob,backfc_field_1d(:,j)) ) )

            call one2three(nlatlonglv,xana_1d,&
                nlong,nlat,nlevels,xana_field(:,:,:,j))

        end do

        !write(121,rec=i+inumber)xana_field
        !write(111,rec=i+inumber)true_mdl_3d
                        !===============================================
            if (i==nassim_times) exit
        end if

        !!--------------integration forward---------------------
        backfc_field_3d = 0.
        backfc_field_1d = 0.
        do j=1,nens
            work_field=0.
            work_field(:,:,:,0)=xana_field(:,:,:,j)
            call integrate(bcycle,bcycle,1,work_field)
            backfc_field_3d(:,:,:,j)=work_field(:,:,:,1)
        !store in backfc_field_1d        
            call three2one(nlong,nlat,nlevels,&
                backfc_field_3d(:,:,:,j),nlatlonglv,backfc_field_1d(:,j))
        end do
    end do

        !do i=1,mm
        !        call mean(nens,x_a(:,i),emean_a(i))
        !end do
        !-----------------------------------------------------------------------
                
        return
        end subroutine
        !---------------------------------------
            
!get the matrix from model space to observation space
subroutine mdl2ob(M_mdl2ob)

include 'sphectra.h'
include 'paramod.h'

integer::id,i,j,ii,jj
real::mdl_lon(nlong),mdl_lat(nlat)
real::ob_lat(nstn),ob_lon(nstn)
real::wht_2d(nlong,nlat),wht_1d(nlatlong)
integer::nxl,nxr,nyb,nyt
real::lonl,lonr,latb,latt
real::w11,w12,w21,w22
!==w12-----w22
!==
!==w11-----w21
M_mdl2ob=0.
!---input the latlong of model grid-----
do i=1,nlong
    mdl_lon(i)=0.+5.625*(i-1)
end do

mdl_lat((nlat/2+1):nlat)=lat_angle(1:nlat/2)
mdl_lat(1:nlat/2)=lat_angle(nlat:(nlat/2+1):-1)

!--input the latlong of observation-----
open (60,file='./data/lonlat_stn.grd',&
     access='direct',form='unformatted',recl=nstn)
read(60,rec=1)ob_lon
read(60,rec=2)ob_lat
!---------------------------------------
do i=1,nstn
    wht_1d=0.
    wht_2d=0.
    call search_nearest_point(ob_lon(i),ob_lat(i),nlong,mdl_lon,&
        nlat,mdl_lat,nxl,nxr,nyb,nyt)

    lonl=mdl_lon(nxl)
    lonr=mdl_lon(nxr)
    latb=mdl_lat(nyb)
    latt=mdl_lat(nyt)

    if (nxr==1) lonr=360.

    call Bilinear_interp_weight(lonl,lonr,&
        latb,latt,ob_lon(i),ob_lat(i),w11,w12,w21,w22)

    !----grid 1,1  nxl,nyb
    wht_2d(nxl,nyb)=w11
    !----grid 1,2  nxl,nyt
    wht_2d(nxl,nyt)=w12
    !----grid 2,1  nxr,nyb
    wht_2d(nxr,nyb)=w21
    !----grid 2,2  nxr,nyt
    wht_2d(nxr,nyt)=w22
    !-----------------------------
    call trans(wht_2d)
    call two2one(nlong,nlat,wht_2d,nlatlong,wht_1d)
    M_mdl2ob(i,1:nlatlong)=wht_1d(:)
end do
    
    do id=2,nlevels
         i=1+(id-1)*nlatlong
         j=id*nlatlong
         ii=1+(id-1)*nstn
         jj=id*nstn
         M_mdl2ob(ii:jj,i:j)=M_mdl2ob(1:nstn,1:nlatlong)
    end do

return
end subroutine
    

!---from hlat,hlon find the nearest point for clon,clat
!---output:
!----nxl,nxr,nyb,nyt
subroutine search_nearest_point(clon,clat,nx,hlon,ny,hlat,nxl,nxr,nyb,nyt)
implicit none
!!Only for grid box with the same grid length
!!Coded by Wang Jincheng, 2010.07.12
real :: clat
real :: clon

integer :: nx
real :: hlon(nx)
integer :: ny
real :: hlat(ny)

integer :: nxl,nxr,nyb,nyt

integer :: ix,iy
real :: hlon_int,lonl(nx)
real :: hlat_int(ny),latl(ny)

!equal distance for lon
hlon_int=hlon(2)-hlon(1)   !get the lon interval
do ix=1,nx
    lonl(ix)=hlon(ix)-clon   !get the distance from each grid to the stn
end do

ix=1
do while (abs(lonl(ix)).GT.hlon_int)
   ix=ix+1
end do
if (ix.EQ.nx) then
  nxl=nx
  nxr=1
else
  nxl=ix
  nxr=ix+1
end if

!not equal distance for lat
do iy=1,ny-1
    hlat_int(iy)=abs(hlat(iy)-hlat(iy+1))
    latl(iy)=abs(hlat(iy)-clat)
    if (latl(iy)<hlat_int(iy)) then
        nyb=iy
        nyt=iy+1
        exit
    end if
end do

end subroutine


subroutine Bilinear_interp_weight(x1,x2,y1,y2,xi,yi,w11,w12,w21,w22)
!!Bilinear interpolation 
!! (x1,y2) ---------------(x2,y2)
!!         |       |     |
!!         |       |     |
!!         |-------|-----|
!!         |    (xi,yi)  |
!!         |       |     |
!!         |       |     |
!! (x1,y1) ---------------(x2,y1)
!!Coded by Wang Jincheng, 2010.07.12
!!----------------------------------------------
implicit none
real ::x1
real ::x2
real ::y1
real ::y2

real ::xi
real ::yi

real ::wx1,wx2
real ::wy1,wy2
real ::w11,w12,w21,w22

wx1=(x2-xi)/(x2-x1)
wx2=(xi-x1)/(x2-x1)
wy1=(y2-yi)/(y2-y1)
wy2=(yi-y1)/(y2-y1)

w11=wy1*wx1
w12=wy2*wx1
w21=wy1*wx2
w22=wy2*wx2

return

end subroutine


!Fifth-order function
!input:
!---nx,ny: number of grid on x and y direction,m=nx*ny
!---d:grid distance
!---c:schur radius
!---BACK:background covariance matrix
subroutine schur(BACK)
include 'sphectra.h'
include 'paramod.h'
    integer::i,j,ix,iy,il1,il2
    real::ra
    real::BACK(nlatlonglv,nlatlonglv)
    real::c_long,c_lat,r_long,r_lat
    real::wht
!----------------------------------------
    do i=1,nlatlong
        call getlonglat(i,c_long,c_lat)  !get the long,lat of the center grid

        do j=i,nlatlong
            call getlonglat(j,r_long,r_lat)   !get the long,lat of the other grid

            call greatcir_distance(RAYTER,c_long,c_lat,r_long,r_lat,ra)

            call fifth_order(c_radius,ra,wht)   !return the localization coefficient

            !if (wht>0.) then
            !    do il1=1,nlevels
            !        ix=i+(il1-1)*nlatlong
            !        iy=j+(il1-1)*nlatlong
            !        BACK(ix,iy)=BACK(ix,iy)*wht
            !    end do
            !else if (wht==0.) then
                do il1=1,nlevels
                    do il2=il1,nlevels
                    ix=i+(il1-1)*nlatlong
                    iy=j+(il2-1)*nlatlong
                    BACK(ix,iy)=BACK(ix,iy)*wht
                    end do
                end do
            !end if

        end do
    end do

    do i=1,nlatlonglv
        do j=i,nlatlonglv
            BACK(j,i)=BACK(i,j)
        end do
    end do
    return
    end subroutine
!=---=---------right-------------
                        
           

!fifth order function
subroutine fifth_order(c_radius,ra,wht)
implicit none
    real::ra,wht
    real::c,c_radius
    c=c_radius/2.
			
    if (ra>=0.and.ra<=c) then
	    wht=(-0.25*(ra/c)**5+&
			0.5*(ra/c)**4+0.625*(ra/c)**3-1.667*(ra/c)**2+1.0)
    else if (ra>=c.and.ra<=c_radius) then
		wht=(0.0833*(ra/c)**5-0.5*(ra/c)**4+&
			0.625*(ra/c)**3+1.667*(ra/c)**2-5.0*(ra/c)+4.0&
			-0.6667*(c/ra))
    else
	    wht=0.0
    end if

    return
    end subroutine
!=====================================right

!-----return the long,lat,level of each grid i
subroutine getlonglat(i,c_long,c_lat)
include 'sphectra.h'

integer::i,ix,iy
real::a,b
real::c_long,c_lat

a=real(i)/real(nlong)
iy=floor(a)+1        !get the lat number

ix=mod(i,nlong)      !get the long number

if (ix==0) then
    ix=nlong
    iy=iy-1
end if

c_long=0.+(ix-1)*long_inter
c_lat=lat_angle(iy)
c_long=PI*c_long/180.0
c_lat=PI*c_lat/180.0
return
end subroutine
!----------right-------------

!calculate the great circle distance between two grids on earth
!each coordinate is c_long,c_lat
!the other is r_long,r_lat
subroutine greatcir_distance(RAYTER,c_long,c_lat,r_long,r_lat,r)
implicit none
real::c_long,c_lat,r_long,r_lat
real::r,RAYTER
real::a,b,c
!------------
a=sin(c_lat)*sin(r_lat)
b=abs(c_long-r_long)
b=cos(b)
b=cos(c_lat)*cos(r_lat)*b
c=a+b
c=acos(c)
r=RAYTER*c
return
end subroutine
!---------right--------------
