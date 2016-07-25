!---------------------------------------------
!!!!generate np bred vector
!--input:
!x0_start: the initial reference state
!p_start: initial random perturbations
!stan: rescaling factor
!ncl: the number of breeding cycles
!tcl: rescaling interval
!--output:
!pstart: output bred vectors
subroutine breeding(x0_start,nper,pstart,stan,ncl,tcl,xout,ptemp)
!subroutine breeding(x0_start,nper,pstart,stan,ncl,tcl,ptemp)
include 'sphectra.h'
include 'paramod.h'

!implicit none
integer::i,j,k
!integer::nx,ny
integer::ncl,tcl,nper
real::x0_start(nlong,nlat,nlevels)
real::pstart(nlong,nlat,nlevels,nper)
real::ptemp(nlong,nlat,nlevels,nper)
real::stan(nlevels)
real::x0(nlong,nlat,nlevels),x0_p(nlong,nlat,nlevels)
real::x1(nlong,nlat,nlevels),x1_p(nlong,nlat,nlevels)
real::xtemp(nlong,nlat,nlevels,0:1)
real::p0(nlong,nlat,nlevels),pp(nlong,nlat,nlevels),pout(nlong,nlat,nlevels)
real::xout(nlong,nlat,nlevels)
real::dl
x0=x0_start
ptemp=pstart
do i=1,ncl
        xtemp=0.
        xtemp(:,:,:,0)=x0(:,:,:)
        !reference state
        call integrate(tcl,tcl,1,xtemp)
        x1(:,:,:)=xtemp(:,:,:,1)  !x1 is the final state at the end of a breeding cycle
        !================================
        do j=1,nper
        p0(:,:,:)=ptemp(:,:,:,j)
        x0_p=x0+p0              !superpose a error vector on the initial state
        xtemp=0.
        !----perturbed states   
        xtemp(:,:,:,0)=x0_p(:,:,:)
        call integrate(tcl,tcl,1,xtemp)
        x1_p(:,:,:)=xtemp(:,:,:,1)
        !================================
        pp=x1_p-x1             !calculate the perturabtions
        do k=1,nlevels
            call field_rescale_2d(nlong,nlat,1,&
                pp(:,:,k),stan(k),pout(:,:,k))  !rescale the perturbations
        end do
        ptemp(:,:,:,j)=pout(:,:,:)
        end do
        x0=x1


!write(132,rec=i)ptemp(:,:,2,:)
!write(41,rec=i)x1
end do
xout=x1    !output the final reference state xout
return
end subroutine
!-----------------------right-----------


!---------------------------------------------
!!!!generate np bred vector and calculate the growth rate
!in each cycle
!--input:
!x0_start: the initial reference state
!p_start: initial random perturbations
!stan: rescaling factor
!ncl: the number of breeding cycles
!tcl: rescaling interval
!--output:
!pstart: output bred vectors
subroutine breeding_errorgrowth(x0_start,nper,pstart,ncl,tcl,egrow)
include 'sphectra.h'
include 'paramod.h'
        
!implicit none
integer::i,j,k
!integer::nx,ny
integer::ncl,tcl,nper
real::x0_start(nlong,nlat,nlevels)
real::pstart(nlong,nlat,nlevels,nper)
real::ptemp(nlong,nlat,nlevels,nper)
real::stan(nlevels)
real::x0(nlong,nlat,nlevels),x0_p(nlong,nlat,nlevels)
real::x1(nlong,nlat,nlevels),x1_p(nlong,nlat,nlevels)
real::xtemp(nlong,nlat,nlevels,0:1)
real::p0(nlong,nlat,nlevels),pp(nlong,nlat,nlevels),pout(nlong,nlat,nlevels)
real::xout(nlong,nlat,nlevels)
real::dl,inisize(nper),endsize(nper)
real::egrow(nlong,nlat,nlevels,nper)
x0=x0_start
ptemp=pstart
!----calculate the size of initial error
!do i=1,nper
!	call weight_rmse(.false.,nlong,nlat,nlevels,&
!	pstart(:,:,:,i),pstart(:,:,:,i),inisize(i))
!end do
!------------
do i=1,ncl
        xtemp=0.
        endsize=0.
        xtemp(:,:,:,0)=x0(:,:,:)
        call integrate(tcl,tcl,1,xtemp)
        x1(:,:,:)=xtemp(:,:,:,1)  !x1 is the final state at the end of a breeding cycle
        !================================
        do j=1,nper
        p0(:,:,:)=ptemp(:,:,:,j)
        x0_p=x0+p0              !superpose a error vector on the initial state
        xtemp=0.
        !----        
        xtemp(:,:,:,0)=x0_p(:,:,:)
        call integrate(tcl,tcl,1,xtemp)
        x1_p(:,:,:)=xtemp(:,:,:,1)
        !================================
        pp=x1_p-x1             !calculate the perturabtions
        !call weight_rmse(.false.,nlong,nlat,nlevels,pp,pp,endsize(j))
        egrow(:,:,:,j)=abs(pp)-abs(p0)
        !--------------------------------
                !do k=1,nlevels
                !call field_rescale_2d(nlong,nlat,1,&
                !pp(:,:,k),stan(k),pout(:,:,k))  !rescale the perturbations
                !end do
        !ptemp(:,:,:,j)=pout(:,:,:)
        end do
        x0=x1

!write(31,rec=i)ptemp
!write(41,rec=i)x1
end do
write(*,*) "X1: ",x1(1,:,2)
pause
xout=x1    !output the final reference state xout
return
end subroutine
!----------------------------------------------

!This subroutine calculates nper NLLVs.
!--input:
!----y0_start:the initial reference state
!----nn:number of field grids.
!----ns:number of error samples 
!----id0:position of initial condition in the data file
!----p0:initial perturbations
!----array:initial sequence from 1 to nn
!----stan:rescaled perturbation size
!----tcl:rescaling interval
!----ncl:number of breeding cycle
!--output:
!----pout:orthogonal vectors after breeding cycles
subroutine reorth(y0_start,nper,pstart,stan,ncl,tcl,pout)
include 'sphectra.h'
include 'paramod.h'
!implicit none
integer::i,j,k,id
integer::ncl,tcl
real::stan(nlevels)
real::per0(nlong,nlat,nlevels),per1(nlong,nlat,nlevels)
real::per1_all(nlong,nlat,nlevels,nper)
real::y0_start(nlong,nlat,nlevels)
real::y0(nlong,nlat,nlevels),y1(nlong,nlat,nlevels)
real::y0_p(nlong,nlat,nlevels),y1_p(nlong,nlat,nlevels)
real::ytemp0(nlong,nlat,nlevels,0:1)
real::pstart(nlong,nlat,nlevels,nper)
real::ptemp(nlong,nlat,nlevels,nper),pout(nlong,nlat,nlevels,nper)
!------------------
ptemp=pstart    !do not change pstart
y0=y0_start     !do not change y0_start(for breeding)
!---------------------
!y0--->y1
!+per0
!y0_p--->y1_p
outer:  do k=1,ncl
                ytemp0=0.
                ytemp0(:,:,:,0)=y0(:,:,:)
                !reference state
                call integrate(tcl,tcl,1,ytemp0)
                y1(:,:,:)=ytemp0(:,:,:,1)    !the integration of reference
                !-------------------------
inner:          do i=1,nper
                        ytemp0=0.
                        per0=0.
                        dist=0.
                        per0(:,:,:)=ptemp(:,:,:,i)
                        y0_p=y0+per0
                        ytemp0(:,:,:,0)=y0_p(:,:,:)
                        call integrate(tcl,tcl,1,ytemp0)
                        y1_p(:,:,:)=ytemp0(:,:,:,1)   !the integration of forecast
                        !!!-------------------------------------
                        per1=y1_p-y1
                        !if(k<=ncl)then  !mod(k,ncl)/=0)then
                        !     do j=1,nlevels
                        !     call field_rescale_2d(nlong,nlat,per1(:,:,j),&
                        !          stan(j),ptemp(:,:,j,i))
                        !     end do
                        !else
                        
                        per1_all(:,:,:,i)=per1(:,:,:)  !if orthogonalizationis necessary
                                                       !deposite all the perturbations
                        !!!---------------------------------------
                        !endif
                end do inner

                !if(k>ncl)then!mod(k,ncl)==0)then
                call gsr_gb(.false.,nlong,nlat,nlevels,nper,nper,&
                    per1_all,stan,ptemp)
                !do j=1,ns
                !        call onetotwo(nn,errorstart(:,j),nx,ny,pstart(:,:,j))
                !end do
                !endif
        y0=y1

        !write(131,rec=(id-1)*ncl+k)ptemp(:,:,2,:)
        end do outer

        pout=ptemp
return
end subroutine
!!!---------------------------------------------------------------

!---------------------------------------------
!!!!generate np bred vector regionally rescaling
!--input:
!x0_start: the initial reference state
!p_start: initial random perturbations
!stan: rescaling factor
!ncl: the number of breeding cycles
!tcl: rescaling interval
!--output:
!pstart: output bred vectors
subroutine breeding_mask(x0_start,nper,pstart,stan,ncl,tcl,ptemp)
include 'sphectra.h'
include 'paramod.h'

!implicit none
integer::i,j,k
!integer::nx,ny
integer::ncl,tcl,nper
real::x0_start(nlong,nlat,nlevels)
real::pstart(nlong,nlat,nlevels,nper)
real::ptemp(nlong,nlat,nlevels,nper)
real::stan(nlong,nlat,nlevels)
real::x0(nlong,nlat,nlevels),x0_p(nlong,nlat,nlevels)
real::x1(nlong,nlat,nlevels),x1_p(nlong,nlat,nlevels)
real::xtemp(nlong,nlat,nlevels,0:1)
real::p0(nlong,nlat,nlevels),pp(nlong,nlat,nlevels),pout(nlong,nlat,nlevels)
real::xout(nlong,nlat,nlevels)
real::dl
x0=x0_start
ptemp=pstart
do i=1,ncl
        xtemp=0.
        xtemp(:,:,:,0)=x0(:,:,:)
        call integrate(tcl,tcl,1,xtemp)
        x1(:,:,:)=xtemp(:,:,:,1)  !x1 is the final state at the end of a breeding cycle
        !================================
        do j=1,nper
        p0(:,:,:)=ptemp(:,:,:,j)
        x0_p=x0+p0              !superpose a error vector on the initial state
        xtemp=0.
        !----        
        xtemp(:,:,:,0)=x0_p(:,:,:)
        call integrate(tcl,tcl,1,xtemp)
        x1_p(:,:,:)=xtemp(:,:,:,1)
        !================================
        pp=x1_p-x1             !calculate the perturabtions
        
        call region_rescale(smth_r,nlong,nlat,nlevels,pp,stan,pout)  
                              !regional rescaling using the mask
        ptemp(:,:,:,j)=pout(:,:,:)
        end do
        x0=x1

!write(31,rec=i)ptemp
!write(41,rec=i)x1
end do
xout=x1    !output the final reference state xout
return
end subroutine
!----------------------------------------------

!This subroutine calculates nper NLLVs regionally rescaled.
!--input:
!----y0_start:the initial reference state
!----nn:number of field grids.
!----ns:number of error samples 
!----id0:position of initial condition in the data file
!----p0:initial perturbations
!----array:initial sequence from 1 to nn
!----stan:rescaled perturbation size
!----tcl:rescaling interval
!----ncl:number of breeding cycle
!--output:
!----pout:orthogonal vectors after breeding cycles
subroutine reorth_mask(y0_start,nper,pstart,stan,ncl,tcl,pout)
include 'sphectra.h'
include 'paramod.h'
!implicit none  
integer::i,j,k  
integer::ncl,tcl
real::stan(nlong,nlat,nlevels)
real::per0(nlong,nlat,nlevels),per1(nlong,nlat,nlevels)
real::per1_all(nlong,nlat,nlevels,nper)
real::y0_start(nlong,nlat,nlevels)
real::y0(nlong,nlat,nlevels),y1(nlong,nlat,nlevels)
real::y0_p(nlong,nlat,nlevels),y1_p(nlong,nlat,nlevels)
real::ytemp0(nlong,nlat,nlevels,0:1)
real::pstart(nlong,nlat,nlevels,nper)
real::ptemp(nlong,nlat,nlevels,nper),pout(nlong,nlat,nlevels,nper)
!------------------
ptemp=pstart    !do not change pstart
y0=y0_start     !do not change y0_start(for breeding)
!---------------------
!y0--->y1
!+per0
!y0_p--->y1_p
outer: do k=1,ncl
                ytemp0=0.
                ytemp0(:,:,:,0)=y0(:,:,:)
                call integrate(tcl,tcl,1,ytemp0)
                y1(:,:,:)=ytemp0(:,:,:,1)    !the integration of reference
                !-------------------------
inner:          do i=1,nper
                        ytemp0=0.
                        per0=0.
                        dist=0.
                        per0(:,:,:)=ptemp(:,:,:,i)
                        y0_p=y0+per0
                        ytemp0(:,:,:,0)=y0_p(:,:,:)
                        call integrate(tcl,tcl,1,ytemp0)
                        y1_p(:,:,:)=ytemp0(:,:,:,1)   !the integration of forecast
                        !!!-------------------------------------
                        per1=y1_p-y1
                        !if(mod(k,ncl)/=0)then
                        !     do j=1,nlevels
                        !     call field_rescale_2d(nlong,nlat,per1(:,:,j),&
                        !          stan(j),ptemp(:,:,j,i))
                        !     end do
                        !else

                        per1_all(:,:,:,i)=per1(:,:,:)  !if orthogonalizationis necessary
                                                       !deposite all the perturbations
                        !!!---------------------------------------
                        !endif
                end do inner

                !if(mod(k,ncl)==0)then
                call gsr_region(.false.,nlong,nlat,nlevels,nper,nper,&
                         per1_all,ptemp)
                do i=1,nper
                call region_rescale(smth_r,nlong,nlat,nlevels,ptemp(:,:,:,i),&
                        stan,pout(:,:,:,i))
                end do
                ptemp = pout
                !do j=1,ns
                !        call onetotwo(nn,errorstart(:,j),nx,ny,pstart(:,:,j))
                !end do
                !endif
        y0=y1
        end do outer

        !call gsr_region(.false.,nlong,nlat,nlevels,nper,nper,&
        !        ptemp,pout)
        !pout=ptemp
return
end subroutine
!!!---------------------------------------------------------------

!---------------------------------------------
!!!!generate np bred vector
!--input:
!x0_start: the initial reference state
!p_start: initial random perturbations
!stan: rescaling factor
!ncl: the number of breeding cycles
!tcl: rescaling interval
!--output:
!pstart: output bred vectors
subroutine breeding_ana(x0_ref,nper,pstart,stan,ncl,tcl,ptemp)
!subroutine breeding(x0_start,nper,pstart,stan,ncl,tcl,ptemp)
include 'sphectra.h'
include 'paramod.h'

!implicit none
integer::i,j,k
!integer::nx,ny
integer::ncl,tcl,nper
real::x0_ref(nlong,nlat,nlevels,0:ncl)
real::pstart(nlong,nlat,nlevels,nper)
real::ptemp(nlong,nlat,nlevels,nper)
real::stan(nlevels)
real::x0(nlong,nlat,nlevels),x0_p(nlong,nlat,nlevels)
real::x1(nlong,nlat,nlevels),x1_p(nlong,nlat,nlevels)
real::xtemp(nlong,nlat,nlevels,0:1)
real::p0(nlong,nlat,nlevels),pp(nlong,nlat,nlevels),pout(nlong,nlat,nlevels)
real::xout(nlong,nlat,nlevels)
real::dl

ptemp=pstart
do i=1,ncl
        !reference state
        x0(:,:,:)=x0_ref(:,:,:,i-1)
        x1(:,:,:)=x0_ref(:,:,:,i)  !x1 is the final state at the end of a breeding cycle
        !================================
        do j=1,nper
        p0(:,:,:)=ptemp(:,:,:,j)
        x0_p=x0+p0              !superpose a error vector on the initial state
        xtemp=0.
        !----perturbed states   
        xtemp(:,:,:,0)=x0_p(:,:,:)
        call integrate(tcl,tcl,1,xtemp)
        x1_p(:,:,:)=xtemp(:,:,:,1)
        !================================
        pp=x1_p-x1             !calculate the perturabtions
        do k=1,nlevels
            call field_rescale_2d(nlong,nlat,1,&
                pp(:,:,k),stan(k),pout(:,:,k))  !rescale the perturbations
        end do
        ptemp(:,:,:,j)=pout(:,:,:)
        end do
end do
!xout=x1    !output the final reference state xout
return
end subroutine
!-----------------------right-----------




!This subroutine calculates nper NLLVs.
!--input:
!----y0_start:the initial reference state
!----nn:number of field grids.
!----ns:number of error samples 
!----id0:position of initial condition in the data file
!----p0:initial perturbations
!----array:initial sequence from 1 to nn
!----stan:rescaled perturbation size
!----tcl:rescaling interval
!----ncl:number of breeding cycle
!--output:
!----pout:orthogonal vectors after breeding cycles
subroutine reorth_ana(y0_ref,nper,pstart,stan,ncl,tcl,pout)
include 'sphectra.h'
include 'paramod.h'
!implicit none
integer::i,j,k,id
integer::ncl,tcl
real::stan(nlevels)
real::per0(nlong,nlat,nlevels),per1(nlong,nlat,nlevels)
real::per1_all(nlong,nlat,nlevels,nper)
real::y0_ref(nlong,nlat,nlevels,0:ncl)
real::y0(nlong,nlat,nlevels),y1(nlong,nlat,nlevels)
real::y0_p(nlong,nlat,nlevels),y1_p(nlong,nlat,nlevels)
real::ytemp0(nlong,nlat,nlevels,0:1)
real::pstart(nlong,nlat,nlevels,nper)
real::ptemp(nlong,nlat,nlevels,nper),pout(nlong,nlat,nlevels,nper)
!------------------
ptemp=pstart    !do not change pstart
!---------------------
!y0--->y1
!+per0
!y0_p--->y1_p
outer:  do k=1,ncl
                !reference state
                y0(:,:,:)=y0_ref(:,:,:,k-1)
                y1(:,:,:)=y0_ref(:,:,:,k)
                !-------------------------
inner:          do i=1,nper
                        ytemp0=0.
                        per0=0.
                        dist=0.
                        per0(:,:,:)=ptemp(:,:,:,i)
                        y0_p=y0+per0
                        ytemp0(:,:,:,0)=y0_p(:,:,:)
                        call integrate(tcl,tcl,1,ytemp0)
                        y1_p(:,:,:)=ytemp0(:,:,:,1)   !the integration of forecast
                        !!!-------------------------------------
                        per1=y1_p-y1
                        !if(k<=ncl)then  !mod(k,ncl)/=0)then
                        !     do j=1,nlevels
                        !     call field_rescale_2d(nlong,nlat,per1(:,:,j),&
                        !          stan(j),ptemp(:,:,j,i))
                        !     end do
                        !else
                        
                        per1_all(:,:,:,i)=per1(:,:,:)  !if orthogonalizationis necessary
                                                       !deposite all the perturbations
                        !!!---------------------------------------
                        !endif
                end do inner

                !if(k>ncl)then!mod(k,ncl)==0)then
                call gsr_gb(.false.,nlong,nlat,nlevels,nper,nper,&
                    per1_all,stan,ptemp)
                !do j=1,ns
                !        call onetotwo(nn,errorstart(:,j),nx,ny,pstart(:,:,j))
                !end do
                !endif

        !write(131,rec=(id-1)*ncl+k)ptemp(:,:,2,:)
        end do outer

        pout=ptemp
return
end subroutine
!!!---------------------------------------------------------------
