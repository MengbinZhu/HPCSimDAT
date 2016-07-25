
!!one sample
!---input:
!----dimen:number of variables
!----yn(dimen):initial state
!----error0:error of observation
!----nstep:steps of prediction
!----npair:pairs of ensemble prediction
!---output:
!----xob:time series of observation
!----xtrue:true trajectory
!----xen:time series of ensemble prediction
subroutine individual(p0,ini_anafield,start_fc,xctl,xtrue,xen,xall_member)
include 'sphectra.h'
include 'paramod.h'
integer::i,j,k
real::p0(nlong,nlat,nlevels,2*npair)      !initial perturbatons
real::xctl(nlong,nlat,nlevels,0:nfc_day)   !control run
real::xtrue(nlong,nlat,nlevels,0:nfc_day) !reference state
real::xen(nlong,nlat,nlevels,0:nfc_day)   !ensemble prediction
real::xsingle(nlong,nlat,nlevels,npair*2)
real::xtemp(nlong,nlat,nlevels,0:nfc_day)
real::xall_member(nlong*nlat,nlevels,npair*2,0:nfc_day)
!--------
xtrue=0.
xctl=0.
xen=0.
xall_member=0.
!!!!form initial true and analysis state
xtrue(:,:,:,0)=start_fc(:,:,:)
xctl(:,:,:,0)=ini_anafield(:,:,:)
!--------------------------------------------------------------------------
call integrate(nfc_hour,nout_inter,nfc_day,xctl)   !form the control forecast
!--------------------------------------------------------------------------
call integrate(nfc_hour,nout_inter,nfc_day,xtrue)  !form the true run
!--------------------------------------------------------------------------
do i=1,2*npair
        xsingle(:,:,:,i)=xctl(:,:,:,0)+p0(:,:,:,i)    !!!form ensemble members-------
end do
!---------integrate all samples
        do j=1,2*npair
                xtemp=0.
                xtemp(:,:,:,0)=xsingle(:,:,:,j)
                call integrate(nfc_hour,nout_inter,nfc_day,xtemp)
                xen=xen+xtemp/real(2*npair)
                do k=0,nfc_day
                    do nl=1,nlevels
                        call two2one(nlong,nlat,xtemp(:,:,nl,k),&
                            nlong*nlat,xall_member(:,nl,j,k))
                    end do
                end do                         
        end do
return
end subroutine
!----------------------------------right-------------
