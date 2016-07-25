     program runtest
!
! Integration test of the SPHERICAL model with dissipation and forcing with mask
!
      include 'sphectra.h'
      include 'paramod.h'
!
      real spsi(nvaria2,nlevels),am,fameq,fammo,famdis
      real sq(nvaria2,nlevels),sq1(nvaria2,nlevels)
      real gpsi(n2long,nlat),gq(n2long,nlat)
      real gpsi_std(nlong,nlat,nlevels)
      real xctl(nlong,nlat,nlevels,0:nfc_day)
      real std(nlevels)
!
      integer::ilong,ilat,nl,it,ii,jj,iflag,nfi
      character*80 fntrue, fnana
      character*80 fnpsi1,fnpsi2,fnpsi3,fstd
      !character*80 fntrue1,fntrue2,fntrue3,fnave1,fnave2,fnave3,fnctl1,fnctl2,fnctl3
      !character*80 fnall1,fnall2,fnall3
!------initial value
      data std /12.,8.13,5.97/
      !std=2.0
      std=std*g/f0
!--------------------
 
      !read *,fnci1
      !read *,fnci2
      !read *,fnci3
      !read *,fnpsi1
      !read *,fnpsi2
      !read *,fnpsi3
!input file
      fnpsi1='./data/outpsi_200.dat'
      fnpsi2='./data/outpsi_500.dat'
      fnpsi3='./data/outpsi_800.dat'
      !fstd='std_ght_mo.dat'

!output file
      fntrue='./data/true_westhem.dat'
      fnana='./data/ana_westhem.dat'
      !fntrue1='outpsi_true_200.dat'
      !fntrue2='outpsi_true_500.dat'
      !fntrue3='outpsi_true_800.dat'

      !fnave1='outpsi_ave_200.dat'
      !fnave2='outpsi_ave_500.dat'
      !fnave3='outpsi_ave_800.dat'

      !fnall1='outpsi_allmember_200.dat'
      !fnall2='outpsi_allmember_500.dat'
      !fnall3='outpsi_allmember_800.dat'

      !fnctl1='outpsi_ctl_200.dat'
      !fnctl2='outpsi_ctl_500.dat'
      !fnctl3='outpsi_ctl_800.dat'

      !read *,iddep
      !read *,nsteps
      !read *,nstok,nsteq
!
!  Ouverture des fichiers
!  input data
    open(21,file=fnpsi1,form='unformatted',access='direct',recl=nlong*nlat)
    open(22,file=fnpsi2,form='unformatted',access='direct',recl=nlong*nlat)
    open(23,file=fnpsi3,form='unformatted',access='direct',recl=nlong*nlat)
    !open(81,file=fstd,form='unformatted',access='direct',recl=nlong*nlat)
!  output data----------------------------------
    !open(111,file=fntrue,form='unformatted',access='direct',recl=nlong*nlat*nlevels)
    !open(121,file=fnana,form='unformatted',access='direct',recl=nlong*nlat*nlevels)

      !open(81,file='topo.dat',&
      !     form='unformatted',access='direct',recl=nlong*nlat)      
      !open(31,file=fntrue1,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)
      !open(32,file=fntrue2,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)
      !open(33,file=fntrue3,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)

      !open(41,file=fnave1,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)
      !open(42,file=fnave2,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)
      !open(43,file=fnave3,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)

      !open(51,file=fnall1,form='unformatted',access='direct',&
      !     recl=(1+nfc_day)*nlong*nlat*2*npair)
      !open(52,file=fnall2,form='unformatted',access='direct',&
      !     recl=(1+nfc_day)*nlong*nlat*2*npair)
      !open(53,file=fnall3,form='unformatted',access='direct',&
      !     recl=(1+nfc_day)*nlong*nlat*2*npair)

      !open(61,file=fnctl1,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)
      !open(62,file=fnctl2,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)
      !open(63,file=fnctl3,form='unformatted',access='direct',recl=(1+nfc_day)*nlong*nlat)
!------------------------------------
!  Initialisations
!
    call initqgmod      !initialize the model
    call mdl2ob(M_mdl2ob)
    call random_seed()
!  read the std field
    !do nl=1,nlevels
    !    read(81,rec=nl)gpsi_std(:,:,nl)
    !    call trans(gpsi_std(:,:,nl))
    !    gpsi_std(:,:,nl)=gpsi_std(:,:,nl)*0.2*g/f0  !STD of streamfunction
    !end do

  
!  
!  Reading in the initial conditions
!
    !iflag=0
        it=401
        !iflag=iflag+1
        print *,"TOTAL NUMBER OF ITERATIONS: ",it
!-------read initial state
        do nl=1,nlevels
            nfi = 20 + nl
            read(nfi,rec=it)start_gpsi(:,:,nl)   !the initial state
            call trans(start_gpsi(:,:,nl))
        end do
        
        !xctl(:,:,:,0)=start_gpsi(:,:,:)
        !call integrate(nfc_hour,nout_inter,nfc_day,xctl)

        call barotropic_ENKF(M_mdl2ob,start_gpsi)
!------output---------
!write(*,*)egrow
!do it = 1,nfc_day
!    do il = 1,nlevels
!        call output_trans(xctl(:,:,il,it))
!    end do
    !write(*,*)xctl(1,:,2,it)
    !pause
!    write(111,rec=it)xctl(:,:,:,it)
!end do
!c
      stop
    end program


