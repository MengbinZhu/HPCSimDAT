!c
!c==================================================================
!c MES PROPRES PARAMETRISATIONS (PAM, voir aussi paramode.h)
!c==================================================================
!c  NBRE DE NIVEAUX
!c
	parameter(NLEVELS=3)
!c++++++++++++++++++++++++++++++++++++++++++++
!c
!c DIMENSION DES CHAMPS
      parameter(NVARIA=NCPLX*NBDEGLIB*NLEVELS        &
               ,NVARIA2=NCPLX*NBDEGLIB)
!c
!c 
!c++++++++++++++++++++++++++++++++++++++++++++
!c  PAS DE TEMPS (secondes)
!c
      parameter(DT=3600)
!!c++++++++++++++++++++++++++++++++++++++++++++
!c RAYONS DE DEFORMATION (m)
!c (taken from Marshall & Molteni : Towards a dynamical
!c   understanding of planetary-scale flow regimes, JAS 93)
!c
!c     couche 200-500 hpa : ray2 
!c            500-800 hpa : ray4
!c
      parameter(ray2=700.e3,ray4=450.e3)
      parameter(ray2n2=1./(ray2*ray2),ray4n2=1./(ray4*ray4))
!c++++++++++++++++++++++++++++++++++++++++++++
!C HAUTEUR D'ECHELLE DE L'ATMOSPHERE
!C H0=9.e3 m
!c
      parameter(h0=9.e3,h0n1=1./h0)
!c
!c+++++++++++++++++++++++++++++++++++++++++++
!c INTERVALLE ENTRE 2 NIVEAUX DE PRESSION SUCCESSIFS (Pascal)
!c
      parameter(DP=300.e2)
!c++++++++++++++++++++++++++++++++++++++++++++
!c PRESSION A LA LIMITE SUP.
!c 
      parameter(PRESTOP=50.e2)
!c++++++++++++++++++++++++++++++++++++++++++++
!c DENSITE DE L'ATMOSPHERE A 950mb (niveau 6)
!c d'apres le US STANDART ATMOSPHERE, 1962
!c (kg.m-3)
!c 
      parameter(RHO6=1.16445)
!c+++++++++++++++++++++++++++++++++++++++++++++
!c VALEUR DU PARMETRE DE CORIOLIS A LA LATITUDE
!c 45 degres nord (s-1)
!c
      parameter(sin45=0.70710678118655)
      parameter(f0=omega2*sin45)
!c VALEUR DU COEF. SPECTRAL DE L'HARMONIQUE
!c n=1 m=0 DU PARAMETRE DE CORIOLIS f(phi)=2*omega*sin(phi)
!c         f1=2*omega/sqrt(3.)
!c
      parameter(sqrt3=1.732050808)
      parameter(f1=omega2/sqrt3)
!c++++++++++++++++++++++++++++++++++++++++++++++++++
!c PARAMETRES UTILISES POUR LE CALCUL DE LA DISSIPATION
!c (taken from Marshall & Molteni : Towards a dynamical
!c   understanding of planetary-scale flow regimes, JAS 92 ou 93)
!c
!c Radiative time scale for temperature relaxation
!c              25 days
!c
      parameter(TORDAYS=25.)
      parameter(TOR=86400.*TORDAYS) !en secondes for sure
      parameter(RAY2N2ONTOR=RAY2N2/TOR                       &
              ,RAY4N2ONTOR=RAY4N2/TOR)
!c
!c Ekman dissipation time scale
!c               3 days
!c
      parameter(TOEDAYS=3.)
      parameter(TOE=TOEDAYS*86400.)
      parameter(ONTOE=1./TOE)
!c
!c Parametre de la dissipation selective
!c
!c
      parameter(TOHDAYS=2.)
      parameter(TOH=TOHDAYS*86400.)
!c+++++++++++++++++++++++++++++++++++++++++++++++++++
!c Rayon de la Terre au carre
!c
      parameter(rayter2=rayter*rayter)
!c++++++++++++++++++++++++++++++++++++++++++++++++++++
!c  Coefficients de la transformation lineaire Q --> PSI
!c
      common /blin/ cq2psi(nlevels,nlevels,nvaria2)
!c++++++++++++++++++++++++++++++++++++++++++++++++++++
!c  Commons de topographie et d'autres
!c
      common /btopo/ stopo(nvaria2),gtopo(n2long,nlat),ftopo(nvaria2)
      common /bekm3/ glmask(n2long,nlat),gek(n2long,nlat),sek(nvaria2)
      common /bforc/ force(nvaria2,nlevels)
      common /bcori/ corio(n2long,nlat)
!c++++++++++++++++++++++++++++++++++++++++++++++++++++
!c  Operateur de dissipation selective
!c
      common /bdise/ oprdissip(nvaria2)
!c++++++++++++++++++++++++++++++++++++++++++++++++++++
!c  Options
!c
      common /biopt/ iopforce
!+++++++++++++++++++++++++++++++++++++++++++++++++++
!   variable for EnKF assimilation
    integer,parameter :: nlatlonglv=nlatlong*nlevels
    integer,parameter :: allnstn=nstn*nlevels
    real :: M_mdl2ob(allnstn,nlatlonglv)
    real :: ini_anafield(nlong,nlat,nlevels)

!+++++++++++++++++++++++++++++++++++++++++++++++++++
!   ensemble prediction variable
    real :: p0_random(nlong,nlat,nlevels,npair*2)  !initial random error
    real :: start_gpsi(nlong,nlat,nlevels)       
    real :: start_fc(nlong,nlat,nlevels)        !initial state for forecasting
    real :: ana_error_gpsi(nlong,nlat,nlevels)
    real :: gpsi_ctl(nlong,nlat,nlevels,0:nfc_day)
    real :: gpsi_true(nlong,nlat,nlevels,0:nfc_day)
    real :: gpsi_allmember(nlong*nlat,nlevels,npair*2,0:nfc_day)

!  different  ensemble prediction methods
    real :: gpsi_bv(nlong,nlat,nlevels,2*npair)
    real :: gpsi_nllv(nlong,nlat,nlevels,2*npair)
    real :: gpsi_EnKF(nlong,nlat,nlevels,2*npair)
    real :: gpsi_mc(nlong,nlat,nlevels,2*npair)
    real :: gpsi_emean_bv(nlong,nlat,nlevels,0:nfc_day)
    real :: gpsi_emean_nllv(nlong,nlat,nlevels,0:nfc_day)
    real :: gpsi_emean_EnKF(nlong,nlat,nlevels,0:nfc_day)
    real :: gpsi_emean_mc(nlong,nlat,nlevels,0:nfc_day)

