!C===============================================================================
!C PARAMETRISATION DE LA TRONCATURE.
!C==============================================================================
!C     PARAMETER ( NTRUNC = 5, NHLAT = 4)
!C     PARAMETER ( NTRUNC = 13, NHLAT = 10)
      PARAMETER ( NTRUNC = 21, NHLAT = 16)
!C     PARAMETER ( NTRUNC = 31, NHLAT = 24)
!C     PARAMETER ( NTRUNC = 42, NHLAT = 32)
!C     PARAMETER ( NTRUNC = 63, NHLAT = 48)
!C     PARAMETER ( NTRUNC = 96, NHLAT = 72)
!C     PARAMETER ( NTRUNC = 170, NHLAT = 128)
!C     PARAMETER ( NTRUNC = 42, NHLAT = 32)
!C===============================================================================
!C PARAMETRISATION DE LA TRONCATURE (AUTOMATISE).
!C==============================================================================
      PARAMETER (NHLONG=NHLAT*2)
      PARAMETER( &
                       NLAT= 2 * NHLAT, NLONG   = 2 * NHLONG,        &
                       NLATLONG= NLAT * NLONG,                       &
                       NH1LONG = NHLONG + 1   ,			     &	
                       N1 = NTRUNC/2, N2 = (NTRUNC+1)/2,             &
                       NP = (N1+1) * (N2+1), NI = (N1+1) * N2,       &
                       NDIM = NP + NI,                               &
       NCPLX   =    2, NBDEGLIB = (NTRUNC + 1) * (NTRUNC + 2) / 2  , &
                       N2LONG  = NLONG + 2, NSIZE = N2LONG * NLAT,   &
                       NC1S = NCPLX)
!C=========================================================================
!C PARAMETRISATION PHYSIQUE ET MATHEMATIQUES (mksa).
!C (NB: Taken from Gill ( Atmosphere-Ocean Dynamics) )
!C=========================================================================
!C MATHEMATICAL CONSTANTS
      PARAMETER ( PI = 3.141592653589793238 , ANGLE = 2.0 * PI / NLONG    &
     ,PI2=2.0*PI) 
!C PHYSICAL CONSTANTS

!C EARTH RADIUS (m)
      PARAMETER(RAYTER=6.371E6)
!C ACCELERATION DUE TO GRAVITATY (m/s2)
      PARAMETER(G=9.81)
!C ROTATION RATE OF EARTH (1/s)
      PARAMETER(OMEGA = 7.292E-5, OMEGA2 = 2.0*OMEGA)
!C GAMMA AND RATIO OF SPECIFIC HEATS FOR DRY AIR
      PARAMETER(GAMAIR = 1.4, RCPAIR = (GAMAIR - 1.0) / GAMAIR)
!C EARTH SOLID ROTATION SPEED AT THE EQUATOR (m/s)
      PARAMETER(EARTHSP=RAYTER*OMEGA)
!C SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR (J/kg/K)
      PARAMETER(CPAIR=1005.)
!C THE GAZ CONSTANT FOR DRY AIR (J/kg/K)
      PARAMETER(RAIR=RCPAIR*CPAIR)
!C SPECIFIC HEAT AT CONSTANT VOLUME FOR DRY AIR (J/kg/K)
      PARAMETER(CVAIR=CPAIR-RAIR)
!C ZERO DEGREE CELSIUS IN KELVIN (K)
      PARAMETER(T0C=273.15)
!C PRESSURE OF 1 BAR (N/m2)
      PARAMETER(PRESEA=1.E5)
!C MASS OF EARTH (1.e18 kg)
      PARAMETER(MASSEAR=5.977E6)
!C MASS OF THE ATMOSPHERE ( 1e18 kg)
      PARAMETER(MASSATM=5.3)
!C MASS OF OCEAN (1e18 kg)
      PARAMETER(MASSOCE=1.4E3)
!C STEFAN'S CONSTANT (W/m2/K4)
      PARAMETER(STEFC=5.67E-8)

!C========================================================================
!C parameters for latitude of the model.
!C========================================================================
      real :: lat_angle(nlat)
      data lat_angle / 2.768903, 8.306703, 13.84448, 19.38223,&
                       24.91993, 30.45755, 35.99508, 41.53246,&
                       47.06964, 52.60653, 58.14296, 63.67863,&
                       69.21297, 74.74454, 80.26878, 85.76059,&
                      -2.768903,-8.306703,-13.84448,-19.38223,&
                      -24.91993,-30.45755,-35.99508,-41.53246,&
                      -47.06964,-52.60653,-58.14296,-63.67863,&
                      -69.21297,-74.74454,-80.26878,-85.76059/
       real,parameter ::long_inter=5.625
!=========================================================================
!  observation parameter.
!=========================================================================
    integer,parameter :: nstn=110 
    integer,parameter :: nens=200
    real,parameter    :: c_radius=6.0E6    !zero correlation (m)
    integer,parameter :: nassim_times=120   !assimilation n times

!C========================================================================
!C parameters for ensemble prediction.
!C========================================================================
      integer,parameter :: nstok = 24
      integer,parameter :: nsteq = 0
      integer,parameter :: npair = 5
      integer,parameter :: nfc_day = 15       !forecast time
      integer,parameter :: hourpday = 24      
      integer,parameter :: nfc_hour = 15*24   !15 is the forecast lead time
      integer,parameter :: nout_inter = 24    !the output interval
      integer,parameter :: nsample = 50    !the number of sample

!C========================================================================
!C parameters for bred vectors
!C========================================================================
      integer,parameter :: bcycle = 12
!      integer,parameter :: ncycle = 21    !for analysis state
      integer,parameter :: ncycle = 100    !for ideal 
      integer,parameter :: smth_r = 3

!C========================================================================
!C COMMON BLOCK POUR LES OPERATIONS DIRECTES ET INVERSES SPECTRALES.
!C========================================================================
      PARAMETER(RAYDEFIV = 0.)
      COMMON/CONTROL /MCTRL,MPRINT,MESSAI,MHEMIS,MINIT,MTLG,IERR,   &
                     MMAX
      COMMON /TRCVERT/MVADEB, MVALG, MVSDEB, MVSLG
      INTEGER MVADEB(0:NTRUNC), MVALG(0:NTRUNC)
      INTEGER MVSDEB(0:NTRUNC), MVSLG(0:NTRUNC)
      COMMON /GAUSS/  PTGAUSS,PXGAUSS,C1S1MX2,GCOS,GSIN            &
     ,THETA
      REAL            PTGAUSS(NHLAT),PXGAUSS(NHLAT),C1S1MX2(NHLAT)  &
     ,GCOS(NLAT),GSIN(NLAT),THETA(NLAT)
      COMMON /COMFFT/ WORK,TRIGS,IFAX
      REAL            WORK (4*NLONG*NLAT)
      REAL            TRIGS(NHLONG * 3)
      INTEGER         IFAX(10)
      COMMON /POLYNLG/PLM,PW,ALM
      REAL            PLM(NBDEGLIB,NHLAT),PW(NBDEGLIB,NHLAT),      &
                     ALM(NBDEGLIB,NHLAT)
      COMMON/PRINTCTR/KDUMPSP ,KFREE01  ,KPLGNVL ,KPRGAUSS,       &
                     KPRBSPL ,KMSSGACT,KECRTACT,KNOWRITE,         &
                     KDUMPACT,KTEMPACT,KFREE03,KFREE04,           &
                     KIOMSSG,KFREE05,KFREE06
      COMMON/LAPLACCM/OPRLAP,OPRLAPIV
      REAL            OPRLAP(NBDEGLIB,2),OPRLAPIV(NBDEGLIB,2)
      COMMON /DERLONG/ OPRM(NBDEGLIB), OPRTRUNC(NBDEGLIB)

      COMMON/GRILLE/ GTPP, GTVP, GTPG, GTVG, GTJA
      REAL            GTPP   (N2LONG,NLAT),                      &
                     GTVP   (N2LONG,NLAT),                      &
                     GTPG   (N2LONG,NLAT),                      &
                     GTVG   (N2LONG,NLAT),                      &
                     GTJA   (N2LONG,NLAT)
      COMMON /SPECPROV/ STPP, STVP
      REAL              STPP   (NCPLX,NBDEGLIB),                &
                       STVP   (NCPLX,NBDEGLIB)
      COMMON/SPECTRAL/ SPPSI, SPVOR, SPJAC
      REAL             SPPSI   (NCPLX,NBDEGLIB),                &
                      SPVOR   (NCPLX,NBDEGLIB),                 &
                      SPJAC   (NCPLX,NBDEGLIB)
      COMMON /WFOUR/ ZFS, ZFA
      REAL           ZFS (NCPLX, 0:NHLONG, NHLAT),              &
                    ZFA (NCPLX, 0:NHLONG, NHLAT)
      PARAMETER (NFOUR = NCPLX * NHLAT * NH1LONG) 
      REAL ZEFS (NFOUR), ZEFA (NFOUR)
      EQUIVALENCE (ZFS,ZEFS), (ZFA,ZEFA)
      COMMON/MESTEMPS/ T1,T2,T3,T4
!C==================================================================
!C LECTURE DES DONNEES METEO-FRANCE
      COMMON /INTERFA/ INTERF
      INTEGER INTERF(NBDEGLIB)
!C==================================================================


