
!C***********************************************************************
!C                                                              MSSG
!C***********************************************************************
!C
      SUBROUTINE MSSG(MSSGSTR) 
!C
!C     R.BUTEL    L.M.D.         6/83
!C
      INCLUDE 'sphectra.h'
      CHARACTER*8 MSSGSTR
!C
      IF(AND(MCTRL,KNOWRITE).NE.0) RETURN
      IF(AND(MCTRL,KMSSGACT).NE.0)THEN
        WRITE(6,9000) MSSGSTR,SECOND()
      ENDIF
      RETURN
 9000 FORMAT(5X,A8,' A',F20.10)
      END
!C
!C
!C***********************************************************************
!C                                                           INIT CTRL
!C***********************************************************************
!C
      SUBROUTINE INITCTRL 
!C
!C     R.BUTEL    L.M.D.         6/83  VERSION 2 (3/10/83)
!C
      INCLUDE 'sphectra.h'
      IF(AND(MCTRL,KMSSGACT).NE.0) CALL MSSG('INITCTRL')
!C   --OPTION DE MISE AU POINT ( PAR DEFAUT )
!C
!C-----------------------------------------------------------------------
!C                                                      INIT CTRL.1
!C--------------------------------------------INITIALISATION DES TRCCTRLS
!C
      MVADEB(0) = 1
      MVSDEB(0) = NI + 1
      MVALG(0)  = 1 + (NTRUNC - 1)/2
      MVSLG(0)  = 1 + NTRUNC/2
!C
      DO 1007 JMOD = 1, NTRUNC
        MVADEB(JMOD) = MVADEB(JMOD-1) + MVALG(JMOD-1)
        MVSDEB(JMOD) = MVSDEB(JMOD-1) + MVSLG(JMOD-1)
        MVALG (JMOD) = 1 + (NTRUNC - JMOD - 1)/2
        MVSLG (JMOD) = 1 + (NTRUNC - JMOD)/2
!C
 1007 CONTINUE
      MVALG(NTRUNC) = 0
!C
!C
!C-----------------------------------------------------------------------
!C                                                       INIT CTRL.3
!C----------------------------------------------INITIALISATION DE CONTROL
      MCTRL = 0
      MPRINT = 0
      NCORIOL = 0
      MINIT   = 0
      MHEMIS  = 2
      MMAX    = 11
!C
!C-----------------------------------------------------------------------
!C                                                  INIT CTRL.4
!C-----------------------------------------------CONTROLE DES IMPRESSIONS
!C
!C   -- INITIALISATION DU COMMON /PRINTCTR/
!C   -- PEUT SERAIT T IL MIEUX DE FAIRE UM BLOC DATA ; VOIR
      KDUMPSP  =    1
      KFREE01   =    2
      KPLGNVL  =    4
      KPRGAUSS =    8
      KPRBSPL  =   16
      KMSSGACT =   32
      KECRTACT =   64
      KNOWRITE =  128
      KDUMPACT =  256
      KTEMPACT =  512
      KFREE03 = 1024
      KFREE04 = 2048
      KIOMSSG  = 4096
      KFREE05 = 8192
      KFREE06 =16384
!C
!C
       RETURN
      END
!C
!C***********************************************************************
!C                                                      INIT TRSF
!C***********************************************************************
!C
      SUBROUTINE INITTRSF(LEGNFIL) 
!C
!C     R.BUTEL    L.M.D.         9/83
!C
      INCLUDE 'sphectra.h'
!C
      IF(AND(MCTRL,KMSSGACT).NE.0) CALL MSSG('INITTRSF')
!C
!C*    1.  Lecture des valeurs des Fonctions de Legendre et coeffts de Gauss
!C     ---------------------------------------------------------------------
!C
      CALL READDSK(LEGNFIL)
!C
!C*    2.  Integration des coefficients dans les tableaux de Legendre
!C     --------------------------------------------------------------
!C
      RAYM1 = 1. / RAYTER
      DO 210 JLAT = 1, NHLAT
        COEF = PXGAUSS (JLAT)
       DO 205 JC = 1, NBDEGLIB
          PW(JC,JLAT) = COEF * PLM(JC,JLAT)
          ALM(JC,JLAT) = RAYM1 * ALM(JC,JLAT)
  205   CONTINUE
  210 CONTINUE
!C
!C*    3.  Initialisation de COMFFT
!C     ----------------------------
!C
      call fftfax(nlong,ifax,trigs)
!C
!C*    4.  Initialisation de OPERA
!C     ---------------------------
!C
      RAYM2 = 1. / ( RAYTER * RAYTER )
      RAYM3 = RAYDEFIV ** 2
      DO 3000 JMOD = 0, NTRUNC
        IS = MVSDEB(JMOD)
        IN = MVSLG (JMOD)
        JD = 0
        IF (JMOD.EQ.0) JD = 1
        DO 3010 JN = JD, IN - 1
          LN = JMOD + 2 * JN
          OPRLAP   (IS + JN, 2) = - FLOAT (LN * (LN+1)) * RAYM2
          OPRLAPIV(IS + JN, 2) = 1. / OPRLAP(IS + JN, 2)
          OPRLAP   (IS + JN, 1) = OPRLAP(IS + JN, 2) - RAYM3
          OPRLAPIV(IS + JN, 1) = 1. / OPRLAP(IS + JN, 1)
          OPRM     (IS + JN) = FLOAT(JMOD) * RAYM1
 3010   CONTINUE
        IS = MVADEB(JMOD)
        IN = MVALG (JMOD)
        IF(IN.EQ.0) GOTO 3025
        JD = 0
        DO 3020 JN = JD, IN - 1
          LN = JMOD + 2 * JN + 1
          OPRLAP   (IS + JN, 2) = - FLOAT (LN * (LN+1)) * RAYM2
          OPRLAPIV(IS + JN, 2) = 1. / OPRLAP(IS + JN, 2)
          OPRLAP   (IS + JN, 1) = OPRLAP(IS + JN, 2) - RAYM3
          OPRLAPIV(IS + JN, 1) = 1. / OPRLAP(IS + JN, 1)
          OPRM     (IS + JN) = FLOAT(JMOD) * RAYM1
 3020   CONTINUE
 3025   CONTINUE
 3000 CONTINUE
      OPRLAP (NI+1,1) = 0.
      OPRLAPIV (NI+1,1) = 0.
      OPRLAP (NI+1,2) = 0.
      OPRLAPIV (NI+1,2) = 0.
      OPRM (NI+1) = 0.
!C
      RETURN
      END
!C
!C
!C***********************************************************************
!C                                                      READ DSK
!C***********************************************************************
!C
      SUBROUTINE READDSK(LEGNFIL) 
!C
!C     R.BUTEL    L.M.D.         9/83
!C
!C
!C     OBJET : LIT LES VALEURS DES POINTS ET POIDS DE GAUSS,
!C             AINSI QUE CELLES DES POLYNOMES DE LEGENDRE ASSOCIES
!C
!C     SPECIFICATIONS :
!C
!C          LA LECTURE A LIEU DANS L'ORDRE SUIVANT
!C
!C            RECORD 1 : N TRUNC, ENTIER DONNANT LA TRONCATURE POUR
!C                               LAQUELLE PT,PX ET PLM ONT ETE CALCULES
!C            LA TRONCATURE UTILISEE EST EN CONTROLE VERTICAL
!C            RECORDS 2 : (PT(*),PX(*),C1S1MX2(*)),*=1,NH LAT
!C
!C            RECORDS 3 : PLM(1,2)
!C
!C            RECORDS 4 : ALM
!C
!C            RECORDS 5 : PW
!C
      character*80 LEGNFIL
      INCLUDE 'sphectra.h'
!C
      IF(AND(MCTRL,KMSSGACT).NE.0) CALL MSSG('READDSK')
!C
!C--LECTURE ET CONTROLE DE LA TRONCATURE
!C
      OPEN (10,FILE=legnfil  & 
          ,FORM='unformatted',STATUS='old')
      REWIND 10
      READ (10) MT
      if(MT.ne.NTRUNC) then
        write(lw,*)'**** ERREUR: coefficients lus pour la troncature ',MT
        write(lw,*) 'ARRET par READ DSK'
        stop
      endif
!c      print*,MT
!C
!C--LECTURE DES POINTS ET POIDS DE GAUSS
!C
      READ (10) PTGAUSS, PXGAUSS, C1S1MX2
!c      print*,'First',PT GAUSS, PX GAUSS, C1S1MX2
!C
!C--LECTURE DES POLYNOMES DE LEGENDRE
!C
        READ (10)     &
      ((PLM(JCDEGLIB,JCLAT),JCDEGLIB=1,NBDEGLIB),JCLAT=1,NHLAT)
!c        print*,'Second',PLM
!C
!C--LECTURE DES PSEUDO DERIVEES DES POLYNOMES DE LEGENDRE
!C
        READ (10)     &
       ((ALM(JCDEGLIB,JCLAT),JCDEGLIB=1,NBDEGLIB),JCLAT=1,NHLAT)
!c        print*,'Third',ALM
!c        stop
!C  
       CLOSE (10)
       RETURN
       END
!C
!!C
!C******************************************************************************
!C                                                           LAPLACE
!C******************************************************************************
!C
      SUBROUTINE LAPLACE (SPOUT,SPIN,IDIR) 
!C
!C**** LAPLACE  laplacien spectral
!C
!C     B. Legras      LMD        25-05-84
!C
!C     Objet: calcule le laplacien ou resout l'equation de Poisson en spectral
!C
!C     Arguments:
!C                SP IN  : champ spectral en entree
!C                SP OUT : champ spectral en sortie, peut etre identique a
!C                         SP IN
!C                IDIR   : choix de l'operation
!C                         IDIR = 2 : laplacien direct
!C                         IDIR =-2 : laplacien inverse(resolution de Poisson)
!C                         IDIR = 1 : laplacien modifie direct
!C                         IDIR =-1 : laplacien modifie inverse
!C
!C     Note: Les valeurs propres du laplacien sont rangees dans un tableau
!C           de dimension NB DEGLIB. En mode hemispherique, on n'utilise
!C           que la premiere partie des tableaux (de longueur NI).
!C
!C******************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL  SPIN (NCPLX,NBDEGLIB)
      REAL  SPOUT(NCPLX,NBDEGLIB)
!C
!C*    1.  Laplacien direct
!C     --------------------
!C
      IF(IDIR.GT.0) THEN
       DO 110 JC = 1, NI + (MHEMIS-1) * NP
         SPOUT(1,JC) = SPIN(1,JC) * OPRLAP(JC, IDIR)
         SPOUT(2,JC) = SPIN(2,JC) * OPRLAP(JC, IDIR)
  110  CONTINUE
       RETURN
!C
!C*    2.  Laplacien inverse
!C     ---------------------
!C
      ELSE
       DO 210 JC = 1, NI + (MHEMIS-1) * NP
         SPOUT(1,JC) = SPIN(1,JC) * OPRLAPIV(JC, -IDIR)
         SPOUT(2,JC) = SPIN(2,JC) * OPRLAPIV(JC, -IDIR)
  210  CONTINUE
      ENDIF
!C
      RETURN
      END
!C
!C*************************************************************************
!C                                                           SP DDL
!C*************************************************************************
!C
      SUBROUTINE SPDDL(SPOUT,SPIN,IDEN) 
!C
!C**** SP DDL   Derivation spectrale en longitude
!C
!C     B. Legras          LMD           25-05-84
!C
!C     Objet:  Calcul de la derivee spectrale en longitude
!C
!C     Arguments :SP IN : Champ spectral en entree et en sortie si IDEN = .TRUE.
!C                SP OUT: Champ spectral en sortie si IDEN = .FALSE.
!C                IDEN  : Variable logique controlant si la transformation a
!C                        lieu en place ou hors place
!C
!C     Note: En mode hemispherique, on utilise la premiere partie des tableaux
!C           (de longueur NI).
!C
!C******************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL SPIN  (NCPLX,NBDEGLIB)
      REAL SPOUT (NCPLX,NBDEGLIB)
      LOGICAL IDEN
!C
!C*    1.  Calcul hors place
!C     ---------------------
!C
      IF (.NOT.IDEN) THEN
        DO 120 JC = 1, NI + (MHEMIS-1) * NP
          SPOUT(2,JC) =   SPIN(1,JC) * OPRM(JC)
          SPOUT(1,JC) = - SPIN(2,JC) * OPRM(JC)
  120   CONTINUE
        RETURN
!C
!C*    2.  Calcul sur place
!C     --------------------
!C
      ELSE
        DO 220 JC = 1, NI + (MHEMIS-1) * NP
          T           =    SPIN(1,JC)
          SPIN(1,JC) =  - SPIN(2,JC) * OPRM(JC)
          SPIN(2,JC) =    T           * OPRM(JC)
  220   CONTINUE
      ENDIF
!C
      RETURN
      END
!C
!C
!C****************************************************************************
!C                                                              JACOB S
!C****************************************************************************
!C
      SUBROUTINE JACOBS (SPJAX,SP1,SP2) 
!C
!C**** JACOB S   Jacobien spectral, version spherique
!C
!C     B. Legras             LMD               Creation: 29-05-84
!C                                             Certification: 16-06-84
!C                                             Performance: 4,412 ms pour
!C                                                          N TRUNC = 21
!C
!C     Objet:  Calcul du jacobien de deux champs scalaires spectraux.
!C             Le resultat est lui meme donne sous forme spectrale.
!C
!C     Arguments:
!C             SP 1  : Premier champ en entree
!C             SP 2  : Deuxieme champ en entree
!C             SP JAX: Champ en sortie
!C
!C     Definition:
!C             SP JAC = 1./(R*R) * D (SP 1./(R*R) * DSP 2) / D (MU,LAMBDA)
!C             MU    : sinus de la latitude
!C             LAMBDA: longitude
!C             R     : rayon de la sphere
!C
!C     Note:  Le facteur 1./(R*R) est contenu pour moitie dans la derivation
!C            en longitude et pour moitie dans le tableau des pseudo-derivees
!C            ALM.
!C
!C******************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL SP1  (NCPLX,NBDEGLIB)
      REAL SP2  (NCPLX,NBDEGLIB)
      REAL SPJAX(NCPLX,NBDEGLIB)
!C
!C*    1. Calcul des derivees sur la grille
!C     ------------------------------------
      CALL TLGINVS (GTPG,SP1,ALM,-1)
      CALL TLGINVS (GTVG,SP2,ALM,-1)
!C
!C     1.2 Derives en longitudes en semi-Fourier
      CALL SPDDL (STPP,SP1,.FALSE.)
      CALL SPDDL (STVP,SP2,.FALSE.)
      CALL TLGINVS (GTPP,STPP,PLM,1)
      CALL TLGINVS (GTVP,STVP,PLM,1)
!C
!C*    1.3 Passage sur la grille
      CALL RFFTMLT(GTPG,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
      CALL RFFTMLT(GTVG,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
      CALL RFFTMLT(GTPP,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
      CALL RFFTMLT(GTVP,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
!C
!C*    2.  Calcul du jacobien sur la grille
!C     ------------------------------------
!C
      DO 200 JC = 1, NLAT * N2LONG
        GTJA(JC,1) = GTPP(JC,1) * GTVG(JC,1)     &
                   - GTVP(JC,1) * GTPG(JC,1)
  200 CONTINUE
!C
!C*    Normalisation due aux pseudo-derivees
      DO 220 JCLAT = 1, NHLAT
        COEF = C1S1MX2(JCLAT)
        DO 220 JC = 1, N2LONG
          GTJA(JC, JCLAT + NHLAT) = COEF * GTJA(JC, JCLAT + NHLAT)
          GTJA(JC, JCLAT         ) = COEF * GTJA(JC, JCLAT         )
  220 CONTINUE
!C
!C*    3.  Calcul du jacobien spectral
!C     -------------------------------
!C
      CALL RFFTMLT(GTJA,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,-1)
      CALL TLGDIRS (SPJAX,GTJA,PW,1)
!C
      RETURN
      END
!CC
!C****************************************************************************
!C                                                              PROGRA S
!C****************************************************************************
!C
      SUBROUTINE PROGRAS (SPJAX,SP1,SP2) 
!C
!C**** PROGRA S   , Produit scalaire de deux gradients (version spherique)
!C
!C     G. Brunet             LMD               Creation: 22-02-93

!C
!C     Objet:  Produit scalaire du  gradient de deux champs scalaires spectraux.
!C             Le resultat est lui meme donne sous forme spectrale.
!C
!C     Arguments:
!C             SP 1  : Premier champ en entree
!C             SP 2  : Deuxieme champ en entree
!C             SP JAX: Champ en sortie
!C
!C     Definition:
!C             SP JAC = 1./(R*R) * D (SP 1./(R*R) * DSP 2) / D (MU,LAMBDA)
!C             MU    : sinus de la latitude
!C             LAMBDA: longitude
!C             R     : rayon de la sphere
!C
!C     Note:  Le facteur 1./(R*R) est contenu pour moitie dans la derivation
!C            en longitude et pour moitie dans le tableau des pseudo-derivees
!C            ALM. Cette routine est une modification de JACOB S (seulement la boucle
!C            200 est modifiee).
!C
!C******************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL SP1  (NCPLX,NBDEGLIB)
      REAL SP2  (NCPLX,NBDEGLIB)
      REAL SPJAX(NCPLX,NBDEGLIB)
!C
!C*    1. Calcul des derivees sur la grille
!C     ------------------------------------
      CALL TLGINVS (GTPG,SP1,ALM,-1)
      CALL TLGINVS (GTVG,SP2,ALM,-1)
!C
!C     1.2 Derives en longitudes en semi-Fourier
      CALL SPDDL (STPP,SP1,.FALSE.)
      CALL SPDDL (STVP,SP2,.FALSE.)
      CALL TLGINVS (GTPP,STPP,PLM,1)
      CALL TLGINVS (GTVP,STVP,PLM,1)
!C
!C*    1.3 Passage sur la grille
      CALL RFFTMLT(GTPG,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
      CALL RFFTMLT(GTVG,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
      CALL RFFTMLT(GTPP,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
      CALL RFFTMLT(GTVP,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,1)
!C
!C*    2.  Calcul du jacobien sur la grille
!C     ------------------------------------
!C
      DO 200 JC = 1, NLAT * N2LONG
        GTJA(JC,1) = GTPP(JC,1) * GTVP(JC,1)    &
                   + GTVG(JC,1) * GTPG(JC,1)
  200 CONTINUE
!C
!C*    Normalisation due aux pseudo-derivees
      DO 220 JCLAT = 1, NHLAT
        COEF = C1S1MX2(JCLAT)
        DO 220 JC = 1, N2LONG
          GTJA(JC, JCLAT + NHLAT) = COEF * GTJA(JC, JCLAT + NHLAT)
          GTJA(JC, JCLAT         ) = COEF * GTJA(JC, JCLAT         )
  220 CONTINUE
!C
!C*    3.  Calcul du jacobien spectral
!C     -------------------------------
!C
      CALL RFFTMLT(GTJA,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NLAT,-1)
      CALL TLGDIRS (SPJAX,GTJA,PW,1)
!C
      RETURN
      END
!C
!C****************************************************************************
!C                                                              JACOB H
!C****************************************************************************
!C
      SUBROUTINE JACOBH (SPJAX,SP1,SP2)
!C
!C**** JACOB H   Jacobien spectral, version hemispherique
!C
!C     B. Legras             LMD               Creation: 29-05-84
!C                                             Certification: 16-06-84
!C                                             Performance: 2,318 ms pour
!C                                                          N TRUNC = 21
!C
!C     Objet:  Calcul du jacobien de deux champs scalaires spectraux.
!C             Le resultat est lui meme donne sous forme spectrale.
!C
!C     Arguments:
!C             SP 1  : Premier champ en entree
!C             SP 2  : Deuxieme champ en entree
!C             SP JAX: Champ en sortie
!C
!C     Definition:
!C             SP JAC = 1./(R*R) * D (SP 1./(R*R) * DSP 2) / D (MU,LAMBDA)
!C             MU    : sinus de la latitude
!C             LAMBDA: longitude
!C             R     : rayon de la sphere
!C
!C     Note:  Le facteur 1./(R*R) est contenu pour moitie dans la derivation
!C            en longitude et pour moitie dans le tableau des pseudo-derivees
!C            ALM.
!C
!C******************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL SP1  (NCPLX,NBDEGLIB)
      REAL SP2  (NCPLX,NBDEGLIB)
      REAL SPJAX(NCPLX,NBDEGLIB)
!C
!C*    1. Calcul des derivees sur la grille
!C     ------------------------------------
      CALL TLGINVH (GTPG,SP1,ALM,-1)
      CALL TLGINVH (GTVG,SP2,ALM,-1)
!C
!C     1.2 Derives en longitudes en semi-Fourier
      CALL SPDDL (STPP,SP1,.FALSE.)
      CALL SPDDL (STVP,SP2,.FALSE.)
      CALL TLGINVH (GTPP,STPP,PLM,1)
      CALL TLGINVH (GTVP,STVP,PLM,1)
!C
!C*    1.3 Passage sur la grille
      CALL RFFTMLT(GTPG,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NHLAT,1)
      CALL RFFTMLT(GTVG,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NHLAT,1)
      CALL RFFTMLT(GTPP,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NHLAT,1)
      CALL RFFTMLT(GTVP,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NHLAT,1)
!C
!C*    2.  Calcul du jacobien sur la grille
!C     ------------------------------------
!C
      DO 200 JC = 1, NHLAT * N2LONG
        GTJA(JC,1) = GTPP(JC,1) * GTVG(JC,1)   &
                   - GTVP(JC,1) * GTPG(JC,1)
  200 CONTINUE
!C
!C*    Normalisation due aux pseudo-derivees
      DO 220 JCLAT = 1, NHLAT
        COEF = C1S1MX2(JCLAT)
        DO 220 JC = 1, N2LONG
          GTJA(JC, JCLAT         ) = COEF * GTJA(JC, JCLAT         )
  220 CONTINUE
!C
!C*    3.  Calcul du jacobien spectral
!C     -------------------------------
!C
      CALL RFFTMLT(GTJA,WORK,TRIGS,IFAX,1,N2LONG,NLONG,NHLAT,-1)
      CALL TLGDIRH(SPJAX,GTJA,PW,1)
!C
!C
      RETURN
      END
!C
!C
!C***********************************************************************
!C                                                        TLG DIR
!C***********************************************************************
!C
      SUBROUTINE TLGDIRT(JRNIV,LEGPOL,PARITE,FCHAMP,SPC) 
!C
!C    JR NIV  : IN NOMBRE DE NIVEAUX VERTICAUX
!C              --
!C    LEGPOL  : IN POLYNOMES UTILISES POUR LA TRANSFORMATION
!C              --
!C    PARITE  : IN PARITE DES POLYNOMES UTILISES
!C              --
!C    F CHAMP : IN CHAMP EN COMPOSANTES FOURIER-LONGITUDE/GRILLE LATITUDE
!C              --
!C    SPC     : OUT CHAMP EN COMPOSANTES SPECTRALES ( HARMONIQUES SPHER.)
!C              ---
!C
!C    FORMULES : SPC(N,M) = [[ F(M)(X) * P(N,M)(X) DX
!C
!C      APPROCHE (ET EXACT SUR LA TRONCATURE ) PAR
!C
!C        SOMME(PT GAUSS)(F(M)(PT)*P(N,M)(PT)*PX(PT)
!C
!C    CALCUL : PT > 0  ET PT < 0 (P(N,M)(-PT)=(-1)^(N+M)*P(N,M)(PT) )
!C
!C    ON EFFECTUE DONC 2 SOMMES : SPN(NORD) ET SPS(SUD)
!C
!C    ET LE RESULTAT EST :
!C
!C      SPC(M + D,M) = SPN(D) + (-1)^D * SPS(D)
!C
!C    MAIS ON PEUT AUSSI DECOMPOSER LE CHAMP INITIAL EN PARTIE SYMETRIQUE
!C    ----
!C         ET ANTISYMETRIQUE ( MULTIPLIEES PAR UN COEFFICIENT 2)
!C
!C         ET LES FORMULES DEVIENNENT :
!C
!C         SP C(M + D,M) = SP S(D)     QUAND M + D EST PAIR
!C                       = SP A(D)     QUAND M + D EST IMPAIR
!C
!C
      INCLUDE 'sphectra.h'
      REAL             SPC    (NCPLX,NBDEGLIB)
      REAL             FCHAMP(NCPLX,NH1LONG,NLAT)
      REAL             LEGPOL(NBDEGLIB,NHLAT)
      INTEGER          PARITE
!C
!C   --M TLG = 1 : VERTICAL HALF
!C             2 : VERTICAL HALF + SCALAR
!C             3 : VERTICAL FULL
!C             4 : VERTICAL FULL + SCALAR
!C             5 : VERTICAL HALF + SCALAR ; REVERSE INDEX
!C             6 : DIAGONAL
!C             7 : version hemispherique antisym de 8
!C             8 : version verticale sym/antisym matr+scalaire eco
!C             9 : VERTICAL FULL + REGROUPEMANT (BLG'S IDEA)
!C            10 : VERTICAL HALF + SCALAR ; REVERSE INDEX ; PACKED DECOMP
!C
      GOTO (1111,2222,3333,4444,5555,6666,7777,8888,9999,11110),MTLG
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.1111
!C----------------------------------------------------------VERTICAL HALF
!C
 1111 CONTINUE
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.2222
!C----------------------------------------------VERTICAL  HALF + SCALAIRE
!C
!C
 2222 CONTINUE
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.3333
!C------------------VERSION MATRICIELLE ; TRONCATURE EN CONTROLE VERTICAL
!C
 3333 CONTINUE
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.4444
!C----------------------------------------------VERTICALE FULL + SCALAIRE
!C
 4444 CONTINUE
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.5555
!C-------------------------------VERTICAL HALF + SCALAIRE; REVERSE INDEX
!C
 5555 CONTINUE
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.6666
!C--------------VERSION NON MATRICIELLE ; TRONCATURE EN CONTROLE DIAGONAL
!C
 6666 CONTINUE
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.7777
!C----------version hemispherique antisymetrique de la transformee 8
!C
 7777 CONTINUE
      MCTRL = KMSSGACT
      CALL MSSG('TD.7 DEB')
      T1 = SECOND()
      CALL TLGDIRH(SPC,FCHAMP,LEGPOL,PARITE)
      T2 = SECOND()
      CALL MSSG('TD.7 FIN')
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.8888
!C--------------Version verticale sym/antisym matr+scalaire eco----------
!C
 8888 CONTINUE
      MCTRL = KMSSGACT
      CALL MSSG('TD.8 DEB')
      T1 = SECOND()
      CALL TLGDIRS(SPC,FCHAMP,LEGPOL,PARITE)
      T2 = SECOND()
      MCTRL = KMSSGACT
      CALL MSSG('TD.8 FIN')
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.9999
!C------------------------------------VERTICAL FULL REGROUPE (BLG'S IDEA)
!C
 9999 CONTINUE
      RETURN
!C
!C-----------------------------------------------------------------------
!C                                                       TLG DIR.11110
!C----------------VERTICAL HALF + SCALAIRE; REVERSE INDEX ; PACKED DECOMP
!C
11110 CONTINUE
      RETURN
      END
!C
!C *****************************************************************************
!C                                                           TLG DIR S
!C *****************************************************************************
!C
      SUBROUTINE TLGDIRS (SPC,FCHAMP,LEGPOL,PARITE)
!C
!C**** TLG DIR S  Transformee de Legendre directe
!C
!C     B. Legras        LMD      creation: 14-03-84
!C                               certification: 16-05-84
!C                               performance: 0.569 ms pour N TRUNC = 21
!C                                            et M MAX = 20
!C
!C     Objet:  Transformee de Legendre directe de la representation
!C             semi-Fourier a la representation spectrale.
!C             Troncature spherique complete.
!C             Description verticale de la troncature.
!C             Troncature parametrisable par le COMMON TRCVERT.
!C             En standard: troncature triangulaire homogene.
!C
!C     Arguments:
!C             LEGPOL: tableau des fonctions de Legendre associees ou de
!C                     leurs pseudo-derivees
!C             PARITE: parite de la fonction de transformation contenue
!C                     dans LEGPOL
!C             F CHAMP:tableau des coefficients de semi-Fourier en entree.
!C                     rangement: l'hemisphere Nord d'abord et dans chaque
!C                                hemisphere, de l'Equateur vers le Pole
!C             SPC   : tableau des coefficients d'harmoniques spheriques
!C                     en sortie
!C
!C     Methode:Les sommations en latitude sont regroupees et calculees
!C             par MXMA pour les premieres colonnes de la troncature.
!C             Les dernieres colonnes sont calculees par SDOT.
!C             On separe en entree les contributions symetriques et
!C             antisymetriques.
!C
!C     Externes: MXMA
!C               SDOT
!C
!C ****************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL    LEGPOL  (NBDEGLIB,NHLAT)
      REAL    SPC     (NCPLX,NBDEGLIB)
      REAL    FCHAMP (NCPLX,NH1LONG,NLAT)
      INTEGER PARITE
!C
!C
!C*    1.  Initialisation des increments pour MXMA
!C     -------------------------------------------
!C
      DATA   INC1P,    INC2P,    INC1F,  INC2F,    INC1S,  INC2S  &
         /      1,NBDEGLIB,  N2LONG,      1,    NC1S,      1/
!C
!C*     2. Separation des contributions symetriques et antisymetriques
!C      --------------------------------------------------------------
!C
  200 CONTINUE
      IPLUS = N2LONG * NHLAT
      IF(PARITE.EQ.1) THEN
      DO 210 JC = 1, N2LONG * NHLAT
        ZEFS(JC) = FCHAMP(JC,1,1)   &
                 + FCHAMP(JC+IPLUS,1,1)
        ZEFA(JC) = FCHAMP(JC,1,1)   &
                 - FCHAMP(JC+IPLUS,1,1)
  210 CONTINUE
      ELSE
      DO 220 JC = 1, N2LONG * NHLAT
        ZEFS(JC) = FCHAMP(JC,1,1)  &
                 - FCHAMP(JC+IPLUS,1,1)
        ZEFA(JC) = FCHAMP(JC,1,1)  &
                 + FCHAMP(JC+IPLUS,1,1)
  220 CONTINUE
      ENDIF
!C
!C*    3.  Calcul matriciel pour les ordres jusqu'a M MAX-1
!C     ----------------------------------------------------
!C
  300 CONTINUE
      DO 350 JMOD = 0, MMAX - 1
!C
!C*    3.1 Definition des parametres de boucle
  310 ISA = MVADEB(JMOD)
      ISS = MVSDEB(JMOD)
      INA = MVALG (JMOD)
      INS = MVSLG (JMOD)
      IK = NCPLX
      IF(JMOD.EQ.0) IK = 1
!C
!C*    3.2 Transformee symetrique
  320   CALL MXMA(       &
              LEGPOL  (   ISS ,      1), INC1P, INC2P, &
              ZFS    (1,      JMOD, 1), INC1F, INC2F,  &
              SPC     (1, ISS         ), INC1S, INC2S, &
              INS, NHLAT, IK)
!C
!C*    3.3 Transformee antisymetrique
  330   IF (INA.NE.0) THEN
        CALL MXMA(   &
              LEGPOL  (   ISA ,      1), INC1P, INC2P, &
              ZFA    (1,      JMOD, 1), INC1F, INC2F,  &
              SPC     (1, ISA         ), INC1S, INC2S, &
              INA, NHLAT, IK)
        ENDIF
!C
  350 CONTINUE
!C
!C*    4.  Calcul scalaire de l'ordre M MAX a N TRUNC
!C     ----------------------------------------------
!C
  400 DO 450 JMOD = MMAX, NTRUNC
        ISA = MVADEB(JMOD)
        ISS = MVSDEB(JMOD)
        INA = MVALG(JMOD)
        INS = MVSLG(JMOD)
!C
!C*    4.1 Partie symetrique
        DO 415 JN = 0, INS - 1
          SPC( 1, ISS+JN ) = SDOT(NHLAT, LEGPOL(ISS+JN , 1), INC2P, &
                                         ZFS  ( 1, JMOD,1), INC1F)
          SPC( 2, ISS+JN ) = SDOT(NHLAT, LEGPOL(ISS+JN , 1), INC2P, &
                                         ZFS  ( 2, JMOD,1), INC1F)
  415   CONTINUE
!C
!C*    4.2 Partie antisymetrique
        IF(INA.NE.0) THEN
        DO 425 JN = 0, INA - 1
          SPC( 1, ISA+JN ) = SDOT(NHLAT, LEGPOL(ISA+JN , 1), INC2P, &
                                         ZFA  ( 1, JMOD,1), INC1F)
          SPC( 2, ISA+JN ) = SDOT(NHLAT, LEGPOL(ISA+JN , 1), INC2P, &
                                         ZFA  ( 2, JMOD,1), INC1F)
  425   CONTINUE
        ENDIF
  450 CONTINUE
!C
!C
      RETURN
      END
!C
!C *****************************************************************************
!C                                                           TLG DIR H
!C *****************************************************************************
!C
      SUBROUTINE TLGDIRH (SPC,FCHAMP,LEGPOL,PARITE)
!C
!C**** TLG DIR H  Transformee de Legendre directe hemispherique antisymetrique
!C
!C     B. Legras        LMD      creation: 17-05-84
!C                               certification: 22-05-84
!C                               performance: 0.283 ms pour N TRUNC = 21
!C                                            et M MAX = 19
!C
!C     Objet:  Transformee de Legendre directe de la representation
!C             semi-Fourier a la representation spectrale.
!C             Transformee hemispherique pour un champ final antisymetrique.
!C             NB: deux possibilites : champ de grille antisym + parite=1
!C                                     champ de grille sym     + parite=-1
!C                 la transformee reste la meme
!C             Description verticale de la troncature.
!C             Troncature parametrisable par le COMMON TRCVERT.
!C             En standard: troncature triangulaire homogene.
!C
!C     Arguments:
!C             LEGPOL: tableau des fonctions de Legendre associees ou de
!C                     leurs pseudo-derivees
!C             PARITE: parite de la fonction de transformation contenue
!C                     dans LEGPOL
!C             F CHAMP:tableau des coefficients de semi-Fourier en entree.
!C                     rangement: ce tableau ne contient qu'un hemisphere
!C                                (l'hemisphere Nord) range de l'Equateur
!C                                vers le pole
!C                                NB: en cas d'equivalence avec un tableau
!C                                spherique complet, tenir compte du fait
!C                                que l'hemisphere Sud est stocke dans la
!C                                premiere moitie de ce dernier.
!C             SPC   : tableau des coefficients d'harmoniques spheriques
!C                     en sortie
!C
!C     Methode:Les sommations en latitude sont regroupees et calculees
!C             par MXMA pour les premieres colonnes de la troncature.
!C             Les dernieres colonnes sont calculees par SDOT.
!C
!C     Externes: MXMA
!C               SDOT
!C
!C     Precautions d'emploi: hormis celles signales plus haut,
!C                           limiter M MAX a N TRUNC car il n'y
!C                           a pas de calcul pour JMOD = N TRUNC
!C
!C ****************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL    LEGPOL  (NBDEGLIB,NHLAT)
      REAL    SPC     (NCPLX,NBDEGLIB)
      REAL    FCHAMP (NCPLX,NH1LONG,NHLAT)
      INTEGER PARITE
!C
!C
!C*    1.  Initialisation des increments pour MXMA
!C     -------------------------------------------
!C
      DATA   INC1P,    INC2P,    INC1F,  INC2F,    INC1S,  INC2S  &
         /      1,NBDEGLIB,  N2LONG,      1,    NC1S,      1/
!C
!C*     2. Contribution totale: multiplication par 2 en entree
!C      ------------------------------------------------------
!C
  200 CONTINUE
      DO 210 JC = 1, N2LONG * NHLAT
        ZEFA(JC) = 2 * FCHAMP(JC,1,1)
  210 CONTINUE
!C
!C*    3.  Calcul matriciel pour les ordres jusqu'a M MAX-1
!C     ----------------------------------------------------
!C
  300 CONTINUE
      DO 350 JMOD = 0, MMAX - 1
!C
!C*    3.1 Definition des parametres de boucle
  310 ISA = MVADEB(JMOD)
      INA = MVALG (JMOD)
      IK = NCPLX
      IF(JMOD.EQ.0) IK = 1
!C
!C*    3.3 Transformee antisymetrique
  330   CALL MXMA(      &
              LEGPOL  (   ISA ,      1), INC1P, INC2P,  &
              ZFA    (1,      JMOD, 1), INC1F, INC2F,  &
              SPC     (1, ISA         ), INC1S, INC2S, &
              INA, NHLAT, IK)
!C
  350 CONTINUE
!C
!C*    4.  Calcul scalaire de l'ordre M MAX a N TRUNC
!C     ----------------------------------------------
!C
  400 DO 450 JMOD = MMAX, NTRUNC - 1
        ISA = MVADEB(JMOD)
        INA = MVALG (JMOD)
!C
!C*    4.2 Partie antisymetrique
        DO 450 JN = 0, INA - 1
          SPC( 1, ISA+JN ) = SDOT(NHLAT, LEGPOL(ISA+JN , 1), INC2P,  &
                                         ZFA  ( 1, JMOD,1), INC1F)
          SPC( 2, ISA+JN ) = SDOT(NHLAT, LEGPOL(ISA+JN , 1), INC2P, &
                                         ZFA  ( 2, JMOD,1), INC1F)
  450 CONTINUE
!C
!C
      RETURN
      END
!C


!C***********************************************************************
!C                                                          TLG INV T
!C***********************************************************************
!C
      SUBROUTINE TLGINVT(JRNIV,LEGPOL,PARITE,SPC,FCHAMP) 
      INCLUDE 'sphectra.h'
      REAL            SPC    (NCPLX,NBDEGLIB,     JRNIV)
      REAL            FCHAMP(NCPLX,NH1LONG,NLAT,JRNIV)
      REAL            LEGPOL (NBDEGLIB,NHLAT)
      INTEGER         PARITE
!C
!C     JR NIV : IN NOMBRE DE NIVEAUX REEL DU CHAMP A TRANSFORMER
!C              --
!C     LEGPOL : IN POLYNOMES UTILISES POUR LA TRANSFORMEE
!C              --
!C     PARITE : IN PARITE DES POLYNOMES UTILISES
!C              --
!C     SPC : IN CHAMP DE COMPOSANTES SPECTRALES
!C           --
!C     F CHAMP : OUT CHAMP COMPOSANTE FOURIER(LONGITUDE)/GRILLE(LATITUDE)
!C               ---
!C
!C     FORMULE :
!C
!C       F(M)(X) = SOMME(N) (P(N,M)(X) * SPC(N,M) )
!C
!C     PAR SYMETRIE :
!C
!C       F(M)(-X) = SOMME(D PAIR) (P(M+D,M)(X) * SPC(M+D,M))
!C                - SOMME(D IMPAIR) (   "    )
!C
!C       ON CALCULE DONC CES 2 SOMMES SEPARAMENT,
!C       LEUR RECOMBINAISON DONNE F(M)(X) ET F(M)(-X) POUR X>0
!C
!C   --M TLG=1 : vertical sym/antisym matricielle reel/imag separes
!C      "   2  : vertical sym/antisym matricielle eco
!C      "   3  : version hemispherique antisym de 2
!C      "   4  : DIAGONAL (NH LAT)
!C      "   5  : DIAGONAL (PROPER) + SCALAR
!C      "   6  : DIAGONAL VLC (NH LAT)
!C      "   7  : DIAGONAL VLC ( PROPER)
!C
!C----------------------------------------------------------------------
!C
      GOTO (1111,2222,3333), MTLG
!C
!C-----------------------------------------------------------------------
!C                                                       TLG INV.1111
!C------------------vertical sym/antisym matricielle reel/imag separes
!C
 1111 CONTINUE
      RETURN
!C
!C------------------------------------------------------------------------
!C                                                        TLG INV.2222
!C------------------vertical sym/antisym matricielle eco
!C
 2222 CONTINUE
      MCTRL = KMSSGACT
      CALL MSSG('TI.2.DEB')
      T1 = SECOND()
      CALL TLGINVS(FCHAMP,SPC,LEGPOL,PARITE)
      T2 = SECOND()
      CALL MSSG('TI.2.FIN')
      RETURN
!C
!C-------------------------------------------------------------------------
!C                                                         TLG INV.3333
!C-------------------version hemispherique antisymetrique de la transformee 2
!C
 3333 CONTINUE
      MCTRL = KMSSGACT
      CALL MSSG('TI.3.DEB')
      T1 = SECOND()
      CALL TLGINVH(FCHAMP,SPC,LEGPOL,PARITE)
      T2 = SECOND()
      CALL MSSG('TI.3.FIN')
      RETURN
!C
      END
!C *****************************************************************************
!C                                                        INTERFACE MF
!C *****************************************************************************
!C
      SUBROUTINE INTERFACE
!C
!C Interface spectrale pour LMD et meteo France 
!C     G. Brunet.
!C

!C ****************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      ntab=0

      do jmod=0,ntrunc
         ISA = MVADEB(JMOD)
         ISS = MVSDEB(JMOD)
         INA = MVALG (JMOD)
         INS = MVSLG (JMOD)  
      do i=1,ins
         interf(ntab+2*I-1)=iss+i-1
      enddo
      do i=1,ina
         interf(ntab+2*i)=isa+i-1
      enddo
      
!c      write(6,*)jmod,ntab,iss,isa,ina,ins
      ntab=ntab+ina+ins
     
      enddo
         
!C     
      RETURN
      END

!CC *****************************************************************************
!C                                                        TLG INV S
!C *****************************************************************************
!C
      SUBROUTINE TLGINVS (FCHAMP,SPC,LEGPOL,PARITE)
!C
!C**** TLG INV S  Transformee de Legendre inverse 
!C
!C     B. Legras        LMD      Creation: 16-04-84
!C                               Certification: 16-05-84
!C                               Performance: 0.5296 ms pour N TRUNC = 21
!C
!C     Objet:  Transformee de Legendre inverse de la representation
!C             spectrale a la representation en semi-Fourier.
!C             Troncature spherique complete.
!C             Description verticale de la troncature parametrisee
!C             par le common TRCVERT.
!C             En standard: troncature triangulaire homogene.
!C
!C     Arguments:
!C             LEGPOL:  Tableau des fonctions de Legendre associees ou de
!C                      leur pseudo-derivees.
!C             PARITE:  Parite de la fonction de transformation contenue
!C                      dans LEGPOL.
!C             SPC   :  Tableau des coefficients d'harmoniques spheriques.
!C                      En entree.
!C             F CHAMP: Tableau des coefficients de semi-Fourier.
!C                      En sortie.
!C                      Rangement: l'hemisphere Nord d'abord et dans chaque
!C                      hemisphere, de l'Equateur vers le Pole.
!C
!C     Methode: Les sommations en l pour les colonnes de la troncature
!C             sont regroupees et calculees ensemble pour les differentes
!C             latitudes.
!C             Les contributions symetriques et antisymetriques sont
!C             calculees separement et combinees en sortie.
!C
!C     Externes: MXMA
!C
!C ****************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL    LEGPOL   (NBDEGLIB,NHLAT)
      REAL    SPC      (NCPLX,NBDEGLIB)
      REAL    FCHAMP  (NCPLX,0:NHLONG,NLAT)
      INTEGER PARITE
!C
!C
!C*    1. Definition des increments de MXMA
!C     ------------------------------------
!C
      DATA      INC1P, INC2P,   INC1S, INC2S,    INC1F,  INC2F  &
         / NBDEGLIB,     1,   NC1S,     1,  N2LONG,      1/
!C
!C*    2. Calcul matriciel de la sommation sur le degre
!C     ------------------------------------------------
!C
  200 CONTINUE
      DO 250 JMOD = 0, NTRUNC
!C
!C*    2.1 Definition des parametres de boucle
  210   ISA = MVADEB(JMOD)
        ISS = MVSDEB(JMOD)
        INA = MVALG (JMOD)
        INS = MVSLG (JMOD)
        IK = NCPLX
        IF(JMOD.EQ.0) IK = 1
!C
!C*    2.2 Transformee symetrique

  220   CALL MXMA(    &
              LEGPOL  (   ISS ,       1), INC1P, INC2P,  &
              SPC     (1, ISS          ), INC1S, INC2S,  &
              ZFS    (1,       JMOD, 1), INC1F, INC2F,  &
              NHLAT, INS, IK)
!C
!C*    2.3 Transformee antisymetrique
  230   IF (INA.NE.0) THEN
        CALL MXMA(     &
              LEGPOL  (   ISA ,       1), INC1P, INC2P,  &
              SPC     (1, ISA          ), INC1S, INC2S,  &
              ZFA    (1,       JMOD, 1), INC1F, INC2F,   &
              NHLAT, INA, IK)
        ELSE
        DO 235 JC = 1, NHLAT
          ZFA(1,JMOD,JC) = 0.
          ZFA(2,JMOD,JC) = 0.
  235   CONTINUE
        ENDIF
!C
  250 CONTINUE
      
!C
!C*    3.  Combinaison des parties symetriques et antisymetriques
!C     ----------------------------------------------------------
!C
  300 CONTINUE
      IF (PARITE.EQ.1) THEN
      DO 305 JC = 1, NHLAT
        FCHAMP(1,0,JC+NHLAT) = ZFS(1,0,JC) - ZFA(1,0,JC)
        FCHAMP(1,0,JC       ) = ZFS(1,0,JC) + ZFA(1,0,JC)
        FCHAMP(2,0,JC+NHLAT) = 0.
        FCHAMP(2,0,JC       ) = 0.
  305 CONTINUE
      DO 310 JMOD = 1, NTRUNC
      DO 310 JC = 1, NHLAT
        FCHAMP(1,JMOD,JC+NHLAT) = ZFS(1,JMOD,JC) - ZFA(1,JMOD,JC)
        FCHAMP(1,JMOD,JC       ) = ZFS(1,JMOD,JC) + ZFA(1,JMOD,JC)
        FCHAMP(2,JMOD,JC+NHLAT) = ZFS(2,JMOD,JC) - ZFA(2,JMOD,JC)
        FCHAMP(2,JMOD,JC       ) = ZFS(2,JMOD,JC) + ZFA(2,JMOD,JC)
  310 CONTINUE
!C
      ELSE
      DO 315 JC = 1, NHLAT
        FCHAMP(1,0,JC+NHLAT) = ZFA(1,0,JC) - ZFS(1,0,JC)
        FCHAMP(1,0,JC       ) = ZFA(1,0,JC) + ZFS(1,0,JC)
        FCHAMP(2,0,JC+NHLAT) = 0.
        FCHAMP(2,0,JC       ) = 0.
  315 CONTINUE
      DO 320 JMOD = 1, NTRUNC
      DO 320 JC = 1, NHLAT
        FCHAMP(1,JMOD,JC+NHLAT) = ZFA(1,JMOD,JC) - ZFS(1,JMOD,JC)
        FCHAMP(1,JMOD,JC       ) = ZFA(1,JMOD,JC) + ZFS(1,JMOD,JC)
        FCHAMP(2,JMOD,JC+NHLAT) = ZFA(2,JMOD,JC) - ZFS(2,JMOD,JC)
        FCHAMP(2,JMOD,JC       ) = ZFA(2,JMOD,JC) + ZFS(2,JMOD,JC)
  320 CONTINUE
      ENDIF
!C
!C*    4.  Mise a zero des modes de Fourier non representes
!C     ----------------------------------------------------
!C
  400 CONTINUE
      DO 410 JMOD = NTRUNC+1, NHLONG
      DO 410 JC = 1, NLAT
        FCHAMP(1,JMOD,JC) = 0.
        FCHAMP(2,JMOD,JC) = 0.
  410 CONTINUE
!C
!C
      RETURN
      END
!C
!C
!C*******************************************************************************
!C                                                   TEST GF
!C*******************************************************************************
!C
      SUBROUTINE TESTGF(GFT,GFR,EPSIL) 
!C     Calcule la difference des deux champs de grille GF T et GF R et verifie
!C     que celle-ci reste partout inferieure a EPSIL.
!C     Impression d'un message diagnostic.
!C     GF R : champ de reference, non modifie.
!C     GF T : champ a tester, contient la difference en sortie.
!C
      INCLUDE 'sphectra.h'
!C
      DIMENSION GFR(N2LONG,NLAT), GFT(N2LONG,NLAT)
!C
      IF(MHEMIS.EQ.2) THEN
        MSIZE = N2LONG * NLAT
        MLAT = NLAT
      ENDIF
      IF(MHEMIS.EQ.1) THEN
        MSIZE = N2LONG * NHLAT
        MLAT = NHLAT
      ENDIF
!C
      DO 100 J = 1, MSIZE
        GFT(J,1) = GFT(J,1) - GFR(J,1)
100   CONTINUE
      DO 110 JC = 1, MLAT
        GFT(NLONG+1,JC) = 0.
        GFT(N2LONG,JC) = 0.
  110 CONTINUE
!C
      GR = 0.
      GD = 0.
      DO 200 J = 1, MSIZE
        GD = AMAX1(GD, ABS(GFT(J,1)))
        GR = AMAX1(GR, ABS(GFR(J,1)))
200   CONTINUE
!C
      IF (ABS(GD/GR).GE.EPSIL) THEN
        WRITE(6,5000) GD, GR
      ELSE
        WRITE(6,5001)
      ENDIF
      RETURN
5000  FORMAT(1X,'CHAMP FINAL SUR LA GRILLE DIFFERENT DU CHAMP INITIAL',&
           /1X,'**ERREUR A RECHERCHER**',&
            /1X,'GD = ',F20.10,&
            /1X,'GR = ',F20.10)
5001  FORMAT(1X,'GR FINAL CORRECT')
      END
!C
!C*******************************************************************************
!C                                                    TEST SP
!C*******************************************************************************
!C
      SUBROUTINE TESTSP(SPT,SPR,EPSIL) 
!C
!C
!C     Calcule la difference de deux champs spectraux et verifie que celle-ci
!C     reste sur tous les modes inferieure a EPSIL.
!C     Imprime un message diagnostic.
!C     SP R : Champ de reference, non modifie.
!C     SP T : Champ a tester, contient la difference en sortie.
!C
      INCLUDE 'sphectra.h'
!C
      DIMENSION SPT(NCPLX, NBDEGLIB)
      DIMENSION SPR(NCPLX, NBDEGLIB)
!C
      IF(MHEMIS.EQ.2) MSIZE = NCPLX * NBDEGLIB
      IF(MHEMIS.EQ.1) MSIZE = NCPLX * NI
!C
      DO 100 J = 1, MSIZE
        SPT(J,1) = SPT(J,1) - SPR(J,1)
100   CONTINUE
!C
      SR = 0.
      SD = 0.
      DO 200 J = 1, MSIZE
        SD = AMAX1(SD,ABS(SPT(J,1)))
        SR = AMAX1(SR,ABS(SPR(J,1)))
200   CONTINUE
!C
      IF (ABS(SD/SR).GE.EPSIL) THEN
        WRITE(6,5000) SD, SR
      ELSE
        WRITE(6,5001)
      ENDIF
      RETURN
5000  FORMAT(1X,'CHAMP FINAL SPECTRAL DIFFERENT DU CHAMP INITIAL',  &
            /1X,'**ERREUR A RECHERCHER**',&
            /1X,'SD = ',F20.10,&
            /1X,'SR = ',F20.10)
5001  FORMAT(1X,'SR FINAL CORRECT')
      END
!C
!C***********************************************************************
!C                                                         ACQUIRDF
!C***********************************************************************
!C
      SUBROUTINE ACQUIRDF(PARITE) 
!C
!C       Generateur d'un champ test sur la grille de collocation
!C       PARITE definit la parite de ce champ
!C
      INCLUDE 'sphectra.h'
      INTEGER PARITE
!C
!C*    1. Definition de quelques fonctions de Legendre associees
!C     ---------------------------------------------------------
!C
!C     1.1 Fonctions antisymetriques (m+l impair)
!C
      Y01(X,Y) =  Y
      Y03(X,Y) = Y * (5.*Y*Y - 3.)
      Y05(X,Y) = Y * (63.*Y**4 - 70.*Y**2 + 15.)
      Y12(X,Y) =  Y            * SQRT(1.0 - Y * Y) * COS(ANGLE * X)
      Y14(X,Y) =  Y*(105.*Y*Y-45.) * SQRT(1.-Y*Y) * SIN(ANGLE*X)
      Y23(X,Y) =  Y * (1. - Y*Y) * (COS(ANGLE*X*2.) + SIN(ANGLE*X*2.))
      Y25(X,Y) =  Y * (1. - Y*Y) * (3.*Y*Y - 1.)  &
                                * (COS(ANGLE*X*2.) - SIN(ANGLE*X*2.))
      Y34(X,Y) =  Y * (1. - Y*Y) * SQRT(1. - Y*Y) * SIN(ANGLE*X*3.)
      Y45(X,Y) =  Y * (1. - Y*Y) ** 2
!C
!C     1.2 Fonctions symetriques (m+l pair)
!C
      Y00(X,Y) = 1.
      Y02(X,Y) = 1. - Y * Y * 3.
      Y11(X,Y) = 1.            * SQRT(1.0 - Y * Y) * COS(ANGLE * X)
      Y22(X,Y) = 1.            *     (1.0 - Y * Y) * COS(ANGLE * X * 2.)
!C
!C--CHAMP SPECTRAL MIS A ZERO
      DO 500 NC = 1, NCPLX * NBDEGLIB
        SPVOR(NC,1) = 0.0
  500 CONTINUE
!C--Champ de grille mis a zero
      DO 510 JC = 1, N2LONG * NLAT
        GTPP(JC,1) = 0.
  510 CONTINUE
!C
!C
!C--LA REPRESENTATION UTILISEE EN LATITUDE EST LA SUIVANTE
!C--   90 [ NORD..........0 [ NORD ; O [ SUD........90 [ SUD
!C--    NH LAT              1    ; NH LAT+1          N LAT
!C--     AFIN D'AVOIR DES INDICES DE BOUCLE // LORS DES CALCULS
!C
      DO 1500 JCLAT = 1, NHLAT
        DO 2500 JCLONG = 1, NLONG
          IF(MHEMIS.EQ.2) THEN
            VAL = Y00(FLOAT(JCLONG - 1),PTGAUSS(JCLAT))  &
               + Y02(FLOAT(JCLONG - 1),PTGAUSS(JCLAT))  &
               + Y11(FLOAT(JCLONG - 1),PTGAUSS(JCLAT))  &
               + Y22(FLOAT(JCLONG - 1),PTGAUSS(JCLAT))
!C
            SIGNE = + 1.0 * PARITE
            GTPP  (JCLONG,JCLAT         ) = VAL
            GTPP  (JCLONG,JCLAT + NHLAT) = VAL * SIGNE
          ENDIF
!C
!C   --ET MAINTENANT LA PARITE ANTISYM. SI NECESSAIRE A UN TEST COMPLET
!C
          VAL = Y01(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y12(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y03(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y05(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y14(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y23(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y25(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y34(FLOAT(JCLONG - 1),PTGAUSS(JCLAT)) &
             + Y45(FLOAT(JCLONG - 1),PTGAUSS(JCLAT))
!C
          SIGNE = - 1.0 * PARITE
            GTPP  (JCLONG,JCLAT + NHLAT) =  &
           GTPP  (JCLONG,JCLAT + NHLAT) + VAL
            GTPP  (JCLONG,JCLAT         ) = &
           GTPP  (JCLONG,JCLAT         ) + VAL * SIGNE
 2500   CONTINUE
 1500 CONTINUE
!C
!C--     PSEUDO LONGITUDES SUPPLEMENTAIRES NECESSAIRES A LA FFT
        DO 2000 JCLAT = 1, NHLAT
!C
          GTPP  (NLONG+1,JCLAT+NHLAT) = 0.0
          GTPP  (NLONG+1,JCLAT       ) = 0.0
          GTPP  (NLONG+2,JCLAT+NHLAT) = 0.0
          GTPP  (NLONG+2,JCLAT       ) = 0.0
 2000   CONTINUE
!C
!C       Recopie de GT PP dans GT VP
        DO 2100 JC = 1, N2LONG*NLAT
          GTVP(JC,1) = GTPP(JC,1)
 2100   CONTINUE
      RETURN
      END
!C******************************************************************************
!C                                                          SP GEN
!C******************************************************************************
!C
      SUBROUTINE SPGEN 
!C
!C     Generation d'un champs initial quelconque en spectral
!C
      INCLUDE 'sphectra.h'
!C
      GERME = 0.5
!C
      DO 100 J = 1, NI*NCPLX
        SPVOR(J,1) = GERME
        SPPSI(J,1) = GERME
        GERME = ALOG(ABS(GERME))
100   CONTINUE
!C
!C     elimination des parties complexes des modes zonaux
!C
      IS = MVADEB(0)
      DO 200 JC = 0, MVALG(0)-1
        SPVOR(2,IS+JC) = 0.
        SPPSI(2,IS+JC) = 0.
  200 CONTINUE
!C
      IF (MHEMIS.EQ.2) THEN
!C
      DO 104 J = NI*NCPLX+1, NBDEGLIB*NCPLX
        SPVOR(J,1) = 0.
        SPPSI(J,1) = 0.
  104 CONTINUE
      RETURN
!C
      ELSE
!C
!C  Generation de la partie symetrique
!C
      DO 105 J = NI*NCPLX+1, NBDEGLIB*NCPLX
        SPVOR(J,1) = GERME
        SPPSI(J,1) = GERME
        GERME = ALOG(ABS(GERME))
  105 CONTINUE
!C
!C  Elimination de la partie complexe des modes zonaux et du mode (0,0)
!C
      IS = MVSDEB(0)
      DO 205 JC = 0, MVSLG(0)-1
        SPVOR(2,IS+JC) = 0.
        SPPSI(2,IS+JC) = 0.
  205 CONTINUE
      SPVOR(1,NP) = 0.
      SPPSI(1,NP) = 0.
      ENDIF
!C
      RETURN
      END
!C
!C*****************************************************************************
!C                                                        SP GE J
!C*****************************************************************************
!C
      SUBROUTINE SPGEJ 
!C
!C     Generation du champ de fonction de courant pour le test comparatif
!C     du jacobien avec la version de reference sur Cyber-CNES.
!C
      INCLUDE 'sphectra.h'
!C
!C     Partie antisymetrique
!     Mode m=0
      IS = MVADEB(0)
      IN = MVALG (0)
      DO 100 JC = 0, IN-1
        SPPSI(1,IS+JC) = 0.01 * (2*JC+1)
        SPPSI(2,IS+JC) = 0.
  100 CONTINUE
!C
!C     Modes m>0
      DO 300 JMOD = 1, NTRUNC
        IS = MVADEB(JMOD)
        IN = MVALG (JMOD)
        IF(IN.NE.0) THEN
        DO 310 JC = 0, IN - 1
          SPPSI(1,IS+JC) = 0.001 * (JMOD +(2*JC+1))
          SPPSI(2,IS+JC) = -0.003 * JMOD
  310   CONTINUE
        ENDIF
  300 CONTINUE
!C
      IF (MHEMIS.EQ.2) THEN
!C
!C     Partie symetrique
!C     Mode m=0
      IS = MVSDEB(0)
      IN = MVSLG (0)
      DO 110 JC = 0, IN-1
        SPPSI(1,IS+JC) = 0.01 * (2*JC)
        SPPSI(2,IS+JC) = 0.
  110 CONTINUE
!C
!C     Modes m>0
      DO 305 JMOD = 1, NTRUNC
        IS = MVSDEB(JMOD)
        IN = MVSLG (JMOD)
        DO 315 JC = 0, IN - 1
          SPPSI(1,IS+JC) = 0.001 * (JMOD +(2*JC))
          SPPSI(2,IS+JC) = -0.003 * JMOD
  315   CONTINUE
  305 CONTINUE
!C
      ELSE
!C     Elimination de la partie symetrique
!C
      DO 400 JC = NI+1,NI+NP
        SPPSI(1,JC) = 0.
        SPPSI(2,JC) = 0.
        SPJAC(1,JC) = 0.
        SPJAC(2,JC) = 0.
  400 CONTINUE
      ENDIF
!C
      RETURN
      END
!C
!C******************************************************************************
!C                                                             DUPLIC
!C******************************************************************************
!C
      SUBROUTINE DUPLIC(FCHAMP,ISYM) 
!C
!C     sous programme de duplication d'un champ de semi Fourier
!C      hemispherique pour usage avec des fonctions spheriques
!!C     ISYM : +1 ou -1 selon que le champ est symetrique ou antisymetrique
!C
      INCLUDE 'sphectra.h'
!C
      REAL FCHAMP(NCPLX,NH1LONG,NLAT)
!C
      IF (ISYM.EQ.1) THEN
        DO 100 JMOD = 1, NH1LONG
        DO 100 JC = 1, NHLAT
          FCHAMP(1,JMOD,JC+NHLAT) = FCHAMP(1,JMOD,JC)
          FCHAMP(2,JMOD,JC+NHLAT) = FCHAMP(2,JMOD,JC)
  100   CONTINUE
      ENDIF
      IF (ISYM.EQ.-1) THEN
        DO 200 JMOD = 1, NH1LONG
        DO 200 JC = 1, NHLAT
          FCHAMP(1,JMOD,JC+NHLAT) = - FCHAMP(1,JMOD,JC)
          FCHAMP(2,JMOD,JC+NHLAT) = - FCHAMP(2,JMOD,JC)
  200   CONTINUE
      ENDIF
!C
      RETURN
      END
!C****************************************************************************
!C                                                          RESOL
!C****************************************************************************
      SUBROUTINE RESOL(U,V,P,K,D,A,B,M,N) 
      include 'sphectra.h'
!C     --------------------------------------------
!C     SI HM=M EN NOTATION NORMALE
!C     M=HM+1
!C     U(1:N)=HU(HM:HMAX)
!C     D(L)=HD(HM,HN) OU HN=M+L-2
!C     SOIT D(L)=HD(HM,HM+L-1)
!C     --------------------------------------------
!C     SI ON SUPPOSE U(0:HM),V(0:HM)  CONNUS
!C     LES RELATIONS A NB D'ONDE ZONAL FIXE HM
!C          UMN= (N-1)DMN PMN-1 + IM KMN -(N-2)DMN+1 PMN+1
!C          VMN=-(N-1)DMN KMN-1 + IM PMN + (N-2)DMN+1 KN+1
!C     FOURNISSENT LE BON NB D'EQ POUR PSI ET KHI
!C     SI ON SUPPOSE PMN+1 ET KMN+1 NULS POUR N=NTRUNC
!C     LA RESOLUTION S'EFFECTUE ALORS PAR DOUBLE BALAYAGE
!C     LA MATRICE DU PB EST TRIDIAGONALE PAR BLOCS DE DIM 2
!C     1.ON LA REND DIAG SUP
!C     2.ELIMINATION
!C     --------------------------------------------------
      COMPLEX U(N),V(N),P(N+1),K(N+1),A(2,2,0:N),B(2,0:N)               
      COMPLEX MAT(2,2),DET,XX,YY,HI                                     
      DIMENSION D(N+1)                                                  
      D(1)=0.                                                           
      DO 1 L=1,N+1                                                      
1     D(L)=SQRT(FLOAT((M+L-2)**2-(M-1)**2)/FLOAT(4*(M+L-2)**2-1))       
      IF(M.GT.1) GO TO 2                                                
!c      PRINT*,(D(I),I=1,5)
      K(1)=0.                                                           
      P(1)=0.                                                           
      K(2)=V(1)/D(2)/2.                                                 
      P(2)=-U(1)/D(2)/2.                                                
      DO 3 L=3,N                                                        
      K(L)=( V(L-1)+FLOAT(L-3)*D(L-1)*K(L-2))/D(L)/FLOAT(L)             
3     P(L)=(-U(L-1)+FLOAT(L-3)*D(L-1)*P(L-2))/D(L)/FLOAT(L)             
      RETURN                                                            
2     HI=(0.,1.)                                                        
      P(N+1)=(0.,0.)                                                    
      K(N+1)=(0.,0.)                                                    
      Z=FLOAT(M-1)                                                      
      DO 4 I=1,2                                                        
      B(I,0)=(0.,0.)                                                    
      DO 4 J=1,2                                                        
4      A(I,J,0)=(0.,0.)                                                 
      DO 5 L=1,N                                                        
      X=FLOAT(L+M)*D(L+1)                                               
      Y=FLOAT(L+M-3)*D(L)                                               
      XX=U(L)-Y*B(1,L-1)                                                
      YY=V(L)+Y*B(2,L-1)                                                
      MAT(1,1)= A(1,1,L-1)*Y                                            
      MAT(2,1)=-A(2,1,L-1)*Y+HI*Z                                       
      MAT(1,2)= A(1,2,L-1)*Y+HI*Z                                       
      MAT(2,2)=-A(2,2,L-1)*Y                                            
      DET=1./(MAT(1,1)*MAT(2,2)-MAT(1,2)*MAT(2,1))                      
      A(1,1,L)= X*DET* MAT(2,2)                                         
      A(2,1,L)=-X*DET* MAT(2,1)                                         
      A(1,2,L)= X*DET*MAT(1,2)                                          
      A(2,2,L)=-X*DET* MAT(1,1)                                         
      B(1,L)=(XX*MAT(2,2)-MAT(1,2)*YY)*DET                              
      B(2,L)=(MAT(1,1)*YY-XX*MAT(2,1))*DET                              
5       CONTINUE                                                        
      DO 6 L=N,1,-1                                                     
      P(L)=A(1,1,L)*P(L+1)+A(1,2,L)*K(L+1)+B(1,L)                       
      K(L)=A(2,1,L)*P(L+1)+A(2,2,L)*K(L+1)+B(2,L)                       
6      CONTINUE                                                         
      RETURN                                                            
      END 
!c=========================================================================
      SUBROUTINE truncspec (SPOUT,SPIN) 
      INCLUDE 'sphectra.h'
      REAL  SPIN (NCPLX,NBDEGLIB)
      REAL  SPOUT(NCPLX,NBDEGLIB)    
       DO 110 JC = 1, NI + (MHEMIS-1) * NP
         SPOUT(1,JC) = SPIN(1,JC) * OPRtrunc(JC)
         SPOUT(2,JC) = SPIN(2,JC) * OPRtrunc(JC)
  110  CONTINUE

      RETURN
      END
!c============================================================
           SUBROUTINE TRUNC(m1) 
           INCLUDE 'sphectra.h'
           DO 3000 JMOD = 0, NTRUNC
              IS = MVSDEB(JMOD)
              IN = MVSLG (JMOD)
              JD = 0
              IF (JMOD.EQ.0) JD = 1
              DO 3010 JN = JD, IN - 1
                 LN = JMOD + 2 * JN
                 if(jmod.le.m1.and.ln.le.m1)then
                 OPRtrunc   (IS + JN) = 1.
                 else
                  OPRtrunc   (IS + JN) = 0. 
                 endif
 3010         CONTINUE
        IS = MVADEB(JMOD)
        IN = MVALG (JMOD)
        IF(IN.EQ.0) GOTO 3025
        JD = 0
        DO 3020 JN = JD, IN - 1
           LN = JMOD + 2 * JN + 1
           if(jmod.le.m1.and.ln.le.m1)then
              OPRtrunc  (IS + JN) = 1.
           else
              OPRtrunc   (IS + JN) = 0.
           endif
 3020   CONTINUE
 3025   CONTINUE
 3000 CONTINUE
      OPRtrunc (NI+1) = 1.
      RETURN
      end
!c********************************************
!C*****************************************************************************
!C                                                                 CALVENT
!C*****************************************************************************
                                                                 
      SUBROUTINE CALVENT(U,V,PSI,KHI,MTRUNC)
      INCLUDE 'sphectra.h'
      COMPLEX U(NI+NP),V(NI+NP),PSI(NI+NP),KHI(NI+NP)
      REAL D(0:NTRUNC+1)
      COMPLEX UM(0:NTRUNC),VM(0:NTRUNC),PM(0:NTRUNC+1),  &
          KM(0:NTRUNC+1)

!C    ----------------------------------------------
!C     CALCULE UCOSTETA ET VCOSTETA PAR
!C     UMN=(N-1)DMNPSIMN-1 - (N-2)DMN+1PSIMN+1 + IMKHIMN
!C     VMN=-(N-1)DMNKHIMN-1 + (N+2)DMN+1KHIMN+1 + IMPSIMN
!C     LE CALCUL DE U A LA TRONCATURE MTRUNC NECESSITE PSI EN MTRUNC+1
!C     DONC PAS D'APPEL AVEC MTRUNC.GE.NTRUNC-1
!C     --------------------------------------------------
      IF(MTRUNC.GE.NTRUNC)THEN
         WRITE(LW,*)'APPEL CALVENT AVEC MTRUNC.GE.NTRUNC'
      STOP
      ENDIF
      DO 10 I=1,NI+NP
      U(I)=0.
 10   V(I)=0.
      DO 20 I=0,NTRUNC+1
      PM(I)=0.
      KM(I)=0.
 20   D(I)=0.
!C
!C     BOUCLE SUR M
!C
      DO 100 M=0,MTRUNC
!C       RECOPIE DANS PM ET KM
!C
      DO 110 N=M,MTRUNC
      IF(MOD(N-M,2).EQ.0)THEN
         PM(N)=PSI(MVSDEB(M)+(N-M)/2)/RAYTER
         KM(N)=KHI(MVSDEB(M)+(N-M)/2)/RAYTER
      ELSE
         PM(N)=PSI(MVADEB(M)+(N-M-1)/2)/RAYTER
         KM(N)=KHI(MVADEB(M)+(N-M-1)/2)/RAYTER
      ENDIF
      
 110  CONTINUE
     
      PM(MTRUNC+1)=0.
      KM(MTRUNC+1)=0.
!C       CALCUL DES DM
      D(M)=0.
      DO 120 N=M+1,MTRUNC+1
 120  D(N)=SQRT((N*N-M*M)/(4.*N*N-1.))
!C
!C     CALCUL DES UM ET VM
      DO 130 N=M,MTRUNC
      UM(N)=(N-1)*D(N)*PM(N-1)-(N+2)*D(N+1)*PM(N+1) &
          +(0.,1.)*M*KM(N)
 130  VM(N)=-(N-1)*D(N)*KM(N-1)+(N+2)*D(N+1)*KM(N+1) &
          +(0.,1.)*M*PM(N)
!C
!C     RECOPIE DES UM ET VM DANS U ET V
      DO 140 N=M,MTRUNC
      IF(MOD(N-M,2).EQ.0)THEN
         U(MVSDEB(M)+(N-M)/2)=UM(N)
         V(MVSDEB(M)+(N-M)/2)=VM(N)
      ELSE
         U(MVADEB(M)+(N-M-1)/2)=UM(N)
         V(MVADEB(M)+(N-M-1)/2)=VM(N)
      ENDIF
    
 140  CONTINUE
!C
!C     FIN DE LA BOUCLE EN M
 100  CONTINUE
      RETURN
      END
!C
!C
!C *****************************************************************************
!C                                                        TLG INV SH
!C *****************************************************************************
!C
      SUBROUTINE TLGINVSH (FCHAMP,SPC,LEGPOL,PARITE)
!C
!C**** TLG INV H  Transformee de Legendre inverse hemispherique symetrique
!C
!C     B. Legras        LMD      Creation: 18-05-84
!C                               Certification: 22-05-84
!C                               Performance: 0.2313 ms pour N TRUNC = 21
!C
!C     Objet:  Transformee de Legendre inverse de la representation
!C             spectrale a la representation en semi-Fourier.
!C             Transformee pour un champ hemispherique symetrique
!C             (c.a.d. ayant des coefficients tq m+l pair).
!C             Description verticale de la troncature parametrisee
!C             par le common TRCVERT.
!C             En standard: troncature triangulaire homogene.
!C
!C     Arguments:
!C             JR NIV:  Nombre de niveaux.
!C                      Cette variable n'est pas utilise; elle est fournie
!C                      pour compatibilite avec les autres versions.
!C             LEGPOL:  Tableau des fonctions de Legendre associees ou de
!C                      leur pseudo-derivees.
!C             PARITE:  Parite de la fonction de transformation contenue
!C                      dans LEGPOL.
!C                      Non pertinent ici.
!C             SPC   :  Tableau des coefficients d'harmoniques spheriques.
!C                      En entree.
!C             F CHAMP: Tableau des coefficients de semi-Fourier.
!C                      En sortie.
!C                      Rangement: l'hemisphere Sud d'abord et dans chaque
!C                      hemisphere, de l'Equateur vers le Pole.
!C
!C     Methode: Les sommations en l pour les colonnes de la troncature
!C             sont regroupees et calculees ensemble pour les differentes
!C             latitudes.
!C             Les contributions symetriques et antisymetriques sont
!C             calculees separement et combinees en sortie.
!C
!C     Externes: MXMA
!C
!C ****************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL    LEGPOL   (NBDEGLIB,NHLAT)
      REAL    SPC      (NCPLX,NBDEGLIB)
      REAL    FCHAMP  (NCPLX,0:NHLONG,NHLAT)
      INTEGER PARITE
!C
!C
!C*    1. Definition des increments de MXMA
!C     ------------------------------------
!C
      DATA      INC1P, INC2P,   INC1S, INC2S,    INC1F,  INC2F &
         / NBDEGLIB,     1,   NC1S,     1,  N2LONG,      1/
!C
!C*    2. Calcul matriciel de la sommation sur le degre
!C     ------------------------------------------------
!C
  200 CONTINUE
      DO 250 JMOD = 0, NTRUNC
!C
!C*    2.1 Definition des parametres de boucle
  210   ISS = MVSDEB(JMOD)
        INS = MVSLG (JMOD)
        IK = NCPLX
        IF(JMOD.EQ.0) IK = 1
!C
!C*    2.3 Transformee symetrique
  230   CALL MXMA(       &
              LEGPOL  (   ISS ,       1), INC1P, INC2P,&
              SPC     (1, ISS          ), INC1S, INC2S,&
              FCHAMP (1,       JMOD, 1), INC1F, INC2F, &
              NHLAT, INS, IK)
!C
  250 CONTINUE
!C
!C*    4.  Mise a zero des modes de Fourier non representes
!C     ----------------------------------------------------
!C
  400 CONTINUE
      DO 405 JC = 1, NHLAT
        FCHAMP(2,0,JC) = 0.
  405 CONTINUE
      DO 410 JMOD = NTRUNC+1, NHLONG
      DO 410 JC = 1, NHLAT
        FCHAMP(1,JMOD,JC) = 0.
        FCHAMP(2,JMOD,JC) = 0.
  410 CONTINUE
!C
!C
      RETURN
      END
!C
!C *****************************************************************************
!C
      SUBROUTINE TLGINVH (FCHAMP,SPC,LEGPOL,PARITE)
!C
!C**** TLG INV H  Transformee de Legendre inverse hemispherique antisymetrique
!C
!C     B. Legras        LMD      Creation: 18-05-84
!C                               Certification: 22-05-84
!C                               Performance: 0.2313 ms pour N TRUNC = 21
!C
!C     Objet:  Transformee de Legendre inverse de la representation
!C             spectrale a la representation en semi-Fourier.
!C             Transformee pour un champ hemispherique antisymetrique
!C             (c.a.d. ayant des coefficients tq m+l impair).
!C             Description verticale de la troncature parametrisee
!C             par le common TRCVERT.
!C             En standard: troncature triangulaire homogene.
!C
!C     Arguments:
!C             JR NIV:  Nombre de niveaux.
!C                      Cette variable n'est pas utilise; elle est fournie
!C                      pour compatibilite avec les autres versions.
!!C             LEGPOL:  Tableau des fonctions de Legendre associees ou de
!C                      leur pseudo-derivees.
!C             PARITE:  Parite de la fonction de transformation contenue
!C                      dans LEGPOL.
!C                      Non pertinent ici.
!C             SPC   :  Tableau des coefficients d'harmoniques spheriques.
!C                      En entree.
!C             F CHAMP: Tableau des coefficients de semi-Fourier.
!C                      En sortie.
!C                      Rangement: l'hemisphere Sud d'abord et dans chaque
!C                      hemisphere, de l'Equateur vers le Pole.
!C
!C     Methode: Les sommations en l pour les colonnes de la troncature
!C             sont regroupees et calculees ensemble pour les differentes
!C             latitudes.
!C             Les contributions symetriques et antisymetriques sont
!C             calculees separement et combinees en sortie.
!C
!C     Externes: MXMA
!C
!C ****************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL    LEGPOL   (NBDEGLIB,NHLAT)
      REAL    SPC      (NCPLX,NBDEGLIB)
      REAL    FCHAMP  (NCPLX,0:NHLONG,NHLAT)
      INTEGER PARITE
!C
!C
!C*    1. Definition des increments de MXMA
!C     ------------------------------------
!C
      DATA      INC1P, INC2P,   INC1S, INC2S,    INC1F,  INC2F  &
         / NBDEGLIB,     1,   NC1S,     1,  N2LONG,      1/
!C
!C*    2. Calcul matriciel de la sommation sur le degre
!C     ------------------------------------------------
!C
  200 CONTINUE
      DO 250 JMOD = 0, NTRUNC - 1
!C
!C*    2.1 Definition des parametres de boucle
  210   ISA = MVADEB(JMOD)
        INA = MVALG (JMOD)
        IK = NCPLX
        IF(JMOD.EQ.0) IK = 1
!C
!C*    2.3 Transformee antisymetrique
  230   CALL MXMA(        &
              LEGPOL  (   ISA ,       1), INC1P, INC2P,  &
              SPC     (1, ISA          ), INC1S, INC2S,  &
              FCHAMP (1,       JMOD, 1), INC1F, INC2F,   &
              NHLAT, INA, IK)
!C
  250 CONTINUE
!C
!C*    4.  Mise a zero des modes de Fourier non representes
!C     ----------------------------------------------------
!C
  400 CONTINUE
      DO 405 JC = 1, NHLAT
        FCHAMP(2,0,JC) = 0.
  405 CONTINUE
      DO 410 JMOD = NTRUNC, NHLONG
      DO 410 JC = 1, NHLAT
        FCHAMP(1,JMOD,JC) = 0.
        FCHAMP(2,JMOD,JC) = 0.
  410 CONTINUE
!C
!C
      RETURN
      END
!C *****************************************************************************
!C                                                           TLG DIR SH
!C *****************************************************************************
!C
      SUBROUTINE TLGDIRSH (SPC,FCHAMP,LEGPOL,PARITE)
!C
!C**** TLG DIR H  Transformee de Legendre directe hemispherique symetrique
!C
!C     B. Legras        LMD      creation: 17-05-84
!C                               certification: 22-05-84
!C                               performance: 0.283 ms pour N TRUNC = 21
!C                                            et M MAX = 19
!C
!C     Objet:  Transformee de Legendre directe de la representation
!C             semi-Fourier a la representation spectrale.
!C             Transformee hemispherique pour un champ final symetrique.
!C             Description verticale de la troncature.
!C             Troncature parametrisable par le COMMON TRCVERT.
!C             En standard: troncature triangulaire homogene.
!C
!C     Arguments:
!C             LEGPOL: tableau des fonctions de Legendre associees ou de
!C                     leurs pseudo-derivees
!C             PARITE: parite de la fonction de transformation contenue
!C                     dans LEGPOL
!C             F CHAMP:tableau des coefficients de semi-Fourier en entree.
!C                     rangement: ce tableau ne contient qu'un hemisphere
!C                                (l'hemisphere Nord) range de l'Equateur
!C                                vers le pole
!C                                NB: en cas d'equivalence avec un tableau
!C                                spherique complet, tenir compte du fait
!C                                que l'hemisphere Sud est stocke dans la
!C                                premiere moitie de ce dernier.
!C             SPC   : tableau des coefficients d'harmoniques spheriques
!C                     en sortie
!C
!C     Methode:Les sommations en latitude sont regroupees et calculees
!C             par MXMA pour les premieres colonnes de la troncature.
!C             Les dernieres colonnes sont calculees par SDOT.
!C
!C     Externes: MXMA
!C               SDOT
!C
!C     Precautions d'emploi: hormis celles signales plus haut,
!C                           limiter M MAX a N TRUNC car il n'y
!C                           a pas de calcul pour JMOD = N TRUNC
!C
!C ****************************************************************************
!C
      INCLUDE 'sphectra.h'
!C
      REAL    LEGPOL  (NBDEGLIB,NHLAT)
      REAL    SPC     (NCPLX,NBDEGLIB)
      REAL    FCHAMP (NCPLX,NH1LONG,NHLAT)
      INTEGER PARITE
!C
!C
!C*    1.  Initialisation des increments pour MXMA
!C     -------------------------------------------
!C
      DATA   INC1P,    INC2P,    INC1F,  INC2F,    INC1S,  INC2S  &
         /      1,NBDEGLIB,  N2LONG,      1,    NC1S,      1/
!C
!C*     2. Contribution totale: multiplication par 2 en entree
!C      ------------------------------------------------------
!C
  200 CONTINUE
      DO 210 JC = 1, N2LONG * NHLAT
        ZEFA(JC) = 2 * FCHAMP(JC,1,1)
  210 CONTINUE
!C
!C*    3.  Calcul matriciel pour les ordres jusqu'a M MAX-1
!C     ----------------------------------------------------
!C
  300 CONTINUE
      DO 350 JMOD = 0, MMAX - 1
!C
!C*    3.1 Definition des parametres de boucle
  310 ISS = MVSDEB(JMOD)
      INS = MVSLG (JMOD)
      IK = NCPLX
      IF(JMOD.EQ.0) IK = 1
!C
!C*    3.3 Transformee symetrique
  330   CALL MXMA(     &
              LEGPOL  (   ISS ,      1), INC1P, INC2P,  &
              ZFA    (1,      JMOD, 1), INC1F, INC2F,   &
              SPC     (1, ISS         ), INC1S, INC2S,  &
              INS, NHLAT, IK)
!C
  350 CONTINUE
!C
!C*    4.  Calcul scalaire de l'ordre M MAX a N TRUNC
!C     ----------------------------------------------
!C
  400 DO 450 JMOD = MMAX, NTRUNC
        ISS = MVSDEB(JMOD)
        INS = MVSLG (JMOD)
!C
!C*    4.2 Partie antisymetrique
        DO 450 JN = 0, INS - 1
          SPC( 1, ISS+JN ) = SDOT(NHLAT, LEGPOL(ISS+JN , 1), INC2P,  &
                                         ZFA  ( 1, JMOD,1), INC1F)
          SPC( 2, ISS+JN ) = SDOT(NHLAT, LEGPOL(ISS+JN , 1), INC2P,  &
                                         ZFA  ( 2, JMOD,1), INC1F)
  450 CONTINUE
!C
!C
      RETURN
      END
!C
!C***************************************************************************
!C                                                        CALPSI
!C***************************************************************************
      SUBROUTINE CALPSI(PSI,KHI,U,V,MTRUNC)
      include 'sphectra.h'

!C     --------------------------------------------------
!C     DONNE UCOSPHI ET VCOSPHI ET RESSORT
!C     PSI ET KHI.
!C     ---------------------------------------------------
!C     LTF ,LTM NB D'ONDES ZONAUX+1
!C     LTN NB D'ONDE TOTAL+1

      COMPLEX U(NI+NP),V(NI+NP),PSI(NI+NP),KHI(NI+NP)
      COMPLEX UM(0:NTRUNC),VM(0:NTRUNC),KM(0:NTRUNC+1),PM(0:NTRUNC+1)        
      COMPLEX A(2,2, 0:NTRUNC+1),B(2,0:NTRUNC+1)                              
      DIMENSION D(0:NTRUNC+1)                                    
!CC     BOUCLE SUR M:ON PRESENTE LES COMPO DE U ET V A M DONNE
!C      AVEC N CROISSANT
      DO 5 I=1,NI+NP
      KHI(I)=0.
5     PSI(I)=0.
!C      PRINT*,'NIP=',NIP
!C      PRINT*,'NSI',NSI
!C      PRINT*,'MVADEB=',MVADEB
!C      CALL PRICHNS(U,'CHAMP U')
      DO 10 M=0,MTRUNC
      DO 11 N=M,MTRUNC
      IF(MOD(N-M,2).EQ.0)THEN
      UM(N-M)=U(MVSDEB(M)+(N-M)/2)*rayter
      VM(N-M)=V(MVSDEB(M)+(N-M)/2)*rayter
      ELSE
      UM(N-M)=U(MVADEB(M)+(N-M-1)/2)*rayter
      VM(N-M)=V(MVADEB(M)+(N-M-1)/2)*rayter
      ENDIF
11    CONTINUE
      NMAX=MTRUNC+1-M
!c      IF(M.LE.2)THEN
!c      PRINT*,'M=',M
!c      PRINT900,UM
!c900   FORMAT('UM=',(4('(',E10.3,'.',E10.3,')')))
!c      ENDIF
      CALL RESOL(UM,VM,PM,KM,D,A,B,M+1,NMAX)                              
!c      IF(M.LE.2)THEN
!c      PRINT*,'M=',M
!c      PRINT901,PM
!c901   FORMAT('PM=',(4('(',E10.3,'.',E10.3,')')))
!c      ENDIF
!C
!C     REMISE EN PLACE
      DO 12 N=M,MTRUNC
      IF(MOD(N-M,2).EQ.0)THEN
      PSI(MVSDEB(M)+(N-M)/2)=PM(N-M)
      KHI(MVSDEB(M)+(N-M)/2)=KM(N-M)
      ELSE
      PSI(MVADEB(M)+(N-M-1)/2)=PM(N-M)
      KHI(MVADEB(M)+(N-M-1)/2)=KM(N-M)
      ENDIF
12    CONTINUE
10      CONTINUE                                                        
!C      CALL PRICHNS(PSI,'COURANT')
      RETURN
        END 
!C*****************************************************************************
!C                                                                 gausstrig
!C*****************************************************************************
        subroutine gausstrig
        include 'sphectra.h'  
!c     Donne les sinus, cosinus et latitude sur la grille de gauss.
        do i=1,nhlat
           gsin(i)=ptgauss(i)
           gcos(i)=sqrt(1.0-gsin(i)**2)
           theta(i)=atan(gsin(i)/gcos(i))
           gsin(i+nhlat)=-gsin(i)
           gcos(i+nhlat)=gcos(i)
           theta(i+nhlat)=-theta(i)
        enddo
        return 
        end
!C*****************************************************************************
!C                                                                 VPASSM
!C*****************************************************************************
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
!C
!C     SUBROUTINE 'VPASSM' - MULTIPLE VERSION OF 'VPASSA'
!C     PERFORMS ONE PASS THROUGH DATA
!C     AS PART OF MULTIPLE COMPLEX FFT ROUTINE
!C     A IS FIRST REAL INPUT VECTOR
!C     B IS FIRST IMAGINARY INPUT VECTOR
!C     C IS FIRST REAL OUTPUT VECTOR
!C     D IS FIRST IMAGINARY OUTPUT VECTOR
!C     TRIGS IS PRECALCULATED TABLE OF SINES ' COSINES
!C     INC1 IS ADDRESSING INCREMENT FOR A AND B
!C     INC2 IS ADDRESSING INCREMENT FOR C AND D
!C     INC3 IS ADDRESSING INCREMENT BETWEEN A'S & B'S
!C     INC4 IS ADDRESSING INCREMENT BETWEEN C'S & D'S
!C     LOT IS THE NUMBER OF VECTORS
!C     N IS LENGTH OF VECTORS
!C     IFAC IS CURRENT FACTOR OF N
!C     LA IS PRODUCT OF PREVIOUS FACTORS
!C
!C     SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
!C
!C     END
      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/, &
          SIN72/0.951056516295154/,COS72/0.309016994374947/, &
          SIN60/0.866025403784437/
!C
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
!C     CHECK FACTORS ARE CORRECT - ENSURE NON-NEGATIVE
      IF (IGO.LE.0) GOTO 998
      IF (IGO.GT.4) GO TO 999
      GO TO (10,50,90,130),IGO
!C
!C     CODING FOR FACTOR 2
!C
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
!C
!C     CODING FOR FACTOR 3
!C
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=    &
         C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))&
        -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=    &
         S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))&
        +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=    &
         C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))&
        -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=    &
         S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))&
        +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
!C
!C     CODING FOR FACTOR 4
!C
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=      &
         C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))&
        -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=      &
         S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))&
        +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=      &
         C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))&
        -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=      &
         S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))&
        +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=      &
         C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))&
        -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=      &
         S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))&
        +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
!C
!C     CODING FOR FACTOR 5
!C
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))&
       -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))&
       +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))&
       +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))&
       -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))&
       -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))&
       +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))&
       +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))&
       -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=       &
         C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))&
           -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))        &
        -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))&
           +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=       &
         S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))&
           -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))        &
        +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))&
           +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=       &
         C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))&
           +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))        &
        -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))&
           -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=       &
         S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))&
           +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))        &
        +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))&
           -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=      &
         C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))&
           -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))        &
        -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))&
           +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=      &
         S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))&
           -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))        &
        +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))&
           +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=      &
         C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))&
           +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))        &
        -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))&
           -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=      &
         S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))&
           +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))        &
        +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))&
           -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
!C  ** ERROR - FACTOR LESS THAN 1  NOT ALLOWED **
998   print *,' VPASSM: FACTORS ARE INCORRECT '
      CALL ABORT
!C  ** ERROR - FACTOR HIGHER THAN 5 NOT ALLOWED **
999   print *,' VPASSM: FACTORS HIGHER THAN 5 ARE NOT SUPPORTED '
      CALL ABORT
      END
      SUBROUTINE FAX(IFAX,N,MODE)
!C***********************************************************************
!C                                                                      *
!C C06-SUMMATIOM
!C C06-SUMMATION OF SERIES                                      B6.1/3  *
!C                                                                      *
!C                                                              FFTRIG  *
!C                                                              FAX     *
!C                                                                      *
!C                                                                      *
!C SUP
!C SUBPROGRAM       SUBROUTINE   FFTRIG                                 *
!C                               FAX                                    *
!C                                                                      *
!C PURPOSE          SETUP ROUTINES FOR FFT PACKAGES                     *
!C                                                                      *
!C                                                                      *
!C VERSION          CYBER                         CRAY-1                *
!C                                                                      *
!C                  JAN 1979 ORIGINAL             JAN 1979 ORIGINAL     *
!C                                                                      *
!C USAGE                                                                *
!C                  CALL FFTRIG(TRIGS,N,3)                              *
!C                  CALL FAX   (IFAX ,N,3)                              *
!C                                                                      *
!C ARGUMENTS        1.DIMENSION                                         *
!C                       TRIGS(DIMENSION 3*N/2 - ADD 1 IF N/2 IS ODD)   *
!C                       IFAX(10)                                       *
!C                                                                      *
!C                  2.INPUT                                             *
!C                      N - THE LENGHT OF THE TRANSFORMS TO BE PERFORMED*
!C                          N MUST BE EVEN.                             *
!C                          THE NUMBER OF WORDS OF IFAX USED INCREASES  *
!C                          LOGARITHMICALLY WITH N.                     *
!C                          IFAX(10) SUFFICES FOR PRACTICAL PURPOSES.   *
!C                          (TRANSFORMS OF LENGHT AT LEAST 10000)       *
!C                                                                      *
!C                  3.OUTPUT                                            *
!C                      TRIGS - FFTRIG RETURNS AN ARRAY OF TRIGONOMETRIC*
!C                              FUNCTION VALUES SUBSEQUENTLY USED BY    *
!C                              FFT ROUTINES.                           *
!C                      IFAX  - FAX FACTORIZES N/2 INTO A PRODUCT OF    *
!C                              4"S AND 2"S AND HIGHER PRIME NUMBERS.   *
!C                              IFAX(1) CONTAINS THE NUMBER OF FACTORS. *
!C                              AND THE FACTORS THEMSELVES ARE STORED   *
!C                              IN ASCENDING ORDER IN IFAX(2),IFAX(3).. *
!C                              IF FAX IS CALLED WITH N ODD ,IFAX(1)    *
!C                              IS SET TO -99(ERROR CONDITION) AND NO   *
!C                              FACTORIZATION IS DONE.                  *
!C                                                                      *
!C WRITE UP         NONE                                                *
!C                                                                      *
!C ENTRY POINTS           FFTRIG,  FAX                                  *
!C                                                                      *
!C COMMON BLOCKS    NONE                                                *
!C I/O              NONE                                                *
!C PRECISION        SINGLE                                              *
!C OTHER ROUTINES   NONE                                                *
!C       REQUIRED                                                       *
!C 7/80                     FFTRIG-1                                    *
!C                                                                      *
!C***********************************************************************
!C                                                                      *
!C CO6-SUMMATION OF SERIES                                       B6.1/3 *
!C                                                                      *
!C                                                              FFTRIG  *
!C                                                              FAX     *
!C                                                                      *
!C ACSSES (OBJECT)  CYBER:                                              *
!C                           ATTACH,ECLIB.                              *
!C                           LDSET(LIB=ECLIB)                           *
!C                  CRAY 1:                                             *
!C                           LDR(LIB=ECLIB...)                          *
!C                                                                      *
!C ACCESS (SOURCE)           ATTACH,OLDPL,ECLIBPL                       *
!C                                                                      *
!C                  CYBER :         %DEFINE CYBER                       *
!C                  CRAY:           %DEFINE CRAY                        *
!C                                  %C   FFTRIG,   FAX                  *
!C                                                                      *
!C LANGUAGE         FORTRAN                                             *
!C                                                                      *
!C SPECIALIST       CLIVE TEMPERTON                                     *
!C                                                                      *
!C HISTORY          WRITTEN BY C.TEMPERTON      JAN     1979            *
!C                                                                      *
!C ALGORITHM                                                            *
!C REFERENCES                                                           *
!C                                                                      *
!C OBJECT SIZE               FFTRIG  FAX  (OCTAL WORDS)                 *
!C                  CYBER:     145   127                                *
!C                  CRAY :     221   157                                *
!C                                                                      *
!C                                                                      *
!C ACCURACY                                                             *
!C                                                                      *
!C TIMING                                                               *
!C                                                                      *
!C PORTABILITY      STANDARD FORTRAN                                    *
!C                                                                      *
!C SYSTEM ROUTINES  NONE                                                *
!C        REQUIRED                                                      *
!C
!C 7/80                      FFTRIG-2                                   *
!C
!C***********************************************************************
!!C     END
      DIMENSION IFAX(10)
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N)  GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
!C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
!C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
!C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40

!C     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
!C     INC ALTERNATIVELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
!C     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
!C     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
!C FROM UPDATE LIBRARY         LMDBIB                CY=2       02/08/85
      SUBROUTINE FFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!************************************************************************
!*                                                                      *
!* C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                      *
!*                                                               FFT99  *
!*                                                               FFT991 *
!*                                                                      *
!*                                                                      *
!* SUBPROGRAM       SUBROUTINE    FFT99                                 *
!*                                FFT991                                *
!*                                                                      *
!* PURPOSE          PERFORM MULTIPLE FAST FOURIER TRANSFORMS            *
!*                                                                      *
!*                                                                      *
!* VERSION          CYBER                         CRAY-1                *
!*                                                                      *
!*                  JAN 1979 ORIGINAL             JAN 1979 ORIGINAL     *
!*                                                                      *
!* USAGE                                                                *
!*                  CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   *
!*                  CALL FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   *
!*                                                                      *
!* ARGUMENTS        1.DIMENSION                                         *
!*                       A(IDIM),WORK((N+1)*M),TRIGS(3*N/2),IFAX(10)    *
!*                       WORK IS A WORK ARRAY                           *
!*                                                                      *
!*                  2.INPUT                                             *
!*                      A - AN ARRAY CONTAINING THE INPUT DATA OR       *
!*                          COEFFICIENT VECTORS.                        *
!*                         THIS ARRAY IS OVERWRITTEN BY THE RESULTS.    *
!*                      TRIGS AND IFAX - ARRAYS SET UP BY FFTRIG AND FAX*
!*                                     - SEE WRITEUP OF FFTRIG AND FAX  *
!*                      INC - THE WORD INCREMENT BETWEEN SUCCESSIVE     *
!*                           ELEMENTS OF EACH DATA OR COEFFICIENT VECTOR*
!*                           E.G. INC=1 FOR CONSECUTIVELY STORED DATA.  *
!*                      JUMP - THE WORD INCREMENT BETWEEN THE FIRST     *
!*                            ELEMENTS OF SUCCESSIVE DATA OR COEFFICIENT*
!*                            VECTORS.                                  *
!*                      N - THE LENGTH OF EACH TRANSFORM. (SEE NOTE X)  *
!*                      M - THE NUMBER OF TRANSFORMS TO BE DONE         *
!*                          SIMULTANEOUSLY.                             *
!*                      ISIGN - +1 FOR A TRANSFORM FROM FOURIER         *
!*                              COEFFICIENTS TO DATA VALUES.            *
!*                              -1 FOR A TRANSFORM FROM DATA VALUES     *
!*                              TO FOURIER COEFFICIENTS.                *
!*                                                                      *
!!*                  3.OUTPUT                                            *
!*                      A - CONTAINS EITHER THE COEFFICIENTS OR THE     *
!*                          DATA VALUES,DEPENDING ON ISIGN.             *
!*                          IN EACH CASE N INDEPENDENT QUANTITIES       *
!*                          OCCUPY N+2 WORDS.   THE COEFFICIENTS ARE    *
!*                          STORED AS SUCCESSIVE PAIRS OF REAL AND      *
!*                          IMAGINARY PARTS -                           *
!*                          A(K),B(K) , K=0,1,...N/2                    *
!*                          B(0) AND B(N/2) ARE STORED ALTHOUGH THEY    *
!*                          MUST BE 0.                                  *
!*                      FOR FFT99 THE DATA IS STORED WITH EXPLICIT      *
!*                          PERIODICITY -                               *
!*                          X(N-1),X(0),X(1),....X(N-1),X(0)            *
!*                      FOR FFT991 THE DATA APPEARS AS -                *
!*                          X(0),X(1),X(2),......X(N-1),0,0             *
!*                                                                      *
!* NOTES            1. ON CRAY-1, ARRANGE DATA SO THAT JUMP IS NOT A    *
!*                     MULTIPLE OF 8 (TO AVOID MEMORY BANK CONFLICTS)   *
!*                                                                      *
!* WRITE UP         COMPUTER BULLETIN B6.6/1                            *
!*                                                                      *
!* ENTRY POINTS        FFT99,FFT991                                     *
!*                                                                      *
!* COMMON BLOCKS    NONE                                                *
!*                                                                      *
!* I/O              NONE                                                *
!*                                                                      *
!* PRECISION        SINGLE                                              *
!*                                                                      *
!* OTHER ROUTINES   FFT99A,FFT99B,VPASSM          (CY)                  *
!*       REQUIRED   CAL99,CPASS                   (CR)                  *
!*                                                                      *
!*                                                                      *
!* 7/80                      FFT99-1                                    *
!*                                                                      *
!************************************************************************
!*                                                                      *
!* C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                      *
!*                                                               FFT99  *
!*                                                               FFT991 *
!*                                                                      *
!* ACCESS (OBJECT)  CYBER:                                              *
!*                           ATTACH,ECLIB.                              *
!*                           LDSET(LIB=ECLIB)                           *
!*                  CRAY 1:                                             *
!*                           LDR(LIB=ECLIB...)                          *
!*                                                                      *
!* ACCESS (SOURCE)           ATTACH,OLDPL,ECLIBPL                       *
!*                                                                      *
!*                  CYBER :         %DEFINE CYBER                       *
!*                  CRAY:           %DEFINE CRAY                        *
!*                                  %C    FFT99,FFT991                  *
!*                                                                      *
!* LANGUAGE         FORTRAN                                             *
!*                  BUT CRAY IMPLEMENTATION OF PASS IS IN CAL           *
!*                                                                      *
!* SPECIALIST       CLIVE TEMPERTON                                     *
!*                                                                      *
!* HISTORY          WRITTEN BY C.TEMPERTON      JAN     1979            *
!*                                                                      *
!* ALGORITHM        THE ALGORITHM IS THE SELF-SORTING (TEMPERTON)       *
!*                  VERSION OF THE FAST FOURIER TRANSFORM               *
!*                                                                      *
!* REFERENCES       ECMWF TECHNICAL REPORT NO.3                         *
!*                  ECMWF INTERNAL REPORT NO.21 -   C.TEMPERTON         *
!*                                                                      *
!* OBJECT SIZE               FFT991  FFT99  (OCTAL WORDS)               *
!*                  CYBER:    2665    2676                              *
!*                  CRAY :    1250    1260                              *
!*                                                                      *
!*                                                                      *
!* ACCURACY                                                             *
!*                                                                      *
!* TIMING           VECTORIZATION IS ON VECTORS OF LENGTH M.      (CR)  *
!*                  HENCE TIMING IS STRONGLY DEPENDENT ON M.            *
!*                  TIME PER TRANSFORM ON CRAY-1 (MICROSECONDS)         *
!*                  N    M=4    M=16    M=64                            *
!*                 64     46      17      10                            *
!*                128     81      33      21                            *
!*                180    150      58      37                            *
!*                192    149      58      36                            *
!*                240    192      76      49                            *
!*                256    191      76      49                            *
!*                288    219      89      58                            *
!*                300    253     102      68                            *
!*                320    248     101      66                            *
!*                360    286     118      79                            *
!*               1024    898     359     238                            *
!*                                                                      *
!* PORTABILITY      STANDARD FORTRAN                                    *
!*                  STANDARD CAL  (CR)                                  *
!*                                                                      *
!* SYSTEM ROUTINES  NONE                                                *
!*        REQUIRED                                                      *
!*                                                                      *
!* 7/80                      FFT99-1                                    *
!*                                                                      *
!************************************************************************
!C
!C     SUBROUTINE 'FFT99' - MULTIPLE FAST REAL PERIODIC TRANSFORM
!C     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
!C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!C     IS GIVEN BY COOLEY, LEWIS ' WELCH (J. SOUND VIB., VOL. 12
!C     (1970), 315-337)
!C
!C     A IS THE ARRAY CONTAINING INPUT ' OUTPUT DATA
!C     WORK IS AN AREA OF SIZE (N+1)*LOT
!C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!C     INC IS THE INCREMENT WITHIN EACH DATA "VECTOR"
!C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!C     N IS THE LENGTH OF THE DATA VECTORS
!C     LOT IS THE NUMBER OF DATA VECTORS
!C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!C
!C     ORDERING OF COEFFICIENTS:
!C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!C
!C     ORDERING OF DATA:
!C         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
!C         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
!C
!C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!C     PARALLEL
!C
!C     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!C
!C     DEFINITION OF TRANSFORMS:
!C     -------------------------
!C
!C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!C
!C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!C
!C
!C     LINE FOLLOWING NEXT IS NOT SUBROUTINE HEADER(ONLY COMMENT)
!C
!C     SUBROUTINE FFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!C
!C     END
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
!C
      NFAX=IFAX(1)
      IF(NFAX.LE.0) GO TO 99
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
!C
!C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=INC+1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
!C
      IGO=60
      GO TO 40
!C
!C     PREPROCESSING (ISIGN=+1)
!C     ------------------------
!C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
!C
!C     COMPLEX TRANSFORM
!C     -----------------
!C
   40 CONTINUE
      IA=INC+1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,  &
        INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,  &
         2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
!C
      IF (ISIGN.EQ.-1) GO TO 130
!C
!C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=IA
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
!C
!C     FILL IN CYCLIC BOUNDARY POINTS
  110 CONTINUE
      IA=1
      IB=N*INC+1
!CDIR$ IVDEP
      DO 120 L=1,LOT
      A(IA)=A(IB)
      A(IB+INC)=A(IA+INC)
      IA=IA+JUMP
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!C
!C     POSTPROCESSING (ISIGN=-1):
!C     --------------------------
!C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!C
  140 CONTINUE
      RETURN
!C  ** ERROR EXIT   IFAX(1) LE 0 **
99    print *,' FFT99 CALLED BUT FACTORS NOT SUPPLIED '
!C99    CALL REMARK(" FFT99 CALLED BUT FACTORS NOT SUPPLIED ")
      CALL ABORT
      END

      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
	  call rfftmlt(a,work,trigs,ifax,inc,jump,n,lot,isign)
	  return
	  end

      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
!C
!C     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
!C     (SPECTRAL TO GRIDPOINT TRANSFORM)
!C
!C     LINE FOLLOWING NEXT IS NOT SUBROUTINE HEADER(ONLY COMMENT)
!C
!C     SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
!C     END
      DIMENSION A(N),WORK(N),TRIGS(N)
      NH=N/2
      NX=N+1
      INK=INC+INC
!C
!C     A(0) ' A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
!CDIR$ IVDEP
      DO 10 L=1,LOT
      WORK(JA)=A(IA)+A(IB)
      WORK(JB)=A(IA)-A(IB)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
!C
!C     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
!C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
!CDIR$ IVDEP
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))-      &
         (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+      &
         (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+  &
         (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-  &
         (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
!C
      IF (IABASE.NE.IBBASE) GO TO 50
!C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
!CDIR$ IVDEP
      DO 40 L=1,LOT
      WORK(JA)=2.0*A(IA)
      WORK(JA+1)=-2.0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
!C
   50 CONTINUE
      RETURN
      END

      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!C
!C     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
!C     (GRIDPOINT TO SPECTRAL TRANSFORM)
!C
!C
!C     LINE FOLLWING NEXT IS NOT SUBROUTINE HEADER(ONLY COMMENT)
!C
!C     SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!C     END
      DIMENSION WORK(N),A(N),TRIGS(N)
!C
      NH=N/2
      NX=N+1
      INK=INC+INC
!C
!C     A(0) ' A(N/2)
      SCALE=1.0/FLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
!CDIR$ IVDEP
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JB)=SCALE*(WORK(IA)-WORK(IB))
      A(JA+INC)=0.0
      A(JB+INC)=0.0
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
!C
!C     REMAINING WAVENUMBERS
      SCALE=0.5*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
!C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
!CDIR$ IVDEP
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))  &
      +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))  &
      -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))&
      +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))&
      -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
!C
      IF (IABASE.NE.IBBASE) GO TO 50
!C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0*SCALE
!CDIR$ IVDEP
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
!C
   50 CONTINUE
      RETURN
      END

 
!************************************************************************
      SUBROUTINE FFTCC(Z,W,EX,IFAX,INC,JUMP,N,NFT,ISIGN)
!* FFT COMPLEX-->COMPLEX USING TEMPERTON'S REAL-->COMPLEX FFT
!* SAME ARGUMENTS AS FFT991
!* N MUST BE EVEN
      IF(ISIGN.EQ.1)THEN
      CALL TRCCC(Z,W,INC,JUMP,N,NFT,ISIGN)
      CALL FFT991(Z,W,EX,IFAX,INC,JUMP,N,NFT,ISIGN)
      RETURN
      ELSE
      CALL FFT991(Z,W,EX,IFAX,INC,JUMP,N,NFT,ISIGN)
      CALL TRCCC(Z,W,INC,JUMP,N,NFT,ISIGN)
      RETURN
      ENDIF
      END
       Subroutine Fftfax (n,ifax,trigs)
!*      Emulation de fftfax par appel a fax et fftrig
	   dimension ifax(10),trigs(n)
	   call fax(ifax,n,3)
	   call fftrig(trigs,n,3)
	   return
	   end
      SUBROUTINE FFTRIG (TRIGS,N,MODE)
      DIMENSION TRIGS(N)
!C    FFTRIG RETURNS AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES
!C    SUBSEQUENTLY USED BY    F F T    ROUTINES
!C    SEE COMMENTS IN ROUTINE    F A X
!C     END
      PI=2.0*ASIN(1.0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/FLOAT(NN)
      L=NN+NN
      DO 10 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(LA+I)=COS(ANGLE)
      TRIGS(LA+I+1)=SIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=2.0*SIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5*DEL
      DO 50 I=2,N
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=SIN(ANGLE)
   50 CONTINUE
      RETURN
      END

      SUBROUTINE rfftmlt(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!C
!C     SUBROUTINE 'FFT991' - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!C     FAST FOURIER TRANSFORM
!!C
!C     The routine has been recalled rfftmlt to agree with new ffts
!C     on Cray computers
!C
!************************************************************************
!*                                                                      *
!* C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                      *
!*                                                               FFT991 *
!*                                                                      *
!*                                                                      *
!* SUBPROGRAM       SUBROUTINE    FFT99                                 *
!*                                FFT991                                *
!*                                                                      *
!* PURPOSE          PERFORM MULTIPLE FAST FOURIER TRANSFORMS            *
!*                                                                      *
!*!                                                                      *
!* VERSION          CYBER                         CRAY-1                *
!*                                                                      *
!*                  JAN 1979 ORIGINAL             JAN 1979 ORIGINAL     *
!*                                                                      *
!* USAGE                                                                *
!*                  CALL FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   *
!*                  CALL RFFTMLT(A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)  *
!*                                                                      *
!* ARGUMENTS        1.DIMENSION                                         *
!*                       A(IDIM),WORK((N+1)*M),TRIGS(3*N/2),IFAX(10)    *
!*                       WORK IS A WORK ARRAY                           *
!*                                                                      *
!*                  2.INPUT                                             *
!*                      A - AN ARRAY CONTAINING THE INPUT DATA OR       *
!*                          COEFFICIENT VECTORS.                        *
!*                         THIS ARRAY IS OVERWRITTEN BY THE RESULTS.    *
!*                      TRIGS AND IFAX - ARRAYS SET UP BY FFTRIG AND FAX*
!*                                     - SEE WRITEUP OF FFTRIG AND FAX  *
!*                      INC - THE WORD INCREMENT BETWEEN SUCCESSIVE     *
!*!                           ELEMENTS OF EACH DATA OR COEFFICIENT VECTOR*
!*                           E.G. INC=1 FOR CONSECUTIVELY STORED DATA.  *
!*                      JUMP - THE WORD INCREMENT BETWEEN THE FIRST     *
!*                            ELEMENTS OF SUCCESSIVE DATA OR COEFFICIENT*
!*                            VECTORS.                                  *
!*                      N - THE LENGTH OF EACH TRANSFORM. (SEE NOTE X)  *
!*                      M - THE NUMBER OF TRANSFORMS TO BE DONE         *
!*                          SIMULTANEOUSLY.                             *
!*                      ISIGN - +1 FOR A TRANSFORM FROM FOURIER         *
!*                              COEFFICIENTS TO DATA VALUES.            *
!*                              -1 FOR A TRANSFORM FROM DATA VALUES     *
!*                              TO FOURIER COEFFICIENTS.                *
!*                                                                      *
!*                  3.OUTPUT                                            *
!*                      A - CONTAINS EITHER THE COEFFICIENTS OR THE     *
!*                          DATA VALUES,DEPENDING ON ISIGN.             *
!*                          IN EACH CASE N INDEPENDENT QUANTITIES       *
!*                          OCCUPY N+2 WORDS.   THE COEFFICIENTS ARE    *
!*                          STORED AS SUCCESSIVE PAIRS OF REAL AND      *
!*                          IMAGINARY PARTS -                           *
!*                          A(K),B(K) , K=0,1,...N/2                    *
!*                          B(0) AND B(N/2) ARE STORED ALTHOUGH THEY    *
!*                          MUST BE 0.                                  *
!*                      FOR FFT99 THE DATA IS STORED WITH EXPLICIT      *
!*                          PERIODICITY -                               *
!*                          X(N-1),X(0),X(1),....X(N-1),X(0)            *
!*                      FOR FFT991 THE DATA APPEARS AS -                *
!*                          X(0),X(1),X(2),......X(N-1),0,0             *
!*                                                                      *
!* NOTES            1. ON CRAY-1, ARRANGE DATA SO THAT JUMP IS NOT A    *
!*                     MULTIPLE OF 8 (TO AVOID MEMORY BANK CONFLICTS)   *
!*                                                                      *
!* WRITE UP         COMPUTER BULLETIN B6.6/1                            *
!*                                                                      *
!* ENTRY POINTS        FFT99,FFT991                                     *
!*                                                                      *
!* COMMON BLOCKS    NONE                                                *
!*                                                                      *
!* I/O              NONE                                                *
!*                                                                      *
!* PRECISION        SINGLE                                              *
!*                                                                      *
!* OTHER ROUTINES   FFT99A,FFT99B,VPASSM          (CY)                  *
!*       REQUIRED   CAL99,CPASS                   (CR)                  *
!*                                                                      *
!*                                                                      *
!* 7/80                      FFT99-1                                    *
!*                                                                      *
!************************************************************************
!*                                                                      *
!* C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                      *
!*                                                               FFT99  *
!*                                                               FFT991 *
!*                                                                      *
!* ACCESS (OBJECT)  CYBER:                                              *
!*                           ATTACH,ECLIB.                              *
!*                           LDSET(LIB=ECLIB)                           *
!*!                  CRAY 1:                                             *
!*                           LDR(LIB=ECLIB...)                          *
!*                                                                      *
!* ACCESS (SOURCE)           ATTACH,OLDPL,ECLIBPL                       *
!*                                                                      *
!*                  CYBER :         %DEFINE CYBER                       *
!*                  CRAY:           %DEFINE CRAY                        *
!*                                  %C    FFT99,FFT991                  *
!*                                                                      *
!* LANGUAGE         FORTRAN                                             *
!*                  BUT CRAY IMPLEMENTATION OF PASS IS IN CAL           *
!*                                                                      *
!* SPECIALIST       CLIVE TEMPERTON                                     *
!*                                                                      *
!* HISTORY          WRITTEN BY C.TEMPERTON      JAN     1979            *
!*                                                                      *
!* ALGORITHM        THE ALGORITHM IS THE SELF-SORTING (TEMPERTON)       *
!*                  VERSION OF THE FAST FOURIER TRANSFORM               *
!*                                                                      *
!* REFERENCES       ECMWF TECHNICAL REPORT NO.3                         *
!*                  ECMWF INTERNAL REPORT NO.21 -   C.TEMPERTON         *
!*                                                                      *
!* OBJECT SIZE               FFT991  FFT99  (OCTAL WORDS)               *
!*                  CYBER:    2665    2676                              *
!*                  CRAY :    1250    1260                              *
!*                                                                      *
!*                                                                      *
!* ACCURACY                                                             *
!*                                                                      *
!* TIMING           VECTORIZATION IS ON VECTORS OF LENGTH M.      (CR)  *
!*                  HENCE TIMING IS STRONGLY DEPENDENT ON M.            *
!*                  TIME PER TRANSFORM ON CRAY-1 (MICROSECONDS)         *
!*                  N    M=4    M=16    M=64                            *
!*                 64     46      17      10                            *
!*                128     81      33      21                            *
!*                180    150      58      37                            *
!*                192    149      58      36                            *
!*                240    192      76      49                            *
!*                256    191      76      49                            *
!*                288    219      89      58                            *
!*                300    253     102      68                            *
!*                320    248     101      66                            *
!*                360    286     118      79                            *
!*               1024    898     359     238                            *
!*                                                                      *
!* PORTABILITY      STANDARD FORTRAN                                    *
!*                  STANDARD CAL  (CR)                                  *
!*                                                                      *
!* SYSTEM ROUTINES  NONE                                                *
!*        REQUIRED                                                      *
!*                                                                      *
!* 7/80                      FFT99-1                                    *
!*                                                                      *
!************************************************************************
!C
!C     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
!C     THAT IN MRFFT2
!C
!C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!C     IS GIVEN BY COOLEY, LEWIS ' WELCH (J. SOUND VIB., VOL. 12
!C     (1970), 315-337)
!C
!C     A IS THE ARRAY CONTAINING INPUT ' OUTPUT DATA
!C     WORK IS AN AREA OF SIZE (N+1)*LOT
!C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!C     INC IS THE INCREMENT WITHIN EACH DATA "VECTOR"
!C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!C     N IS THE LENGTH OF THE DATA VECTORS
!C     LOT IS THE NUMBER OF DATA VECTORS
!C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!C
!C     ORDERING OF COEFFICIENTS:
!C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!C
!C     ORDERING OF DATA:
!C         X(0),X(1),X(2),...,X(N-1)
!C
!C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!C     PARALLEL
!C
!C     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!C
!C     DEFINITION OF TRANSFORMS:
!C     -------------------------
!C
!C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!C
!C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!C
!C     SUBROUTINE RFFTMLT(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!C     
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
!C
      NFAX=IFAX(1)
      IF(NFAX.LE.0) GO TO 99
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
!C
!C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
!C
      IGO=60
      GO TO 40
!C
!C     PREPROCESSING (ISIGN=+1)
!C     ------------------------
!C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
!C
!C     COMPLEX TRANSFORM
!C     -----------------
!C
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS, &
        INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,  &
         2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
!C
      IF (ISIGN.EQ.-1) GO TO 130
!C
!C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
!CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
!C
!C     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1
!CDIR$ IVDEP
      DO 120 L=1,LOT
      A(IB)=0.0
      A(IB+INC)=0.0
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!C
!C     POSTPROCESSING (ISIGN=-1):
!C     --------------------------
!C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!C
  140 CONTINUE
      RETURN
!C   **  ERROR     IFAX(1) LE 0  **
99    print *,' RFFTMLT CALLED BUT FACTORS NOT SUPPLIED '
      CALL ABORT
      END
!C*****************************************************************************
!C                                                                 second
!C***************************************************************************** 
      real function second()
!c      real tuser, tsyst
!c     Commente par G. Brunet
!c      total = etime(tuser,tsyst)
      second = 9999999.
      return 
      end
!C*****************************************************************************
!C                                                                 mxma
!C*****************************************************************************
      subroutine mxma (a,na,iad,b,nb,ibd,c,nc,icd,nac,nab,nbc)
      dimension a(1),b(1),c(1)
!*     Multiplication matricielle
!*     c = a.b
!*
!* 	a: 	premiere matrice du produit
!*	na: 	espacement entre les elements de colonne de a
!* 	iad:	espacement entre les elements de ligne de a
!* 	b: 	deuxieme matrice du produit
!*	nb: 	espacement entre les elements de colonne de b
!* 	ibd:	espacement entre les elements de ligne de b
!* 	c: 	matrice du produit
!*	nc: 	espacement entre les elements de colonne de c
!* 	icd:	espacement entre les elements de ligne de c
!*	nac:	nombre de lignes dans a et c
!*	nab:	nombre de colonnes dans a et de lignes dans b
!*	nbc:	nombre de colonnes dans b et c
!*
!*       On calcule pour i = 0, nac-1
!*                    et j = 0, nbc-1
!*       c(1+i*nc+j*icd) = sum ( k = 0, nab-1) a(1+i*na+k*iad) * b(1+k*nb+j*ibd)
!*  
        do 1000 j = 0, nbc-1
        do 1000 i = 0, nac-1
           sum = 0.
           do 500 k = 0, nab-1
              sum = sum + a(1+i*na+k*iad) * b(1+k*nb+j*ibd)
 500       continue
 1000   c(1+i*nc+j*icd) = sum
        return 
        end 
!************************************************************************
      SUBROUTINE TRCCC(F,W,INC,JUMP,N,NFT,ISIGN)
!*
!* THIS ROUTINE IS USED TO COMPUTE COMPLEX TO COMPLEX TRANSFORMS
!*       USING TEMPERTON'S ROUTINE FFT991
!*
!*
!*
!* AFTER FT ALONG X1, ONE HAS
!*
!*    SUM   F(X)COS(K1X1)      SUM   F(X)SIN(K1X1)
!*     X1                       X1
!*
!*--------------------------------------------------------------
!* AFTER FT ALONG X2 USING FFT991 ONE HAS
!*
!*     SUM  F(X)COS(K1X1)SIN(K2X2)   SUM F(X)SIN(K1X1)SIN(K2X2)
!*      X                             X
!*
!*     SUM  F(X)COS(K1X1)COS(K2X2)   SUM F(X)SIN(K1X1)COS(K2X2)
!*      X                             X
!*
!*---------------------------------------------------------------
!* WHAT WE WANT IS
!*
!*     SUM  F(X)( COS(K1X1)COS(K2X2)-SIN(K1X1)SIN(K2X2) )
!*      X
!*
!*     SUM  F(X)( COS(K1X1)SIN(K2X2)+SIN(K1X1)COS(K2X2) )
!*      X
!*
!*---------------------------------------------------------------
      REAL F(1),W(2,0:1)
      NH=N/2
      NFTH=NFT/2
      I1=1
      INCD=INC*2
      IF(ISIGN.EQ.-1)THEN
        DO 1 IFT=1,NFTH
        I2=I1+1
        I1P=I1+INC
        I2P=I2+INC
        DO 3 J=0,NH
        W(1,J)=F(I1+INCD*J)-F(I2P+INCD*J)
        W(2,J)=F(I2+INCD*J)+F(I1P+INCD*J)
        W(1,N-J)=F(I1+INCD*J)+F(I2P+INCD*J)
        W(2,N-J)=F(I2+INCD*J)-F(I1P+INCD*J)
  3     CONTINUE
        DO 5 J=0,N-1
        F(I1+INC*J)=W(1,J)
        F(I2+INC*J)=W(2,J)
  5     CONTINUE
        I1=I1+JUMP*2
  1     CONTINUE
      ELSE
        DO 2 IFT=1,NFTH
        I2=I1+1
        I1N=I1+N*INC
        I2N=I2+N*INC
        DO 4 J=1,NH-1
        W(1,2*J)=.5*(F(I1+INC*J)+F(I1N-INC*J))
        W(2,2*J)=.5*(F(I2+INC*J)+F(I2N-INC*J))
        W(1,2*J+1)=.5*(F(I2+INC*J)-F(I2N-INC*J))
        W(2,2*J+1)=.5*(F(I1N-INC*J)-F(I1+INC*J))
  4     CONTINUE
        W(1,0)=F(I1)
        W(2,0)=F(I2)
        W(1,1)=0.
        W(2,1)=0.
        W(1,N)=F(I1N)
        W(2,N)=F(I2N)
        W(1,N+1)=0.
        W(2,N+1)=0.
        DO 6 J=0,N+1
        F(I1+INC*J)=W(1,J)
        F(I2+INC*J)=W(2,J)
  6     CONTINUE
        I1=I1+JUMP*2
  2     CONTINUE
      ENDIF
      RETURN
      END
