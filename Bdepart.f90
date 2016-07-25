
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dialogue
!c
!c  Interface utilisateur (boutons du modele
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      character*80 fntopo,fnforce,fnlmask,fndissip
!c
!c  Reading of the names of files topography (grid), lmask (% of land), force (constant)
!c
      fntopo = '/scratch3/BMC/nim/Brian.Etherton/QGmodel1/topo.t21_le'
      fnlmask = '/scratch3/BMC/nim/Brian.Etherton/QGmodel1/lmask.sus_le'
      fnforce = '/scratch3/BMC/nim/Brian.Etherton/QGmodel1/force1.rec_le'
!c
!c  Foring options:
!c!  0) Forcing equal to 0
!c!  1) Constant forcing read from the fnforce file
!c!  2)
!c!
      iopforce = 1
!c!
      open(71,file=fntopo,form='unformatted')
      open(72,file=fnlmask,form='unformatted')
      open(73,file=fnforce,form='unformatted')
!c!
      print *,'Dialogue ends...'
      return
      end
!c!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine initread
!c!
!c!  Initialisation des parametres lus sur fichier
!c!
      include 'sphectra.h'
      include 'paramod.h'
!c!
      real wrk(nvaria2),gwrk(n2long,nlat),gtopo1(n2long,nlat)
!c!
!c!  Initialisations
!c!
      call sset(n2long*nlat,0.,gtopo,1) 
      call sset(nvaria2,0.,wrk,1)
      call sset(n2long*nlat,0.,gek,1)
!c!
!c  Parametre de Coriolis
!c!
      wrk(1) = f1
      call spectogrid(corio,wrk)
!c!
!c!  Topographie: definition de ftopo,stopo,gtopo
!c!
      read(71)ii,((gtopo(i,j),i=1,nlong),j=1,nlat)
      !gtopo1(1:nlong,1:nlat/2)=gtopo(1:nlong,nlat:nlat/2+1:-1)
      !gtopo1(1:nlong,nlat/2+1:nlat)=gtopo(1:nlong,1:nlat/2)
      !write(81,rec=1)gtopo1(1:nlong,:)
      call scopy(n2long*nlat,gtopo,1,gtopo1,1)
      call gridtospec(stopo,gtopo1)
      call shprod(n2long*nlat,corio,1,gtopo,1,gwrk,1)
      call sscal(n2long*nlat,h0n1,gwrk,1)
      call gridtospec(ftopo,gwrk)
!c!
!c!  Coefficient de friction d'Ekman: calcul de gek et sek.
!c!
      read(72)ii,((glmask(i,j),i=1,nlong),j=1,nlat)
      do i=1,nlong
         do j=1,nlat
            fh = 1. - exp(-gtopo(i,j)/1000.)
            gek(i,j) = ontoe*(1. + 0.5*(glmask(i,j) + fh))
         enddo
      enddo
      call scopy(n2long*nlat,gek,1,gwrk,1)
      call gridtospec(sek,gwrk)
!c!
!c  Forcage
!c
      if(iopforce.eq.0) then
        call sset(nvaria,0.,force,1)
      elseif(iopforce.eq.1) then
        do nl=1,nlevels
           read(73)ii,((gwrk(i,j),i=1,nlong),j=1,nlat)
           call gridtospec(force(1,nl),gwrk)
        enddo
      endif
!c
      print *,'Initread ends...'
!c
      return
      end
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine initdissip
!c
!c  Initialisation des parametres de la dissipation
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
!c  Coeffecients de la dissipation selective
!c
       do nv=1,nvaria2
         index = (nv - 1)/2 + 1
         oplap = oprlap(index,2)*rayter2
         oprdissip(nv) = (oplap/float(ntrunc*(ntrunc+1)))**4/toh
      enddo
!c
!c  Tabulation de la dissipation totale (fonction de PSI)
!c
      print *,'Initdissip ends...'
      return
      end
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine initcq2psi
!c
!c  Calcul des coefficients de la transformation q-->psi dans le cas d'un modele
!  QG a 3 niveaux. NLEVELS = 3 obligatoirement.
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real xm(nlevels,nlevels)
!c
!c  On separe le cas du mode constant
!c
      nvcons1 = 2*mvsdeb(0) - 1
      nvcons2 = nvcons1 + 1
      do nv=1,nvaria2
!c
!c  Adresse du mode complexe associe au mode reel nv
!c
         nvc = (nv - 1)/2 + 1
!c
!c  Coefficient du laplacien associe
!c
         oplap = oprlap(nvc,2)
!c
!c  Matrice des coefficients de la transformation directe
!c
         xm(1,1) = oplap - ray2n2
         xm(1,2) = ray2n2
         xm(1,3) = 0.
         xm(2,1) = ray2n2
         xm(2,2) = oplap -ray2n2 - ray4n2
         xm(2,3) = ray4n2
         xm(3,1) = 0.
         xm(3,2) = ray4n2
         xm(3,3) = oplap - ray4n2
!c
!c  Inversion de la matrice dans le cas non constant
!c  et mise a zero sinon
!c
         if((nv.eq.nvcons1).or.(nv.eq.nvcons2)) then
           call sset(nlevels**2,0.,cq2psi(1,1,nv),1)
         else
           call linrg(3,xm,3,cq2psi(1,1,nv),3)
!*          print *,nv
!*          print *,((xm(i,j),i=1,3),j=1,3)
!*          print *,((cq2psi(i,j,nv),i=1,3),j=1,3)
         endif
!*        if(nv.eq.1) then
!*          do i=1,3
!*             do j=1,3
!*                print *,i,j,xm(i,j),cq2psi(i,j,nv)
!*             enddo
!*          enddo
!*        endif
      enddo
!c
      print *,'Initcq2psi ends...'
      return
      end
!c!
      subroutine initqgmod
!c!
!c!  Initialisations de bases du modele QG3 .
!c!
      call initsph     !Bsimplif.f90
      call dialogue
      call initread
      call initdissip
      call initcq2psi
!cc
      return
      end
