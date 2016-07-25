!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c SOUS-PROGRAMMES UTILISES POUR LE SCHEMA NUMERIQUE
!c DU MODELE SPHERIQUE et les differents termes des equations
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine onestep(y,x)
!c
!c  Performs one step forward of predictor-corrector scheme
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real x(nvaria),y(nvaria)
      real wk1(nvaria),wk2(nvaria)
!c
!c  Etape de prediction
!c
      call fun(wk1,x)
      do nv=1,nvaria
         wk1(nv) = x(nv) + dt*wk1(nv)
      enddo
!c
!c  Etape de correction
!c
      call fun(wk2,wk1)
      do nv=1,nvaria
         wk2(nv) = x(nv) + dt*wk2(nv)
         y(nv) = 0.5*(wk1(nv) + wk2(nv))
      enddo
!c
!c
      return
      end
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fun(y,x)
!c
!c Calcule le second membre f(x) defini par dx/dt=f(x)
!c Ici, y=f(x)=-J(psi,x)-D(psi,x)+Forcage
!c J etant le jacobien et D les phenomenes de dissipation
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real x(nvaria2,nlevels),y(nvaria2,nlevels)
      real psi(nvaria2,nlevels),tseta(nvaria2,nlevels)
      real sdisp(nvaria2,nlevels)
!c
!c Calcul de la fonction de courant
!c
!c      print*,x
      call q2psi(psi,x)
!c
!c Calcul du tourbillon relatif
!c
      do nl=1,nlevels
         call laplace(tseta(1,nl),psi(1,nl),2)
      enddo
!c
!c Calcul de la dissipation
!c
      call dissip(sdisp,x,psi,tseta)
!c
!c  Calcul du forcage (le resultat est en common)
!c
      call forcing(psi,x,tseta)
!c
!c Calcul pour chaque niveau 
!c
      do nl=1,nlevels
        call jacobs(y(1,nl),x(1,nl),psi(1,nl))
        call saxpy(nvaria2,-1.,sdisp(1,nl),1,y(1,nl),1)
        call saxpy(nvaria2,1.,force(1,nl),1,y(1,nl),1)
      enddo
!c
	return
	end
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine forcing(psi,q,tseta)
!c
!c  Calcul du forcage selon iopforce
!c
       include 'sphectra.h'
       include 'paramod.h'
!c
       real q(nvaria2,nlevels),psi(nvaria2,nlevels),tseta(nvaria2,nlevels)
!c
       return
       end
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dissip(sdisp,q,psi,tseta)
!c
!c Calcul le terme de dissipation
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real q(nvaria2,nlevels),psi(nvaria2,nlevels),tseta(nvaria2,nlevels)
      real sdisp(nvaria2,nlevels)
      real sdisp1(nvaria2),sdisp2(nvaria2,nlevels),sdisp3(nvaria2,nlevels)
!c
!c  Dissipation d'Ekman'
!c
      call dekman(sdisp1,tseta(1,nlevels),psi(1,nlevels))
!*     call sset(nvaria2,0.,sdisp1,1)
!c
!c  Dissipation selective
!c
      call dissel(sdisp2,q)
!*     call sset(nvaria,0.,sdisp2,1)
!c
!c  Diffusion thermique
!c
      call dtherm(sdisp3,psi)
!*     call sset(nvaria,0.,sdisp3,1)
!c
!c  Somme (Termes D1,D2,D3 (+))
!c
      do nv=1,nvaria2
         do nl=1,nlevels-1
            sdisp(nv,nl) = sdisp3(nv,nl) + sdisp2(nv,nl)
         enddo
         sdisp(nv,nlevels) = sdisp3(nv,nlevels) + sdisp2(nv,nlevels) + sdisp1(nv)
      enddo
!c
      return
      end
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dekman(sdisp,tseta,psi)
!c
!c  Calcul de la dissipation d'Ekman (Version Marshall & Molteni, JAS 93  Il s'
!c  agit du terme +EK3 de leur appendix
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real sdisp(nvaria2),tseta(nvaria2),psi(nvaria2)
      real swrk(nvaria2),gtseta(n2long,nlat),gdisp(n2long,nlat)
!c
!c  Calcu du premier terme k.lap(psi)
!c
      call spectogrid(gtseta,tseta)
      call shprod(n2long*nlat,gtseta,1,gek,1,gdisp,1)
      call gridtospec(sdisp,gdisp)
!c
!c  Calcul du 2eme terme grad(k).grad(psi)
!c
      call progras(swrk,sek,psi)
!c
!c  Addition
!c
      call saxpy(nvaria2,1.,swrk,1,sdisp,1)
!c
      return
      end
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dissel(sdisp,q)
!c
!c  Dissipation selective (Termes +Hi dans Marshall & Molteni JAS 93)
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real q(nvaria2,nlevels),sdisp(nvaria2,nlevels)
      real wrk(nvaria2)
!c
      do nl=1,nlevels
         call scopy(nvaria2,q(1,nl),1,wrk,1)
         wrk(1) = wrk(1) - f1
         if(nl.eq.nlevels) then
           call saxpy(nvaria2,-1.,ftopo,1,wrk,1)
         endif
         do nv=1,nvaria2
            sdisp(nv,nl) = wrk(nv)*oprdissip(nv)
         enddo
      enddo
!c
      return
      end
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dtherm(sdisp,psi)
!c
!c  Diffusion thermique
!c  Signe + dans MM 93
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real sdisp(nvaria2,nlevels),psi(nvaria2,nlevels)
      real tr12(nvaria2),tr23(nvaria2)
!c
!c  Calcul de TR12 et TR23
!c
      do nv=1,nvaria2
         tr12(nv) = ray2n2ontor*(psi(nv,1) - psi(nv,2))
         tr23(nv) = ray4n2ontor*(psi(nv,2) - psi(nv,3))
      enddo
!c
!c  Niveau par niveau
!c
      do nv=1,nvaria2
         sdisp(nv,1) = -tr12(nv)
         sdisp(nv,2) = tr12(nv) - tr23(nv)
         sdisp(nv,3) = tr23(nv)
      enddo
!c
      return
      end
