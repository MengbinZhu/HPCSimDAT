      Subroutine q2psi(psi,q)
!c
!c  Transformation Q-->PSI
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real psi(nvaria2,nlevels),q(nvaria2,nlevels)
      real trv(nvaria2,nlevels)
      real trv1(nlevels),trv2(nlevels)
!c
!c  Translation de q pour obtenir le TP relatif
!c  C est Q - f - fh/Ho
!c
      call scopy(nvaria2*nlevels,q,1,trv,1)
      call saxpy(nvaria2,-1.,ftopo,1,trv(1,3),1)
      do nl=1,nlevels
         trv(1,nl) = trv(1,nl) - f1
      enddo
!c
!c  Partie lineaire de la transformation
!c
      do nv=1,nvaria2
         do nl=1,nlevels
            trv1(nl) = trv(nv,nl)
         enddo
         call murrv(nlevels,nlevels,cq2psi(1,1,nv),nlevels,nlevels,trv1,1,nlevels,trv2)
         do nl=1,nlevels
            psi(nv,nl) = trv2(nl)
         enddo
      enddo
!c
      return
      end
!c
      subroutine psi2q(q,psi)
!c
!c Programme passant de psi, la fonction de courant (geost.), a q,
!c la vorticite potentielle.
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real psi(nvaria2,nlevels),q(nvaria2,nlevels)
!c
      call sset(nvaria2*nlevels,0.,q,1)
!c     f1=omega2/sqrt(3.)
!c     namax=2*mvadeb(ntrunc-1)
!c
!c
      do nv2=1,nvaria2
         oplap = oprlap((nv2-1)/2+1,2)
         q(nv2,1)=oplap*psi(nv2,1)-(psi(nv2,1)-psi(nv2,2))*ray2n2
         q(nv2,2)=oplap*psi(nv2,2)        &
                 -(psi(nv2,2)-psi(nv2,3))*ray4n2    &
                 +(psi(nv2,1)-psi(nv2,2))*ray2n2
         q(nv2,3)=oplap*psi(nv2,3)        &
                 +(psi(nv2,2)-psi(nv2,3))*ray4n2   &
                 +ftopo(nv2)
      enddo
!c
!c  Addition de La vorticite planetaire
!c
      do nl=1,nlevels
         q(1,nl)=q(1,nl)+f1
      enddo
!c
      return
      end
!c
      subroutine dq2dpsi(psi,q)
!c
!c  Transformation Q-->PSI. Ne sont retenus que les termes lineaires
!c  Les termes constants sont retires.
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real psi(nvaria2,nlevels),q(nvaria2,nlevels)
      real trv(nvaria2,nlevels)
      real trv1(nlevels),trv2(nlevels)
!c
      call scopy(nvaria2*nlevels,q,1,trv,1)
!c
!c  Partie lineaire de la transformation
!c
      do nv=1,nvaria2
         do nl=1,nlevels
            trv1(nl) = trv(nv,nl)
         enddo
         call murrv(nlevels,nlevels,cq2psi(1,1,nv),nlevels,nlevels,trv1,1,nlevels,trv2)
         do nl=1,nlevels
            psi(nv,nl) = trv2(nl)
         enddo
      enddo
!c
      return
      end
!c
      subroutine dpsi2dq(q,psi)
!c
!c Programme passant de psi, la fonction de courant (geost.), a q,
!c la vorticite potentielle. N est retenu ici que la partie lineaire
!c de la transformation
!c
      include 'sphectra.h'
      include 'paramod.h'
!c
      real psi(nvaria2,nlevels),q(nvaria2,nlevels)
!c
      call sset(nvaria2*nlevels,0.,q,1)
!c     f1=omega2/sqrt(3.)
!c     namax=2*mvadeb(ntrunc-1)
!c
!c
      do nv2=1,nvaria2
         oplap = oprlap((nv2-1)/2+1,2)
         q(nv2,1)=oplap*psi(nv2,1)-(psi(nv2,1)-psi(nv2,2))*ray2n2
         q(nv2,2)=oplap*psi(nv2,2)      &
                 -(psi(nv2,2)-psi(nv2,3))*ray4n2  &
                 +(psi(nv2,1)-psi(nv2,2))*ray2n2
         q(nv2,3)=oplap*psi(nv2,3)      &
                 +(psi(nv2,2)-psi(nv2,3))*ray4n2
      enddo
!c
      return
      end
!c
        subroutine zetap2q(q,zetap)
!c       On a q=lap(psi)+R psi + b = zetap +b. Le vecteur b est le vecteur constant
!c       contennant la vorticite planetaire et la topographie.

      include 'sphectra.h'
      include 'paramod.h'
      real zetap(nvaria2,nlevels),q(nvaria2,nlevels)

      do nl=1,nlevels
         call scopy(nvaria2,zetap(1,nl),1,q(1,nl),1)
         q(1,nl) = q(1,nl) + f1
         if(nl.eq.nlevels) then
           call saxpy(nvaria2,1.,ftopo,1,q(1,nl),1)
         endif
      enddo
      return
      end
!c
!c
      subroutine q2zetap(zetap,q)
!c    On a q=lap(psi)+R psi + b = zetap +b. Le vecteur b est le vecteur constant
!c    contennant la vorticite planetaire et la topographie.

      include 'sphectra.h'
      include 'paramod.h'
      real zetap(nvaria2,nlevels),q(nvaria2,nlevels)

      do nl=1,nlevels
         call scopy(nvaria2,q(1,nl),1,zetap(1,nl),1)
         zetap(1,nl) = zetap(1,nl) - f1
         if(nl.eq.nlevels) then
           call saxpy(nvaria2,-1.,ftopo,1,zetap(1,nl),1)
         endif
      enddo
      return
      end

