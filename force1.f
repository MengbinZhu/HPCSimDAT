      program forcage
c
c Programme du calcul du forcage S, obtenu par annulation
c de la moyenne temporelle des dq/dt pour chaque niveau.
c Cf Marshall & Molteni JAS 93.
c S=J([psi],[q])+[J(psi',q')]+D([psi])   [ ]=moy. temp.
c Dissipation AVEC Masque
c
      include 'sphectra.h'
      include 'paramod.h'
c
      real spsi(nvaria2,nlevels),sq(nvaria2,nlevels)
      real tendance(nvaria),tendancem(nvaria2,nlevels)
      real ptend(nvaria2,nlevels),gwrk(n2long,nlat),gpsi(n2long,nlat)
      integer mallo(12)
      character*80 fnp1,fnp2,fnp3,fndates,fnforce
c
c  Dialogue
c
      read *,fnp1
      read *,fnp2
      read *,fnp3
      read *,fndates
      read *,fnforce
      read *,nmallo
      read *,(mallo(nm),nm=1,nmallo)
c
c  Ouverture des fichiers
c
      open(11,file=fnp1,form='unformatted')
      open(12,file=fnp2,form='unformatted')
      open(13,file=fnp3,form='unformatted')
      open(14,file=fndates,form='unformatted')
      open(21,file=fnforce,form='unformatted')
c
c  Initialisation du modele
c
      call initqgmod
c
c  Autres initialisations
c
      call sset(nvaria2*nlevels,0.,tendancem,1)
      call sset(n2long*nlat,0.,gpsi,1)
      cnt = 0.
c
c  Lecture et calcul
c
      do n=1,100000
         read(14,end=1001)idate,igap
         call ddate(idate,iy,im,ij,ih)
         if((ibelong(im,nmallo,mallo).eq.1).and.(igap.ne.0)) then
           cnt = cnt + 1.
           print *,n,cnt,idate
           do nl=1,nlevels
              nfi = 10 + nl
              read(nfi)id,ig,((gpsi(i,j),i=1,nlong),j=1,nlat)
              call gridtospec(spsi(1,nl),gpsi)
           enddo
           call psi2q(sq,spsi)
           call fun(tendance,sq)
           call saxpy(nvaria,1.,tendance,1,tendancem,1)
         else
           do nl=1,nlevels
              nfi = 10 + nl
              read(nfi)id
           enddo
         endif
      enddo
1001  continue
c
c  Traitement des sorties
c
      print *,'Nombre de champs dans la moyenne:',cnt
      fac = -1./cnt
      call sscal(nvaria,fac,tendancem,1)
      do nl=1,nlevels
         call spectogrid(gwrk,tendancem(1,nl))
         write(21)ii,((gwrk(i,j),i=1,nlong),j=1,nlat)
      enddo
      call dq2dpsi(ptend,tendancem)
      do nl=1,nlevels
         call spectogrid(gwrk,ptend(1,nl))
         write(21)ii,((gwrk(i,j),i=1,nlong),j=1,nlat)
      enddo
c
      stop
      end
