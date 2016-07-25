      subroutine initsph
!c
!c  Initialisations
!c     
      character*80 qegnfil
!c
      qegnfil='./PLYLG21'
      call initctrl
      call inittrsf(qegnfil)
      call gausstrig
!c
      return
      end
!c
      subroutine spectogridh(grf,spf)
!c
!c  Passage en points de grille hemispherique antisymetrique
!c
      include 'sphectra.h'
!c
      real grf(n2long,nlat)
      real spf(ncplx,nbdeglib)
!c
      call tlginvh (grf,spf,plm,1)
      call fft991(grf,work,trigs,ifax,1,n2long,nlong,nhlat,1)
!c
      return
      end
!c
      subroutine gridtospech(spf,grf)
!c
!c  Passage en spectral hemispherique antisymetrique
!c
      include 'sphectra.h'
!c
      real grf(n2long,nlat),gwork(n2long,nlat)
      real spf(ncplx,nbdeglib)
!c
      call scopy(n2long*nlat,grf,1,gwork,1)
      call fft991(gwork,work,trigs,ifax,1,n2long,nlong,nhlat,-1)
      call tlgdirh (spf,gwork,pw,1)
!c
      return
      end
!c
      subroutine spectogridsh(grf,spf)
!c
!c  Passage en points de grille hemispherique symetrique
!c
      include 'sphectra.h'
!c
      real grf(n2long,nlat)
      real spf(ncplx,nbdeglib)
!c
      call tlginvsh (grf,spf,plm,1)
      call fft991(grf,work,trigs,ifax,1,n2long,nlong,nhlat,1)
!c
      return
      end
!c
      subroutine gridtospecsh(spf,grf)
!c
!c  Passage en spectral hemispherique symetrique
!c
      include 'sphectra.h'
!c
      real grf(n2long,nlat),gwork(n2long,nlat)
      real spf(ncplx,nbdeglib)
!c
      call scopy(n2long*nlat,grf,1,gwork,1)
      call fft991(gwork,work,trigs,ifax,1,n2long,nlong,nhlat,-1)
      call tlgdirsh (spf,gwork,pw,1)
!c
      return
      end
!c
      subroutine spectogrid(grf,spf)
!c
!c  Passage en points de grille spherique
!c
      include 'sphectra.h'
!c
      real grf(n2long,nlat)
      real spf(ncplx,nbdeglib)
!c
      call tlginvs (grf,spf,plm,1)
      call fft991(grf,work,trigs,ifax,1,n2long,nlong,nlat,1)
!c
      return
      end
!c!
      subroutine gridtospec(spf,grf)
!c!
!c!  Passage en spectral spherique
!c!
      include 'sphectra.h'
!c!
      real grf(n2long,nlat),gwork(n2long,nlat)
      real spf(ncplx,nbdeglib)
!c!
      call scopy(n2long*nlat,grf,1,gwork,1)
      call fft991(gwork,work,trigs,ifax,1,n2long,nlong,nlat,-1)
      call tlgdirs (spf,gwork,pw,1)
!c!
      return
      end subroutine
!c
      real function specscal(a,b)
!c
!c  calcule le produit scalaire entre a et b pour a,b en spectrale
!c  Le produit scalaire est la MOYENNE du produit des champs sur la
!c  sphere.
!c
      include 'sphectra.h'
!c
      real a(ncplx,nbdeglib),b(ncplx,nbdeglib)
!c
      sc=2.*sdot(ncplx*nbdeglib,a,1,b,1)
      p=0.
!c
      do k=1,mvalg(0)
         l = mvadeb(0) + k - 1
         p = a(1,l)*b(1,l) + p
      enddo
!c
      do k=1,mvslg(0)
         l = mvsdeb(0) + k - 1
         p = a(1,l)*b(1,l) + p
      enddo
!c
      specscal = sc - p
!c
      return
      end
