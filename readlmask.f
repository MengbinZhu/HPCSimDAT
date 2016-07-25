      PROGRAM EXTFOR 
C
c  Lecture du Masque pour la proportion Terre/ocean
c
      include 'sphectra.h'
      PARAMETER (NGP=NLONg*nlat)
C
      REAL WG1(nlat,nlong),wg2(nlong,nlat)
c
      open(11,file='/d12/qgmod/mask_ggrid')
      open(21,file='/d12/qgmod/lmask.sus',form='unformatted')
c
      DO 10 L=1,1
         read(11,100) WG1
         do i=1,nlong
            do j=1,nhlat
               wg2(i,j) = wg1(j+nhlat,i)
               wg2(i,j+nhlat) = wg1(nhlat-j+1,i)
            enddo
         enddo
         ii = 0
         write(21)ii,wg2
10    continue
100   FORMAT (10f8.6)     
c
      stop
      end
