      subroutine estatlatm(dum,x)
!c
      real x(66,32),dum(22,9)
!c
      do i=51,64
         do j=5,13
            dum(i-50,j-4)=x(i,j)
         enddo
      enddo
      do i=1,8
         do j=5,13
            dum(i+14,j-4)=x(i,j)
         enddo
      enddo
!c
      return
      end
!c
!c
      subroutine estpacatm(dum,x)
!c
      real x(66,32),dum(25,9)
!c
      do i=22,46
         do j=5,13
            dum(i-21,j-4)=x(i,j)
         enddo
      enddo
!c
      return
      end      
!c

