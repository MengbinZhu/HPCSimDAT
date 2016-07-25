!C     IMSL routines needed by the QG 3-layer model, rewritten to match IMSL API.

!C      SCOPY: Copy a vector X to a vector Y, both single-precision
!C      N:  Length of vectors X and Y
!C      SX: Real vector (input)
!C      INCX: Displacement between desired elements of SX
!C      SY: Real vector (output)
!C      INCY: displacement between desired elements of SY

      subroutine scopy (n,sx,incx,sy,incy)
      integer n,incx,incy
      real sx(1),sy(1)
      integer iy,ix,i

      if (n.le.0) return
      if ((incx.eq.incy).and.(incx.eq.1)) go to 10
      ix=1
      iy=1
      if (incx.lt.0) ix=1+(1-N)*incx
      if (incy.lt.0) iy=1+(1-N)*incy
      do 5 i=1,n
         sy(iy) = sx(ix)
         ix=ix+incx
         iy=iy+incy
 5    continue
      return
!c     Special case for incx=incy=1
 10   do 15 i=1,n
         sy(i) = sx(i)
 15   continue
      return
      end

!C     SSET: Set the components of a vector to a scalar, all single precision
!C     INCX must be greater than zero.
!C     N: Length of vector X
!C     SA: Real scalar
!C     SX: Real vecotr of length (N*incx) (input/output)
!C     INCX: Displacement between elements of SX

      subroutine sset(n,sa,sx,incx)
      integer n,incx
      real sa,sx(1)
      integer i

      do 10 i=1,n*incx,incx
         sx(i)=sa
 10   continue
      return
      end

!C      SSCAL: Multiply a vector by a scalar, single-precision: y = ay
!C      N:  Length of vector Y
!C      SA: Real scalar
!C      SX: Real vector (input/output)
!C      INCX: Displacement between desired elements of SX
!C      
      subroutine sscal (n,sa,sx,incx)
      integer n,incx
      real sa,sx(1)
      integer ix,i

      if (n.le.0) return
      if (incx.eq.1) go to 10
      ix=1
      if (incx.lt.0) ix=1+(1-n)*incx
      do 5 i=1,n
         sx(ix)=sa*sx(ix)
         ix=ix+incx
 5    continue
      return
 10   do 15 i=1,n
         sx(i)=sa*sx(i)
 15   continue
      return
      end

!C     SAXPY: Compute scalar times vector plus vector, y = ax+y,
!C     all single precision.
!C     N:    Length of X and Y (input)
!C     SA:   Real scalar (input)
!C     SX:   Real vector (input)
!C     INCX: Displacement between desired elements of SX
!C     SY:   Real vector (input/output)
!C     INCY: Displacement between desired elements of SY

      subroutine saxpy (n,sa,sx,incx,sy,incy)
      integer n,incx,incy
      real sa,sx(1),sy(1)
      integer iy,ix,i

      if ((n.le.0).or.(sa.eq.0.e0)) return
      if ((incx.eq.incy).and.(incx.eq.1)) go to 10
      ix=1
      iy=1
      if (incx.lt.0) ix=1+(1-N)*incx
      if (incy.lt.0) iy=1+(1-N)*incy
      do 5 i=1,n
         sy(iy) = sa*sx(ix)+sy(iy)
         ix=ix+incx
         iy=iy+incy
 5    continue
      return
!c     Special case for incx=incy=1
 10   do 15 i=1,n
         sy(i) = sa*sx(i)+sy(i)
 15   continue
      return
      end

!C  SDOT:  Single-precision dot product x*y
!C  N:  length of vectors X and Y (input)
!C  sx: Real vector (input)
!C  incx: Displacement between successive x(i) in sx
!C  sy: Real vector (input)
!C  incy: Displacement between successive y(i) in sy

      real function sdot (n,sx,incx,sy,incy)
      integer n,incx,incy
      real sx(1),sy(1)
      integer ix,iy,i

      sdot = 0.0e0
      if (n.le.0) return
      if ((incx.eq.incy).and.(incx.eq.1)) go to 10
      ix=1
      iy=1
      if (incx.lt.0) ix=1+(1-N)*incx
      if (incy.lt.0) iy=1+(1-N)*incy
      do 5 i=1,n
         sdot = sdot + sx(ix)*sy(iy)
         ix=ix+incx
         iy=iy+incy
 5    continue
      return
!c     special case for incx=incy=1
 10   do 15 i=1,n
         sdot = sdot + sx(i)*sy(i)
 15   continue
      return
      end

      
!C      SHPROD: Compute the Hadamard product of two single-prec vectors
!C              ( z(i) = x(i)*y(i) )
!C      N:  Length of vectors X, Y, and z
!C      SX: Real vector (input)
!C      INCX: Displacement between desired elements X(i) within SX (input)
!C      SY: Real vector (input)
!C      INCY: displacement between desired elements X(i) within SY (input)
!C      SZ: Real vector (output)
!C      INCZ: displacement between desired elements X(i) within SY (input)

      subroutine shprod (n,sx,incx,sy,incy,sz,incz)
      integer n,incx,incy,incz
      real sx(1),sy(1),sz(1)
      integer iy,ix,iz,i

      if (n.le.0) return
      if ((incx.eq.1).and.(incx.eq.1).and.(incz.eq.1)) go to 10
      ix=1
      iy=1
      iz=1
      if (incx.lt.0) ix=1+(1-N)*incx
      if (incy.lt.0) iy=1+(1-N)*incy
      if (incz.lt.0) iz=1+(1-N)*incz
      do 5 i=1,n
         sz(iz) = sx(ix)*sy(iy)
         ix=ix+incx
         iy=iy+incy
         iz=iz+incz
 5    continue
      return
!c     Special case for incx=incy=incz=1
 10   do 15 i=1,n
         sz(i) = sx(i)*sy(i)
 15   continue
      return
      end

!C MURRV: Multiply a real rectangular matrix by a vector.  Y = A*X
!C ****WARNING****  Don't try using this code as a replacement for IMSL's MURRV
!C outside the QG3L code.  This code *ignores* LDA and assumes the dimensions of
!C A are nra-by-nca, and assumes nx = nca and ny=nra, without checking them.  
!C It assumes IPATH=1 (that is, it won't compute Y=A'*X).
!C
!C NRA: Number of rows in A
!C NCA: Number of cols in A
!C A:   Real array
!C LDA: Ignored
!C NX:  Number of elements in X
!C X:   Real vector
!C NY:  Number of elements in Y
!C Y:   Real vector (output)
      subroutine murrv(nra,nca,a,lda,nx,x,ipath,ny,y)
      integer nra,nca,lda,nx,ny,ipath
      real a(nra,nca),x(nx),y(ny)
      integer icol,irow

      do 10 irow=1,nra
         y(irow)=0.0e0
         do 15 icol=1,nca
            y(irow)=y(irow)+a(irow,icol)*x(icol)
 15      continue
 10   continue
      end

!C LINRG: Compute the inverse of a real general matrix.
!C This is done using Gauss-Jordan elimination with full pivoting, as
!C described in "Numerical Recipes in Fortran", page 28-29.  I have removed
!C the parts of the code which solve Ax=b for a set of b vectors, keeping
!C only the inverse-computation parts.  In addition, the Recipes code
!C inverts a matrix in place, destroying the original matrix: thus, I use 
!C scopy to make a copy of the matrix before inverting.  As in MURRV, lda is
!C unused.
!
!C     N: Order of the matrix A.
!C     A: N by N matrix to be inverted
!C     LDA: Ignored
!C     AI: N by N matrix containing inverse(A) (output)
!C     LDAI: Ignored

      subroutine linrg(n,a,lda,ai,ldai)
      integer n,lda,ldai
      real a(n,n), ai(n,n)
      
      integer i,j,k,l,ll,nmax,ipiv,indxr,indxc,irow,icol
      real big,dum,pivinv
!C     See if this can be changed to a runtime variable:
      parameter (nmax=50)
      dimension ipiv(nmax),indxr(nmax),indxc(nmax)
      
!C     Copy the matrix before inverting it:
      call scopy(n*n,a,1,ai,1)

      do 11 j=1,n
         ipiv(j)=0
 11   continue
!C     This is the main loop over the columns to be reduced
      do 22 i=1,n
         big = 0.
!C     This is the outer loop of the search for a pivot element
         do 13 j=1,n
            if (ipiv(j).ne.1) then
               do 12 k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(ai(j,k)).ge.big) then
                        big=abs(ai(j,k))
                        irow=j
                        icol=k
                     endif
                  else if (ipiv(k).gt.1) then
                     pause 'Singular matrix'
                  endif
 12            continue
            endif
 13      continue
         ipiv(icol)=ipiv(icol)+1
!C     We now have the pivot element, so we interchange rows if needed to put
!C     the pivot element on the diagonal.  The columns are not physically
!C     interchanged, only relabeled: indx(i), the column of the ith pivot
!C     element, is the ith column that is reduced, while indxr(i) is the row
!C     in which that pivot element was originally located.  If indxr(i) !=
!C     indxc(i) there is an implied column interchnage.  This means the inverse
!C     matrix will be scrambled by columns.
         if (irow.ne.icol) then
            do 14 l=1,n
               dum=ai(irow,l)
               ai(irow,l) = ai(icol,l)
               ai(icol,l) = dum
 14         continue
         endif
         indxr(i)=icol
         indxc(i)=irow
!C     We are now ready to divide the pivot row by the pivot element, located
!C     at IROW and ICOL
         if (ai(icol,icol).eq.0) pause 'Singular matrix.'
         pivinv=1./ai(icol,icol)
         ai(icol,icol)=1.
         do 16 l=1,n
            ai(icol,l) = ai(icol,l)*pivinv
 16      continue
!C     Next we reduce the rows...
         do 21 ll=1,n
!C     ...except for the pivot one, of course.
            if (ll.ne.icol) then
               dum=ai(ll,icol)
               ai(ll,icol)=0.
               do 18 l=1,n
                  ai(ll,l)=ai(ll,l)-ai(icol,l)*dum
 18            continue
            endif
 21      continue
 22   continue
!C     This is the end of the main loop over columns of the reduction.
!C     It only remains to unscramble the slution in view of the column inter-
!C     changes.  We do this by interchanging pairs of columns in the reverse
!C     order that the permutation was built up.
      do 24 l=n,1,-1
         if (indxr(l).ne.indxc(l))then
            do 23 k=1,n
               dum=ai(k,indxr(l))
               ai(k,indxr(l))=ai(k,indxc(l))
               ai(k,indxc(l))=dum
 23         continue
         endif
 24   continue
!C     And we are done.
      return
      end
