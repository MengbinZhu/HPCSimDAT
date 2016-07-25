

!----------------------------------------------------
        !!!!form random error of normal distribution
        !!!s is the standard deviation
        !each member is from a normal distribution of mean is 0; standard deviation is s.
        subroutine error_N(s,n,x)
        include 'sphectra.h' 
   
        real::u1,u2
        integer::n,i
        real::x1,x(n)
        real::s
        do i=1,n
                call random_number(u1)
                call random_number(u2)
                x1=sqrt(-2*log(u1))*cos(2*pi*u2)
                x(i)=x1*s
        end do
        
        return
        end subroutine
        !-----------------------------------right



!!return random_mtrx(1:m,1:n), m is the order of Cov (or variance)
!!supply Cov if the random vector have correlationship with each other; 
!!otherwise supply variance is OK
subroutine genVGaussRandom(nl,nall,mean,Cov,ne,random_mtrx)
implicit none
    integer :: nl,ne,nall,nlatlong,i,j,k,p
    real :: Cov(nl,nl)
    real :: mean
    real :: random_1grid(nl)
    real :: random_mtrx(nall,ne)
    real,allocatable :: random_gs01(:)
    integer :: info
    real,allocatable :: Lmtrx(:,:),EVal(:)
    nlatlong=nall/nl

        allocate(Lmtrx(nl,nl),EVal(nl)) 
        Lmtrx=Cov
        
        !call syevd(Lmtrx,EVal,jobz='V',uplo='L',info=info)
        call eof(nl,Lmtrx,Eval)
        p= count(EVal>1e-5)
      !  print*,Lmtrx
      !  print*,'-------------'
      !  print*,EVal
        do i=1,p    !m,m-p+1,-1
            Lmtrx(:,i)=Lmtrx(:,i)*sqrt(EVal(i))
        end do
        allocate(random_gs01(p))
       
        do j=1,nlatlong
            do i=1,ne
                call error_N(1.,p,random_gs01)
                random_1grid(:)=matmul(Lmtrx(:,1:p),random_gs01)+mean
                
                do k=1,nl
                random_mtrx(j+(k-1)*nlatlong,i)=random_1grid(k)
                end do

            end do
        end do
        deallocate(random_gs01)
        deallocate(Lmtrx,EVal)
        return
end subroutine
