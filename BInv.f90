! =======  ==========    =========
!2011/10/19     Chen Lei  Original code
!calculate the inverse of a matrix   A(original matrix)=>A(inverse matrix of A)
      Module my_inverse
      Implicit none
      Integer, Parameter, private :: dbl=Selected_Real_Kind(p=4)
      Contains 
      Subroutine inverse(A,flag)
      Implicit none
      Integer :: i , j, ipivot, jpivot, m, k  !cycle standard
      Integer :: istat  ! Status, istat=0 for success
      integer, Intent(out) :: flag   ! flag=1,failure
      Real(Kind=DBL), Intent(inout), Dimension(:,:) :: A !
      Real, Allocatable, Dimension(:,:) :: P, Q 
      Real (Kind=DBL), Dimension(:), Allocatable :: Temp
      Real (Kind=DBL) :: pivot
      pivot =0.0_dbl
      m=Size(A,1)
      flag=0
                If  (m==1.or. m/=Size(A,2)) then
                     write(*,*)'A is not a square matrix or is a scalar'
                     flag=1                     
                     return
                End If              
                Allocate (P(m,m), Q(m,m),Temp(m),Stat=istat)
      P=0.0_DBL
      Q=0.0_DBL
      Forall ( i=1:m:1) 
        P(i,i)=1.0_DBL
        Q(i,i)=1.0_DBL
      End Forall                              
      main: Do k=1,m-1,1
                   pivot =0._DBL                    
                   row:   Do i=k,m,1   
                       column:  Do j=k,m,1                                                 
                                        If (abs(A(i,j))>abs(pivot)) then                   
                                                pivot=A(i,j)            
                                                ipivot=i             
                                                jpivot=j            
                                         End If                                                                    
                                       End do column                     
                                End do row                 
                       If (abs(pivot)<=1.0E-10)  then                                            
                           !write(*,*) 'A is a degraded matirix or at least a very ill-conditional'
                           flag=1
                           return
                       End If                                                                                                  
                       Call shift 
                       Call update 
                End do main                               
                If (abs(A(m,m))>1.0E-10) Then                               
                   p(m,:)=p(m,:)/A(m,m)  
                   A=Matmul(Q,P) 
                Else                               
                   !write(*,*) 'A is a degraded matirix or at least a very ill-conditional' 
                   flag=1
                   return
                End If 
                Deallocate(P,Q,Stat=istat)                                    

            
            Contains                                     
                Subroutine shift
                Implicit none
                   If (k/=ipivot) Then
                       Temp=A(k,:)                                 
                       A(k,:)=A(ipivot,:)                                   
                       A(ipivot,:)=Temp
                       Temp=P(k,:)
                       P(k,:)=P(ipivot,:)                               
                       P(ipivot,:)=Temp 
                   End If
                   If (k/=jpivot) Then
                      Temp=A(:,k)
                      A(:,k)=A(:,jpivot)
                      A(:,jpivot)=Temp
                      Temp=Q(:,k)                                         
                      Q(:,k)=Q(:,jpivot)
                      Q(:,jpivot)=Temp
                   End if   
                End Subroutine shift 
            
               Subroutine update
               Implicit none
                   Real (Kind=DBL) :: factor
                   Integer :: iupdate                                       
                   factor=1.0_DBL/pivot                                     
                   Forall ( iupdate=k+1:m:1)                                       
                   A(iupdate,k+1:m:1)=A(iupdate,k+1:m:1)-&
                   factor*A(iupdate,k)*A(k,k+1:m:1)                                  
                   P(iupdate,:)=P(iupdate,:)-factor*A(iupdate,k)*P(k,:)                                
                   Q(:,iupdate)=Q(:,iupdate)-factor*A(k,iupdate)*Q(:,k)                                
                   End Forall                                      
                   p(k,:)=factor*P(k,:)                                     
            End subroutine update                                       
        End Subroutine inverse
            
            
       End Module  My_inverse       
