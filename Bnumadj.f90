        subroutine adjekman(sdisp,dq)
        include 'sphectra.h'
        include 'paramod.h'
        real sdisp(nvaria2,nlevels),sdisp1(nvaria2,nlevels),swrk(nvaria2)
        real gdisp(n2long,nlat),gwrk(n2long,nlat),dq(nvaria2,nlevels)
       real gwrk2(n2long,nlat),swrk2(nvaria2)
        call sset(nvaria,0.,sdisp,1)
       call sset(nvaria,0.,sdisp1,1)
!c
!ccc     call laplace(swrk,dq(1,nlevels),2)
!ccc       call spectogrid(gwrk,swrk)
!ccc     call shprod(n2long*nlat,gwrk,1,gek,1,gdisp,1)
!ccc     call gridtospec(sdisp1(1,nlevels),gdisp)
!cccc
!ccc       call sset(nvaria2,0.,swrk,1)
!ccc     call progras(swrk,sek,dq(1,nlevels))
!ccc       call saxpy(nvaria2,1.,swrk,1,sdisp1(1,nlevels),1)
!cc
!c adjoint of k lap(dpsi)
       call spectogrid(gwrk,dq(1,nlevels))
       call shprod(n2long*nlat,gwrk,1,gek,1,gdisp,1)
       call gridtospec(swrk,gdisp)
       call laplace(sdisp1(1,nlevels),swrk,2)

!c adjoint of grad(k).grad(dpsi)
       call laplace(swrk,sek,2)
       call spectogrid(gwrk2,swrk)
       call shprod(n2long*nlat,gwrk,1,gwrk2,1,gdisp,1)
       call gridtospec(swrk2,gdisp)
!c
       call progras(swrk,sek,dq(1,nlevels))
!c
       call saxpy(nvaria2,1.,swrk2,1,swrk,1)
!c
!c sum
       call saxpy(nvaria2,-1.,swrk,1,sdisp1(1,nlevels),1)
!cc
        call dq2dpsi(sdisp,sdisp1)
!cc
        return
        end

!c Routine qui calcule le modele lineaire tangent d (dq) / dt = F (q, dq );  q est
!c       la vorticite potentielle (en spectrale), dq la perturbation de q et F vient
!c       sur dqt

        subroutine linfun(q,dq,dqt)
        include 'sphectra.h'
        include 'paramod.h'
        real q(nvaria2,nlevels),dq(nvaria2,nlevels),dqt(nvaria2,nlevels),psi(nvaria2,nlevels)
        real dpsi(nvaria2,nlevels),dsdisp1(nvaria2,nlevels)
       real dsdisp2(nvaria2,nlevels),dsdisp3(nvaria2,nlevels)
        real dtzeta(nvaria2),zetap(nvaria2,nlevels),yj(nvaria2,nlevels)

        call q2psi(psi,q)
        call dq2dpsi(dpsi,dq) 
        call laplace(dtzeta,dpsi(1,nlevels),2)
        call zetap2q(zetap,dq)
        call sset(nvar,0.,dqt,1)

        do nl=1,nlevels
                call jacobs(dqt(1,nl),q(1,nl),dpsi(1,nl))
                call jacobs(yj(1,nl),dq(1,nl),psi(1,nl))
        enddo

        call dtherm(dsdisp1,dpsi)
        call dissel(dsdisp2,zetap)

       call dekman(dsdisp3(1,nlevels),dtzeta,dpsi(1,nlevels))

        call saxpy(nvaria,1.,yj,1,dqt,1)

        call saxpy(nvaria,-1.,dsdisp1,1,dqt,1)
        call saxpy(nvaria,-1.,dsdisp2,1,dqt,1)
        call saxpy(nvaria,-1.,dsdisp3,1,dqt,1)

        return
        end



!c Routine qui calcule le modele lineaire tangent adjoint d (dq) / -dt = F*(q, dq );  q est
!c       la vorticite potentielle (en spectrale), dq la perturbation de q et F vient
!c       sur dqt

        subroutine adjfun(q,dq,dqt)
        include 'sphectra.h'
        include 'paramod.h'
        real q(nvaria2,nlevels),dq(nvaria2,nlevels),dqt(nvaria2,nlevels),psi(nvaria2,nlevels)
        real dpsi(nvaria2,nlevels),dsdisp1(nvaria2,nlevels)
        real dsdisp2(nvaria2,nlevels),dsdisp3(nvaria2,nlevels)
        real zetap(nvaria2,nlevels),yj(nvaria2,nlevels),yj1(nvaria2,nlevels)
        real swrk(nvaria2,nlevels),swrk1(nvaria2)

        call q2psi(psi,q)
        call dq2dpsi(dpsi,dq)
        call zetap2q(zetap,dq)
        call sset(nvar,0.,dqt,1)

        call sset(nvar,0.,swrk,1)
        call laplace(swrk1,dq(1,nlevels),2)

        do nl=1,nlevels
                call jacobs(yj1(1,nl),dq(1,nl),q(1,nl))
                call jacobs(yj(1,nl),psi(1,nl),dq(1,nl))
        enddo
        call dq2dpsi(dqt,yj1)  

        call dtherm(dsdisp1,dpsi)
        call dissel(dsdisp2,zetap)

!c       call dekman(swrk(1,nlevels),swrk1,dq(1,nlevels))
!c       call dq2dpsi(dsdisp3,swrk) 
       call adjekman(dsdisp3,dq)
!c
        call saxpy(nvaria,1.,yj,1,dqt,1)
        call saxpy(nvaria,-1.,dsdisp1,1,dqt,1)
        call saxpy(nvaria,-1.,dsdisp2,1,dqt,1)
        call saxpy(nvaria,-1.,dsdisp3,1,dqt,1)

!c       call sscal(nvar,-1.,dqt,1)
        return
        end


