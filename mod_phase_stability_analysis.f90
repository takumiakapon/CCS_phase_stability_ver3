module mod_phase_stability_analysis
    use mod_autodiff
    use mod_condition
    use mod_input
    use mod_fugacity
    use mod_cubic_eq_d
    
    implicit none
    
    integer,private::i,j
    
    contains
    
    subroutine phase_stability_liquid(alpha0,P0,z0,g)
        !!設定用
        implicit none
        real(8),intent(inout),dimension(com_2phase)::alpha0
        real(8),intent(in)::P0,z0(com_2phase+com_ion)
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::alpha(:)
        real(8),allocatable,dimension(:)::x0
        !!計算用
        type(diffs),dimension(com_2phase)::x,y,w,lnfai_V,lnfai_L
        type(diffs)::wt
        
        real(8),allocatable,dimension(:)::kakuninn
        
        
        allocate(x0(4))
        
        !!自動微分の下準備
        do i=1,com_2phase
            x0(i) = alpha0(i)
        end do
        
        call diffsset1(x0,xd)
        call sizeset(x0,g)
        
        alpha => xd(1:com_2phase)
        call outxs(alpha,kakuninn)
        !write(*,*) kakuninn
        
        do i=1,com_2phase
            w(i) = (alpha(i) / 2.0d0) ** 2.0d0
        end do
        
        
        call outxs(w,kakuninn)
        !write(*,*) kakuninn
        !!気相のモル分率
        call residualvectorset3(com_2phase,wt)
        do i=1,com_2phase
            wt = wt + w(i)
        end do
         call outxs(alpha,kakuninn)
        !write(*,*) kakuninn
        do i=1,com_2phase
            y(i) = w(i) /wt
        end do
        
        !!液相のモル分率
        do i=1,com_2phase
            call residualvectorset4(z0(i),com_2phase,x(i))
        end do
       
        call vapor_fugacity(y,lnfai_V,P0)
        
        call liquid_fugacity(lnfai_L,P0)
        call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn
        
        
        
        do i=1,com_2phase
            g(i) = sqrt(w(i)) * (log(w(i)) + lnfai_V(i) -log(x(i)) - lnfai_L(i))
        end do
        
    end subroutine
    
    
    subroutine phase_stability_vapor(alpha0,P0,z0,g)
        implicit none
        !!設定用
        real(8),intent(inout),dimension(com_2phase)::alpha0
        real(8),intent(in)::P0,z0(com_2phase+com_ion)
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::alpha(:)
        real(8),allocatable,dimension(:)::x0
        
        !!計算用
        type(diffs),dimension(com_2phase)::x,y,w,lnfai_L,lnfai_V
        type(diffs)::wt
        real(8),allocatable,dimension(:)::kakuninn
        
        allocate(x0(com_2phase))
        
        !!自動微分の下準備
        do i=1,com_2phase
            x0(i) = alpha0(i)
        end do
        
        call diffsset1(x0,xd)
        call sizeset(x0,g)
        
        alpha => xd(1:com_2phase)
        
        
        do i=1,com_2phase
            if (alpha0(i) == 0.0d0) then
                call residualvectorset3(com_2phase,w(i))
            else
                w(i) = (alpha(i) / 2.0d0) ** 2.0d0
            end if
        end do
        
        call residualvectorset3(com_2phase,wt)
        do i=1,com_2phase
            wt = wt + w(i)
        end do
        !!液相のモル分率
        do i =1,com_2phase
            x(i) = w(i) / wt
        end do
        
        !気相のモル分率
        do i=1,com_2phase
            call residualvectorset4(z0(i),com_2phase,y(i))
        end do
        
        call vapor_fugacity(y,lnfai_V,P0)
        call liquid_fugacity(lnfai_L,P0)
        
        call outxs(lnfai_L,kakuninn)
        !write(*,*) kakuninn
        
        
        
        do i=1,com_2phase
            g(i) = sqrt(w(i)) * (log(w(i)) + lnfai_L(i) -log(y(i)) - lnfai_V(i))
        end do
        
        
    end subroutine
    
    
end module