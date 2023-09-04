module mod_ini_flash
    use mod_autodiff
    use mod_condition
    use mod_input
    use mod_fugacity
    use mod_cubic_eq_d

    implicit none
    integer,private::i,j

    contains

    subroutine ini_fla(P,z0,lnk0,v0,z_factor0,g)
        implicit none
        real(8),intent(inout)::P,V0,z_factor0
        real(8),intent(inout),dimension(com_2phase)::lnk0
        real(8),intent(inout),dimension(com_2phase+com_ion)::z0
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::lnk(:),V
        real(8),allocatable,dimension(:)::x0
        
        
        type(diffs)::fra
        type(diffs),dimension(com_2phase)::k,y,lnfai_V,lnfai_L
        type(diffs),dimension(com_2phase+com_ion)::z,x
        
        real(8),allocatable,dimension(:)::kakuninn
        real(8)::kaku
        
        
        allocate(x0(com_2phase+1))
        
        !!自動微分の下準備----------------------------
        do i=1,com_2phase
            x0(i)=lnk0(i)
        end do
        x0(com_2phase+1) = v0

        call diffsset1(x0,xd)
        call sizeset(x0,g)

        lnk => xd(1:com_2phase)
        V => xd(com_2phase+1)
        !?-----------------------------------------
        call outxs(lnk,kakuninn)
        !write(*,*) kakuninn
        
        
        !!rachford-rice
        do i=1,com_2phase+com_ion
            call residualvectorset4(z0(i),com_2phase+1,z(i))
        end do
        call residualvectorset3(com_2phase+1,fra)
        do i=1,com_2phase
            fra=fra+(1.0d0-exp(lnk(i)))*z(i)/(1.0d0-V+V*exp(lnk(i)))
            call out_diffsx(fra,kaku)
            !write(*,*) kaku
        end do
        !do i=com_2phase+1,com_2phase+com_ion
        !    fra = fra+z(i)/(1.0d0-V)   
        !end do
        
        
            
        do i=1,com_2phase
            k(i)=exp(lnk(i))
            x(i)=z(i)/(1.d0-V+V*k(i))
            y(i)=x(i)*k(i)
        end do
        do i=com_2phase+1,com_2phase+com_ion
            x(i)=z(i)/(1.0d0-V)
        end do
        
        call vapor_fugacity_ini(y,lnfai_V,P,z_factor0)
        
        call liquid_fugacity_ini(lnfai_L,P)
        
        
        do i=1,com_2phase
            g(i) = lnk(i)+lnfai_V(i)-lnfai_L(i)
        end do 
        g(com_2phase+1) = fra
        call out_diffsx(fra,kaku)
        !write(*,*) kaku
        
        !!相安定解析と初期フラッシュでは主要変数の数が違うから、それを「mod_fugacity」では考慮できていない→主要変数の数も引数にしちゃう？
        
        
        
        
    end subroutine
    
    
    end module
    