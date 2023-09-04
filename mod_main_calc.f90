module mod_main_calc
    use mod_autodiff
    use mod_condition
    use mod_cubic_eq_d
    use mod_input
    use mod_fugacity

    implicit none

    integer,private::i,j

    contains
    subroutine main_calc(V0,lnk0,Nc0,Nc0old,Nm0,Nm0old,P0,P0old,Pb0,q_judge,phase_judge,g)
        implicit none
        integer,intent(inout)::q_judge,phase_judge(n)
        real(8),intent(inout)::Pb0
        real(8),intent(inout),dimension(n)::V0,P0,P0old
        real(8),intent(inout),dimension(com_2phase,n)::lnk0
        real(8),intent(inout),dimension(com_2phase+com_ion,n)::Nc0,Nc0old
        real(8),intent(inout),dimension(com_mine,n)::Nm0,Nm0old
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::P(:),V(:),lnk1(:),lnk2(:),lnk3(:),lnk4(:),Nc1(:),Nc2(:),Nc3(:),Nc4(:),Nc5(:),Nc6(:),Nc7(:)
        type(diffs),pointer::Nc8(:),Nc9(:),Nc10(:),Nc11(:),Nc12(:),Nc13(:),Nc14(:),Nm1(:),Nm2(:),Nm3(:),Nm4(:),Nm5(:),Pb
        real(8),allocatable,dimension(:)::x0,kakuninn
        real(8)::kaku

        !!
        real(8)::z_factor0(n),v_L2(n),v_L3(n),v_L4(n)
        type(diffs)::lnk(com_2phase,n),Nc(com_2phase+com_ion,n),Nm(com_mine,n),Nt(n),z(com_2phase+com_ion,n),rach(n)
        type(diffs)::kakuninnnnnn(com_2phase+com_ion),x(com_2phase+com_ion,n),y(com_2phase,n),k(com_2phase,n)
        type(diffs)::lnfai_L(com_2phase,n),lnfai_V(com_2phase,n),z_factor(n),v_L1(n),MV_V(n),MD_V(n)
        type(diffs)::v_L(com_2phase+com_ion,n),MV_L(n),MD_L(n),Sw(n),Sg(n),L(n)
        

        allocate(x0(n*eq+q_judge))

        !!自動微分の下準備
        do i=1,n
            !?相平衡
            x0(i*eq-24)=lnk0(1,i)
            x0(i*eq-23)=lnk0(2,i)
            x0(i*eq-22)=lnk0(3,i)
            x0(i*eq-21)=lnk0(4,i)
            !?物質収支
            x0(i*eq-20)=Nc0(1,i)
            x0(i*eq-19)=Nc0(2,i)
            x0(i*eq-18)=Nc0(3,i)
            x0(i*eq-17)=Nc0(4,i)
            x0(i*eq-16)=Nc0(5,i)
            x0(i*eq-15)=Nc0(6,i)
            x0(i*eq-14)=Nc0(7,i)
            x0(i*eq-13)=Nc0(8,i)
            x0(i*eq-12)=Nc0(9,i)
            x0(i*eq-11)=Nc0(10,i)
            x0(i*eq-10)=Nc0(11,i)
            x0(i*eq-9)=Nc0(12,i)
            x0(i*eq-8)=Nc0(13,i)
            x0(i*eq-7)=Nc0(14,i)
            x0(i*eq-6)=Nm0(1,i)
            x0(i*eq-5)=Nm0(2,i)
            x0(i*eq-4)=Nm0(3,i)
            x0(i*eq-3)=Nm0(4,i)
            x0(i*eq-2)=Nm0(5,i)
            !?飽和率の制約式
            x0(i*eq-1)=P0(i)
            !?rachford-rice
            x0(i*eq)=V0(i)
        end do
        if (q_judge == 1) then
            x0(n*eq+q_judge)=Pb0
        end if

        call diffsset1(x0,xd)
        call sizeset(x0,g)

        lnk1 => xd(eq-24:n*eq-24:eq)
        lnk2 => xd(eq-23:n*eq-23:eq)
        lnk3 => xd(eq-22:n*eq-22:eq)
        lnk4 => xd(eq-21:n*eq-21:eq)
        Nc1 => xd(eq-20:n*eq-20:eq)
        Nc2 => xd(eq-19:n*eq-19:eq)
        Nc3 => xd(eq-18:n*eq-18:eq)
        Nc4 => xd(eq-17:n*eq-17:eq)
        Nc5 => xd(eq-16:n*eq-16:eq)
        Nc6 => xd(eq-15:n*eq-15:eq)
        Nc7 => xd(eq-14:n*eq-14:eq)
        Nc8 => xd(eq-13:n*eq-13:eq)
        Nc9 => xd(eq-12:n*eq-12:eq)
        Nc10 => xd(eq-11:n*eq-11:eq)
        Nc11 => xd(eq-10:n*eq-10:eq)
        Nc12 => xd(eq-9:n*eq-9:eq)
        Nc13 => xd(eq-8:n*eq-8:eq)
        Nc14 => xd(eq-7:n*eq-7:eq)
        Nm1 => xd(eq-6:n*eq-6:eq)
        Nm2 => xd(eq-5:n*eq-5:eq)
        Nm3 => xd(eq-4:n*eq-4:eq)
        Nm4 => xd(eq-3:n*eq-3:eq)
        Nm5 => xd(eq-2:n*eq-2:eq)
        P => xd(eq-1:n*eq-1:eq)
        V => xd(eq:n*eq:eq)
        if (q_judge == 1) then
            Pb => xd(n*eq+q_judge)
        end if
        !?-------------------------------

        do i=1,n
            lnk(1,i)=lnk1(i)
            lnk(2,i)=lnk2(i)
            lnk(3,i)=lnk3(i)
            lnk(4,i)=lnk4(i)
            Nc(1,i)=Nc1(i)
            Nc(2,i)=Nc2(i)
            Nc(3,i)=Nc3(i)
            Nc(4,i)=Nc4(i)
            Nc(5,i)=Nc5(i)
            Nc(6,i)=Nc6(i)
            Nc(7,i)=Nc7(i)
            Nc(8,i)=Nc8(i)
            Nc(9,i)=Nc9(i)
            Nc(10,i)=Nc10(i)
            Nc(11,i)=Nc11(i)
            Nc(12,i)=Nc12(i)
            Nc(13,i)=Nc13(i)
            Nc(14,i)=Nc14(i)
            Nm(1,i)=Nm1(i)
            Nm(2,i)=Nm2(i)
            Nm(3,i)=Nm3(i)
            Nm(4,i)=Nm4(i)
            Nm(5,i)=Nm5(i)
            do j=1,com_2phase
                k(j,i) = exp(lnk(j,i))
            end do
        end do

        !!rachford-rice
        !#TODO相の数で判断するか否か
        do i=1,n
            call residualvectorset3(n*eq+q_judge,Nt(i))
            do j=1,com_2phase+com_ion
                Nt(i) = Nt(i) + Nc(j,i)
            end do
            do j=1,com_2phase+com_ion
                z(j,i) = Nc(j,i) /Nt(i)
            end do
            call residualvectorset3(n*eq+q_judge,rach(i))
            do j=1,com_2phase
                rach(i) = rach(i) +(1.0d0-k(j,i))*z(j,i)/(1.0d0-V(i)+V(i)*k(j,i))
            end do
        end do
        !do j=1,com_2phase+com_ion
        !    kakuninnnnnn(j) = z(j,1)
        !end do
        !call outxs(kakuninnnnnn,kakuninn)
        !write(*,*) kakuninn

        !!モル分率----------------------
        do i=1,n
            if (phase_judge(i) == 2) then
                do j=1,com_2phase
                    x(j,i) = z(j,i)/(1.0d0-V(i)+V(i)*k(j,i))
                    y(j,i) = x(j,i)*k(j,i)
                end do
                do j=com_2phase+1,com_2phase+com_ion
                    x(j,i) = z(j,i)/(1.0d0-V(i))
                end do
            else
                do j=1,com_2phase
                    call residualvectorset3(n*eq+q_judge,y(j,i))
                end do
                do j=1,com_2phase+com_ion
                    x(j,i) = z(j,i)
                end do
            end if
        end do
        !do j=1,com_2phase+com_ion
        !    kakuninnnnnn(j) = x(j,1)
        !end do
        !call outxs(kakuninnnnnn,kakuninn)
        !write(*,*) kakuninn

        !!フガシティ
        call vapor_fugacity_main(y,lnfai_V,P,z_factor0,z_factor,q_judge,phase_judge)
        call liquid_fugacity_main(lnfai_L,P,v_L1,v_L2,v_L3,v_L4)


        
        !!モル密度、体積
        !?気相
        do i=1,n
            if (phase_judge(i) == 2) then
                MD_V(i) = P(i)/z_factor(i)/R/temp
                MV_V(i) = 1.0d0/MD_V(i)
            else
                call residualvectorset3(n*eq+q_judge,MD_V(i))
                call residualvectorset3(n*eq+q_judge,MV_V(i))
            end if
        end do
        !?液相
        do i=1,n
            v_L(1,i) = v_L1(i)
            call residualvectorset4(v_L2(i),n*eq+q_judge,v_L(2,i))
            call residualvectorset4(v_L3(i),n*eq+q_judge,v_L(3,i))
            call residualvectorset4(v_L4(i),n*eq+q_judge,v_L(4,i)) 
            do j=com_2phase+1,com_2phase+com_ion
                call residualvectorset3(n*eq+q_judge,v_L(j,i))!!要検討====================
            end do
            call residualvectorset3(n*eq+q_judge,MV_L(i))
            do j=1,com_2phase+com_ion
                MV_L(i) = MV_L(i) + x(j,i)*v_L(j,i)
            end do
            MD_L(i) = 1.0d0/MV_L(i)
        end do

        !call outxs(MD_V,kakuninn)
        !write(*,*) kakuninn


        !!飽和率の算出
        do i=1,n
            L(i)=1.0d0-V(i)
            Sw(i)=Nt(i)*L(i)*MV_L(i)
            if(phase_judge(i) == 2) then
                Sg(i)=Nt(i)*V(i)*MV_V(i)
            else
                call residualvectorset3(n*eq+q_judge,Sg(i))
            end if
        end do
        














        
        do i=1,n
            do j=1,com_2phase
                g(i*eq-eq+j)=lnk(j,i)+lnfai_V(j,i)-lnfai_L(j,i)
            end do
            do j=1,com_2phase+com_ion
                g(i*eq-eq+j+com_2phase)=Nc(j,i)
            end do
            do j=1,com_mine
                g(i*eq-eq+j+com_2phase+com_2phase+com_ion)=Nm(j,i)
            end do

            
            g(i*eq-1)=P(i)!Sw(i)+Sg(i)-1.0d0
            g(i*eq)=rach(i)
        end do

        !!流量制御
        if (q_judge == 1) then
            g(n*eq+q_judge) = Pb
        end if

        call outxs(rach,kakuninn)
        write(*,*) kakuninn
        
    end subroutine

end module
