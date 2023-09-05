module mod_main_calc
    use mod_autodiff
    use mod_condition
    use mod_cubic_eq_d
    use mod_input
    use mod_fugacity
    use mod_injection_data

    implicit none

    integer,private::i,j

    contains
    subroutine main_calc(V0,lnk0,Nc0,Nc0old,Nm0,Nm0old,Nmini,P0,P0old,Pb0,fai0,q_judge,phase_judge,phase,Swd,krgd,krwd,g)
        implicit none
        integer,intent(inout)::q_judge,phase_judge(n),phase(n)
        real(8),intent(inout)::Pb0,Nmini(com_mine)
        real(8),intent(inout),dimension(n)::V0,P0,P0old,fai0
        real(8),intent(inout),dimension(com_2phase,n)::lnk0
        real(8),intent(inout),dimension(com_2phase+com_ion,n)::Nc0,Nc0old
        real(8),intent(inout),dimension(com_mine,n)::Nm0,Nm0old
        real(8),intent(in),dimension(21)::Swd,krgd,krwd
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::P(:),V(:),lnk1(:),lnk2(:),lnk3(:),lnk4(:),Nc1(:),Nc2(:),Nc3(:),Nc4(:),Nc5(:),Nc6(:),Nc7(:)
        type(diffs),pointer::Nc8(:),Nc9(:),Nc10(:),Nc11(:),Nc12(:),Nc13(:),Nc14(:),Nm1(:),Nm2(:),Nm3(:),Nm4(:),Nm5(:),Pb
        real(8),allocatable,dimension(:)::x0,kakuninn
        real(8)::kaku

        !!
        real(8)::z_factor0(n),v_L2(n),v_L3(n),v_L4(n),Mw(com_2phase+com_ion)
        type(diffs)::lnk(com_2phase,n),Nc(com_2phase+com_ion,n),Nm(com_mine,n),Nt(n),z(com_2phase+com_ion,n),rach(n)
        type(diffs)::kakuninnnnnn(com_2phase+com_ion),x(com_2phase+com_ion,n),y(com_2phase,n),k(com_2phase,n)
        type(diffs)::lnfai_L(com_2phase,n),lnfai_V(com_2phase,n),z_factor(n),v_L1(n),MV_V(n),MD_V(n)
        type(diffs)::v_L(com_2phase+com_ion,n),MV_L(n),MD_L(n),Sw(n),Sg(n),L(n),Mw_ave_L(n),Mw_ave_V(n)
        type(diffs)::phase_d_L(n),phase_d_V(n)

        real(8)::beta_v,myu_H2O_20,sa,sisuu,myu_L_normal,Vc(com_2phase),Tc(com_2phase),Pc(com_2phase),myu_c(com_2phase)
        real(8)::av_0,av_1,av_2,av_3,av_4
        type(diffs)::myu_L(n),Pc_ave_V(n),Tc_ave_V(n),Vc_ave_V(n),zeta(n),ro_r_V(n),myu_V_SC(n),myu_V(n),krg(n),krw(n)
        real(8),allocatable,dimension(:)::Sw0
        type(diffs)::faimine(n),fai(n),kk(n),ramda
        real(8)::faimineold(n),faiold(n),NmMD(com_mine)
        
        real(8)::re,W(n)
        type(diffs)::q,MD_injection_d,T_L(com_2phase+com_ion,n),T_V(com_2phase,n)

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
        call outxs(Sw,Sw0)
        
        
        !!相質量密度[kg/m^3]
        Mw(1)=18.01528d0 !モル質量[g/mol]
        Mw(2)=44.01d0
        Mw(3)=16.04d0
        Mw(4)=32.082d0
        Mw(5)=1.0079d0
        Mw(6)=61.0171d0
        Mw(7)=60.00092d0
        Mw(8)=17.0073d0
        Mw(9)=40.078d0
        Mw(10)=24.0d0
        Mw(11)=60.0d0
        Mw(12)=28.9799d0
        Mw(13)=33.0735d0
        Mw(14)=32.07d0

        do i=1,n
            call residualvectorset3(n*eq+q_judge,Mw_ave_L(i))
            call residualvectorset3(n*eq+q_judge,Mw_ave_V(i))
            do j=1,com_2phase+com_ion
                Mw_ave_L(i)=Mw_ave_L(i)+Mw(i)*x(j,i)
            end do
            do j=1,com_2phase
                Mw_ave_V(i)=Mw_ave_V(i)+Mw(i)+y(j,i)
            end do
            phase_d_L(i)=Mw_ave_L(i)*MD_L(i)*10.0d0**(-3.0d0) ![kg/m^3]
            phase_d_V(i)=Mw_ave_V(i)*MD_V(i)*10.0d0**(-3.0d0)
        end do

        !!相粘度
        !?液相
        beta_v=-1.297d0+(0.574*10.0d0**(-1.0d0))*(temp-273.15)-(0.697d0*10.0d0**(-3.0d0))*(temp-273.15)**2.0d0+&
            (0.447d0*10.0d0**(-5.0d0))*(temp-273.15)**3.0d0-(0.105d0*10.0d0**(-7.0d0))*(temp-273.15)**4.0d0 ![1/GPa]
        myu_H2O_20=1002.d0 !20[℃] ![μPa・s]
        sa=20.0d0+273.15d0-temp
        sisuu=(sa*(1.2378d0-(1.303d0*10.0d0**(-3.0d0))*sa+&
            (3.06d0*10.0d0**(-6.0d0))*sa**2.0d0+&
            (2.55d0*10.0d0**(-8.0d0))*sa**3.0d0)/(96.0d0+temp-273.15d0))
        myu_L_normal=myu_H2O_20*10.0d0**sisuu ![μPa・s]
        do i=1,n
            myu_L(i)=(myu_L_normal*10.0d0**(-6.0d0))*(1.0d0+beta_v*(P(i)*10.0d0**(-9.0d0))) ![Pa・s]
        end do

        !?気相
        !vc臨界体積[m^3/mol],1:H2O
        Vc(1)=56.3d0*10.0d0**(-6.0d0)
        Vc(2)=94.0d0*10.0d0**(-6.0d0)
        Vc(3)=99.0d0*10.0d0**(-6.0d0)
        Vc(4)=98.5d0*10.0d0**(-6.0d0)
        Tc(1) = 647.30d0
        Tc(2) = 304.2d0
        Tc(3) = 190.6d0
        Tc(4) = 373.2d0
        Pc(1) = 217.6d0 * 101325.d0
        Pc(2) = 72.8d0 * 101325.d0
        Pc(3) = 45.4d0 * 101325.d0
        Pc(4) = 88.2d0 * 101325.d0
        
        do i=1,n
            call residualvectorset3(n*eq+q_judge,Pc_ave_V(i))
            call residualvectorset3(n*eq+q_judge,Tc_ave_V(i))
            call residualvectorset3(n*eq+q_judge,Vc_ave_V(i))
            if (phase_judge(i) == 2) then
                do j=1,com_2phase
                    Pc_ave_V(i)=Pc_ave_V(i)+Pc(j)/101325.0d0*y(j,i) ![Pa→atm]
                    Tc_ave_V(i)=Tc_ave_V(i)+Tc(j)*y(j,i) ![K]
                    Vc_ave_V(i)=Vc_ave_V(i)+Vc(j)*y(j,i) ![m^3/mol]
                end do
                zeta(i)=Tc_ave_V(i)**(1.0d0/6.0d0)/((Mw_ave_V(i)**(1.0d0/2.0d0))*(Pc_ave_V(i)**(2.0d0/3.0d0)))
                ro_r_V(i)=MD_V(i)*Vc_ave_V(i)
            else
                call residualvectorset3(n*eq+q_judge,zeta(i))
                call residualvectorset3(n*eq+q_judge,ro_r_V(i))
            end if
        end do

        do i=1,com_2phase
            myu_c(i)=sqrt(Mw(i))*(Pc(i)/101325.d0)**(2.0d0/3.0d0)/(Tc(i)**(1.0d0/6.0d0))*(4.610d0*(temp/Tc(i))**0.618d0-&
                    2.04d0*exp(-0.449d0*(temp/Tc(i)))+1.94d0*exp(-4.058d0*(temp/Tc(i)))+0.1d0)*10.0d0**(-4.0d0) ![cP]
        end do
        !write(*,*) myu_c

        do i=1,n
            if (phase_judge(i) == 2) then
                !myu_V_SC(i)=(y(1,i)*myu_c(1)*sqrt(Mw(1))+y(2,i)*myu_c(2)*sqrt(Mw(2)))/&
                !            (y(1,i)*sqrt(Mw(1))+y(2,i)*sqrt(Mw(2))) ![cP]
                myu_V_SC(i)=(y(1,i)*myu_c(1)*sqrt(Mw(1))+y(2,i)*myu_c(2)*sqrt(Mw(2))&
                            +y(3,i)*myu_c(3)*sqrt(Mw(3))+y(4,i)*myu_c(4)*sqrt(Mw(4)))/&
                            (y(1,i)*sqrt(Mw(1))+y(2,i)*sqrt(Mw(2))+y(3,i)*sqrt(Mw(3))+y(4,i)*sqrt(Mw(4))) ![cP]
            else
                call residualvectorset3(n*eq+q_judge,myu_V_SC(i))
            end if
        end do

        av_0=1.0230d0*10.0d0**(-1.0d0)
        av_1=2.3364d0*10.0d0**(-2.0d0)
        av_2=5.8533d0*10.0d0**(-2.0d0)
        av_3=-4.0758d0*10.0d0**(-2.0d0)
        av_4=9.3324d0*10.0d0**(-3.0d0)

        do i=1,n
            if (phase_judge(i) == 2) then
                myu_V(i)=(myu_V_SC(i)+(((av_0+av_1*ro_r_V(i)+av_2*ro_r_V(i)**2.0d0+av_3*ro_r_V(i)**3.0d0+&
                         av_4*ro_r_V(i)**4.0d0)**4.0d0-10.0d0**(-4.0d0))/zeta(i)))*10.0d0**(-3.0d0) !cP→Pa・s
            else
                call residualvectorset3(n*eq+q_judge,myu_V(i))
            end if
        end do
        
        call outxs(myu_L,kakuninn)
        write(*,*) kakuninn!!粘度よさそう

        !!相対浸透率について
        do i=1,n
            do j=1,20
                if (Swd(j) < Sw0(i) .and. Sw0(i) <= Swd(j+1)) then
                    krg(i)=(krgd(j+1)-krgd(j))/(Swd(j+1)-Swd(j))*(Sw(i)-Swd(j))+krgd(j)
                    krw(i)=(krwd(j+1)-krwd(j))/(Swd(j+1)-Swd(j))*(Sw(i)-Swd(j))+krwd(j)
                end if
            end do
        end do

        !!孔隙率について
        NmMD(1)=Nm1_MD
        NmMD(2)=Nm2_MD
        NmMD(3)=Nm3_MD
        NmMD(4)=Nm4_MD
        NmMD(5)=Nm5_MD
        do i=1,n
            call residualvectorset4(fai0(i),n*eq+q_judge,faimine(i))
            faimineold(i)=fai0(i)
            do j=1,com_mine
                faimine(i)=faimine(i)-(Nm(j,i)-Nmini(j))/NmMD(j)
                faimineold(i)=faimineold(i)-(Nm0old(j,i)-Nmini(j))/NmMD(j)
            end do
            fai(i)=faimine(i)*(1.0d0+Cr*(P(i)-iniPressure))
            faiold(i)=faimineold(i)*(1.0d0+Cr*(P0old(i)-iniPressure))            
        end do
        

        !!絶対浸透率
        do i=1,n
            !kk(i)=k_ini*(fai(i)/faiini)**3.0d0*((1.0d0-faiini)/(1.0d0-fai(i)))**2.0d0
            call residualvectorset4(k_ini,n*eq+q_judge,kk(i))
        end do


        !!injection
        if (phase_judge(1) == 2) then !?2相の時
            ramda = krg(1)/myu_V(1) + krw(1)/myu_L(1)
        elseif (phase(1) == 1) then !?液相のみ
            ramda = 1.0d0/myu_L(1)
        else !?気相のみ 
            ramda = 1.0d0/myu_V(1)
        end if
        re=0.14d0*sqrt(dx**2.0d0+dy**2.0d0)
        if (q_judge == 1) then
            q=2.0d0*pi*kk(1)*dz*ramda*(Pb-P(1))/(log(re/rw)+skin)/dx/dy/dz !?流量制御
        else
            !?q=2.0d0*pi*kk(1)*dz*ramda*(Pbh0-P(1))/(log(re/rw)+skin)/dx/dy/dz !?圧力制御
            !?call residualvectorset3(n*eq+q_judge,q) !?圧力制御 !!ここ覚えていない、上の方が正しそうだけど、なんだっけ？
        end if

        call injection_data_d(P(1),MD_injection_d)

        !!weighting factor------
        do i=1,n-1
            if (P0(i+1) >= P0(i)) then
                W(i)=0.0d0
            else 
                w(i)=1.0d0
            end if
        end do
        !write(1,*)w
        !?----------------------

        !!トランスミッシビリティー-----
        do i=1,n
            if (phase_judge(i) == 2) then !?2相の時
                do j=1,com_2phase+com_ion
                    T_L(j,i)=kk(i)*krw(i)*MD_L(i)*x(j,i)/myu_L(i)
                end do
                do j=1,com_2phase
                    T_V(j,i)=kk(i)*krg(i)*MD_V(i)*y(j,i)/myu_V(i)
                end do
            elseif (phase(i) == 1) then !?液相だけの時
                do j=1,com_2phase+com_ion
                    T_L(j,i)=kk(i)*MD_L(i)*x(j,i)/myu_L(i)
                end do
                do j=1,com_2phase
                    call residualvectorset3(n*eq+q_judge,T_V(j,i))
                end do
            else !?気相だけの時
                do j=1,com_2phase+com_ion
                    call residualvectorset3(n*eq+q_judge,T_L(j,i))
                end do
                do j=1,com_2phase
                    T_V(j,i)=kk(i)*MD_V(i)*y(j,i)/myu_V(i)
                end do
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
