!変わる条件

module mod_input
    use mod_condition
    implicit none
    
    !grid数
    integer,parameter::n=5
    
    !gridの長さ
    real(8),parameter::dx=100.d0 ![m]
    real(8),parameter::dy=100.d0 ![m]
    real(8),parameter::dz=100.d0 ![m]
    
    !支配方程式の数
    integer,parameter::eq=25 !4+19+1+1
    
    !化学反応式の数
    integer,parameter::chemi=5
    
    !鉱物反応式の数
    integer,parameter::mine=5
    
    !成分の数
    integer,parameter::com_2phase=4 !2相に存在
    integer,parameter::com_ion=10 !イオン
    integer,parameter::com_mine=5 !鉱物
    integer,parameter::com_all=19 !全部
    
    
    !!化学反応の平衡定数
    real(8),parameter::ke1=10.0d0**(-6.54924d0+0.00900174d0*temp_C-0.0001021115d0*temp_C**2.0d0+2.76188d0*10.0d0&
                    **(-7.0d0)*temp_C**3.0d0&
                    -3.56142d0*10.0d0**(-10.0d0)*temp_C**4.0d0)!exp(-6.735d0)!4.47d0*10.0d0**(-7.0d0)
    real(8),parameter::ke2=10.0d0**(10.608d0-0.0127676d0*temp_C+0.000120258d0*temp_C**2.0d0-3.01731d0*10.0d0&
                **(-7.0d0)*temp_C**3.0d0&
                    +2.69372d0*10.0d0**(-10.0d0)*temp_C**4.0d0)!exp(-10.56d0)!4.67d0*10.0d0**(-8.0d0)
    real(8),parameter::ke3=10.0d0**(14.9282d0-0.0418762d0*temp_C+0.000197367d0*temp_C**2.0d0-5.54951d0*10.0d0&
                **(-7.0d0)*temp_C**3.0d0&
                    -9.15378d0*10.0d0**(-11.0d0)*temp_C**4.0d0)!exp(-13.56d0)!1.8d0*10.0d0**(-13.0d0)
    real(8),parameter::ke4=10.0d0**(-7.61239d0-0.0292971d0*temp_C-0.000267437d0*temp_C**2.0d0-8.99916d0*10.0d0&
                **(-7.0d0)*temp_C**3.0d0&
                    -1.21531d0*10.0d0**(-9.0d0)*temp_C**4.0d0)!exp(-10.56d0)!4.67d0*10.0d0**(-8.0d0)
    real(8),parameter::ke5=10.0d0**(14.6855d0-0.0325462d0*temp_C+0.000063246d0*temp_C**2.0d0-8.87466d0*10.0d0&
                **(-8.0d0)*temp_C**3.0d0&
                    -3.41634d0*10.0d0**(-11.0d0)*temp_C**4.0d0)!exp(-13.56d0)!1.8d0*10.0d0**(-13.0d0)
    
    !!鉱物反応の平衡定数
    real(8),parameter::km1=10.0d0**(31.7457d0-0.201254d0*temp_C+0.00059589d0*temp_C**2.0d0-9.04116d0*10.0d0&
                **(-7.0d0)*temp_C**3.0d0&
                    +9.15378d0*10.0d0**(-11.0d0)*temp_C**4.0d0)!exp(1.806d0)!2.29d0*10.0d0**4.0d0
    real(8),parameter::km2=10.0d0**(12.826d0-0.0581923d0*temp_C+0.000151291d0*temp_C**2.0d0-1.72044d0*10.0d0&
                **(-7.0d0)*temp_C**3.0d0&
                    -6.81994d0*10.0d0**(-11.0d0)*temp_C**4.0d0)!exp(1.806d0)!2.29d0*10.0d0**4.0d0
    real(8),parameter::km3=10.0d0**(31.4728d0-0.143258d0*temp_C+0.000424441d0*temp_C**2.0d0-6.99719d0*10.0d0&
                **(-7.0d0)*temp_C**3.0d0&
                    +2.90931d0*10.0d0**(-10.0d0)*temp_C**4.0d0)!exp(1.806d0)!2.29d0*10.0d0**4.0d0
    real(8),parameter::km4=10.0d0**(2.06889d0-0.0142668d0*temp_C-6.06096d0*10.0d0**(-6.0d0)*temp_C**2.0d0+1.45921*10.0d0&
                **(-7.0d0)*temp_C**3.0d0&
                    -4.18928d0*10.0d0**(-10.0d0)*temp_C**4.0d0)!exp(1.806d0)!2.29d0*10.0d0**4.0d0
    real(8),parameter::km5=10.0d0**(3.11813d0-0.0283728d0*temp_C+4.14154d0*10.0d0**(-5.0d0)*temp_C**2.0d0&
                    +4.75659d0*10.0d0**(-8.0d0)*temp_C**3.0d0-3.48756d0*10.0d0**(-10.0d0)*temp_C**4.0d0)
                    !exp(10.28d0)!3.0d0*10.0d0**4.0d0
    
    !!初期の比表面積[m^2/m^3]
    real(8),parameter::AA1=330 !carbfix site
    real(8),parameter::AA2=330
    real(8),parameter::AA3=330 !carbfix site
    real(8),parameter::AA4=330
    real(8),parameter::AA5=330 !carbfix site
    


    
    
    !岩石密度
    real(8),parameter::Rd=3.0d0*10.0d0**6.0d0 ![g/cm^3]→[g/m^3]
    
    !!!!!!!!CaCO3!!!!!!!
    !鉱物の初期割合
    real(8),parameter::Nm1_per=0.01
    
    !鉱物のMw
    real(8),parameter::Mwm1=100.089
    
    !初期の鉱物のモル数
    real(8),parameter::Nm1_ini=Rd*Nm1_per/Mwm1
    
    !鉱物のモル密度[mol/m^3]
    real(8),parameter::Nm1_MD=2.70995d0*10.0d0**6.0d0/Mwm1
    
    !!!!!!!!MgCO3!!!!!!!
    !鉱物の初期割合
    real(8),parameter::Nm2_per=0.01
    
    !鉱物のMw
    real(8),parameter::Mwm2=84.3142
    
    !初期の鉱物のモル数
    real(8),parameter::Nm2_ini=Rd*Nm2_per/Mwm2
    
    !鉱物のモル密度[mol/m^3]
    real(8),parameter::Nm2_MD=3.00929d0*10.0d0**6.0d0/Mwm2
    
    
    !!!!!!!!CaAl2Si2O8!!!!!!!
    !鉱物の初期割合
    real(8),parameter::Nm3_per=0.1
    
    !鉱物のMw
    real(8),parameter::Mwm3=278.209
    
    !初期の鉱物のモル数
    real(8),parameter::Nm3_ini=Rd*Nm3_per/Mwm3
    
    !鉱物のモル密度[mol/m^3]
    real(8),parameter::Nm3_MD=2.76029d0*10.0d0**6.0d0/Mwm3
    
    !!!!!!!!MgSiO3!!!!!!!
    !鉱物の初期割合
    real(8),parameter::Nm4_per=0.1
    
    !鉱物のMw
    real(8),parameter::Mwm4=100.389
    
    !初期の鉱物のモル数
    real(8),parameter::Nm4_ini=Rd*Nm4_per/Mwm4
    
    !鉱物のモル密度[mol/m^3]
    real(8),parameter::Nm4_MD=3.20977d0*10.0d0**6.0d0/Mwm4
    
    !!!!!!!!Mg2SiO4!!!!!!!
    !鉱物の初期割合
    real(8),parameter::Nm5_per=0.1
    
    !鉱物のMw
    real(8),parameter::Mwm5=140.693
    
    !初期の鉱物のモル数
    real(8),parameter::Nm5_ini=Rd*Nm5_per/Mwm5
    
    !鉱物のモル密度[mol/m^3]
    real(8),parameter::Nm5_MD=3.2129d0*10.0d0**6.0d0/Mwm5
    
    !!孔隙率
    real(8),parameter::faiini=0.085
    
    !岩石圧縮率
    real(8),parameter::Cr=8.073d0*10.0d0**(-10.0d0) !苫小牧
    
    !!timestep
    real(8),parameter::dt=1.0d0*60.0d0*60.0d0!*24.0d0
    
    
    
    
    
end module