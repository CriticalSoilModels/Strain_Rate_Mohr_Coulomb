module mod_SRMC_NAMC
    implicit none
    
contains
    
!_____________________________________________________________________________________________
!##    ##    ###    ##     ##  ######     ##     ##  ######  ########
!###   ##   ## ##   ###   ### ##    ##    ##     ## ##    ## ##     ##
!####  ##  ##   ##  #### #### ##          ##     ## ##       ##     ##
!## ## ## ##     ## ## ### ## ##          #########  ######  ########
!##  #### ######### ##     ## ##          ##     ##       ## ##   ##
!##   ### ##     ## ##     ## ##    ##    ##     ## ##    ## ##    ##
!##    ## ##     ## ##     ##  ######     ##     ##  ######  ##     ##

subroutine NAMC_HSR(NOEL, G_0, nu, M_tc, N, D_min, h, alpha_G, alpha_K, alpha_D, D_part, G_s,&
    switch_smooth, RefERate, N_S, switch_original, &
    G, K, eta_y, Dp, EpsP, I_coeff, switch_yield, N_i, SUM_rate,&
    dEps, Sig_0, Sig, Erate, DTIME,&
    Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, &
    Error_Yield_max, FTOL, max_stress_iters)

    !*********************************************************************************
    !
    ! Elasto-plastic constitutive model with strain softening and strain rate effects,
    ! Non associated MOHR-COULOMB
    ! Explicit MODIFIED EULER INTEGRATION SCHEME with automatic error control.
    ! Final correction of the yield surface drift (END OF STEP CORRECTION).
    !
    !*********************************************************************************
    implicit none

    !input variables
    integer, intent(in) :: NOEL !Global ID of Gauss point or particle
    integer, intent(in) :: N_S, max_stress_iters
    double precision, intent(in)::  G_0, nu, M_tc, N, D_min, alpha_G, alpha_K, &
       alpha_D, D_part, G_s, &
       RefERate, DTIME, FTOL
    double precision, dimension(6), intent(in)::  dEps
    double precision, dimension(6), intent(inout):: Sig_0

    logical, intent(in):: switch_smooth, switch_original

    !output variables
    integer, intent(inout):: N_i

    double precision, intent(inout):: h, G, K, eta_y, Dp, I_coeff, &
       Error_yield_1, Error_yield_2, Error_Euler_max,&
       Error_Yield_last, Error_Yield_max, SUM_rate
    double precision, dimension(6), intent(inout):: Sig, Erate, EpsP

    logical, intent(inout):: switch_yield

    !local variables
    integer:: counter, MAXITER, SubStepping_MaxIter
    double precision:: p, q, theta, M, I_act, I_0, Gu, Ku, eta_yu, Du, dI, &
       G1, K1, eta_y1, Dp1, G2, K2, eta_y2, Dp2, &
       p_t, q_t, dI_t, dI_TT, I_TT, dD1, dD2
    double precision:: epsq_rate, epsq_p, eps_v
    double precision:: F0, FT, alpha
    double precision:: STOL, DTmin , LTOL, R_TT, qR, RTOL
    double precision:: dummyvar(3), D1, D2
    double precision:: DE(6,6), dSig_el(6), Sig_t(6), dEps_t(6), dEps_TT(6), &
       Sig1(6), Sig2(6), dEpsp1(6), dEpsp2(6), dEpsp(6), &
       dSig1(6), dSig2(6), Epspt(6)
    double precision:: T, DT

    logical:: ApplyStrainRateUpdates=.false., &!Check if the strain rate path crosses the reference line
       IsUnloading, & !If true material unloads (not elastic unloading)
       Failed=.false. !If rel residual is too high failed is true
    !print *, NOEL
    !______________________________________________________________________________
    ! Error tolerances
    !FTOL=1.0e-8		!Yield surface tolerance
    STOL=1.0e-3		!Relative error tolerance
    DTmin=1.0e-9	!Minimum pseudo-time increment, originally Dtmin = 1.0e-9
    LTOL=0.01d0		!Tolerance for elastic unloading
    !MAXITER=20		!Max. number of iterations
    RTOL = STOL * 1.0e-1
    !______________________________________________________________________________
    !______________________________________________________________________________
    ! Initialization of error trackers
    Error_yield_1   =0.0	!Drift after first approximation
    Error_yield_2   =0.0	!Drift after second approximation
    Error_Euler_max =0.0    !Max relative residual
    Error_Yield_max =0.0    !max drift after averaging
    Error_Yield_last=0.0    !max abs drift after drift correction
    !______________________________________________________________________________
    !Initialization of state variables
    call Get_invariants(Sig_0, p, q, theta)
    call Get_strain_invariants(EpsP, eps_v, epsq_p)!plastic deviatoric strain
    I_0=RefERate
    !call Get_I_coeff(D_part, G_s, -100.0, RefERate, I_0)!Reference inertial coefficient
    call Get_M(M_tc, theta, M)!Get M
    if (G==0.0d0) then
       G=G_0
    endif
    if (K==0.0d0) then
       K=2*G_0*(1+nu)/(3*(1-2*nu))
    endif
    !if (I_coeff==0.0d0) then
    !    call Get_I_coeff(D_part, G_s, p, RefERate, I_coeff)
    !endif

    ! Dp testing: This should be the first time that the dilatancy can be updated
    if (Dp==0.0d0) then ! Since plastic deviatoric strain epsq_p will be zero in elastic Dp will be be zero
       call Get_Dp(h, D_min, I_coeff, I_0, epsq_p, alpha_D, ApplyStrainRateUpdates, Dp)
    endif

    !print *, eta_y
    if (eta_y==0.0d0) then
       call Get_invariants(dEps, dummyvar(1), dummyvar(2), theta)
       call Get_M(M_tc, theta, M)
       eta_y=M-Dp*(1.0-N)
    endif
    !print *, eta_y

    !_____________________________________________________________________________
    !Evaluate yield function at initial stress-state variable

    call YieldFunction(q, p, eta_y, F0)

    !_____________________________________________________________________________
    !Compute the current deviatoric strain rate

    call TwoNormTensor_strain(Erate,6,dummyvar(1))!Norm of the strain rate

    !Compute a smoothed strain rate if switch_smooth is true
    if (switch_smooth) then
       N_i=N_i+1
       !if (Noel ==1) then
       !    print *, N_i
       !endif
       if (N_i<N_s) then !not enough values
          Sum_rate=sum_rate+dummyvar(1) !accumulate the strain rate
          dummyvar(2)=Sum_rate/n_i !takes the average
       else !enough values
          Sum_rate=Sum_rate*(1.0-1.0/N_s) !approximate sum without one term
          Sum_rate=Sum_rate+dummyvar(1)
          dummyvar(2)=Sum_rate/N_s !averaged strain rate
       endif
       if (dummyvar(1)==0.0d0) then !in case in first step this value is zero
          Erate=0.0d0
       else
          Erate=(dummyvar(2)/dummyvar(1))*Erate !corrected strain rate tensor
       endif
       call TwoNormTensor_strain(Erate,6,dummyvar(1)) ! recalculate the two norm so the updated erate value can be tracked
    endif

    !________________________________________________________________________________

    !Update state variables due to HSR
    !Store state variables
    Gu=G
    Ku=K
    eta_yu=eta_y
    Du=Dp

    call Get_strain_invariants(Erate, dummyvar(1), epsq_rate)! deviatoric strain rate
    call Get_I_coeff(D_part, G_s, p, epsq_rate, I_act)!actual inertial coefficient

    !I_act = 0 ! Test to see if this is the reason the models get different results for no updates and D_part = 0
    dI=I_act-I_coeff !change in inertial coefficient
    call check4crossing(I_coeff, I_act, dI, I_0, ApplyStrainRateUpdates) !Check if update needed

    ! Dp Testing: Second time that dilatancy can be updated.
    if (applystrainrateupdates) then !update
       !h=h*(i_act/i_0)**alpha_g
       call update_gk(g_0, nu, i_act, i_0, alpha_g, alpha_k, gu, ku) !updates the elastic properties
       call get_dp(h, d_min, i_act, i_0, epsq_p, alpha_d, applystrainrateupdates, du) !updates the dilation
       eta_yu=m-du*(1.0-n) !updates eta
    endif
    !_________________________________________________________________________________

    ! Fill elastic material matrix
    D1  = Ku+(4*Gu/3)
    D2  = Ku-(2*Gu/3)
    DE  = 0.0
    DE(1:3,1:3) = D2
    DE(1,1) = D1
    DE(2,2) = D1
    DE(3,3) = D1
    DE(4,4) = Gu
    DE(5,5) = Gu
    DE(6,6) = Gu

    !Get elastic predictor stress
    call MatVec(DE,6,dEps,6,dSig_el)
    Sig_t = Sig_0 + dSig_el

    !Get new invariant stresses
    call Get_invariants(Sig_t, p_t, q_t, theta)

    !Evaluate yield function
    call YieldFunction(q_t, p_t, eta_yu, FT)

    !___________________________________________________________________________
    !Now check elastic loading, unloading
    if (FT<-FTOL) then !Elastic behavior
       !Update state parameters
       G=Gu
       K=Ku
       eta_y=eta_yu
       Dp=Du
       I_coeff=I_act

       !Update stress
       Sig=Sig_t
       switch_yield=.false.

       !___________________________________________________________________________
       !*************************  Plastic behavior  ******************************
       !***************************************************************************
       !***************************************************************************

       !___________________________________________________________________________
       !***********Checking surface bisection and viscoplastic unloading***********
       !===========================================================================
    else !plastic behavior
       !if (F0<-FTOL) then !Elasto-plastic transition
       !    call Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc, N, D_min, h, D_part, &
       !G_s, epsq_p, I_coeff, I_act, I_0, alpha_K, alpha_G, alpha_D,&
       !Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL,&
       !F0, Sig_0, alpha)
       !    !Track Deformation category
       !    DeformCateg = 200.0
       !
       !
       !else!pure plastic deformations or unloading
       !    call Check_Unloading(M_tc, eta_y, eta_yu, dI, Sig_0, dSig_el, LTOL, IsUnloading)  !Determines if is unloading path
       !    if (IsUnloading) then !Find new alpha
       !        call Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc,N, D_min, h, D_part, &
       !	G_s, epsq_p, I_coeff, I_act, I_0, alpha_K, alpha_G, alpha_D,&
       !	Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL,&
       !	F0, Sig_0, alpha)
       !        !Track Deformation category
       !        DeformCateg = 300.0
       !
       !    else !Pure plasticity
       !        alpha=0.0d0
       !        !Track Deformation category
       !        DeformCateg = 400.0
       !    endif
       !endif

       ! Now that the updated state parameters have been found plug integrate the
       ! Need to make dEpsP a state variable (Done)
       call Ortiz_Simo_Integration(G_0, nu, M_tc, M, N, D_min, h, G, K, eta_y, Dp, &
          I_0, I_coeff, dI, alpha_G, alpha_K, alpha_D, Sig_0, EpsP, dEps, &
          FTOL, NOEL, max_stress_iters)
       Sig = Sig_0 ! Update the stresses

    endif
 end subroutine NAMC_HSR
 
end module mod_SRMC_NAMC