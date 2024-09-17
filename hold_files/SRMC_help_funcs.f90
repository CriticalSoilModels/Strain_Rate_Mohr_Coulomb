module mod_SRMC_funcs
    use kind_precision_module, only: real_type => dp
    implicit none
    
contains
    


!______________________________________________________________________________
!##     ##    ###    ########  ########
!##     ##   ## ##   ##     ## ##     ##
!##     ##  ##   ##  ##     ## ##     ##
!######### ##     ## ########  ##     ##
!##     ## ######### ##   ##   ##     ##
!##     ## ##     ## ##    ##  ##     ##
!##     ## ##     ## ##     ## ########
! ######   #######  ######## ########       ###    ##    ## ########
!##    ## ##     ## ##          ##         ## ##   ###   ## ##     ##
!##       ##     ## ##          ##        ##   ##  ####  ## ##     ##
! ######  ##     ## ######      ##       ##     ## ## ## ## ##     ##
!      ## ##     ## ##          ##       ######### ##  #### ##     ##
!##    ## ##     ## ##          ##       ##     ## ##   ### ##     ##
! ######   #######  ##          ##       ##     ## ##    ## ########
!########  ######## ########  #### ##     ##    ###    ######## #### ##     ##
!##     ## ##       ##     ##  ##  ##     ##   ## ##      ##     ##  ##     ##
!##     ## ##       ##     ##  ##  ##     ##  ##   ##     ##     ##  ##     ##
!##     ## ######   ########   ##  ##     ## ##     ##    ##     ##  ##     ##
!##     ## ##       ##   ##    ##   ##   ##  #########    ##     ##   ##   ##
!##     ## ##       ##    ##   ##    ## ##   ##     ##    ##     ##    ## ##
!########  ######## ##     ## ####    ###    ##     ##    ##    ####    ###
!########  ######
!##       ##    ##
!##       ##
!######    ######
!##             ##
!##       ##    ##
!########  ######
!Derivatives inclosed here
   subroutine Get_Dp(h, D_min, I, I_0, eps_q, k, ApplyRateUpdating, D)
    !*********************************************************************
    ! Returns the dilation for current inertial coefficient and dev.     *
    ! strain															 *
    !*********************************************************************
    implicit none
    logical, intent(in):: ApplyRateUpdating
    double precision, intent(in):: h, D_min, I, I_0, eps_q, k
    !out
    double precision, intent(out):: D
    !local variables
    double precision:: D_mm
    if (ApplyRateUpdating) then
       D_mm=D_min*(I/I_0)**k !strain/rate hardening
    else
       D_mm=D_min
    endif

    D=h*D_mm*eps_q*exp(1.0-h*eps_q) !hardening rule
 end subroutine Get_Dp

 subroutine Update_GK(G_0, nu, I, I_0, k_G, k_K, G, K)
    !*********************************************************************
    ! Returns updated elastic modulus                                    *
    !																	 *
    !*********************************************************************
    implicit none
    !input
    double precision, intent(in):: G_0, nu, I, I_0, k_G, k_K
    !output
    double precision, intent(out):: G, K
    !local variables
    double precision:: K_0
    G=G_0*(I/I_0)**k_G! updated modulus

    K_0=2*G_0*(1+nu)/(3*(1-2*nu))!bulk modulus
    K=K_0*(I/I_0)**k_K! updated modulus
 end subroutine Update_GK

 subroutine Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec)
    !************************************************************************
    ! Returns the derivative of the yield function with respect to the		*
    ! stress tensor 														*
    ! n=dF/dSigma =dF/dp*dp/dSigma+ dF/dq*dq/dSigma +dF/dtheta*dtheta/dSigma*
    ! n is a (1X6) vector													*
    !************************************************************************
    implicit none
    !input
    double precision, intent(in):: M_tc, eta_y, Sig(6)
    !output
    double precision, dimension(6):: n_vec
    !local variables
    double precision:: p, q, theta, pi=2.0*acos(0.0d0), &
       J2, J3, dJ3dsig(6), dfdtheta, &
       dpdsig(6), dqdsig(6), dev(6), dev2(6), &
       TrS2, II(6), dthetadSig(6), COS_3THETA
    !Get the invariants
    call Get_invariants(Sig, p, q, theta)
    !Get dF/dp=eta_y and dF/dq=1
    !Get dF/dtheta
    dfdtheta=0.45*p*M_tc*((cos(1.5d0*theta+0.25d0*pi))**0.2)*sin(1.5d0*theta+0.25d0*pi)
    !___________________________________________________________________________
    !1) Get dp/dSig=1/3 Imat
    dpdsig=0.0d0
    dpdsig(1)=1.0/3.0
    dpdsig(2)=1.0/3.0
    dpdsig(3)=1.0/3.0

    !2) Get dq/dsig= 2 *dev/3*q
    dev=Sig
    dev(1)=dev(1)-p
    dev(2)=dev(2)-p
    dev(3)=dev(3)-p

    dqdSig=( 3.0 / ( 2.0*q ) ) * dev

    !3) Get dtheta/dSigma= (1/3cos3theta) d/dsigma((J3/2) * (3/J2)^1.5)
    J2=(q**2)/3.0

    J3=dev(1)*dev(2)*dev(3) - dev(1)*dev(6)**2 - dev(2)*dev(4)**2 - dev(3)*dev(5)**2 + 2.0*dev(4)*dev(5)*dev(6)

    !Fill S.S
    dev2(1)=dev(1)**2 + dev(4)**2 + dev(5)**2
    dev2(2)=dev(2)**2 + dev(4)**2 + dev(6)**2
    dev2(3)=dev(3)**2 + dev(5)**2 + dev(6)**2
    dev2(4)=dev(4) * ( dev(1) + dev(2) ) + dev(5)*dev(6)
    dev2(5)=dev(5) * ( dev(1) + dev(3) ) + dev(4)*dev(6)
    dev2(6)=dev(6) * ( dev(2) + dev(3) ) + dev(4)*dev(5)

    !Compute dJ3dSig
    TrS2 = dev2(1) + dev2(2) + dev2(3)

    II=0.0d0!Identity tensor
    II(1)=1.0
    II(2)=1.0
    II(3)=1.0

    dJ3dsig = dev2 - ( TrS2*II / 3.0d0 )

    !Compute dtheta/dsig

    dthetadSig = dJ3dsig - ( 1.5*J3 / J2 ) *dev
    COS_3THETA = cos( 3.0*theta )
    dthetadSig = ( sqrt(3.0) / ( 2.0*COS_3THETA*J2**1.5 ) ) * dthetadSig

    !__________________________________________________________________
    !Get n_vec=dF/dSig
    n_vec=(eta_y*dpdsig)+dqdSig+(dfdtheta*dthetadSig) !n_vec=dF/dSig
 end subroutine Get_dF_to_dSigma

 subroutine Get_dD_to_dI(D_min, h, I_0, kD, eps_q, I, b)
    !************************************************************************
    ! Returns the derivative of the Dilation with respect to the inertial	*
    ! coefficient 															*
    ! b=dD/dI																*
    ! b is a scalar															*
    !************************************************************************
    implicit none
    !input
    double precision, intent(in):: D_min, h, I_0, kD, eps_q, I
    !output
    double precision, intent(out)::b
    !local variables

    b=h*D_min*eps_q*exp(1-h*eps_q)*kD*((I/I_0)**(kD-1.0))/I_0

 end subroutine Get_dD_to_dI

 subroutine Get_dP_to_dSigma(D, Sig, m_vec)
    !************************************************************************
    ! Returns the derivative of the plastic potential function with respect *
    ! to the stress tensor													*
    ! m=dP/dSigma =dP/dp*dp/dSigma+ dP/dq*dq/dSigma							*
    ! m is a (1X6) vector													*
    !************************************************************************
    implicit none
    !input
    double precision, intent(in):: D, Sig(6)
    !output
    double precision, dimension(6):: m_vec
    !local variables
    double precision:: p, q, theta, pi=2.0*acos(0.0d0), &
       J2, J3, dJ3dsig(6), dfdtheta, &
       dpdsig(6), dqdsig(6), dev(6), dev2(6), &
       TrS2, II(6), dthetadSig(6), COS_3THETA
    !Get the invariants
    call Get_invariants(Sig, p, q, theta)
    !Get dP/dp=-D and dF/dq=1
    !___________________________________________________________________________
    !1) Get dp/dSig=1/3 Imat
    dpdsig=0.0d0
    dpdsig(1)=1.0/3.0
    dpdsig(2)=1.0/3.0
    dpdsig(3)=1.0/3.0

    !2) Get dq/dsig= 2 *dev/3*q
    dev=Sig
    dev(1)=dev(1)-p
    dev(2)=dev(2)-p
    dev(3)=dev(3)-p

    dqdSig=(3.0/(2.0*q))*dev

    !__________________________________________________________________
    !Get m_vec=dP/dSig
    m_vec=(-D*dpdsig)+dqdSig !m_vec=dP/dSig
 end subroutine Get_dP_to_dSigma

 subroutine Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
    EpsP, I, ApplyStrainRateUpdate, a)
    !************************************************************************
    ! Returns the derivative of the Dilation with respect to the plastic    *
    ! strain																*
    ! a=dXs/dEpsp= dD/dEpsq * dEPsq/dEpsp									*
    ! a is a (1X6) vector													*
    !************************************************************************
    implicit none
    !input
    logical, intent(in):: ApplyStrainRateUpdate
    double precision, intent(in):: D_min, h, I_0, k_D, epsq_p, epsv_p, &
       EpsP(6), I
    !output
    double precision, intent(out):: a(6)
    !local variables
    double precision:: D, dDdEpsq_p, dev(6),dEpsq_pdEpsp(6)

    !________________________________________________________________________
    !Get dD/dEpsq_p
    if (ApplyStrainRateUpdate) then
       D=D_min*(I/I_0)**k_D
    else
       D=D_min
    end if

    dDdEpsq_p=h*D*exp(1.0-h*epsq_p)*(1.0-h*epsq_p)
    !_______________________________________________________________________

    !_______________________________________________________________________
    !Get dEpsp_Q/dEpsp=
    dev=EpsP
    dev(1)=dev(1)-(epsv_p/3.0)
    dev(2)=dev(2)-(epsv_p/3.0)
    dev(3)=dev(3)-(epsv_p/3.0) !deviatoric stress tensor

    if (epsq_p>0.0d0) then !in case of zero plastic strain
       dEpsq_pdEpsp=(2.0/(3.0*epsq_p))*dev
    else
       dEpsq_pdEpsp=0.0d0
    endif
    !_______________________________________________________________________

    !_______________________________________________________________________
    !Get a=dXs/dEpsp

    a=dDdEpsq_p*dEpsq_pdEpsp
    !______________________________________________________________________
 end subroutine Get_dD_to_dEpsP

 subroutine Get_dEpsq_to_dEps(Epsq, Eps, dEqdEpsq)
    !************************************************************************
    ! Returns the derivative of the deviatoric strain with respect to the   *
    ! deviatoric strain	tensor						     					*
    ! dEqdEpsq is a (1X6) vector											*
    !************************************************************************
    implicit none
    !input
    double precision, intent(in):: Epsq, Eps(6)
    !output
    double precision, intent(out):: dEqdEpsq(6)
    !local variables
    double precision:: evol, dev(6)

    evol=Eps(1)+Eps(2)+Eps(3)!vol strain

    dev=Eps
    dev(1)=dev(1)-evol/3.0
    dev(2)=dev(2)-evol/3.0
    dev(3)=dev(3)-evol/3.0 !deviatoric strain tensor

    if (Epsq>0.0d0) then !in case of zero plastic strain
       dEqdEpsq=(2.0/(3.0*Epsq))*dev
    else
       dEqdEpsq=0.0d0
    endif
 end subroutine Get_dEpsq_to_dEps
!**********************************************************************************


!_________________________________________________________________________________
!######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
!##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
!##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
!######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
!##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
!##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
!##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######


 Subroutine Get_strain_invariants(Eps, Eps_v, Eps_q)
    !*********************************************************************
    ! Takes the strain tensor and returns deviatoric and vol. strain     *
    !																	 *
    !*********************************************************************
    implicit none
    !input
    double precision, dimension(6), intent(in):: Eps
    !output
    double precision, intent(out):: Eps_v, Eps_q
    !local variables
    double precision:: dev(6)
    Eps_v=Eps(1)+Eps(2)+Eps(3)! vol strain

    dev=Eps
    dev(1)=dev(1)-(Eps_v/3.0)
    dev(2)=dev(2)-(Eps_v/3.0)
    dev(3)=dev(3)-(Eps_v/3.0)!deviatoric strain tensor

    call TwoNormTensor_strain(dev, 6, Eps_q)
    Eps_q=Eps_q*sqrt(2.0/3.0) ! dev strain
 end subroutine  Get_strain_invariants


 subroutine Get_invariants(Sig, p, q, theta)
    !*********************************************************************
    ! Takes the stress tensor Sig and return invariants p, q, and theta  *
    !																	 *
    !*********************************************************************
    implicit none
    !input variables
    double precision, dimension(6), intent(in):: Sig
    !output variables
    double precision, intent(out)::p, q, theta
    !local variables
    double precision:: dev(6), J2, J3, sin3theta

    p=(Sig(1)+Sig(2)+Sig(3))/3.0 !mean stress
    dev=Sig
    dev(1)=dev(1)-p !computes deviatoric stress tensor
    dev(2)=dev(2)-p
    dev(3)=dev(3)-p

    call TwoNormTensor(dev, 6, J2)
    J2=(J2**2)/2.0 !J_2 invariant
    q=sqrt(3*J2) ! deviatoric stress

    !J3 stress invariant
    J3 = dev(1)*dev(2)*dev(3) - dev(1)*dev(6)**2 - dev(2)*dev(4)**2 - dev(3)*dev(5)**2 + 2.0*dev(4)*dev(5)*dev(6)

    !sin3theta
    if (J2>0.0d0) then
       sin3theta=0.5*J3*(3.0/J2)**(1.5d0)
    else !Assume triaxial compression
       sin3theta=-1.0d0
    endif
    if (sin3theta<-1.0) sin3theta=-1.0d0
    if (sin3theta>1.0) sin3theta=1.0d0


    theta=-asin(sin3theta)/3.0d0 !Lode's angle

 end subroutine Get_invariants


 subroutine Check_Unloading(M_tc, eta_y, eta_yu, dI, Sig, dSig,&
    LTOL, IsUnloading)
    !*********************************************************************
    ! Returns true if stress path is viscoplastic unloading              *
    !																	 *
    !*********************************************************************
    !input
    implicit none
    double precision, intent(in):: M_tc, eta_y, eta_yu, dI, &
       Sig(6), dSig(6), LTOL
    !output
    logical, intent(out):: IsUnloading
    !local variables
    double precision:: deta, n_vec(6), n_norm, Sig_norm,&
       dSIg_inner_n, beta, phi
    IsUnloading=.false.

    deta=eta_yu-eta_y!change in state parameter
    call Get_dF_to_dSigma(M_tc, eta_yu, Sig, n_vec)!Normal to surface
    call TwoNormTensor(n_vec, 6, n_norm) !norm of n_vec
    call TwoNormTensor(dSig, 6, Sig_norm) !norm of dSig
    call TensorInnerProduct(dSig, n_vec, 6,dSIg_inner_n) !inner product

    beta=acos(deta/(n_norm*Sig_norm))!conical aperture is a plane for inviscid mat.
    phi=acos(dSIg_inner_n/(n_norm*Sig_norm))!angle between stress rate and normal

    if (phi-beta>LTOL) IsUnloading=.true. !condition for unloading

 end subroutine Check_Unloading


 subroutine Get_I_coeff(D_part, G_s, p, eps_rate, I)
    !*********************************************************************
    ! Returns the inertial coefficient                                   *
    !																	 *
    !*********************************************************************
    implicit none
    !input
    double precision, intent(in):: D_part, G_s, p, eps_rate
    !output
    double precision, intent(out):: I
    !local variables
    I=D_part*eps_rate*sqrt(G_s/abs(p))
 end subroutine Get_I_coeff

 subroutine Get_M(M_tc, theta, M)
    !*********************************************************************
    ! Returns M															 *
    !																	 *
    !*********************************************************************
    implicit none
    !in
    double precision, intent(in):: M_tc, theta
    !out
    double precision, intent(out):: M
    !local
    double precision:: COS_VAL, pi=2*acos(0.0d0)
    COS_VAL=cos(1.5*theta+0.25*pi)
    M=M_tc*(1+0.25*COS_VAL**1.2)
 end subroutine Get_M


 subroutine YieldFunction(q, p, eta_y, F)
    !*********************************************************************
    ! Returns the value of the yield function evaluated at q, p , eta    *
    !																	 *
    !*********************************************************************
    implicit none
    !in
    double precision, intent(in):: q, p, eta_y
    !out
    double precision, intent(out):: F
    !local variables

    F=q+eta_y*p !sign is due to compression being negative in UMAT
 end subroutine YieldFunction
!***********************************************************************************************

!_______________________________________________________________________________________________
!##     ##    ###    ######## ##     ##
!###   ###   ## ##      ##    ##     ##
!#### ####  ##   ##     ##    ##     ##
!## ### ## ##     ##    ##    #########
!##     ## #########    ##    ##     ##
!##     ## ##     ##    ##    ##     ##
!##     ## ##     ##    ##    ##     ##

 Subroutine TensorInnerProduct(TensorA, TensorB, N, Re)
    !***********************************************************************
    !
    !     Calculate 2NormTensor = sqrt(A:A)
    !
    ! I   Tensor  : (Square or vector of dimension N)
    ! I   N     :   Number of elements
    ! O   2Norm : Resulting norm
    !
    !***********************************************************************
    implicit none

    real(real_type),intent(in) :: TensorA(N), TensorB(N)
    real(real_type),intent(out) :: Re
    integer, intent(in) :: N
    !***********************************************************************

    ! Local variables
    integer :: X, I

    X=N/2
    Re=0.0d0
    Do I=1,X
       Re=Re+TensorA(I)*TensorB(I)
    end Do
    Do I=X+1,N
       Re=Re+2*(TensorA(I)*TensorB(I))
    end do
 end subroutine TensorInnerProduct


 Subroutine TwoNormTensor(Tensor, N, TwoNorm)
    !***********************************************************************
    !
    !     Calculate 2NormTensor = sqrt(A:A)
    !
    ! I   Tensor  : (Square or vector of dimension N)
    ! I   N     :   Number of elements
    ! O   2Norm : Resulting norm
    !
    !***********************************************************************
    implicit none
    real(real_type) :: Tensor(N)
    integer :: N
    real(real_type) :: TwoNorm
    !***********************************************************************
    ! Local variables
    integer :: X, I
    
    X=N/2
    TwoNorm=0.0d0
    Do I=1,X
       TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
    end Do
    Do I=X+1,N
       TwoNorm=TwoNorm+2*(Tensor(I)*Tensor(I))
    end do
    TwoNorm=sqrt(TwoNorm)

 end subroutine TwoNormTensor

 Subroutine TwoNormTensor_strain(Tensor, N, TwoNorm)
    !***********************************************************************
    !
    !     Calculate 2NormTensor = sqrt(A:A)
    !
    ! I   Tensor  : (Square or vector of dimension N)
    ! I   N     :   Number of elements
    ! O   2Norm : Resulting norm
    !
    !***********************************************************************
    implicit none  
    real(real_type) :: Tensor(N), TwoNorm
    integer :: N
    !***********************************************************************
    ! Local variables
    integer :: X, I

    X=N/2
    TwoNorm=0.0d0
    Do I=1,X
       TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
    end Do
    Do I=X+1,N
       TwoNorm=TwoNorm+0.5*(Tensor(I)*Tensor(I))!The convention in UMAT is to use engineering shear strains
    end do
    TwoNorm=sqrt(TwoNorm)

 end subroutine TwoNormTensor_strain

 Subroutine MatVec(xMat,IM,Vec,N,VecR)
    !***********************************************************************
    !
    !     Calculate VecR = xMat*Vec
    !
    ! I   xMat  : (Square) Matrix (IM,*)
    ! I   Vec   : Vector
    ! I   N     : Number of rows/colums
    ! O   VecR  : Resulting vector
    !
    !***********************************************************************
    implicit none
    real(real_type), intent(in)  :: xMat(N, N), Vec(N)
    integer, intent(in)          :: IM, N
    real(real_type), intent(out) :: VecR(N)

    !***********************************************************************
    ! Local variables
    integer :: I, J
    real(real_type) :: X

    Do I=1,N
       X=0
       Do J=1,N
          X=X+xMat(I,J)*Vec(J)
       End Do
       VecR(I)=X
    End Do
    Return
 End Subroutine MatVec

 Subroutine DotProduct_2(VecA, VecB,N, Dp)
    !***********************************************************************
    !
    !     Calculate the dot product of A(Nx1) and B(1xN)
    !
    ! I   VecA VecB  : Vectors
    ! I   N     :   Dimension
    ! O   Dp : Dot product
    !
    !***********************************************************************
    implicit none
    real(real_type), intent(in)  :: VecA(N), VecB(N)
    integer, intent(in)          :: N
    real(real_type), intent(out) :: Dp

    !***********************************************************************
    ! Local variables
    integer :: I
    Dp=0.0d0
    Do I=1,N
       Dp=Dp+VecA(I)*VecB(I)
    end do

 end subroutine DotProduct_2

 subroutine dbltobool(A,B)
    !******************************************************************
    ! Takes a double which values are either 1.0 or 0.0 and returns a *
    ! Boolean
    !******************************************************************
    implicit none
    double precision, intent(in):: A
    logical, intent(out):: B
    if (A<1.0) then
       B=.false.
    else
       B=.true.
    endif
 end subroutine dbltobool

 real(real_type) function logic2dbl(a)
    logical, intent(in) :: a

    if (a) then
       logic2dbl = 1.d0
    else
       logic2dbl = 0.d0
    end if
 end function logic2dbl

 Subroutine check4crossing(IErate0I, IErateI, dErate_eff,RateRef, Apply)
    !******************************************************************
    ! determines if strain rate updating must occur                   *
    ! Boolean                                                         *
    !******************************************************************
    ! IErate0I: The previous (initial reference strain rate. State parameter pulled from the previous time step
    ! IErateI: The new inertial coefficient for this strain rate
    ! dErate_eff: The increment of strain rate change (This is calculated in this subroutine)
    ! RateRef: Reference strain rate other variables are compared to
    ! Apply: Boolean keeping track to determine if strain rate updates should be applied
    implicit none
    double precision, intent(inout):: IErate0I, IErateI, dErate_eff
    double precision, intent(in)   :: RateRef
    logical:: cond1, cond2, cond3
    logical, intent(out)::Apply
    Apply=.false.
    ! If the rate from the last time step is less than or equalt ot the reference rate, update the previous time step value to be the reference rate
    if(IErate0I<=RateRef) IErate0I=RateRef

    ! If the current rate is less than the reference rate than update the current rate to be the reference rate
    if (IErateI<=RateRef) IErateI=RateRef

    ! Cond1 - Checks if the rate has moved from slower than reference to faster than reference on this time step (Rate increased)
    cond1=(IErate0I==RateRef).and.(IErateI>RateRef)

    ! Cond2 - Checks if the rate has moved from greater than the reference to slower than the reference (Rate slowed)
    cond2=(IErate0I>RateRef).and.(IErateI==RateRef)

    ! Calculate the rate increment
    dErate_eff=IErateI-IErate0I

    ! Cond3 - Check if the current and previous value is greater than the reference rate. if they are that means that rate affects should be applied
    cond3=(IErate0I>RateRef).and.(IErateI>RateRef)

    ! Check if any of the conditions are true, if so strain rate effects need to be applied
    if (cond1.or.cond2.or.cond3) Apply=.true.
 end Subroutine check4crossing

end module mod_SRMC_funcs