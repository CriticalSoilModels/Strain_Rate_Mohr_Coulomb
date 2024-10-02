module mod_yield_function
   use kind_precision_module, only: real_type => dp
   use mod_stress_invariants, only : calc_inc_driver_J3_invariant, Get_invariants, calc_dev_stess, calc_J2_invariant
   use mod_voigt_functions  , only : square_inc_driver_voigt_vector

   implicit none

contains

   pure function calc_dF_to_dtheta(M_tc, p, theta) result(dfdtheta)
      real(kind = real_type), intent(in) :: M_tc, p, theta
      real(kind = real_type) :: dfdtheta

      ! Local variables
      real(kind = real_type), parameter :: PI = 4.0d0 * atan(1.0)

      !Get dF/dtheta
      dfdtheta = 0.45 * p * M_tc * ( ( cos(1.5d0 * theta + 0.25d0 * PI) )**0.2) * sin( 1.5d0 * theta + 0.25d0 * PI)

   end function calc_dF_to_dtheta

   pure function calc_dJ3_to_dSigma(dev) result(dJ3dSig)
      real(kind = real_type), intent(in) :: dev(6)
      real(kind = real_type) :: dJ3dSig(6)

      ! Local variables
      real(kind=  real_type) :: II(6), dev2(6), TrS2

      !Fill S.S
      dev2 = square_inc_driver_voigt_vector( dev )

      !Compute dJ3dSig
      TrS2 = dev2(1) + dev2(2) + dev2(3)

      II=0.0d0!Identity tensor
      II(1)=1.0
      II(2)=1.0
      II(3)=1.0

      ! This is equaivalent to s^{2} - 2/3 J_{2} \matr{1}
      ! J_{2}(\matr{s}) = 2 * I_{1}(\matr{ s^{2} })
      ! See Appendix B. Invariant Notes Moore, Jonathan Thesis for more details
      dJ3dsig = dev2 - ( TrS2*II / 3.0d0 )

   end function calc_dJ3_to_dSigma

   subroutine Get_dF_to_dSigma_2(M_tc, eta_y, Sig, n_vec)
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
      double precision:: p, q, theta, &
         J2, J3, dJ3dsig(6), dfdtheta, &
         dpdsig(6), dqdsig(6), dev(6), dthetadSig(6)

      !Get the invariants
      call Get_invariants(Sig, p, q, theta)

      !Get dF/dp=eta_y and dF/dq=1
      !Get dF/dtheta
      dfdtheta= calc_dF_to_dtheta(M_tc, p, theta)

      !___________________________________________________________________________
      !1) Get dp/dSig=1/3 Imat
      dpdsig=0.0d0
      dpdsig(1)=1.0/3.0
      dpdsig(2)=1.0/3.0
      dpdsig(3)=1.0/3.0

      !2) Get dq/dsig= 2 *dev/3*q
      dev= calc_dev_stess(Sig, p)

      ! TODO: The shear terms aren't being doubled. They need to be doubled because of the symmetry of the matrix
      ! dqdSig=( 3.0 / ( 2.0*q ) ) * dev
      dqdSig = calc_dq_to_dSigma(dev, q)

      !3) Get dtheta/dSigma= (1/3cos3theta) d/dsigma((J3/2) * (3/J2)^1.5)
      J2 = calc_J2_invariant(dev)

      J3 = calc_inc_driver_J3_invariant( dev )

      dJ3dsig = calc_dJ3_to_dSigma(dev)

      !Compute dtheta/dsig
      dthetadSig = calc_Zam_dtheta_to_dSigma(dJ3dsig, dev, J3, J2, theta)

      print *, "------ Lode angle derivative test -------"
      print *, "Lode angle: ", theta
      print *, "J2 invariant", J2
      print *, "det (s)    :", J3
      print *, "dJ / dsigma check: ", 1.0_real_type / (2.0_real_type * sqrt(J2)) * calc_dJ2_to_dSigma(dev)
      print *, "dJ3 / dsigma: ", dJ3dsig
      print *, "My dtheta / dsigma", calc_dtheta_to_dSigma(dJ3dsig, dev, j3, j2, theta)
      print *, "dtheta / dsigma", dthetadSig

      print *, "------ End Lode angle Derivative test ---"
      !__________________________________________________________________
      !Get n_vec=dF/dSig
      n_vec = ( eta_y * dpdsig ) + dqdSig + ( dfdtheta * dthetadSig) !n_vec=dF/dSig
   end subroutine Get_dF_to_dSigma_2

   pure function calc_dq_to_dSigma(dev, q) result(dq_dSigma)
      ! Calc dq/dSigma
      real(kind = real_type), intent(in) :: dev(6), q
      real(kind = real_type) :: dq_dSigma(6)

      ! Local variables
      real(kind = real_type) :: dJ2_dsigma(6)

      ! Calc dJ2/dsigma
      dJ2_dSigma = calc_dJ2_to_dSigma(dev)

      ! Calc dq/dSigma
      dq_dSigma = 3.0_real_type/ (2.0_real_type * q) * dJ2_dsigma

   end function calc_dq_to_dSigma

   pure function calc_dJ2_to_dSigma(dev) result(dJ2_dSigma)
      real(kind = real_type), intent(in) :: dev(6)
      real(kind = real_type) :: dJ2_dSigma(6)

      dJ2_dSigma = dev

      ! Double the shear terms
      dJ2_dSigma(4:6) = 2.0_real_type * dJ2_dSigma(4:6)

   end function calc_dJ2_to_dSigma


   function calc_Zam_dtheta_to_dSigma(dJ3dsig, dev, J3, J2, theta) result(dthetadSig)
      ! Function to calcute dtheta / dSigma following the formula Zambrano originally had
      ! This formula is differnt than what Pott's has
      real(kind = real_type) :: dJ3dSig(6), dev(6), J3, J2, theta
      real(kind = real_type) :: dthetadSig(6)

      ! Local variables
      real(kind = real_type) :: COS_3THETA

      dthetadSig = dJ3dsig - ( 1.5*J3 / J2 ) * dev
      COS_3THETA = cos( 3.0*theta )
      dthetadSig = ( sqrt(3.0) / ( 2.0*COS_3THETA*J2**1.5 ) ) * dthetadSig
   end function calc_Zam_dtheta_to_dSigma

   function calc_dtheta_to_dSigma(dJ3_dSigma, dev, J3, J2, theta) result(dtheta_dSigma)
      real(kind = real_type) :: dJ3_dSigma(6), dev(6)
      real(kind = real_type) :: J3, J2, theta
      real(kind = real_type) :: dtheta_dSigma(6)

      ! Local variables
      real(kind = real_type) :: cos_term, outside_term, inside_term_1(6), dJ2_dSigma(6)
      real(kind = real_type), parameter :: tolerance = 1e-12_real_type
      real(kind = real_type), parameter :: THREE = 3.0_real_type, &
                                           TWO   = 2.0_real_type, &
                                           ZERO  = 0.0_real_type
      ! Calc cos(3 \theta)
      cos_term = cos(3 * theta)

      ! If cos term is zero( Trx compression or tension) set to tiny value
      if ( abs(cos_term) <= tolerance) then
         cos_term = tolerance
      end if

      ! Calc the fraction before the parenthesis
      outside_term = sqrt(THREE) / (TWO * cos_term * J2**1.5_real_type)

      ! Calc the first term inside the parenthesis
      dJ2_dSigma = calc_dJ2_to_dSigma(dev)

      inside_term_1 = THREE * J3 / (2 * J2) * dJ2_dSigma

      ! Calc the full term
      dtheta_dSigma = outside_term * (inside_term_1 - dJ3_dSigma)

   end function calc_dtheta_to_dSigma

end module mod_yield_function
