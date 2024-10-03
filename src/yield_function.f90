module mod_yield_function
   use kind_precision_module, only: real_type => dp
   use mod_stress_invariants, only : calc_inc_driver_J3_invariant, Get_invariants, calc_dev_stess, calc_J2_invariant
   use mod_stress_invar_deriv, only: calc_dp_to_dSigma, calc_dq_to_dSigma, calc_dJ3_to_dSigma, calc_Zam_dtheta_to_dSigma
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
      dpdsig = calc_dp_to_dSigma()

      dev= calc_dev_stess(Sig, p)
      
      !2) Get dq/dsig
      dqdSig = calc_dq_to_dSigma(dev, q)

      !3) Get dtheta/dSigma
      J2 = calc_J2_invariant(dev)

      J3 = calc_inc_driver_J3_invariant( dev )

      dJ3dsig = calc_dJ3_to_dSigma(dev)

      !Compute dtheta/dsig
      dthetadSig = calc_Zam_dtheta_to_dSigma(dJ3dsig, dev, J3, J2, theta)

      ! TODO: Remove this
      ! print *, "------ Lode angle derivative test -------"
      ! print *, "Lode angle: ", theta
      ! print *, "J2 invariant", J2
      ! print *, "det (s)    :", J3
      ! print *, "dJ / dsigma check: ", 1.0_real_type / (2.0_real_type * sqrt(J2)) * calc_dJ2_to_dSigma(dev)
      ! print *, "dJ3 / dsigma: ", dJ3dsig
      ! print *, "My dtheta / dsigma", calc_dtheta_to_dSigma(dJ3dsig, dev, j3, j2, theta)
      ! print *, "dtheta / dsigma", dthetadSig

      ! print *, "------ End Lode angle Derivative test ---"
      !__________________________________________________________________
      !Get n_vec=dF/dSig
      n_vec = ( eta_y * dpdsig ) + dqdSig + ( dfdtheta * dthetadSig) !n_vec=dF/dSig
   end subroutine Get_dF_to_dSigma_2

end module mod_yield_function
