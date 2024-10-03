! Module for holding the derivatives of the strain invariants

module mod_strain_invar_deriv
   use kind_precision_module, only: real_type => dp

   implicit none

contains
   subroutine Get_dEpsq_to_dEps(Epsq, Eps, dEqdEpsq)
      !************************************************************************
      ! Returns the derivative of the deviatoric strain with respect to the   *
      ! deviatoric strain	tensor						     					*
      ! dEqdEpsq is a (1X6) vector											*
      !************************************************************************
      implicit none
      !input
      real(kind = real_type), intent(in):: Epsq, Eps(6)
      !output
      real(kind = real_type), intent(out):: dEqdEpsq(6)
      !local variables
      real(kind = real_type):: evol, dev(6)

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
end module mod_strain_invar_deriv
