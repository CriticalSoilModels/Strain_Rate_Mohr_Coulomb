module test_invariants_suite
   ! Local imports
   use kind_precision_module, only : dp, i32
   use mod_SRMC_funcs       , only : Get_strain_invariants
   use mod_stress_invariants, only : Get_invariants, calc_theta_s, calc_J2_invariant, &
                                     calc_inc_driver_J3_invariant, calc_dev_stess, calc_mean_stress
   
   ! Testdrive imports
   use testdrive, only : new_unittest, unittest_type, error_type, check

   implicit none

   private
   public :: collect_invariants_suite
   ! Note the convention used in incremental driver is
contains


   subroutine collect_invariants_suite(testsuite)
      ! Collection of tests
      ! Inidividual tests are stored in unitest_type

      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      ! Make the test suite using the subroutines defined below
      testsuite = [ &
         new_unittest("strain_invariant_check", test_strain_invariant), &
         new_unittest("lode_angle_check"      , test_lode_angle )     , &
         new_unittest("stress_invariant_check", test_stress_invariant)  &

         ]

   end subroutine collect_invariants_suite

   subroutine test_strain_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local varaibles
      real(kind = dp) :: strain(6), eps_v, eps_q
      real(kind = dp), parameter :: exp_eps_v = 9.0_dp               , &
         exp_eps_q = 10.878112977260887_dp, &
         tol       = epsilon(eps_v)

      strain = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp]

      ! Calc the strain invariant for the predifed strain (in voigt notation)
      call Get_strain_invariants(strain, eps_v, eps_q)

      ! Check if the strain invariants meet the expected values
      call check(error, eps_v, exp_eps_v)
      if(allocated(error)) return

      call check(error, eps_q, exp_eps_q)
      if(allocated(error)) return

   end subroutine test_strain_invariant

   subroutine test_lode_angle(error)
      ! Testing the Lode angle calculation for different stress conditions
      type(error_type), allocatable, intent(out) :: error
   
      ! Define stress states and the value of PI
      real(kind = dp), parameter :: &
         trx_compression(6) = [77.0_dp, 16.0_dp, 16.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         trx_extension(6)   = [5.0_dp, 5.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         shear(6)           = [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         PI = 4.0_dp * atan(1.0_dp)   ! Value of PI
   
      real(kind = dp) :: compress_lode, exten_lode, shear_lode
      real(kind = dp) :: dev(6), mean_stress
   
      ! Test for triaxial compression (expected Lode angle is -pi/6)
      mean_stress = calc_mean_stress(trx_compression)
      dev = calc_dev_stess(trx_compression, mean_stress)
      compress_lode = calc_theta_s(calc_J2_invariant(dev), calc_inc_driver_J3_invariant(dev))

      call check(error, compress_lode, -PI / 6.0_dp, more = "compression lode angle")
      if (allocated(error)) return
   
      ! Test for triaxial extension (expected Lode angle is pi/6)
      mean_stress = calc_mean_stress(trx_extension)
      dev = calc_dev_stess(trx_extension, mean_stress)
      exten_lode = calc_theta_s(calc_J2_invariant(dev), calc_inc_driver_J3_invariant(dev))

      call check(error, exten_lode, PI / 6.0_dp, more = "extension lode angle")
      if (allocated(error)) return
   
      ! Test for shear stress (expected Lode angle is 0)
      mean_stress = calc_mean_stress(shear)
      dev = calc_dev_stess(shear, mean_stress)
      shear_lode = calc_theta_s(calc_J2_invariant(dev), calc_inc_driver_J3_invariant(dev))

      call check(error, shear_lode, 0.0_dp, more = "shear lode angle")
      if (allocated(error)) return
   
   end subroutine test_lode_angle
   
   subroutine test_stress_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local varaibles
      real(kind = dp) :: stress(6), p, q, theta
      real(kind = dp), parameter :: &
         exp_p     = 3.0_dp , &
         exp_q     = 32.07802986469088_dp

      real(kind = dp) :: exp_theta, dev(6), J3, J2

      stress = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp]

      ! Calc the stress invariants
      call Get_invariants(stress, p, q, theta)

      ! Check the mean stress
      call check(error, p, exp_p, more = "Mean Stress Test^")
      if (allocated(error)) return

      ! Check the dev. stress, (q, von mises stress)
      call check(error, q, exp_q, more = "Dev. Stress Test^")
      if(allocated(error)) return

      ! Check the lode angle
      ! See function calc_theta_s for more information about the lode angle used here
   
      dev = calc_dev_stess(stress, p)
      J3 = calc_inc_driver_J3_invariant(dev)
      J2 = calc_J2_invariant(dev)
      ! Formula from Potts and Zdravković
      exp_theta = - asin(3.0_dp * sqrt(3.0_dp) / 2.0_dp * J3 / J2**1.5_dp) / 3.0_dp
      
      call check(error, theta, exp_theta)
      if(allocated(error)) return
   end subroutine test_stress_invariant

end module test_invariants_suite
