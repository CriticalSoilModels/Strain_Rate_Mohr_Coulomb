module mod_test_J3_invariant
   ! Local imports
   use kind_precision_module, only: dp, i32
   use mod_stress_invariants, only : calc_eig_J3, calc_inc_driver_J3_invariant, calc_Zam_J3_invariant, &
                                     calc_dev_stess, calc_mean_stress
   use stdlib_linalg, only: eigh, det

   ! Testdrive imports
   use testdrive, only : new_unittest, unittest_type, error_type, check

   implicit none

   private
   public :: collect_J3_invariant
   ! Note the convention used in incremental driver is
contains


    subroutine collect_J3_invariant(testsuite)
        ! Collection of tests
        ! Inidividual tests are stored in unitest_type

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        ! Make the test suite using the subroutines defined below
        testsuite = [ &
        new_unittest("J3_eig_invariant_check", test_J3_eig_invariant)  &

        ]

    end subroutine collect_J3_invariant

   ! Test calcualting J3 using the eigenvalues
   subroutine test_J3_eig_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local variables
      real(kind = dp) :: stress(6), dev(6), dev_matrix(3,3)
      real(kind = dp) :: mean_stress, J3_eig, J3_inc_driver, dev_eig_vals(3)
      integer :: i

      ! Incremental driver Voigt order
      !={
      !    11 (xx),
      !    22 (yy),
      !    33 (zz),
      !    12 (xy),
      !    13 (xz),
      !    23 (yz)
      !}

      stress = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp]
     
     mean_stress = calc_mean_stress(stress)

     dev = calc_dev_stess(stress, mean_stress)
     
     do i = 1, 3
        dev_matrix(i,i) = dev(i)
     end do

      dev_matrix(1, 2) = dev(4)
      dev_matrix(2, 1) = dev(4)
      dev_matrix(1, 3) = dev(5)
      dev_matrix(3, 1) = dev(5)
      dev_matrix(2, 3) = dev(6)
      dev_matrix(3, 2) = dev(6)

      ! Find the eigen values
      call eigh( dev_matrix, dev_eig_vals )

      ! Find the J3_eig value
      J3_eig = calc_eig_J3(dev_eig_vals)

      J3_inc_driver = calc_inc_driver_J3_invariant(dev)

      ! Compare the results
      call check(error, J3_inc_driver, J3_eig, thr = 1e-8_dp)
      if(allocated(error)) return

   end subroutine test_J3_eig_invariant

end module mod_test_J3_invariant
