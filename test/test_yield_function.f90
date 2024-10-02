module mod_test_yield_function
    ! Local imports
    use kind_precision_module, only: real_type => dp
    use kind_precision_module, only: i32
    use mod_yield_function, only : Get_dF_to_dSigma_2
    use mod_SRMC_funcs    , only : Get_dF_to_dSigma
    use mod_shape_checker , only : check_matrix_shape
    use mod_tensor_value_checker, only: check_tensor_values
    ! Testdrive imports
    use testdrive, only : new_unittest, unittest_type, error_type, check
 
    implicit none
 
    private
    public :: collect_yield_function_deriv
    ! Note the convention used in incremental driver is
 contains
 
 
     subroutine collect_yield_function_deriv(testsuite)
         ! Collection of tests
         ! Inidividual tests are stored in unitest_type
 
         type(unittest_type), allocatable, intent(out) :: testsuite(:)
 
         ! Make the test suite using the subroutines defined below
         testsuite = [ &
         new_unittest("dF_to_dSigma_check", test_dF_to_dSigma)  &
 
         ]
 
     end subroutine collect_yield_function_deriv
    
    subroutine test_dF_to_dSigma(error)
        ! Check the derivatives of the yield function
        type(error_type), allocatable, intent(out) :: error
        
        ! Local variables
        real(kind = real_type) :: &
            M_tc  = 1.0  , &
            eta_y = 1.5  , &
            Sig(6)       , &
            n_vec(6), &
            tol = 1e-9
        real(kind = real_type) :: exp_n_vec(6)
        logical :: passed

        Sig = [1.0, 3.0, 5.0, 7.0, 11.0, 13.0]

        ! The expected n_vector
        
        ! Call the original df_dsigma
        call Get_dF_to_dSigma(M_tc, eta_y, Sig, exp_n_vec)
        
        ! Call the updated df_dsigma
        call Get_dF_to_dSigma_2(M_tc, eta_y, Sig, n_vec)

        ! Compare the values
        call check_tensor_values(n_vec, exp_n_vec, tol, passed)
        
        ! Check for error (logical)
        call check(error, passed,.True.)
        if ( allocated(error) ) return 
        
            
    end subroutine test_dF_to_dSigma

 end module mod_test_yield_function
 