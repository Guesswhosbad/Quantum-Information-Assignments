!To compile with:
!gfortran Ex1-ColomboMassimo.f90 -o Ex1.out -llapack
    
    !**********************************************************
    !**********************    INTRO    ***********************
    !**********************************************************
    !In this module variables get defined. The use of each variable is made explicit.
    
  
    module intro 

        implicit none                       !To prevent implicit variable definition. 
       
        !Loop variables: 
        integer*4 ii, jj                        
        
        !Variable used to check existence of a file with (INQUIRE):
        logical                                 :: exist, converged = .false.
        
        !Input variables:
        integer*4                               :: local_dim, Ndiag, size, max_iterations
        real*8                                  :: lambda
        real*8                                  :: eps

        !complex matrix
        complex*16, dimension(:,:), allocatable :: Ising_ham, AA, BB
        real*8, dimension(:), allocatable       :: eigenvalues
        real*8                                  :: groundstate
        complex*16, dimension(:,:), allocatable :: projector_n, Ising_2n 

        !Math constants.
        complex*16 :: i_math = (0d0, 1d0)
        real*8     :: PI = 4.D0*DATAN(1.D0)
    end module intro


    !**********************************************************
    !*****************   COMMAND LINE ARGS   ******************
    !**********************************************************

    module command_line_args
        implicit none
        character(len=32)           :: arg
        integer*4                   :: n_arg
        integer*4 :: debug_value = 0
        logical :: script = .false.
        
        contains

        !this subroutine read commandline INPUT
        !It returns the value of variables as: SCRIPT(logical), DEBUG_level(integer)
        subroutine check_commandline_args()
            do n_arg = 1, command_argument_count()
                call get_command_argument(n_arg, arg)
        
                select case (arg)
                    case ('-h', '--help')
                        call print_help()
                        stop

                    case ('-s', '--script')
                        script = .true.

                    case ('-d0', '--debug:0')
                        debug_value = 0

                    case ('-d1', '--debug:1')
                        debug_value = 1

                    case ('-d2', '--debug:2')
                        debug_value = 2
                    
                    case ('-d3', '--debug:3')
                        debug_value = 3
        
                    case default
                        print '(2a, /)', 'unrecognised command-line option: ', arg
                        call print_help()
                        stop
                end select
            end do
        end subroutine
    
        !in the previous routine, this subroutine can be called.
        !INPUT = command line ("-h",  "--help")
        !"OUTPUT = STDOUT"
        subroutine print_help()
            print '(a, /)', 'command-line options:'
            print '(a)', '  -h,  --help       print usage information and exit'
            print '(a)',    '  -d0, --debug:0    debug mode OFF'
            print '(a)',    '  -d1, --debug:1    debug mode ON: prints main variables'
            print '(a)',    '  -d2, --debug:2    debug mode ON: prints all variables'
            print '(a)',    '  -d3, --debug:3    debug mode ON: prints all variables and checkpoints'
            print '(a, /)',    '  -s,  --script     script mode ON: remove user interface. (use it together with -d0)'
        end subroutine

    end module

    !**********************************************************
    !****************   DEBUG and PRINTING   ******************
    !**********************************************************
    
    module debug_and_printing
        implicit none
        character(len=*), parameter :: fmt_complex = "(F15.8,a3,F15.8,a)"
        character(len=*), parameter :: fmt_C_mat = "(25(f10.6,a3,f9.5,a4))"
        character(len=*), parameter :: fmt_REAL_mat = "(25(f15.8))"
        integer*4                   :: iid, jjd

        contains
        
        !this subroutine prints out a STDOUT statement when debug variable has a value different from 0
        !INPUT = debug_value(integer)
        !OUTPUT = STDOUT.
        subroutine debug_mode(debug_value) !Prints debug title.
            integer*4 :: debug_value
            IF (debug_value .GT. 1) THEN 
                PRINT*, new_line('a')
                PRINT*, "-------------------------------------", new_line('a')
                IF (debug_value .EQ. 1) THEN
                    PRINT*, "--- DEBUG ON: only main variables ---", new_line('a')
                ELSE IF (debug_value .EQ. 2) THEN
                    PRINT*, "------ DEBUG ON: all variables ------", new_line('a')
                ELSE IF (debug_value .EQ. 3) THEN
                    PRINT*, "------ DEBUG ON: vars & checks ------", new_line('a')
                END IF  
                PRINT*, "-------------------------------------", new_line('a')
            END IF
        end subroutine

        !this subroutine prints out a STDOUT statement when level is lower than debug_value.
        !INPUT = statement (CHARACTER), var (class(*)) , level and debug_value (integers)
        !OUTPUT = STDOUT. the statement followed by the var (variable you want to print). 
        subroutine debug_print(statement, var, level, debug_value)
            integer*4 :: debug_value
            class(*) var(..)
            character(*) :: statement
            integer*4 :: level

            IF (level .LE. debug_value) THEN
                WRITE(*, '(a, i0)', advance = 'yes') statement
                CALL print_var(var)
                PRINT*, new_line('a')
            END IF
        end subroutine

        !This subroutine allows to print any type of datain a specified format.
        subroutine print_var(var)
            class(*) :: var(..)

            SELECT RANK (var)
                RANK(0)
                SELECT TYPE (var)
                    TYPE IS (character(len = *))
                        SELECT CASE (var)
                            CASE ('none')
                                return
                            CASE DEFAULT
                                PRINT*, var
                        END SELECT
                    
                    TYPE IS (integer(kind = 1))
                        PRINT*, var
                    TYPE IS (integer(kind = 2))
                        PRINT*, var
                    TYPE IS (integer(kind = 4))
                        PRINT*, var
                    TYPE IS (integer(kind = 8))
                        PRINT*, var
                    TYPE IS (integer(kind = 16))
                        PRINT*, var
                    
                    TYPE IS (logical(kind = 1))
                        PRINT*, var
                    TYPE IS (logical(kind = 2))
                        PRINT*, var
                    TYPE IS (logical(kind = 4))
                        PRINT*, var
                    TYPE IS (logical(kind = 8))
                        PRINT*, var
                    TYPE IS (logical(kind = 16))
                        PRINT*, var

                    TYPE IS (real(kind = 4))
                        PRINT"(f0.6)", var
                    TYPE IS (real(kind = 8))
                        PRINT"(f0.10)", var
                    TYPE IS (real(kind = 10))
                        PRINT"(f0.12)", var
                    TYPE IS (real(kind = 16))
                        PRINT"(f0.15)", var

                    TYPE IS (complex(kind = 4))
                        CALL print_C_number(cmplx(realpart(var), imagpart(var) ,8))
                    TYPE IS (complex(kind = 8))
                        CALL print_C_number(var)
                    TYPE IS (complex(kind = 16))
                        CALL print_C_number(cmplx(realpart(var), imagpart(var), 8))
        
                END SELECT

                RANK(1)
                SELECT TYPE(var)

                    TYPE IS (integer(kind = 2))
                        WRITE(*, "(a, *(i0, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"
                    TYPE IS (integer(kind = 4))
                        WRITE(*, "(a, *(i0, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"
                    TYPE IS (integer(kind = 8))
                        WRITE(*, "(a, *(i0, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"
                    TYPE IS (integer(kind = 16))
                        WRITE(*, "(a, *(i0, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"

                    
                    TYPE IS (real(kind = 4))
                        WRITE(*, "(a, *(f0.5, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"
                    TYPE IS (real(kind = 8))
                        WRITE(*, "(a, *(f0.8, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"
                    TYPE IS (real(kind = 10))
                        WRITE(*, "(a, *(f0.9, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"
                    TYPE IS (real(kind = 16))
                        WRITE(*, "(a, *(f0.10, :,  ', '))", advance='no') "( ", var(:)
                        PRINT *, ")"

                    TYPE IS (complex(kind = 4))
                        DO iid = 1, size(var)
                            CALL print_C_number(cmplx(realpart(var(iid)), imagpart(var(iid)) ,8))
                        ENDDO
                    TYPE IS (complex(kind = 8))
                        DO iid = 1, size(var)
                            CALL print_C_number(var(iid))
                        ENDDO
                    TYPE IS (complex(kind = 16))
                        DO iid = 1, size(var)
                            CALL print_C_number(cmplx(realpart(var(iid)), imagpart(var(iid)) ,8))
                        ENDDO
                END SELECT

                RANK(2)
                SELECT TYPE(var)
                    TYPE IS (real(kind = 4))
                        CALL print_REAL_matrix(real(var, 8), size(var, 1), size(var, 2))
                        
                    TYPE IS (real(kind = 8))
                        CALL print_REAL_matrix(real(var, 8), size(var, 1), size(var, 2))
                        
                    TYPE IS (real(kind = 10))
                        CALL print_REAL_matrix(real(var, 8), size(var, 1), size(var, 2))
                        
                    TYPE IS (real(kind = 16))
                        CALL print_REAL_matrix(real(var, 8), size(var, 1), size(var, 2))   
                        
                    TYPE IS (complex(kind = 4))
                        CALL print_C_matrix(cmplx(realpart(var), imagpart(var) ,8), size(var, 1), size(var, 2))
                    TYPE IS (complex(kind = 8))
                        CALL print_C_matrix(var, size(var, 1), size(var, 2))
                    TYPE IS (complex(kind = 16))
                        CALL print_C_matrix(cmplx(realpart(var), imagpart(var) ,8), size(var, 1), size(var, 2))
                END SELECT
                
                RANK DEFAULT
                    PRINT*, "Printing a (n>2)-dimensional array has not be implemented."
                

            END SELECT

        end subroutine

        subroutine print_C_number_into_file(zz, unit)
            complex*16 zz
            integer*4 unit
            IF  (imagpart(zz) >= 0) THEN
                WRITE(unit, "(F30.15,a3,F30.15,a)", advance = 'no') realpart(zz), " +", abs(imagpart(zz)), "i "
            ELSE
                WRITE(unit, "(F30.15,a3,F30.15,a)", advance = 'no') realpart(zz), " -", abs(imagpart(zz)), "i "
            END IF
        end subroutine

        subroutine print_C_number(zz)
            complex*16 zz
            IF  (imagpart(zz) >= 0) THEN
                WRITE(*, fmt_complex) realpart(zz), " +", abs(imagpart(zz)), "i "
            ELSE
                WRITE(*, fmt_complex) realpart(zz), " -", abs(imagpart(zz)), "i "
            END IF
        end subroutine

        subroutine print_C_matrix(AA, nn, mm)
            integer*4 nn, mm
            complex*16, dimension(nn, mm) :: AA
            character(len=4), DIMENSION(size(AA, 1), size(AA, 2)) :: sign  !put the "right" sign in front of the imaginary part. 
        
            WHERE(AIMAG(AA)>=0.) sign = '  +' 
            WHERE(AIMAG(AA)<0.) sign = '  -'
            
            DO iid = 1, size(AA, 1)
                WRITE(*, fmt_C_mat) (REALPART(AA(iid,jjd)), sign(iid,jjd), &
                                    ABS(IMAGPART(AA(iid,jjd))), "i   ", jjd= 1, size(AA, 2))    
            END DO
            
        end subroutine

        subroutine print_REAL_matrix(AA, nn, mm)
            integer*4 nn, mm
            real*8, dimension(nn, mm) :: AA 

            DO iid = 1, size(AA, 1)
                WRITE(*, fmt_REAL_mat) (AA(iid,jjd), jjd = 1, size(AA, 2))    
            END DO
            
        end subroutine


        !This subroutine's goal is to read and check a  defined type STDIN. 
        !It keeps on receiving input until the correct one is inserted or the maximum number of wrong ones is reached.
        !INPUT = . a specified type variable
        !OUTPUT = it associates a correct value to the previous variable.
        subroutine read_value(input)
            class(*) input
            integer*4 ierror, ii
            ii = 0
            DO
                SELECT TYPE(input)
                    TYPE IS (real(kind = 4))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (real(kind = 8))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (real(kind = 16))
                    READ(*,*, iostat = ierror) input

                    TYPE IS (integer(kind = 1))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (integer(kind = 2))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (integer(kind = 4))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (integer(kind = 8))
                    READ(*,*, iostat = ierror) input

                    TYPE IS (logical(kind = 1))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (logical(kind = 2))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (logical(kind = 4))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (logical(kind = 8))
                    READ(*,*, iostat = ierror) input
                    TYPE IS (logical(kind = 16))
                    READ(*,*, iostat = ierror) input
            
                END SELECT
                
                IF (ierror == 0) THEN    !IOSTAT = 0 -> READ was executed flawlessly, all variables have received their input values.
                    RETURN
                ELSE IF (ierror .GT. 0) THEN
                    PRINT*, "!IOSTAT > 0, the previous READ has encountered some problem."
                    PRINT*, "A common problem is illegal data: supplying a real number to an integer variable." 
                ELSE IF (ierror .LT. 0) THEN    
                    PRINT*, "IOSTAT < 0, the end of the input has reached."
                    PRINT*, "Under this circumstance, some or all of the variables in the READ may not receive input values."
                END IF

                ii = ii + 1            !ii it's used to exit the do-loop after a certain amount of wrong entries.
                IF (ii > 10) THEN
                   PRINT*, 'Too many wrong entries. exiting the program.'
                   STOP
                END IF
            END DO !crtl + c = CANCEL to force the exit. 

        end subroutine    

        subroutine assign_logical(logic, debug_value)
            logical                     :: logic
            character(50)               :: yesno
            integer*4                   :: ii, debug_value
            DO ii = 1, 10
                READ(*,*) yesno
                yesno = trim(yesno)
                SELECT CASE (yesno)
                    case("Y", "Yes", "y", "yes", "absolutely")
                        logic = .true.
                        CALL debug_print("logical state set to: ", logic, 2, debug_value)
                        exit
                    case("N","No","n","no","nope")
                        logic = .false.
                        CALL debug_print("logical state set to: ", logic, 2, debug_value)
                        exit
                    case default
                        PRINT*, "only y/yes or n/no answer allowed"
                END SELECT
                IF (ii == 10) THEN
                    PRINT*, 'Too many wrong entries. exiting the program.'
                    STOP 
                END IF

            END DO !crtl + c = CANCEL to force the exit.

        end subroutine
        
        subroutine condition_print(statement, condition)
            character(len = *)  :: statement
            logical             :: condition
            IF (condition .eqv. .true.) THEN
                PRINT*, statement
            END IF
        end subroutine

        subroutine checkpoint(debug_value) !Prints checkpoint statement.
            integer*4 :: debug_value
                IF (debug_value .GE. 3) THEN
                    PRINT*, new_line('a'), "--------- CHECKPOINT ---------", new_line('a')
                END IF
        end subroutine

        subroutine inquire_open_file(namefile, index, append)
            character(len=*) :: namefile
            integer*4        :: index
            logical          :: exist, append
            INQUIRE(file=namefile, exist=exist)
            IF (exist) THEN
                IF (append) THEN
                    OPEN(index, file=namefile, status="old", action="write", position="append")
                ELSE 
                    OPEN(index, file=namefile, status="old", action="write")
                END IF
            ELSE
            OPEN(index, file=namefile, status="new", action="write")
            END IF
        end subroutine
   
        subroutine write_in_txt_file(index, format, var, advance) ! for now only for real*8 (to extend)
            character(len=*) :: format, advance
            integer*4        :: index
            class(*)         :: var(..)
            
            SELECT RANK (var)
                RANK(0)
                SELECT TYPE (var)
                    TYPE IS (character(len = *))
                        SELECT CASE (var)
                            CASE ('none')
                                return
                            CASE DEFAULT
                                WRITE(index, format, advance = advance)  var
                        END SELECT
                    
                    TYPE IS (integer(kind = 1))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (integer(kind = 2))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (integer(kind = 4))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (integer(kind = 8))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (integer(kind = 16))
                        WRITE(index, format, advance = advance)  var
                    
                    TYPE IS (logical(kind = 1))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (logical(kind = 2))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (logical(kind = 4))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (logical(kind = 8))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (logical(kind = 16))
                        WRITE(index, format, advance = advance)  var

                    TYPE IS (real(kind = 4))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (real(kind = 8))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (real(kind = 10))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (real(kind = 16))
                        WRITE(index, format, advance = advance)  var

                    TYPE IS (complex(kind = 4))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (complex(kind = 8))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (complex(kind = 16))
                        WRITE(index, format, advance = advance)  var
        
                END SELECT

                RANK(1)
                SELECT TYPE(var)

                    TYPE IS (integer(kind = 2))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (integer(kind = 4))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (integer(kind = 8))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (integer(kind = 16))
                        WRITE(index, format, advance = advance)  var

                    
                    TYPE IS (real(kind = 4))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (real(kind = 8))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (real(kind = 10))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (real(kind = 16))
                        WRITE(index, format, advance = advance)  var

                    TYPE IS (complex(kind = 4))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (complex(kind = 8))
                        WRITE(index, format, advance = advance)  var
                    TYPE IS (complex(kind = 16))
                        WRITE(index, format, advance = advance)  var
                END SELECT

                RANK(2)
                SELECT TYPE(var)
                    TYPE IS (real(kind = 4))
                            !print_C_number_into_file(var, index)
                        
                    TYPE IS (real(kind = 8))
                        !
                        
                    TYPE IS (real(kind = 10))
                        !
                        
                    TYPE IS (real(kind = 16))
                        !   
                        
                    TYPE IS (complex(kind = 4))
                        !
                    TYPE IS (complex(kind = 8))
                        !
                    TYPE IS (complex(kind = 16))
                        !
                END SELECT
                
                RANK DEFAULT
                    PRINT*, "Printing a (n>2)-dimensional into a file has not be implemented yet XD."
            END SELECT

        end subroutine

    end module

    !**********************************************************
    !*****************   TENSOR PROD & MAT  *******************
    !**********************************************************


    module tensor_product_matrix

        interface operator(.tens.)
            module procedure matA_tens_matB
        end interface

        contains
        function matA_tens_matB(AA, BB) result(CC)
            complex*16, dimension(:,:), intent(IN)  :: AA, BB
            complex*16, allocatable, dimension(:,:) :: CC
            integer*4, dimension(2)                 :: dimAA, dimBB, dimCC
            integer*4                               :: iiA, jjA, iiB, jjB, iiC, jjC
            
            dimAA(:) = (/size(AA, 1), size(AA, 2)/)
            dimBB(:) = (/size(BB, 1), size(BB, 2)/)
            dimCC(:) = (/dimAA(1)*dimBB(1), dimAA(2)*dimBB(2)/)

            ALLOCATE(CC(dimCC(1), dimCC(2)))

            DO jjB = 0, dimBB(2) -1
                DO iiB = 0, dimBB(1) -1
                    DO jjA = 1, dimAA(2)
                        DO iiA = 1, dimAA(1)
                        iiC = iiA + iiB * dimAA(1)
                        jjC = jjA + jjB * dimAA(2)
                        CC(iiC, jjC) = AA(iiA, jjA) * BB(iiB +1, jjB +1) 
                        END DO
                    END DO
                END DO
            END DO

        end function


        !Initialize a d_to_the_N identity matrix
        function d_power_N_Id(local_dim, Nbodies) result(big_identity)
            integer*4                               :: local_dim, Nbodies, size, ii
            complex*16, allocatable, dimension(:,:) :: big_identity

            size = local_dim**Nbodies
            ALLOCATE(big_identity(size, size))
            big_identity = 0d0
            DO ii = 1, size
                big_identity(ii, ii) = 1d0
            END DO
        
        end function
       
    end module

    module real_space_RG
        use tensor_product_matrix
    
        complex*16, dimension(2, 2), parameter :: pauliZ = reshape((/ 1d0, 0d0, 0d0, -1d0/), (/ 2,2 /))
        complex*16, dimension(2, 2), parameter :: pauliX = reshape((/ 0d0, 1d0, 1d0,  0d0/), (/ 2,2 /))

        contains

        function get_Ising_Hamiltonian(Nbodies, Lambda) result(Ising_ham)
            integer*4                               :: Nbodies, size

            complex*16, allocatable, dimension(:,:) :: NI_ham, I_ham, Ising_ham
            real*8                                  :: Lambda
            integer*4                               :: NN
            
            size = 2**Nbodies
            ALLOCATE(NI_ham(size, size), I_ham(size, size), Ising_ham(size, size))

            NI_ham = 0d0
            I_ham = 0d0

            DO NN = 1, Nbodies
                NI_ham = NI_ham + (d_power_N_Id(2, NN-1).tens.(pauliZ)).tens.d_power_N_Id(2, Nbodies - NN)
            END DO

            DO NN = 1, Nbodies -1
                I_ham = I_ham + ((d_power_N_Id(2, NN-1).tens.(pauliX)).tens.pauliX).tens.d_power_N_Id(2, Nbodies - NN -1)
            END DO

            Ising_ham = Lambda * NI_ham + I_ham

            DEALLOCATE(NI_ham, I_ham)

        end function

        function get_RG_ham_2n_order(ham_n, AA, BB, local_dim, NN) result (ham_2n)
            complex*16, dimension(:,:)              :: ham_n, AA, BB
            complex*16, allocatable, dimension(:,:) :: ham_2n, term1, term2, term3
            integer*4                               :: dim, local_dim, NN, dim2n

            dim = size(ham_n, 1)
            dim2n = dim ** 2 

            ALLOCATE(ham_2n(dim2n, dim2n), term1(dim2n, dim2n), term2(dim2n, dim2n), term3(dim2n, dim2n))

            term1 = ham_n.tens.d_power_N_Id(local_dim, NN)
            term2 = d_power_N_Id(local_dim, NN).tens.ham_n 
            term3 = AA.tens.BB

            ham_2n = 0
            ham_2n = term1 + term2 + term3

            DEALLOCATE(term1, term2, term3)
        end function

    end module

    !**********************************************************
    !****************   HERMITIAN MATRICES  *******************
    !**********************************************************

    module hermitian_matrices

        contains

        function diagonalize_herm_mat(hermitian_mat) result(eigenvalues)
            complex*16, dimension(:,:)              :: hermitian_mat 
            complex*16, allocatable, dimension(:,:) :: temp_mat
            complex*16, allocatable, dimension(:)   ::  WORK, pre_WORK
            real*8, allocatable, dimension(:)       :: eigenvalues, RWORK
            integer :: INFO, LWORK, dim

            dim = size(hermitian_mat, 1)
            ALLOCATE(eigenvalues(dim), pre_WORK(1), RWORK(max(1, 3*dim-2)), temp_mat(dim, dim))
            temp_mat = hermitian_mat

            LWORK = -1  !With LWORK = -1, it optimizes the dimension of WORK array. returning pre_WORK(1) as that value.
            CALL zheev("N","U", dim, temp_mat, dim, eigenvalues, pre_WORK, LWORK, RWORK, INFO)
            
            LWORK = int(pre_WORK(1))    !After allocating WORK with optimal dimension it can compute the eigenvalues:
            ALLOCATE(WORK(MAX(1,LWORK)))
            CALL zheev("N","U", dim, temp_mat, dim, eigenvalues, WORK, LWORK, RWORK, INFO)

            DEALLOCATE(WORK, RWORK, pre_WORK)

        end function

        function get_eigenvectors(hermitian_mat) result(eigvects)
            complex*16, dimension(:,:)            :: hermitian_mat
            complex*16, allocatable, dimension(:,:)   :: eigvects, temp_mat
            integer*4                             :: dim
            real*8, allocatable, dimension(:)     :: eigenvals
            complex*16, allocatable, dimension(:) :: WORK, preWORK
            real*8, allocatable, dimension(:)     :: RWORK
            integer*4                             :: INFO, LWORK

            dim = size(hermitian_mat, 1)
            ALLOCATE(eigenvals(dim), temp_mat(dim, dim), preWORK(1), RWORK(max(1, 3*dim-2)))
            temp_mat = hermitian_mat 
            !calculate optimal values for diagonalization.
            CALL zheev("V","U", dim, temp_mat, dim, eigenvals, preWORK, -1, RWORK, INFO)
            
            !It diagonalizes.
            LWORK = int(preWORK(1))
            ALLOCATE(WORK(MAX(1,LWORK)))
            CALL zheev("V","U", dim, temp_mat, dim, eigenvals, WORK, LWORK, RWORK, INFO)

            !Save eigenvectors:
            ALLOCATE(eigvects(dim, dim))
            eigvects = temp_mat

            DEALLOCATE(eigenvals, WORK, preWORK, RWORK)

        end function

        function get_first_k_eigenvecs_instable(hermitian_mat, kk, abstol) result(eigenvecs) !Unstable, idk why. Cannot be use in iteration.
            complex*16, dimension(:,:), intent(IN)  :: hermitian_mat
            integer*4, intent(IN)                   :: kk
            real*8, intent(IN)                      :: abstol
            complex*16, dimension(:,:), allocatable :: eigenvecs, temp_mat
            
            integer*4                               :: MM, dim
            real*8, allocatable, dimension(:)       :: eigenvals
            integer*4                               :: INFO
            complex*16, allocatable, dimension(:)   :: WORK, preWORK
            integer*4, allocatable, dimension(:)    :: IWORK, preIWORK, ISUPPZ
            real*8, allocatable, dimension(:)       :: RWORK, preRWORK
        
            dim = size(hermitian_mat, 1)
            ALLOCATE(temp_mat(dim, dim))
            temp_mat = hermitian_mat
        
            ALLOCATE(eigenvals(dim), eigenvecs(dim, kk), ISUPPZ(2*max(1, kk)))
            !LWORK = max(1, 2*dim)
            !LRWORK = max(1, 24*dim)
            !LIWORK = max(1, 10*N)
            ALLOCATE(preWORK(1), preRWORK(1), preIWORK(1))
        
            CALL zheevr('V', 'I', 'U', dim, temp_mat, dim, 0d0, 0d0, 1, kk, abstol, MM, &
                         & eigenvals, eigenvecs, dim, ISUPPZ, preWORK, -1, preRWORK, -1, preIWORK, -1, INFO)
            
            ALLOCATE(WORK(int(preWORK(1))), RWORK(int(preRWORK(1))), IWORK(preIWORK(1)))
        
            CALL zheevr('V', 'I', 'U', dim, temp_mat, dim, 0d0, 0d0, 1, kk, abstol, MM, eigenvals, &
            & eigenvecs, dim, ISUPPZ, WORK, int(preWORK(1)), RWORK, int(preRWORK(1)), IWORK, int(preIWORK(1)), INFO)
        
            DEALLOCATE(WORK, preWORK, ISUPPZ, IWORK, preIWORK, RWORK, preRWORK)
        
        end function

        function project_C_mat(Cmat, projector) result(proj_mat)
            complex*16, dimension(:,:)              :: Cmat, projector
            complex*16, allocatable, dimension(:,:) :: proj_mat
            integer*4, dimension(2)                 :: dimC, dimP

            dimC = shape(Cmat)
            dimP = shape(projector)

            !Checks
            IF (dimC(2) /= dimP(1)) THEN
                PRINT*, "Projector and matrix have inconsistent dimensions"
                return
            END IF
            IF (dimP(2) .GT. dimC(2)) THEN
                PRINT*, "Projector have dimension greater than the matrix"
                stop
            END IF

            ALLOCATE(proj_mat(dimP(2), dimP(2)))

            proj_mat = 0d0
            !proj_mat =   matmul(transpose(conjg(projector)), matmul(Cmat, projector))
            proj_mat =   matmul(matmul(transpose(conjg(projector)), Cmat), projector)

        end function

    end module

    !**********************************************************
    !******************    MAIN PROGRAM    ********************
    !**********************************************************
    !This program's goal is to compute Matrix multiplication of 2 matrices of given (STDIN) dimensions
    !with three different algorithms and to measure the Cpu time required by each of them.
   
    program Ex1 
    
        use intro
        use command_line_args
        use debug_and_printing
        use real_space_RG
        use hermitian_matrices
        implicit none

        CALL check_commandline_args() 

    !dimension input
        local_dim = 2

        !INPUT: Nbodies
        CALL condition_print("Enter the number of particles you wanna start with&
            & (integer positive value only)", .NOT. script)
        CALL read_value(Ndiag) 
        IF (Ndiag .LE. 0 ) THEN
            CALL condition_print("Only positive values accepted", .NOT. script)
            stop
        END IF

        !Debug statement
        CALL debug_mode(debug_value)

    !Ising hamiltonian construction.
        size = 2**Ndiag 
        ALLOCATE(Ising_ham(size, size), AA(size, size), BB(size, size), eigenvalues(size))

        CALL condition_print("Enter the Lambda factor you want to multiply&
        & the Non-interacting Hamiltonian with (integer positive value only)", .NOT. script)
        CALL read_value(lambda)
        !lambda = 0

        Ising_ham =  get_Ising_Hamiltonian(Ndiag, lambda)
        CALL debug_print("Ising Hamiltonian:     ", Ising_ham, 2, debug_value)
       
    !REAL SPACE RENORMALIZATION GROUP
        !A and B operator declaration:
        AA = 0 !probably not needed.
        BB = 0
        AA = d_power_N_Id(local_dim, Ndiag -1).tens.pauliX
        BB = pauliX.tens.d_power_N_Id(local_dim, Ndiag -1)
        CALL debug_print("AA Operator:     ", AA, 2, debug_value)
        CALL debug_print("BB Operator:     ", BB, 2, debug_value)
        
    !Iterative algorithm:
        ii = 0
        max_iterations = 1000 
        eps = 1d-15!eps = Absolute error that is accepted to stop the converging iteration.

        CALL inquire_open_file("groundstate.txt", 1, .true.)
        CALL write_in_txt_file(1, "(f0.10, ',')", lambda, "no")

        DO WHILE (.not.converged)

            !Update(Ising_ham)
            CALL RG_algorithm(Ising_ham, AA, BB, local_dim, Ndiag)

            !Get the Eigenvalues.
            eigenvalues = diagonalize_herm_mat(Ising_ham) / Ndiag 

            !Check convergence 
            IF (ii /= 0) THEN
            converged = abs(groundstate - eigenvalues(1)) .LT. eps
            END IF
            groundstate = eigenvalues(1)
            DEALLOCATE(eigenvalues)

            !Increase count of iteration
            ii = ii +1 

            !Write groundstate of each order n of the iteration in a file  //(Not helpful for the script to read data)
            !CALL write_in_txt_file(1, "(f0.10, ',')", groundstate, "no")

            !Stop the DO WHILE after reaching "max_iteration" iterations:
            IF (ii == max_iterations) THEN
                CALL condition_print("Max number or iterations reached", .NOT. script)
                return
            END IF
        END DO

        !Save converged Eigenvalue.
        CALL write_in_txt_file(1, "(f0.10, ',')", groundstate, "no")
        CALL write_in_txt_file(1, "(i0)", ii, "no")

        CLOSE(1)

        CALL debug_print("Number of iteration:     ", ii, 1, debug_value)
        DEALLOCATE(Ising_ham, AA, BB)
    
    
        contains 

        subroutine RG_algorithm(Ising_ham, AA, BB, local_dim, Ndiag)
            complex*16, dimension(:,:), allocatable :: Ising_ham, Ising_ham1, AA, BB, AA1, BB1
            complex*16, dimension(:,:), allocatable :: projector_n, Ising_2n, eigenvectors
            integer*4                               :: local_dim, Ndiag, size, size2

            size = local_dim ** Ndiag
            size2 = local_dim ** (2 * Ndiag)
            ALLOCATE(projector_n(size2, size)) 
                !2N ISING HAMILTONIAN of order N:
                Ising_2n = get_RG_ham_2n_order(Ising_ham, AA, BB, local_dim, Ndiag)
    
                !PROJECTORS: Diagonalize it to get first "N" eigenvectors:
                eigenvectors = get_eigenvectors(Ising_2n)   !projector_n = get_first_k_eigenvecs(Ising_2n, size, 1d-10) --> with zheevr function it doesnt work :(. To be fixed.
                projector_n(:, :)  = eigenvectors(:, :size)
    
                CALL debug_print("D^2N hamiltonian:     ", Ising_2n, 3, debug_value)
                CALL debug_print("Projector:     ", projector_n, 2, debug_value)
    
                !Projecting ham_2N, AA and BB into N components to get new ham_n, AA', BB'
                Ising_ham1 = project_C_mat(Ising_2n, projector_n)
                AA1        = project_C_mat(d_power_N_Id(local_dim, Ndiag).tens.AA, projector_n)
                BB1        = project_C_mat(BB.tens.d_power_N_Id(local_dim, Ndiag), projector_n)
    
                Ising_ham = Ising_ham1 /2d0
                AA = AA1 /sqrt(2d0)
                BB = BB1 /sqrt(2d0)
    
                CALL debug_print("Ising Ham o(N):     ", Ising_ham, 2, debug_value)
                CALL debug_print("AA' Operator:     ", AA1, 2, debug_value)
                CALL debug_print("BB' Operator:     ", BB1, 2, debug_value)
                
                DEALLOCATE(AA1, BB1, Ising_ham1, Ising_2n, projector_n, eigenvectors)

        end subroutine

    end program Ex1
