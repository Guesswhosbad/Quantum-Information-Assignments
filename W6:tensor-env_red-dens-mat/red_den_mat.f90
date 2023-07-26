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
        logical exist
        
        !Input variables:
        integer*4                               :: inner_dim, Nbodies, system_index
        logical                                 :: sep, entanglement, choice_state

        !complex matrix
        complex*16, dimension(:,:), allocatable :: Dmat, Dmat_sq, red_Dmat1, red_Dmat2, red_Dmat1_sq, red_Dmat2_sq
        complex*16  trace_den, trace_den_sq, trace1, trace2, trace3, trace4 
        real*8 entropy1, entropy2

        !Cpu_time variables
        real*8 start, finish

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
        character(len=*), parameter :: fmt_C_mat = "(25(f10.5,a3,f10.5,a4))"
        integer*4                   :: iid, jjd

        contains
        
        !this subroutine prints out a STDOUT statement when debug variable has a value different from 0
        !INPUT = debug_value(integer)
        !OUTPUT = STDOUT.
        subroutine debug_mode(debug_value) !Prints debug title.
            integer*4 :: debug_value
            IF (debug_value /= 0) THEN 
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
    end module

    !**********************************************************
    !*******************   Tensor State   *********************
    !**********************************************************

    module tensor_state
        implicit none

        type MB_wave
        complex*16, allocatable, dimension(:) :: comp
        real*8                                :: norm
        integer*4                             :: dim
        logical                               :: sep        !if set to true state is separable.
        end type

        interface operator(.getnorm.)  !normalize a wave function.         Here i would like to implement an interface with 2 inputs.
        module procedure getnorm
        end interface

        interface operator(.normalize.)  !normalize a wave function.
        module procedure normalize
        end interface

        interface operator(.trace.)  !normalize a wave function.
        module procedure get_trace
        end interface

        contains

        subroutine init_mb_state(state, inner_dim, Nbodies, sep)
            type(MB_wave)   :: state
            integer*4       :: inner_dim, Nbodies
            logical         :: sep
            !random_vec variables 
            INTEGER :: iseed(4)
            real*8, dimension(3) :: uu

            CALL random_number(uu)
            !seed(:) = FLOOR((4096)*uu(:)) 
            iseed = (/FLOOR((4096)*uu(:)),1/) !Initializes random ISEED
            
            state%sep = sep

            IF (sep) THEN
                state%dim = inner_dim * Nbodies
            ELSE
                state%dim = inner_dim ** Nbodies
            END IF
            ALLOCATE(state%comp(state%dim))
            CALL zlarnv(2, iseed, state%dim, state%comp)
            state = .normalize.state
        
        end subroutine

        function getnorm(wave) result(norm) 
            implicit none
            type(MB_wave), intent(in) :: wave 

            real*8 :: normsq, norm
            integer*4 :: ii

            normsq = 0 
            DO ii = 1, size(wave%comp)
                normsq = normsq + (realpart(wave%comp(ii))**2 +imagpart(wave%comp(ii))**2)
            END DO
            norm = sqrt(normsq)

        end function

        function normalize(wave) result(norm_wave) 
            type(MB_wave), intent(in) :: wave
            type(MB_wave) :: norm_wave

            norm_wave%dim = size(wave%comp)
            norm_wave%sep = wave%sep
            ALLOCATE(norm_wave%comp(size(wave%comp)))
            norm_wave%comp(:) = wave%comp(:) / .getnorm.wave
            norm_wave%norm = .getnorm.norm_wave

        end function

        function get_density_matrix(wave) result(density_matrix)
            type(MB_wave), intent(in) :: wave
            complex*16, dimension(:,:), allocatable :: density_matrix
            complex*16, dimension(:,:), allocatable :: dual, vec
            
            allocate(density_matrix(wave%dim, wave%dim), dual(1, wave%dim), vec(wave%dim,1))
            dual(1,:) = conjg(wave%comp)
            vec(:,1) = wave%comp
            density_matrix = matmul(vec,dual)

        end function

        function get_red_density_matrix(density_matrix, inner_dim, Nbodies, system_index) result(red_density_matrix)
            complex*16, dimension(:,:), allocatable :: red_density_matrix
            complex*16, dimension(:,:)              :: density_matrix
            complex*16                              :: trace
            integer*4 inner_dim, Nbodies, system_index, AA, BB, CC, DD
            integer*4 i_v, i_d, j_d, j_v, tt !i_v = i_vector, i_d = i_dual, same for j.
            
            ALLOCATE(red_density_matrix(inner_dim**(Nbodies-1), inner_dim**(Nbodies-1)))

            AA = inner_dim ** (Nbodies - system_index)
            BB = system_index - 1
            CC = inner_dim ** BB
            DD = inner_dim ** system_index

            DO j_v = 0, AA - 1
                DO j_d = 0, AA - 1
                    DO i_v = 1, CC
                        DO i_d = 1, CC
                            trace = 0
                            DO tt = 0, inner_dim -1
                            trace = trace + density_matrix(i_v + j_v * DD + tt * CC, i_d + j_d * DD + tt * CC)
                            END DO
                            red_density_matrix(i_v + j_v * CC,  i_d + j_d * CC) = trace
                        END DO
                    END DO
                END DO
            END DO
        
        end function

        function diagonalize_herm_mat(hermitian_mat) result(eigenvalues)
            complex*16, dimension(:,:)            :: hermitian_mat
            complex*16, allocatable, dimension(:) ::  WORK, RWORK, pre_WORK
            real*8, allocatable, dimension(:)     :: eigenvalues
            integer :: INFO, LWORK, dim

            dim = size(hermitian_mat, 1)
            ALLOCATE(eigenvalues(dim), pre_WORK(1), RWORK(max(1, 3*dim-2)))

            LWORK = -1  !With LWORK = -1, it optimizes the dimension of WORK array. returning pre_WORK(1) as that value.
            CALL zheev("N","U", dim, hermitian_mat, dim, eigenvalues, pre_WORK, LWORK, RWORK, INFO)
            
            LWORK = int(pre_WORK(1))    !After allocating WORK with optimal dimension it can compute the eigenvalues:
            ALLOCATE(WORK(MAX(1,LWORK)))
            CALL zheev("N","U", dim, hermitian_mat, dim, eigenvalues, WORK, LWORK, RWORK, INFO)

        end function

        function get_Neumann_entropy(density_matrix) result(entropy) 
            complex*16, dimension(:,:)              :: density_matrix
            real*8, dimension(:), allocatable         :: eigenval
            real*8                                  :: entropy
            integer*4                               :: dim, ii

            dim = size(density_matrix, 1)
            eigenval = diagonalize_herm_mat(density_matrix)
            entropy = 0
            
            DO ii = 1, dim
                IF (eigenval(ii) .GT. 0) THEN
                entropy = entropy - real(eigenval(ii) * log(eigenval(ii)))
                END IF
            END DO

        end function

        function get_trace(Cmatrix) result(trace)
            complex*16, dimension(:,:), intent(IN)  :: Cmatrix
            complex*16                              :: trace
            integer*4                               :: dim1, dim2, ii
            dim1 = size(Cmatrix, 1)
            dim2 = size(Cmatrix, 2)
            trace = 0
            IF (dim1 == dim2) THEN
                DO ii = 1, dim1
                    trace = trace + Cmatrix(ii, ii)
                END DO
            END IF

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
        use tensor_state

        implicit none

        type(MB_wave) :: state
        
        CALL check_commandline_args() 

        !INPUT: space.
        CALL condition_print ("Enter Inner space dimension of the smallest subsystem (integer positive value only): ", .NOT. script)
        CALL read_value(inner_dim)
        IF (inner_dim .LE. 0 ) THEN
            CALL condition_print("Only positive values accepted", .NOT. script)
            stop
        END IF

        CALL condition_print("Enter the amount of bodies in the system&
            & (integer positive value only)", .NOT. script)
        CALL read_value(Nbodies) 
        IF (Nbodies .LE. 0 ) THEN
            CALL condition_print("Only positive values accepted", .NOT. script)
            stop
        END IF

        CALL condition_print("Do you want to consider a separable state? (Y/N)", .NOT. script)

        CALL assign_logical(sep, debug_value)
        !Debug statement
        CALL debug_mode(debug_value)

        !LATTICE definition (space and time)

        CALL cpu_time(start)
        CALL init_mb_state(state, inner_dim, Nbodies, sep)
        CALL cpu_time(finish)

        CALL debug_print("Inner dimension of the state is: ", inner_dim, 2, debug_value)
        CALL debug_print("The number of bodies in the system are: ", Nbodies, 2, debug_value)
        CALL debug_print("The dimension of the state results: ", state%dim, 1, debug_value)
        CALL debug_print("State components are: ", state%comp, 1, debug_value)
        CALL debug_print("State norm is: ", state%norm, 2, debug_value)
        CALL debug_print("the cpu_time [s] measured is: ", finish - start, 1, debug_value)

        IF (state%sep .eqv. .true.) THEN
            INQUIRE(file="separable_cputime.txt", exist=exist)
                IF (exist) THEN
                    OPEN(1, file="separable_cputime.txt", status="old", action="write", position="append")
                ELSE
                    OPEN(1, file="separable_cputime.txt", status="new", action="write")
                END IF
                WRITE(1, *) "# Cpu time to initialize state with inner dimension d = ", inner_dim, "and N bodies = ", Nbodies
                WRITE(1, *) finish - start
                CLOSE(1)
        ELSE
            INQUIRE(file="pure_cputime.txt", exist=exist)
                IF (exist) THEN
                    OPEN(2, file="pure_cputime.txt", status="old", action="write", position="append")
                ELSE
                    OPEN(2, file="pure_cputime.txt", status="new", action="write")
                END IF
                WRITE(2, *) "# Cpu time to initialize state with inner dimension d = ", inner_dim, "and N bodies = ", Nbodies
                WRITE(2, *) finish - start
                CLOSE(2)
        END IF 

        IF (script .eqv. .true.) THEN
            stop
        END IF 

        IF (sep) THEN
            stop
        END IF
    
    !Test qubits
        IF (inner_dim == 2 .and. Nbodies == 2 .and. (sep .eqv. .false.)) THEN
            CALL condition_print("Do you want to study a specific entangled state? (Y/N)", .NOT. script)
            CALL assign_logical(entanglement, debug_value)
            IF (entanglement) THEN 
                CALL condition_print("Do you want to study a maximally entangled&
                                    & state (Y) or a separable state (N)? (Y/N)", .NOT. script)
                CALL assign_logical(choice_state, debug_value)   
                IF (choice_state) THEN                 
                    state%comp(1) = 1
                    state%comp(2) = 0
                    state%comp(3) = 0
                    state%comp(4) = 1 
                ELSE
                    state%comp(1) = 1
                    state%comp(2) = i_math 
                    state%comp(3) = 0
                    state%comp(4) = 0 
                END IF
            END IF
            state = .normalize.state
        END IF

    !Density matrix section
        ALLOCATE(Dmat(state%dim, state%dim), Dmat_sq(state%dim, state%dim))
        Dmat = get_density_matrix(state)
        Dmat_sq = matmul(Dmat, Dmat)
        
        CALL debug_print("density matrix of the state results: ", Dmat, 0, debug_value)
        CALL debug_print("Sq. dens. mat. of the state results: ", Dmat, 1, debug_value)

        trace_den = .trace.Dmat
        trace_den_sq = .trace.Dmat_sq

        CALL debug_print("Density matrix trace = ", trace_den, 1, debug_value)
        CALL debug_print("Sq. Dens. mat. trace = ", trace_den_sq, 1, debug_value)

    !Reduced density matrix section.
        ALLOCATE(red_Dmat1(state%dim-1, state%dim-1), red_Dmat2(state%dim-1, state%dim-1))
        ALLOCATE(red_Dmat1_sq(state%dim-1, state%dim-1), red_Dmat2_sq(state%dim-1, state%dim-1))

        !If you want to choose which subsystem calculate the reduced matrix of:
        !CALL condition_print ("Respect which subsystem you want to calculate the&
        !                    &reduce density matrix to? (only integer positive value)", .NOT. script)
        !CALL read_value(system_index)

        red_Dmat1 = get_red_density_matrix(Dmat, inner_dim, Nbodies, 1)
        red_Dmat2 = get_red_density_matrix(Dmat, inner_dim, Nbodies, 2)
        red_Dmat1_sq = matmul(red_Dmat1, red_Dmat1)
        red_Dmat2_sq = matmul(red_Dmat2, red_Dmat2)
    
        CALL debug_print("Reduced density matrix for state 1 results: ", red_Dmat1, 0, debug_value)
        CALL debug_print("Squared Red. dens. mat for state 1 results: ", red_Dmat1_sq, 0, debug_value)
        CALL debug_print("Reduced density matrix for state 2 results: ", red_Dmat2, 0, debug_value)
        CALL debug_print("Squared Red. dens. mat for state 2 results: ", red_Dmat2_sq, 0, debug_value)

        trace1 = .trace.red_Dmat1
        trace2 = .trace.red_Dmat1_sq
        trace3 = .trace.red_Dmat2
        trace4 = .trace.red_Dmat2_sq

        CALL debug_print("1st reduced mat trace = ", trace1, 1, debug_value)
        CALL debug_print("Sq. 1nd red mat trace = ", trace2, 1, debug_value)
        CALL debug_print("2nd reduced mat trace = ", trace3, 1, debug_value)
        CALL debug_print("Sq. 2nd red mat trace = ", trace4, 1, debug_value)

        entropy1 = get_Neumann_entropy(red_Dmat1)
        entropy2 = get_Neumann_entropy(red_Dmat2)

        CALL debug_print("Entropy of 1st state = ", entropy1, 1, debug_value)
        CALL debug_print("Entropy of 2nd state = ", entropy2, 1, debug_value)

        
    !Deallocate matrices.
        DEALLOCATE(state%comp, Dmat, Dmat_sq, red_Dmat1, red_Dmat1_sq, red_Dmat2, red_Dmat2_sq)

    end program Ex1
