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
        integer*4, allocatable, dimension(:)    :: index_arr   
        
        !Variable used to check existence of a file with (INQUIRE):
        logical                                 :: exist, converged = .false.
        
        !Input variables:
        integer*4                               :: local_dim, Nblock, block_dim, max_iterations
        integer*4                               :: md_dim, Nparticles
        real*8                                  :: lambda, eps, groundstate

        !Big hamiltonian 4 blocks variables
        complex*16, dimension(:,:), allocatable :: H1, H2, H3, H4, H12, H23, H34, Htot
        !Reduced density matrix variables
        real*8, dimension(:), allocatable       :: eigenvalues
        
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

        subroutine comparing_H_matrix(AA, BB, eps, level, debug_value)
            integer*4           :: ii, jj, err_counts, level, debug_value
            complex*16          :: AA(:,:), BB(:,:)
            real*8              :: eps
            err_counts = 0

            IF (size(AA, 1) /= size(BB, 1)) THEN
                PRINT*, "Matrices have different dimensions."
                return
            END IF

            DO jj = 1, size(AA, 1)
                DO ii = 1, size(AA, 1)             
                   IF(ABS(AA(ii,jj) - BB(ii,jj)) > eps) THEN
                        PRINT*, 'Entry (',ii, ',', jj,')', 'of the matrices are different'
                        err_counts = err_counts + 1
                    ENDIF
                ENDDO
            ENDDO
            IF (level .LE. debug_value) THEN
            IF(err_counts /= 0) THEN
                PRINT*, new_line('a'), 'There are', err_counts, 'entries in the two matrices that differ more than a value of', eps
            ELSE
                PRINT*, 'Matrices result the same up to a value', eps
            ENDIF
    
            END IF
        endsubroutine 

    end module

    !**********************************************************
    !**********************   SORTING  ************************
    !********************************************************** . !quick sort, merge sort

    module sorting

        contains
        
        !real*8 for now
        subroutine sort_ascending(array, index_arr)
            real*8, dimension(:)                    :: array
            integer*4, allocatable, dimension(:)    :: index_arr
            real*8                                  :: tt
            integer*4                               :: dim, ii, jj, ind

            dim = size(array)
            ALLOCATE(index_arr(dim))
            index_arr = (/(i, i=1,dim, 1)/)

            DO ii=1, dim
                DO jj= 1, dim -1
                    IF (array(ii) .LT. array(jj)) THEN
                    tt = array(ii) 
                    array(ii) = array(jj)
                    array(jj) = tt

                    ind = index_arr(ii)
                    index_arr(ii) = index_arr(jj)
                    index_arr(jj) = ind
                    END IF
                END DO
            END DO

        end subroutine

        subroutine sort_descending(array, index_arr)
            real*8, dimension(:)                    :: array
            integer*4, allocatable, dimension(:)    :: index_arr
            real*8                                  :: tt
            integer*4                               :: dim, ii, jj, ind

            dim = size(array)
            !ALLOCATE(index_arr(dim))
            index_arr = (/(i, i=1,dim, 1)/)

            DO ii=1, dim
                DO jj= 1, dim
                    IF (array(ii) .GT. array(jj)) THEN
                    tt = array(ii) 
                    array(ii) = array(jj)
                    array(jj) = tt

                    ind = index_arr(ii)
                    index_arr(ii) = index_arr(jj)
                    index_arr(jj) = ind
                    END IF
                END DO
            END DO

        end subroutine

        subroutine sort_with_index_arr(array, index_arr)
            real*8, dimension(:)                    :: array
            real*8, dimension(:), allocatable       :: temp_array
            integer*4, dimension(:)                 :: index_arr
            integer*4                               :: dim, ii

            dim = size(array)
            IF (dim /= size(index_arr)) THEN
                Print*, "ERROR: Arrays with different size encountered"
                return
            END IF
            ALLOCATE(temp_array(dim))
            temp_array = array
            DO ii = 1, dim
                array(ii) = temp_array(index_arr(ii))
            END DO

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

    !**********************************************************
    !***************   INFINITE DMRG - ISING  *****************
    !**********************************************************

    module infDMRG
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
    !**************   WAVEFUNCTION & DENS MAT  ****************
    !**********************************************************

    module wavefunction
        use hermitian_matrices

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

        function get_red_density_matrix_tracing_out_Ith_system(density_matrix, inner_dim, Nbodies, system_index)&
            & result(red_density_matrix)
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

        function get_red_density_matrix_first_kk_subsystems(density_matrix, local_dim, Ntotalbodies, kk) &
                    & result(red_density_matrix)
            complex*16, dimension(:,:), allocatable :: red_density_matrix, temp_red1, temp_red2
            complex*16, dimension(:,:)              :: density_matrix
            integer*4       local_dim, Ntotalbodies, kk, dim_temp1, dim_temp2, ii

            ALLOCATE(red_density_matrix(local_dim**kk, local_dim**kk))
            
            DO ii = 1, Ntotalbodies - kk !fare subroutine
                dim_temp1 = (local_dim ** Ntotalbodies) / (local_dim**(ii -1))
                dim_temp2 = (local_dim ** Ntotalbodies) / (local_dim**ii)
                ALLOCATE(temp_red1(dim_temp1, dim_temp1))
    
                IF (ii == 1) THEN 
                    temp_red1 = density_matrix
                ELSE 
                    temp_red1 = temp_red2
                    DEALLOCATE(temp_red2)
                END IF
    
                temp_red2 = get_red_density_matrix_tracing_out_Ith_system(temp_red1, &
                            & local_dim, Ntotalbodies - (ii -1), Ntotalbodies - (ii -1))
                
                DEALLOCATE(temp_red1)
                IF (ii == Ntotalbodies - kk) THEN 
                    IF (size(temp_red2, 1) /= size(red_density_matrix, 1)) THEN
                        Exit
                    END IF
                    red_density_matrix = temp_red2
                    DEALLOCATE(temp_red2)
                END IF
    
            END DO
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
   
    program InfiniteDMRG
        
        use intro
        use command_line_args
        use debug_and_printing
        use sorting
        use infDMRG
        use wavefunction
        
        implicit none

        CALL check_commandline_args() 

    !dimension input
        local_dim = 2

        !INPUT: Nbodies, Lambda.
        CALL condition_print("Enter the number of particles of the 1st block&
            & (integer positive value only)", .NOT. script)
        CALL read_value(Nblock) 
        IF (Nblock .LE. 0 ) THEN
            CALL condition_print("Only positive values accepted", .NOT. script)
            stop
        END IF

        CALL condition_print("Enter the Lambda factor you want to multiply&
        & the Non-interacting Hamiltonian with (integer positive value only)", .NOT. script)
        CALL read_value(lambda)
        !lambda = 5d0

        !Debug statement
        CALL debug_mode(debug_value)

    !Ising hamiltonian allocation and declaration (Split in 4 blocks).
        block_dim = local_dim**Nblock
        md_dim = block_dim * local_dim
        CALL debug_print("Dimension of 1st block:     ", block_dim, 1, debug_value)

        ALLOCATE(H1(block_dim, block_dim), H4(block_dim, block_dim))!, eigenvalues(sysdim))
        ALLOCATE(H2(local_dim, local_dim), H3(local_dim, local_dim))
        ALLOCATE(H12(block_dim, block_dim), H23(local_dim*local_dim, local_dim*local_dim), H34(block_dim, block_dim))
        ALLOCATE(Htot(md_dim**2, md_dim**2))

        H1 =  get_Ising_Hamiltonian(Nblock, lambda)                          !dim = Nblock
        H2 =  get_Ising_Hamiltonian(1, lambda)                               !dim = 1
        H3 =  H2                                                             !dim = 1
        H4 =  H1                                                             !dim = Nblock
        H12 = (d_power_N_Id(local_dim, (Nblock -1)).tens.pauliX) !.tens.pauliX !dim = Nblock
        H23 = pauliX.tens.pauliX                                             !dim = 2
        H34 = pauliX.tens.d_power_N_Id(local_dim, (Nblock -1))      !(pauliX.tens.pauliX).tens.d_power_N_Id(local_dim, (Nblock -1)) !dim = Nblock
        Htot = 0

        CALL debug_print("Ham 1:     ", H1, 2, debug_value)
        CALL debug_print("Ham 2:     ", H2, 2, debug_value)
        CALL debug_print("Ham 3:     ", H3, 2, debug_value)
        CALL debug_print("Ham 4:     ", H4, 2, debug_value)
        CALL debug_print("Ham 12:     ", H12, 2, debug_value)
        CALL debug_print("Ham 23:     ", H23, 2, debug_value)
        CALL debug_print("Ham 34:     ", H34, 2, debug_value)

    !Iterative algorithm:
        ii = 0
        CALL condition_print("Enter the maximum number of iterations&
        & (integer positive value only)", .NOT. script)
        CALL read_value(max_iterations)
        !max_iterations = 1
        eps = 1d-7 !eps = Absolute error that is accepted to stop converging iteration.
        Nparticles = (Nblock+1) * 2

        CALL inquire_open_file("groundstate_DMRG_test.txt", 2, .true.)
        CALL write_in_txt_file(2, "(f0.10, ',')", lambda, "no")

    !ITERATION
       DO WHILE (.not.converged)
            CALL debug_print("ii = ", ii, 1, debug_value)
            CALL infiniteDMRG_algorithm(H1, H2, H3, H4, H12, H23, H34, Htot, local_dim, Nblock, Nparticles)

            !SAVING iteration of order n Eigenvalues.
            eigenvalues = diagonalize_herm_mat(Htot) / Nparticles
            !Check convergence
            IF (ii .gt. 10) THEN
                converged = abs(groundstate - eigenvalues(1)) .LT. eps
            ELSE IF (ii == 0) THEN
                CALL debug_print("Check 1st loop Htot ", "none", 1, debug_value)
                CALL comparing_H_matrix(Htot, get_Ising_Hamiltonian(Nparticles, lambda), eps, 1, debug_value)
                CALL debug_print("Htot ", Htot, 4, debug_value)
                CALL debug_print("correct Htot ", get_Ising_Hamiltonian(Nparticles, lambda), 4, debug_value)

            END IF
            groundstate = eigenvalues(1)
            CALL debug_print("Groundstate value:     ", groundstate, 1, debug_value)
            DEALLOCATE(eigenvalues)
                ii = ii +1
            Nparticles = Nparticles +2

            !Stop the iteration after max_iteration times
            IF (ii == max_iterations) THEN
                CALL condition_print("Max number or iterations reached", .NOT. script)
                converged = .true.
            END IF
        END DO

        !Write groundstate of order n in a file
        CALL write_in_txt_file(2, "(f0.10, ',')", groundstate, "no")
        CALL write_in_txt_file(2, "(i0)", ii, "no")
        

        CLOSE(2)

        CALL debug_print("Number of iteration:     ", ii, 1, debug_value)
        DEALLOCATE(H1, H2, H3, H4, H12, H23, H34, Htot)



        contains 

        subroutine infiniteDMRG_algorithm(H1, H2, H3, H4, H12, H23, H34, Htot, local_dim, Nblock, Nparticles)
            type(MB_wave) :: groundstate
            complex*16, dimension(:,:)              :: H1, H2, H3, H4, H12, H23, H34, Htot
            complex*16, dimension(:,:), allocatable :: eig_vec_4blocks
            complex*16, dimension(:,:), allocatable :: eig_vec_rdm, dens_mat, red_dens_mat, projector, temp_eigvec
            complex*16, dimension(:,:), allocatable :: H1proj, H12proj!, H34proj, H4proj, temp_eigvec2, projector2
            real*8, dimension(:), allocatable       :: red_eigenvalues, Htot_eigenvalues
            integer*4                               :: local_dim, Nblock, block_dim, Nparticles, ind
    

            block_dim = local_dim**Nblock

            CALL debug_print("Nparticles:     ", Nparticles, 2, debug_value)

            CALL debug_print("Ham 1:     ", H1, 2, debug_value)
            CALL debug_print("Ham 2:     ", H2, 2, debug_value)
            CALL debug_print("Ham 3:     ", H3, 2, debug_value)
            CALL debug_print("Ham 4:     ", H4, 2, debug_value)
            CALL debug_print("Ham 12:     ", H12, 2, debug_value)
            CALL debug_print("Ham 23:     ", H23, 2, debug_value)
            CALL debug_print("Ham 34:     ", H34, 2, debug_value)
           
        !Constructing Htot from Hinteraction and H non interactiong (4 blocks structure)
            Htot = 0
            Htot = H1.tens.d_power_N_Id(local_dim, 2 + Nblock)                                                  !H1
            Htot = Htot + (d_power_N_Id(local_dim, Nblock).tens.H2.tens.d_power_N_Id(local_dim, Nblock +1))     !H2
            Htot = Htot + ((d_power_N_Id(local_dim, Nblock +1).tens.H3).tens.d_power_N_Id(local_dim, Nblock))   !H3
            Htot = Htot + ((d_power_N_Id(local_dim, 2 + Nblock)).tens.H4)                                       !H4
            Htot = Htot + ((H12.tens.pauliX).tens.(d_power_N_Id(local_dim, Nblock +1)))                         !H12
            Htot = Htot + (d_power_N_Id(local_dim, Nblock).tens.H23.tens.d_power_N_Id(local_dim, Nblock))       !H23
            Htot = Htot + (d_power_N_Id(local_dim, Nblock + 1).tens.pauliX.tens.H34)                            !H34
            !Htot = Htot / Nparticles

            !CALL debug_print("Htot:     ", Htot, 1, debug_value)
            IF (debug_value .GE. 3) THEN
                PRINT*, "Dim Htot = (m * d)^2 = ", (block_dim*local_dim)**2 
            END IF
            CALL debug_print("Dimension Htot:     ", shape(Htot), 2, debug_value)
            
        !Density matrix from GroundState
            !GroundState:
            ALLOCATE(groundstate%comp((block_dim*local_dim)**2))

            eig_vec_4blocks = get_eigenvectors(Htot)
            Htot_eigenvalues = diagonalize_herm_mat(Htot)
            groundstate%comp = eig_vec_4blocks(:, 1)
            groundstate = .normalize.groundstate !This step is needed because .normalize. set wave%dim, %norm, %ecc that could be use afterward.

            CALL debug_print("Htot Eigenvalues:     ", Htot_eigenvalues, 2, debug_value)
            CALL debug_print("Groundstate (after .normalize.):     ", groundstate%comp, 3, debug_value)
            CALL debug_print("Norm of the Groundstate:     ", .getnorm.groundstate, 2, debug_value)
        
            !Density matrix:
            !ALLOCATE(dens_mat((block_dim*local_dim)**2,(block_dim*local_dim)**2)) !It doesnt need to allocate because the subroutine does it already.
            dens_mat = get_density_matrix(groundstate)

            !CALL debug_print("Density matrix:    ", dens_mat, 3, debug_value)
            CALL debug_print("Density matrix shape:     ", shape(dens_mat), 3, debug_value)
            CALL debug_print("Density matrix vs squared Dens mat:    ", "none", 2, debug_value)
            CALL comparing_H_matrix(dens_mat, matmul(dens_mat, dens_mat), 1d-15, 2, debug_value)

        !Reduced density matrix for 1st block +1 subsystem
            red_dens_mat = get_red_density_matrix_first_kk_subsystems(dens_mat, local_dim, 2*(Nblock +1), Nblock +1)

            CALL debug_print("red_dens_mat:    ", red_dens_mat, 2, debug_value)
            CALL debug_print("TRACE of Reduced Density Matrix:    ", .trace.red_dens_mat, 2, debug_value)

        !Truncation
            ALLOCATE(projector(block_dim * local_dim, block_dim), red_eigenvalues(block_dim * local_dim))
            ALLOCATE(temp_eigvec(block_dim * local_dim, block_dim * local_dim))
            !ALLOCATE(temp_eigvec2(block_dim * local_dim, block_dim), projector2(block_dim * local_dim, block_dim))
            
            eig_vec_rdm = get_eigenvectors(red_dens_mat)
            !projector2 = eig_vec_rdm(:, (block_dim + 1):)
            red_eigenvalues = diagonalize_herm_mat(red_dens_mat)
            CALL debug_print("Eigenvalues of red density matrix:    ", red_eigenvalues, 3, debug_value)

            CALL sort_descending(red_eigenvalues, index_arr)
            
            CALL debug_print("Eigenvectors reduced density matrix:    ", eig_vec_rdm, 3, debug_value)
            temp_eigvec = eig_vec_rdm
            !temp_eigvec2 = projector2

            !DO ind = 1, block_dim * local_dim
            !    projector2 (ind, :) = temp_eigvec2(index_arr(ind), :)
            !END DO

            DO ind = 1, block_dim * local_dim
                eig_vec_rdm (:, ind) = temp_eigvec(:, index_arr(ind))
            END DO

            projector = eig_vec_rdm(:, :block_dim)

            CALL debug_print("Sorted Eigenvectors red. density matrix:    ", eig_vec_rdm, 2, debug_value)
            CALL debug_print("Ordered Eigenvalues of red density matrix:    ", red_eigenvalues, 2, debug_value)
            CALL debug_print("index array:    ", index_arr, 2, debug_value)
            CALL debug_print("Projector on block dimension:    ", projector, 2, debug_value)
            !CALL debug_print("Projector2 on block dimension:    ", projector2, 2, debug_value)

        !Projecting new H1', H2', ecc... adding 1 subsystem.
            H1proj  = (H1.tens.d_power_N_Id(local_dim, 1)) + (H12.tens.pauliX) + (d_power_N_Id(local_dim, Nblock).tens.H2)
            !H12proj = d_power_N_Id(local_dim, 1).tens.H12
            H12proj = (H12.tens.d_power_N_Id(local_dim, 1))
           

            H1 = project_C_mat(H1proj, projector)
            !H2 = H2
            !H3 = H3
            H4 = H1
            !H4 = project_C_mat(H4proj, projector) !for simmetry reasons
            H12 = project_C_mat(H12proj, projector)
            !H23 = H23
            !H34 = project_C_mat(H34proj, projector2)
            H34 = H12

            !H1 = H1 / Nparticles
            !H2 = H2 / Nparticles
            !H3 = H3 / Nparticles
            !H4 = H4 / Nparticles
            !H12 = H12 !/ Nparticles
            !H23 = H23 !/ Nparticles
            !H34 = H34 !/ Nparticles
            
            DEALLOCATE(eig_vec_4blocks, eig_vec_rdm, dens_mat, red_dens_mat, temp_eigvec)
            DEALLOCATE(projector, H1proj, H12proj, red_eigenvalues, Htot_eigenvalues) !, H34proj, H4proj
        end subroutine

    end program InfiniteDMRG

