   
!In this Program lapack routines have been used, so to compile:
!gfortran name.f90 -0 name.out -llapack
    
    
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
        
        !Lattice characteristics:
        integer*4                               :: NN   !number of points (-1)
        real*8                                  :: x_min, x_max, dx, omega
        real*8, allocatable, dimension(:)       :: ax_grid
       
        !Matrices definition:
        complex*16, allocatable, dimension(:,:) :: AA

        !Variables to get eigenvalues and eigenvectors
        real*8, allocatable, dimension(:)   :: eigenvalues, WORK, diagonal, subdiagonal
        real*8, allocatable, dimension(:,:) :: eigenvectors, norm_eigenvectors
        integer                             :: INFO

        !Trace function variables.
        real*8                              :: sum_eigen = 0, trace = 0
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
                WRITE(*, '(a, i0)', advance = 'no') statement
                CALL print_var(var)
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
                        CALL print_C_mumber(cmplx(realpart(var), imagpart(var) ,8))
                    TYPE IS (complex(kind = 8))
                        CALL print_C_mumber(var)
                    TYPE IS (complex(kind = 16))
                        CALL print_C_mumber(cmplx(realpart(var), imagpart(var), 8))
        
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
                            CALL print_C_mumber(cmplx(realpart(var(iid)), imagpart(var(iid)) ,8))
                        ENDDO
                    TYPE IS (complex(kind = 8))
                        DO iid = 1, size(var)
                            CALL print_C_mumber(var(iid))
                        ENDDO
                    TYPE IS (complex(kind = 16))
                        DO iid = 1, size(var)
                            CALL print_C_mumber(cmplx(realpart(var(iid)), imagpart(var(iid)) ,8))
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

        subroutine print_C_mumber(zz)
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
    
        subroutine condition_print(statement, condition)
            character(len=*) :: statement
            logical          :: condition
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
    !******************    MAIN PROGRAM    ********************
    !**********************************************************
    !This program's goal is to compute Matrix multiplication of 2 matrices of given (STDIN) dimensions
    !with three different algorithms and to measure the Cpu time required by each of them.
   
    program Ex1 
    
        use intro
        use command_line_args
        use debug_and_printing
        implicit none
        omega = 1.0
        
        CALL check_commandline_args() 

        CALL condition_print("Enter min(x) and max(x) of the interval you want to discretize (real values only): ", .NOT. script)
        CALL read_value(x_min)
        CALL read_value(x_max)
        IF (x_min .GE. x_max) THEN 
            PRINT*, "max(x) has to be greater then min(x)."
            STOP
        END IF

        CALL condition_print("Enter the number of points you want to discretize the axis&
            &interval with(integer positive value only)", .NOT. script)
        CALL read_value(NN)
        IF (NN .LE. 0) THEN 
            PRINT*, "Only POSITIVE values allowed."
            STOP
        END IF

        !Debug statement
        CALL debug_mode(debug_value)

        !lattice definition
        dx = (x_max - x_min) / (NN - 1)
        CALL debug_print("discrete step: dx = ", dx, 2, debug_value)

        ALLOCATE(ax_grid(NN))
        DO ii = 1, NN
            ax_grid(ii) = x_min + (ii - 1) * dx  
        END DO    
        CALL debug_print("x-grid = ", ax_grid, 1, debug_value)

        !define tridiagonal matrix
        ALLOCATE(diagonal(NN), subdiagonal(NN-1), WORK(max(1,2*NN-2)), eigenvalues(NN))
        ALLOCATE(eigenvectors(NN, NN), norm_eigenvectors(NN,NN))
        DO ii = 1, NN
            diagonal(ii) = 2/dx**2 + (omega * ax_grid(ii))**2
            trace = diagonal(ii) + trace
            IF (ii /= NN) THEN
            subdiagonal(ii) = -1/dx**2 
            END IF
        END DO
        CALL debug_print("trace(H) = ", trace, 2, debug_value)

        CALL dsteqr("I", NN, diagonal, subdiagonal, eigenvectors, NN, WORK, INFO)
        CALL CHECKPOINT(debug_value)

        IF (INFO .EQ. 0) THEN
            CALL debug_print("INFO = 0, successful exit.", "none", 2, debug_value)
            CALL debug_print("", "none", 2, debug_value)
        ELSE IF (INFO .LT. 0) THEN
            CALL debug_print("if INFO = -i, the i-th argument had an illegal value, INFO = ", 0, INFO, debug_value)
        ELSE
            CALL debug_print("ERROR in the dsteqr algorithm: INFO =", INFO, 0, debug_value)
        END IF

        eigenvalues = diagonal
        sum_eigen = sum(eigenvalues(:))
        CALL debug_print("trace(H) after diagonalization = ", sum_eigen, 2, debug_value)

        IF (abs(trace- sum_eigen) .GT. 1d10-6) THEN
            CALL debug_print("trace(H) not preserved, difference = ", trace - sum_eigen, 0, debug_value)
            STOP
        END IF
        
        CALL debug_print("Eigenvalues = ", eigenvalues, 1, debug_value)

        
    !Write into a file the results
            INQUIRE(file="Eigenvalues.txt", exist=exist)
            IF (exist) THEN
                OPEN(1, file="EigenValues.txt", status="old", action="write")
            ELSE
                OPEN(1, file="EigenValues.txt", status="new", action="write")
            END IF
            WRITE(1, *) "# Eigenvalues:"
            WRITE(1, "(*(f0.15, :,  ', '))") eigenvalues(:)
            CLOSE(1)

            INQUIRE(file="EigenVectors.txt", exist=exist)
            IF (exist) THEN
                OPEN(2, file="EigenVectors.txt", status="old", action="write")
            ELSE
                OPEN(2, file="EigenVectors.txt", status="new", action="write")
            END IF
            WRITE(2, *) "# Axis grid:"
            WRITE(2, "(*(f0.15, :,  ', '))") ax_grid(:)
            DO ii = 1, NN
                norm_eigenvectors(:,ii) = eigenvectors(:, ii) / sqrt(sum(eigenvectors(:,ii)**2))
            WRITE(2, *) "# Normalized Eigenvector(", ii, ") =" 
            WRITE(2, "(*(f0.15, :,  ', '))")  norm_eigenvectors(:, ii)
            END DO
            CLOSE(2)

    !Deallocate matrices.
        DEALLOCATE(ax_grid, diagonal, subdiagonal, WORK, eigenvalues, eigenvectors, norm_eigenvectors)

    end program Ex1


   


