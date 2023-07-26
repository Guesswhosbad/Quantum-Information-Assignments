  !When executing try both square and non-square matrices.

    !**********************************************************
    !*********************    MODULES    **********************
    !**********************************************************
module Complex_Matrices
    
    implicit none

    !---------------------- DERIVED TYPES ----------------------!
    type C_matrix
        complex*16, allocatable, dimension(:,:) :: elements     !complex*16 = doublecomplex
        integer*4, dimension(2) :: dims = (/0,0/)               !I dont know why, but in order to work this has to be initialized.
    end type C_matrix

    type, extends(C_matrix) :: square_C_matrix
        complex*16 :: trace = (0d0, 0d0)                        !Same reason here.
    end type square_C_matrix

    character(len=*), parameter :: fmt = "(25(f6.2,a4,f5.2,a3))"!Format variable that will be used to print complex matrices.

    !----------------------- INTERFACES ------------------------!
    interface operator(.adj.)  !Adjoint.
        module procedure mat_adjoint
    end interface

    interface operator(.trace.) !Trace.
        module procedure mat_trace
    end interface
    
    contains 

    !----------------------- SUBROUTINES -----------------------!
    !This subroutine initializes a complex matrix (square complex matrix if n_rows = n_columns)
    subroutine init_C_matrix(mat, n_rows, n_columns)
        class(C_matrix), allocatable :: mat                     !Defining mat as a class, it generalizes the initialization also to extended type.
                                                                !(To be honest there is a more general and proper way, but i didnt already have understood it correctly) 
        integer*4 n_rows, n_columns         
        real re_part(n_rows,n_columns)                          !In order to inizialize a random Complex matrix, it first initialize 2 random real matrices.
        real im_part(n_rows,n_columns)                          !One will be the real part, the other the imaginary one.

        !It associates the TYPE given the "correct" condition.
        IF (n_rows == n_columns) THEN                           
            mat = square_C_matrix()
        ELSE
            mat = C_matrix()
        END IF

        !From here on, it initializes matrix's dimension, elements and trace.
        mat%dims(1) = n_rows                
        mat%dims(2) = n_columns

        ALLOCATE(mat%elements(n_rows, n_columns))
        CALL random_number(re_part)
        CALL random_number(im_part)
        mat%elements = cmplx(re_part, im_part)

        SELECT TYPE(mat)                                        !Only for square matrices it initializes the trace as its correct value.
            TYPE IS(square_C_matrix)
                mat%trace = .trace.mat
        END SELECT
     
    end subroutine

        
    !This subroutine prints out a Complex Matrix with named "name" iff debug mode is on.
    subroutine print_C_matrix_debug(AA, name, debug)
        type(C_matrix) AA 
        character(*) name
        logical debug
        integer*4 kk,jj
        !put the "right" sign in front of the imaginary part. 
        character(len=4), DIMENSION(AA%dims(1),AA%dims(2)) :: sign
            WHERE(AIMAG(AA%elements)>=0.) sign = '  + ' 
            WHERE(AIMAG(AA%elements)<0.) sign = '  - '

        IF (debug .eqv. .true.) then
            PRINT*, name, ' =' 
                DO kk = 1, AA%dims(1)
                    WRITE(*, fmt) (REALPART(AA%elements(kk,jj)), sign(kk,jj), &
                                   ABS(IMAGPART(AA%elements(kk,jj))), "i  ", jj = 1, AA%dims(2))    
                END DO
            PRINT*, new_line('a')
        END IF
    end subroutine


    !This subroutine prints a complex matrix "AA" on a file associated to the unit value. 
    subroutine print_C_matrix_intofile(AA, name, unit)
        type(C_matrix) AA 
        character(*) :: name
        integer*4 unit
        integer*4 kk,jj
        character(len=4), DIMENSION(AA%dims(1),AA%dims(2)) :: sign
            WHERE(AIMAG(AA%elements)>=0.) sign = '  + '
            WHERE(AIMAG(AA%elements)<0.) sign = '  - '
    
            WRITE(unit, *) name, ' dimensions:   Rows = ', AA%dims(1), 'Columns = ', AA%dims(2), new_line('a')

            WRITE(unit, *) name, ' elements =' 
            DO kk = 1, AA%dims(1)
                WRITE(unit, fmt) (REALPART(AA%elements(kk,jj)), sign(kk,jj), &
                                ABS(IMAGPART(AA%elements(kk,jj))), "i  ", jj = 1, AA%dims(2))    
            END DO

    end subroutine


    subroutine print_named_complex_number(re_part, im_part, name)
        real*8 re_part, im_part
        character(*) name
        IF  (im_part >= 0) THEN
            WRITE(*, '(a)', advance='no') name
            WRITE(*, fmt) re_part, " +", im_part, "i "
        ELSE 
            WRITE(*, '(a)', advance='no') name
            WRITE(*, fmt) re_part, " -", abs(im_part), "i "
        END IF
    end subroutine


    !------------------------ FUNCTIONS ------------------------!
    !Evaluate trace of a square matrix.
    function mat_trace(s_mat) result(trace)
        class(square_C_matrix), intent(IN) :: s_mat
        complex*16 trace
        integer*4 ii
        trace = 0  
        DO ii = 1, s_mat%dims(1)                !The percentage symbol % is used to access the members of a derived type.
            trace = s_mat%elements(ii,ii) + trace
        ENDDO
    end function

    !Returns the adjoint matrix. (If input matrix is square, it also returns adjoint trace).
    function mat_adjoint(mat) result(adj_mat)
        class(C_matrix), intent(IN) :: mat
        class(C_matrix), allocatable :: adj_mat

        SELECT TYPE (mat)
        TYPE IS (C_matrix)
        adj_mat = C_matrix()
        TYPE IS (square_C_matrix)
        adj_mat = square_C_matrix()
        END SELECT

        adj_mat%dims(1) = mat%dims(2)
        adj_mat%dims(2) = mat%dims(1)
        adj_mat%elements = (conjg(transpose(mat%elements)))

        SELECT TYPE (mat)
        TYPE IS (square_C_matrix)
        SELECT TYPE (adj_mat)
        TYPE IS (square_C_matrix)
        adj_mat%trace = conjg(mat%trace)
        END SELECT
        END SELECT
    end function

end module Complex_Matrices

module intro
    use Complex_Matrices
    integer*4 n_rows, n_columns
    class(C_matrix), allocatable :: AA, AA_adj

end module intro

    !**********************************************************
    !******************    MAIN PROGRAM    ********************
    !**********************************************************
  
  program main
    
    use intro
    implicit none

    !Enter the dimension of the matrix.
    PRINT *, 'Enter the number of rows and columns of the matrix: (#rows, #columns)'
    CALL read_and_check_mat_dim(n_rows, n_columns) 

    !Randomly initialize the matrix, also initializing the trace in case of square matrix.
    CALL init_C_matrix(AA, n_rows, n_columns)
    
    !Print matrix 
    CALL print_C_matrix_debug(AA, "Random matrix", .true.)

    AA_adj = .adj.AA 

    CALL print_C_matrix_debug(AA_adj, "Adjoint random matrix", .true.)

    SELECT TYPE(AA)
            TYPE IS(square_C_matrix)
                PRINT*, "Traces: "
                CALL print_named_complex_number(REALPART(AA%trace), IMAGPART(AA%trace), "Random matrix trace = ")
    END SELECT

    SELECT TYPE(AA_adj)
            TYPE IS(square_C_matrix)
                CALL print_named_complex_number(REALPART(AA_adj%trace), IMAGPART(AA_adj%trace), "Adjoint random matrix trace = ")
                PRINT*, new_line('a')
    END SELECT

    !Writes input matrix on a file, specifing dimensions and elements.
    OPEN(1, file = 'complex_matrix.txt')   
    CALL print_C_matrix_intofile(AA, "Matrix", 1)
     
end program main

    !**********************************************************
    !*******************    SUBROUTINES    ********************
    !**********************************************************  

    !This subroutine's goal is to read as STDIN the matrices dimension. 
    !It accepts only positive valued integers and informs about possible mistakes.
    !It keeps on receiving input until the correct one is inserted or the maximum number of wrong ones is reached.
    subroutine read_and_check_mat_dim(n_rows, n_columns)
        integer*4 ierror, n_rows, n_columns, ii
        ii = 0
        DO
            READ(*,*, iostat = ierror) n_rows, n_columns !Read inputs should be integers, "ierror value" takes sign of this condition.
        
            IF (ierror == 0) THEN    
            !IOSTAT = 0 means the READ was executed flawlessly and all variables have received their input values.
                IF (n_rows < 0 .or. n_columns < 0) THEN                     !No negative input allowed.
                    PRINT*, "Matrix dimensions accept only POSITIVE integers" 
                ELSE IF (int(n_rows) * int(n_columns) > 10**7) THEN         !If the matrix dimensions are too high the algorithms could take ages to compute the product.
                    PRINT*, "Are you sure the computer can handle this?"
                ELSE
                    RETURN                                                  !If the dimension are plausible it accepts the input.
                END IF
            ELSE IF (ierror /= 0) THEN
            !IOSTAT > 0, the previous READ has encountered some problem. A common problem is illegal data: supplying a real number to an integer variable.
            !IOSTAT < 0, it means the end of the input has reached. Under this circumstance, some or all of the variables in the READ may not receive input values.
                PRINT*, "Matrix dimensions accept only integers"
            END IF
            ii = ii + 1            !ii it's used to exit the do-loop after a certain amount of wrong entries.
            IF (ii > 10) THEN
                PRINT*, 'Too many wrong entries. exiting the program.'
                STOP
            END IF
        END DO !crtl + c = CANCEL to force the exit. 
    end subroutine






