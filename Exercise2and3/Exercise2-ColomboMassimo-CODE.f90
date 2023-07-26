    !In this Program lapack routines have been used so to compile:
    !gfortran Exercise2-ColomboMassimo-CODE.f90 -o Exercise2 -llapack
    
    
    !**********************************************************
    !**********************    INTRO    ***********************
    !**********************************************************
    !In this module variables get defined. The use of each variable is made explicit.
    
    module intro 

        !To prevent implicit variable definition. 
        implicit none

        !Logical variables:
        !- DEBUG variable - allows to run the program in "debugging mode":
        !CHOOSE ONE:
        !logical :: debug = .true. !<--> Debug mode ON;
        logical :: debug = .false. !<--> Debug mode Off.
    
        !- SCRIPT variable - is used together with the python script:
        !CHOOSE ONE:
        logical :: script = .true. !<--> Removes all the STDOUT-sentences except outputs;
        !logical :: script = .false. !<--> Keeps all the STDOUTs that explain the results.

        !Dimensions of the matrices: 
        !   "n_columns_**" = Number of columns of matrix *.
        integer*4 nn, ii, jj !nn = matrix dimension.

        logical :: exist
        logical :: diagonal = .true. !if .true. the program utilize as eigen values n random real values in the range [0,1]
        !logical :: diagonal = .false.

        !Matrices type definition.
        complex*16, allocatable, dimension(:,:) :: AA

        !Variables to get eigenvalues
        complex*16, allocatable, dimension(:) ::  WORK, RWORK, pre_WORK
        real*8, allocatable, dimension(:) :: eigenvalues, spacing, norm_spacing
        integer :: INFO, LWORK

        !Trace function variables.
        real*8 :: sum_eigen = 0, avg_eigen, avg_spacing, TT
        complex*16 :: trace_AA
        complex*16, external :: trace_C_mat
    end module intro


    !**********************************************************
    !******************    MAIN PROGRAM    ********************
    !**********************************************************
    !This program's goal is to compute Matrix multiplication of 2 matrices of given (STDIN) dimensions
    !with three different algorithms and to measure the Cpu time required by each of them.
   
    program Ex1 
    
        use intro
        implicit none

    !Enter the matrices dimension (by keyboard).
        IF (script .eqv. .false.) then
            PRINT *, 'Enter the dimension of the square matrix: '
        END IF
        CALL read_and_check_mat_dim(nn)  !Subroutine's explanation is presented in the subroutine section.
 
    !Debug statement
        CALL debug_mode(debug)   

    !Initializing random hermitian matrix.
        ALLOCATE(AA(nn, nn)) !ALLOCATE() Assign the input dimensions to the matrix. 
        CALL random_hermitian(AA, nn, debug)
        IF (debug .eqv. .true.) THEN
            PRINT*, "-------- RANDOM MATRIX OUTPUT --------", new_line('a')
        CALL print_C_matrix_debug(AA, nn, nn, "Hermitian matrix AA", debug)
        END IF   
    
    !Sequence that finds eigenvalues, ("Diagonalizing").
        ALLOCATE(eigenvalues(nn), pre_WORK(1), RWORK(max(1, 3*nn-2)))

        LWORK = -1  !With LWORK = -1, it optimizes the dimension of WORK array. returning pre_WORK(1) as that value.
        CALL zheev("N","U", nn, AA, nn, eigenvalues, pre_WORK, LWORK, RWORK, INFO)
        
        LWORK = int(pre_WORK(1))    !After allocating WORK with optimal dimension it can compute the eigenvalues:
        ALLOCATE(WORK(MAX(1,LWORK)))
        CALL zheev("N","U", nn, AA, nn, eigenvalues, WORK, LWORK, RWORK, INFO)

    !Check that the trace (sum of all eigenvalues) has been preserved by diagonalization.
        sum_eigen = sum(eigenvalues(:))
        trace_AA = trace_C_mat(AA, nn)
        IF (trace_AA%re - sum_eigen > 1d-10 .or. trace_AA%im /= 0) THEN         !Max error set to 1d-10
            PRINT*, "ERROR: Trace has not been preserved by diagonalization."
            STOP
        END IF
    
    !Prints different infromations (debug mode))
        IF (debug .eqv. .true.) THEN
            WRITE(*,*) "-------- EIGENVALUES OUTPUT --------", new_line('a')
            WRITE(*, '(a, I0)', advance = 'no') "INFO: ", INFO
            WRITE(*, '(a, I0)', advance = 'no') "       Best dimension= ", size(WORK), new_line('a')
            WRITE(*,'(a, *(f10.5,", "), i0)') "Eigenvalues are:", eigenvalues(:)
            WRITE(*,'(a, i0)', advance = 'no') new_line('a')
            WRITE(*,'(a, i0)', advance = 'no') "Sum of Eigenvalues            = "
            PRINT'(f20.10)', sum_eigen
            CALL print_named_complex_number(trace_AA, "Trace of the Hermitian matrix")
            PRINT*, new_line('a')
        END IF
   

    !Normalized Spacing section.
        IF (diagonal .eqv. .true.) THEN
            CALL Random_number(eigenvalues)
            DO ii=1, nn
                DO jj= 1, nn-1
                    IF (eigenvalues(ii) .LT. eigenvalues(jj)) THEN
                    TT = eigenvalues(ii) 
                    eigenvalues(ii) = eigenvalues(jj)
                    eigenvalues(jj) = TT 
                    END IF
                END DO
            END DO
        END IF
    
        ALLOCATE(spacing(nn-1), norm_spacing(nn-1))
        DO ii = 1, nn - 1
            spacing(ii) = (eigenvalues(ii + 1) - eigenvalues(ii))
        END DO
        avg_spacing = sum(spacing(:)) / size(spacing)

        DO ii = 1, nn - 1
            norm_spacing(ii) = spacing(ii)/avg_spacing
        END DO
        

        IF (debug .eqv. .true.) THEN 
        WRITE(*,'(a, *(f10.5,", "), i0)') "Normalized spacing values = ", norm_spacing(:)
        PRINT*, new_line('a'), "Average spacing = ", avg_spacing
        END IF
        
    !Write into a file the results
        IF (diagonal .eqv. .true.) THEN
            INQUIRE(file="Normalized_Spacing_random.txt", exist=exist)
            IF (exist) THEN
                OPEN(1, file="Normalized_Spacing_random.txt", status="old", position="append", action="write")
            ELSE
                OPEN(1, file="Normalized_Spacing_random.txt", status="new", action="write")
            END IF
            WRITE(1, '(f30.15)') norm_spacing(:) 
            CLOSE(1)
        ELSE     
            INQUIRE(file="Normalized_Spacing.txt", exist=exist)
            IF (exist) THEN
                OPEN(2, file="Normalized_Spacing.txt", status="old", position="append", action="write")
            ELSE
                OPEN(2, file="Normalized_Spacing.txt", status="new", action="write")
            END IF
            WRITE(2, '(f30.15)') norm_spacing(:) 
            CLOSE(2)
        END IF

    !Deallocate matrices.
        DEALLOCATE(AA, pre_WORK, WORK, RWORK, norm_spacing, spacing)

    end program Ex1


    !**********************************************************
    !*******************    SUBROUTINES    ********************
    !**********************************************************    
    !subroutines used in the main program. 

    !This subroutine's goal is to read as STDIN the matrices dimension. 
    !It accepts only positive valued integers and informs about possible mistakes.
    !It keeps on receiving input until the correct one is inserted or the maximum number of wrong ones is reached.
    subroutine read_and_check_mat_dim(nn)
        integer*4 ierror, nn, ii
        ii = 0
        DO
            READ(*,*, iostat = ierror) nn !Read inputs should be integers, "ierror value" takes sign of this condition.
        
            IF (ierror == 0) THEN    
            !IOSTAT = 0 means the READ was executed flawlessly and all variables have received their input values.
                IF (nn < 0) THEN                     !No negative input allowed.
                    PRINT*, "Matrix dimensions accept only POSITIVE integers" 
                ELSE IF (int(nn) > 10**6) THEN         !If the matrix dimensions are too high the algorithms could take ages to compute the product.
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

    !This subroutine initializes a square complex nn-matrix as a Hermitian matrix
    !AGGIUNGERE DEBUG OPTION
    subroutine random_hermitian(C_mat, dim, debug)
        integer*4 dim, ii, jj
        logical debug
        complex*16, dimension(dim, dim) :: C_mat
        complex*16, allocatable, dimension(:) :: rand_vec
        INTEGER :: ISEED(4)
        real*8, dimension(3) :: uu
        !integer, dimension(3) :: seed


        ALLOCATE (rand_vec(dim*(dim + 1)/2))
        
        CALL random_number(uu)
        !seed(:) = FLOOR((4096)*uu(:)) 
        ISEED = (/FLOOR((4096)*uu(:)),1/) !Initializes random ISEED
        CALL zlarnv(2, iseed, dim*(dim + 1)/2, rand_vec)
        IF (debug .eqv. .true.) THEN
            WRITe(*,'(a)', advance = "yes") "Random vector for random hermitian matrix:"
            CALL print_complex_vector(rand_vec, dim*(dim + 1)/2)
        END IF
    
        DO ii = 1, dim
            DO jj = ii, dim
                C_mat(ii, jj) = rand_vec(dim * (ii - 1) - (ii - 2) * (ii - 1)/2  + jj - ii + 1)
                IF (ii /= jj) THEN
                    C_mat(jj, ii) = conjg(C_mat(ii,jj))
                END IF
                IF (ii == jj) THEN
                C_mat(ii, jj)%im = 0 
                END IF
            END DO
        END DO

        DEALLOCATE(rand_vec)


    end subroutine

    !This subroutine prints a complex vector in a readable form:
    subroutine print_complex_vector(vec, len)
        complex*16 vec(len)
        character(len=*), parameter :: fmt = "(F20.10,a3,F20.10,a)"
        DO ii = 1, size(vec)
            PRINT fmt, realpart(vec(ii)), " + ", imagpart(vec(ii)), "i" 
        ENDDO
        PRINT*, new_line("a")
    end subroutine

    !This subroutine prints out a Complex Matrix with named "name" iff debug mode is on.
    subroutine print_C_matrix_debug(AA, nn, mm, name, debug)
        integer*4 kk, jj, nn, mm
        complex*16, dimension(nn, mm) :: AA 
        character(*) name
        logical debug
        character(len=*), parameter :: fmt = "(25(f10.5,a3,f10.5,a4))"
        !put the "right" sign in front of the imaginary part. 
        character(len=4), DIMENSION(size(AA, 1), size(AA, 2)) :: sign
            WHERE(AIMAG(AA)>=0.) sign = '  +' 
            WHERE(AIMAG(AA)<0.) sign = '  -'

        IF (debug .eqv. .true.) then
            WRITE(*,'(a, i0)', advance = 'no') name
            PRINT*,' =' 
                DO kk = 1, size(AA, 1)
                    WRITE(*, fmt) (REALPART(AA(kk,jj)), sign(kk,jj), &
                                   ABS(IMAGPART(AA(kk,jj))), "i   ", jj = 1, size(AA, 2))    
                END DO
            PRINT*, new_line('a')
        END IF
    end subroutine

    !This function calculate the trace of a complex matrix
    function trace_C_mat(C_mat, nn) result(sum)
        implicit none
        complex*16 sum
        integer*4 ii, nn
        complex*16, dimension(nn, nn) :: C_mat
        
        sum = 0
        DO ii = 1 , nn
            sum = C_mat(ii, ii) + sum
        END DO
    end function
    
    !This subroutine prints in a readable form a "named" complex number.
    subroutine print_named_complex_number(zz, name)
        complex*16 zz
        real*8 re_part, im_part
        character(*) name
        character(len=*), parameter :: fmt = "(F20.10,a3,F20.10,a)"
        re_part = realpart(zz)
        im_part = imagpart(zz)
        IF  (im_part >= 0) THEN
            WRITE(*,'(a, i0)', advance = 'no') name
            WRITE(*, '(a)', advance='no') " = "
            WRITE(*, fmt) re_part, " +", im_part, "i "
        ELSE 
            WRITE(*,'(a, i0)', advance = 'no') name
            WRITE(*, '(a)', advance='no') " = "
            WRITE(*, fmt) re_part, " -", abs(im_part), "i "
        END IF
    end subroutine

    !This subroutine checks if two matrices are the same up to a given error (eps).
    !it makes the difference between corresponding entries of same dimensions matrices, 
    !then if the ij-difference has an absolute value higher than eps, the matrices are considered different.
    subroutine comparing_matrix(AA, BB, nn, mm, eps)
        integer*4 nn, mm, ii, jj, err_counts
        real*4 AA(nn, mm), BB(nn, mm), eps
        err_counts = 0

        DO ii = 1, nn
            DO jj = 1, mm
                IF(ABS(AA(ii,jj) - BB(ii,jj)) > eps) THEN
                    PRINT*, 'Entry (',ii, ',', jj,')', 'of the matrices are different' !It returns which matrix entry is different.
                    err_counts = err_counts + 1        !it counts how many entries are different.
                ENDIF
            ENDDO
        ENDDO

        IF(err_counts /= 0) THEN
            PRINT*, new_line('a'), 'There are', err_counts, 'entries in the two matrices that differ more than a value of', eps
        ELSE
            PRINT*, '    Matrices result the same up to a value', eps
        ENDIF

    end subroutine 

    !Prints debug title.
    subroutine debug_mode(debug)
        logical :: debug
        IF (debug .eqv. .true.) then
            PRINT*, new_line('a')
            PRINT*, "-------------------------------------", new_line('a')
            PRINT*, "--------- DEBUGGING MODE ON ---------", new_line('a')
            PRINT*, "-------------------------------------", new_line('a')
        END IF
    end subroutine

    !Prints checkpoint statement.
    subroutine checkpoint(debug)
        logical :: debug
        IF (debug .eqv. .true.) then
            PRINT*, "--------- CHECKPOINT ---------", new_line('a')
        END IF
    end subroutine
