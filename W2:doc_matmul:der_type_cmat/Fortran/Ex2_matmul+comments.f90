    !**********************************************************
    !**********************    INTRO    ***********************
    !**********************************************************
    !In this module variables get defined. The use of each variable is made explicit.
    
    module intro 

        !To prevent implicit variable definition. 
        implicit none

        !For-loops variables.
        integer*4 ii, jj, kk

        !Logical variables:
        !- DEBUG variable - allows to run the program in "debugging mode":
        !CHOOSE ONE:
        logical :: debug = .true. !<--> Debug mode ON;
        !logical :: debug = .false. !<--> Debug mode Off.
    
        !- SCRIPT variable - is used together with the python script:
        !CHOOSE ONE:
        !logical :: script = .true. !<--> Removes all the STDOUT-sentences except the cpu_times;
        logical :: script = .false. !<--> Keeps all the STDOUTs that explain the results.

        !Dimensions of the matrices: 
        !   "n_rows_**" = Number of rows of matrix *;
        !   "n_columns_**" = Number of columns of matrix *.
        integer*4 n_rows_AA, n_rows_BB, n_columns_AA, n_columns_BB 
        
        !Cpu_time variables.
        real*4 start, finish

        !Max error used in the subroutine: "comparing_matrix()".
        real*4 max_error_allowed

        !Matrices type definition.
        real*4, allocatable, dimension(:,:) :: AA, BB, CC_1, CC_2, CC_matmul

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
            PRINT *, 'Enter the number of rows and columns of the first matrix: (#rows, #columns)'
        END IF
        CALL read_and_check_mat_dim(n_rows_AA, n_columns_AA)  !Subroutine's explanation is presented in the subroutine section.
 

        IF (script .eqv. .false.) then
        PRINT *, 'Enter the number of rows and columns of the second matrix: (#rows, #columns)'
        END IF
        CALL read_and_check_mat_dim(n_rows_BB, n_columns_BB)

    !Check if the matrices are compatibles for the product (1st matrix's #columns has to be equal to 2nd matrix's rows). Otherwise, stop.
        IF (n_columns_AA /= n_rows_BB) THEN
            IF (script .eqv. .false.) then
            PRINT*, 'In order to do the matrix product the #columns of A has to be equal to #rows of B'
            END IF
            STOP
        END IF

        ALLOCATE(AA(n_rows_AA, n_columns_AA), BB(n_rows_BB, n_columns_BB)) !ALLOCATE() Assign the input dimensions to the matrices.

    !Debugging mode statement.
        IF (debug .eqv. .true.) then
            PRINT*, new_line('a')
            PRINT*, "-------------------------------------", new_line('a')
            PRINT*, "--------- DEBUGGING MODE ON ---------", new_line('a')
            PRINT*, "-------------------------------------", new_line('a')
        END IF

    
    !CHOICE OF FILLING THE MATRIX ENTRIES:
    !ADD CHOICE ON KEYBOARDING INPUT OR RANDOM FILLINGS (OR JUST REMOVE) 
    
    !1) Keyboard input:

        !PRINT *, 'Enter the elements of the first matrix'            
        !READ(*,*) ((AA(ii,jj),jj=1,n_columns_AA),ii=1,n_rows_AA) 
      
        !PRINT *, 'Enter the elements of the second matrix'
        !READ(*,*) ((BB(ii,jj),jj=1,n_columns_BB),ii=1,n_rows_BB)
        
    !2) Random filling:

        CALL random_number(AA) !Randomly fills matrix's entries with value from 0 to 1.
        CALL random_number(BB)
        
    !If Debug Mode ON: Print matrices used in the product.

        CALL print_matrix_debug(AA,n_rows_AA,n_columns_AA, 'Matrix A', debug) !Subroutine's explanation is presented in the subroutine section.
        CALL print_matrix_debug(BB,n_rows_BB,n_columns_BB, 'Matrix B', debug)

    
    !MATRIX MULTIPLICATION IN THREE WAYS (with cpu_time measurement):
    !1) 3 for-loops, Usual Matrix Product: 
        
        ALLOCATE(CC_1(n_rows_AA, n_columns_BB)) !CC will be the name of the matrix obtained by A*B product.
        CC_1 = 0.0                              !CC_* --> * will refer to the algorithm used to do the product.

        CALL cpu_time(start) !Cpu time() marks start and finish time of the 3for-loop.
        !The next 3for-loop's order of computing corresponds to the handmade matrix product.
        DO ii = 1, n_rows_AA 
            DO jj = 1, n_columns_BB
                DO kk = 1, n_columns_AA 
                    CC_1(ii,jj) = CC_1(ii,jj) + AA(ii, kk) * BB(kk, jj)
                ENDDO
            ENDDO
        ENDDO
        CALL cpu_time(finish)

        IF (script .eqv. .false.) then
            PRINT*, 'CPU TIME RESULTS:'
            PRINT*,'The "Usual Matrix Product" cpu_time measured is:', finish - start, 's'
        ELSE
            PRINT*, finish - start !difference of the marked times by cpu_time = time of loop's execution.
        END IF
            
    !2) 3 for-loops with inverted indices:
        
        ALLOCATE(CC_2(n_rows_AA, n_columns_BB))
        CC_2 = 0.0
         
        !3 for-loops algorithm with inverted indeces
        CALL cpu_time(start)
            !Here the order of the 3for-loop is inverted so that the calculation is the same but it's allocating memory in a different way.
        DO kk = 1, n_columns_AA
            DO ii = 1, n_rows_AA
                DO jj = 1, n_columns_BB
                    CC_2(ii,jj) = CC_2(ii,jj) + AA(ii, kk) * BB(kk, jj)
                ENDDO
            ENDDO
        ENDDO
        CALL cpu_time(finish)

        IF (script .eqv. .false.) then
            PRINT*,'The "Inverted Indices Mat. Prod." cpu_time measured is:', finish - start, 's'
        ELSE
            PRINT*, finish - start
        END IF

    !3) MATMUL product:

        ALLOCATE(CC_matmul(n_rows_AA, n_columns_BB))
        CC_matmul = 0.0

        CALL cpu_time(start)
        !Matmul is the matrix multiplication's algorithm provided by BLAS.
        CC_matmul = matmul(AA,BB)
        CALL cpu_time(finish)

        IF (script .eqv. .false.) then
            PRINT*,'The "Blas Matmul" cpu_time measured is:', finish - start, 's'
            PRINT*, new_line('a')
        ELSE
            PRINT*, finish - start
        END IF
        
    !Check if the results are the same.
        IF (debug .eqv. .true.) then 
            max_error_allowed = 0.001 !max absolute value of the difference between same 2 entries that we allow.
            !CC_1(1,2) = CC_1(1,2) + 1 !Modifying by hand the mat. prod. to check the comparing_matrix subroutine.
            
            PRINT*, 'MATRIX COMPARISON:'
            PRINT*, '* Comparing - Usual mat. prod. - and - Inverted Indices mat. prod. - matrices:'
            CALL comparing_matrix(CC_1, CC_2, size(CC_1,1), size(CC_1,2), max_error_allowed)
            PRINT*, '* Comparing - Usual mat. prod. - and - Blas Matmul - matrices:'
            CALL comparing_matrix(CC_1, CC_matmul, size(CC_1,1), size(CC_1,2), max_error_allowed)
            PRINT*, '* Comparing - Inverted Indices mat. prod. - and - Blas Matmul - matrices:'
            CALL comparing_matrix(CC_2, CC_matmul, size(CC_1,1), size(CC_1,2), max_error_allowed)
            PRINT*, new_line('a')
        END IF
    
    !If Debug Mode ON: Print the different algorithms matrix product.
        CALL print_matrix_debug(CC_1, size(CC_1, 1), size(CC_1, 2), 'The "Usual way" matrix product is A*B', debug)
        CALL print_matrix_debug(CC_2, size(CC_2, 1), size(CC_2, 2), 'The "Inverted indices" matrix product is A*B', debug)
        CALL print_matrix_debug(CC_matmul, size(CC_matmul, 1), size(CC_matmul, 2), 'The Blas solution for A*B', debug)
       
    !Deallocate matrices.
        DEALLOCATE(CC_1, CC_2)
        DEALLOCATE(AA, BB, CC_matmul)

    end program Ex1


    !**********************************************************
    !*******************    SUBROUTINES    ********************
    !**********************************************************    
    !subroutines used in the main program. 

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

    !This subroutine prints out a Matrix AA with nn-rows, mm-columns, named "name" iff debug mode is on.
    subroutine print_matrix_debug(AA, nn, mm, name, debug)
        integer*4 nn, mm, kk
        real*4 AA(nn,mm) !nn = #rows, mm = #columns
        character(*) name
        logical debug
        
        IF (debug .eqv. .true.) then
            PRINT*, name, ' =' 
            DO kk=1,nn
                PRINT '(20f6.2)', AA(kk,1:mm)
            ENDDO
            PRINT*, new_line('a')
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