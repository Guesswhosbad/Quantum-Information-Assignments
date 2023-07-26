    module intro 
        implicit none

        !for-loops variables.
        integer*4 ii, jj, kk

        !dimensions of the matrices. 
        integer*4 n_rows_AA, n_rows_BB, n_columns_AA, n_columns_BB !numbers of rows and columns of the 2 matrices
        
        !cpu_time variables.
        real*4 start, finish
        real*4 max_error_allowed

        !matrices
        real*4, allocatable, dimension(:,:) :: AA, BB, CC_1, CC_2, CC_matmul

    end module intro

    program Ex3 !matrix multiplication
    
        use intro
        implicit none
        
    !enter the matrices dimension (by keyboard).
        !PRINT *, 'Enter the number of rows and columns of the first matrix: (#rows, #columns)'
        READ(*,*) n_rows_AA, n_columns_AA
        !PRINT *, 'Enter the number of rows and columns of the second matrix: (#rows, #columns)'
        READ(*,*) n_rows_BB, n_columns_BB
        ALLOCATE(AA(n_rows_AA, n_columns_AA), BB(n_rows_BB, n_columns_BB))

    !check if the matrices are compatibles for the product.
        IF (n_columns_AA /= n_rows_BB) THEN
            PRINT*, 'In order to do the matrix product the #columns of A has to be equal to #rows of B'
            return
        END IF
       
        
    !CHOICE OF FILLING THE MATRIX ENTRIES:

    !1) Keyboard input:

        !PRINT *, 'Enter the elements of the first matrix'            
        !READ(*,*) ((AA(ii,jj),jj=1,n_columns_AA),ii=1,n_rows_AA) 
      
        !PRINT *, 'Enter the elements of the second matrix'
        !READ(*,*) ((BB(ii,jj),jj=1,n_columns_BB),ii=1,n_rows_BB)
        
    !2) Random filling:

        CALL random_number(AA) !type output should depend on the variable's type's input
        CALL random_number(BB)
        
    !Print both matrices.

        !PRINT *, new_line('a'), 'Matrix A ='
        !CALL print_matrix(AA,n_rows_AA,n_columns_AA)
        !PRINT*, new_line('a')
        !PRINT *, 'Matrix B ='
        !CALL print_matrix(BB,n_rows_BB,n_columns_BB)
        !PRINT*, new_line('a')

    !MATRIX MULTIPLICATION IN THREE WAYS (with cpu_time measurement):
    
    !1) 3 for-loops, Usual Matrix Product: 
        
        ALLOCATE(CC_1(n_rows_AA, n_columns_BB))
        CC_1 = 0.0

        !PRINT *, 'the matrix product through 3-for-loops (Usual way) is A*B ='
            
        CALL cpu_time(start)

            DO ii = 1, n_rows_AA
                DO jj = 1, n_columns_BB
                    DO kk = 1, n_columns_AA 
                        CC_1(ii,jj) = CC_1(ii,jj) + AA(ii, kk) * BB(kk, jj)
                    ENDDO
                ENDDO
            ENDDO

        CALL cpu_time(finish)

        !CALL print_matrix(CC_1, size(CC_1, 1), size(CC_1, 2))
        !PRINT*,'the cpu_time measured is:', finish - start, 's'
        PRINT*, finish - start
        !PRINT*, new_line('a')

            
    !2) 3 for-loops with inverted indices:
        
        ALLOCATE(CC_2(n_rows_AA, n_columns_BB))
        CC_2 = 0.0
            
        !PRINT *, 'the matrix product through 3-for-loops with inverted indices is A*B ='
            
        CALL cpu_time(start)

        !3 for-loops algorithm with inverted indeces
            DO kk = 1, n_columns_AA
                DO ii = 1, n_rows_AA
                    DO jj = 1, n_columns_BB
                        CC_2(ii,jj) = CC_2(ii,jj) + AA(ii, kk) * BB(kk, jj)
                    ENDDO
                ENDDO
            ENDDO

        CALL cpu_time(finish)

        !CALL print_matrix(CC_2, size(CC_2, 1), size(CC_2, 2))
        !PRINT*,'the cpu_time measured is:',
        PRINT*, finish - start
        !PRINT*, new_line('a')

    !3) MATMUL product:

        ALLOCATE(CC_matmul(n_rows_AA, n_columns_BB))
        CC_matmul = 0.0

        !PRINT *, 'The Blas solution for A*B is:'

        CALL cpu_time(start)
        CC_matmul = matmul(AA,BB)
        CALL cpu_time(finish)

        !CALL print_matrix(CC_matmul, size(CC_matmul, 1), size(CC_matmul, 2))
        !PRINT*,'the cpu_time measured is:', finish - start, 's'
        PRINT*, finish - start
        !PRINT*, new_line('a')

    !Check if the result is the same. 
       ! max_error_allowed = 0.001 !max absolute value of the difference between same 2 entries the we allow.
       ! CALL comparing_matrix(CC_1, CC_2, size(CC_1,1), size(CC_1,2), max_error_allowed)
       ! CALL comparing_matrix(CC_1, CC_matmul, size(CC_1,1), size(CC_1,2), max_error_allowed)
       ! CALL comparing_matrix(CC_2, CC_matmul, size(CC_1,1), size(CC_1,2), max_error_allowed)
       
    !Deallocate matrices
        DEALLOCATE(CC_1, CC_2)
        DEALLOCATE(AA, BB, CC_matmul)

    end program Ex3

    subroutine print_matrix(AA,nn,mm)
        integer*4 nn, mm, kk
        real*4 AA(nn,mm) !nn = #rows, mm = #columns
        DO kk=1,nn
        PRINT '(20f6.2)', AA(kk,1:mm)
        ENDDO
    endsubroutine   

   subroutine comparing_matrix(AA, BB, nn, mm, eps)
        integer*4 nn, mm, ii, jj, err_counts
        real*4 AA(nn, mm), BB(nn, mm), eps
        err_counts = 0

        DO ii = 1, nn
            DO jj = 1, mm
                IF(ABS(AA(ii,jj) - BB(ii,jj)) > eps) THEN
                    PRINT*, 'Entry (',ii, ',', jj,')', 'of the matrices are different'
                    err_counts = err_counts + 1
                ENDIF
            ENDDO
        ENDDO

        IF(err_counts /= 0) THEN
            PRINT*, new_line('a'), 'There are', err_counts, 'entries in the two matrices that differ more than a value of', eps
        ELSE
            PRINT*, 'Matrices result the same up to a value', eps
        ENDIF

    endsubroutine 