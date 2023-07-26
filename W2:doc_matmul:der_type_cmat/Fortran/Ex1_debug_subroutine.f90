program Ex1
    logical :: debug = .false. !Changing this to .true., it will enter debug mode.
    integer*4 :: n_rows= 4, n_columns=5 !Random dimensions user can choose.
    real*4 :: mat(4, 5) = 0
    character(*), parameter :: name  = "Exercise's matrix"

    CALL print_matrix_debug(mat, n_rows, n_columns, name, debug)
    IF (debug .eqv. .false.) THEN
        PRINT*, "Debug mode is OFF"
    END IF

end program Ex1

!"Checkpoint" subroutine.
subroutine print_matrix_debug(AA, nn, mm, mat_name, debug)
    integer*4 nn, mm, kk
    real*4 AA(nn,mm) !nn = #rows, mm = #columns
    character(*) mat_name
    logical debug

    
    IF (debug .eqv. .true.) then
        PRINT*, mat_name, ' =' 
        DO kk=1,nn
            PRINT '(20f6.2)', AA(kk,1:mm)
        ENDDO
        PRINT*, new_line('a')
    END IF
end subroutine