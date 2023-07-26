    module intro !declare variables
        implicit none
        integer*4 int4_1, int4_2 !i declare only the highest-value/precision data type's variable.
        real*8 pi_8, sqrt2_8

    end module intro

    program Ex2
    use intro
    implicit none
       
        !integer*2 :: int2_1 = 2000000 !Arithmetic overflow converting INTEGER(4) to INTEGER(2) at (1). 
                                        !This check can be disabled with the option '-fno-range-check'
        
        int4_1 = 2000000
        int4_2 = 1

!        print*, new_line('a'), 'Insert two Integers:' !these 2 lines are useful to understand overflowing sums.
!        read*, int4_1, int4_2                         !you can enter different values.
        
        pi_8 = 4.D0*DATAN(1.D0) * 1d1**32 !scientific notation
        sqrt2_8 =  sqrt(2d0) * 1d1**21

        print*, new_line('a'), 'Using Integer*2:'
        call sum2int2(int4_1, int4_2)
        print*, new_line('a'), 'Using Integer*4:'
        call sum2int4(int4_1, int4_2)
        print*, new_line('a'), 'Using Real*4:'
        call sum2re4(pi_8, sqrt2_8)
        print*, new_line('a'), 'Using Real*8:'
        call sum2re8(pi_8, sqrt2_8)
        print*, new_line('a')

    end program Ex2

    subroutine sum2int2(x4, y4)
        integer*4 x4, y4
        integer*2 x2, y2, sum2
        x2 = int(x4, 2)         !transform the integer*4 into int*2.
        y2 = int(y4, 2)
        sum2 = x2 + y2
        print*, 'The sum of', x4, 'and', y4, 'is', sum2
    end subroutine sum2int2

    subroutine sum2int4(x4, y4)
        integer*4 x4, y4, sum4
        sum4 = x4 + y4
        print*, 'The sum of', x4, 'and', y4, 'is', sum4
    end subroutine sum2int4

    subroutine sum2re4(x, y)
        real*8 x, y
        real*4 re4_1, re4_2, sum
        re4_1 = real(x, 4)
        re4_2 = real(y, 4)
        sum = re4_1 + re4_2
        print*, 'The sum of', re4_1, 'and', re4_2, 'is', sum 
    end subroutine sum2re4

    subroutine sum2re8(x, y)
        real*8 x, y, sum
        sum = x + y
        print*, 'The sum of', x, 'and', y, 'is', sum 
    end subroutine sum2re8