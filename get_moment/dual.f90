        subroutine dual_2d(a1,a2,dual_a1,dual_a2)
            real*8,dimension(2),intent(in):: a1,a2
            real*8,dimension(2),intent(out):: dual_a1,dual_a2
            real*8::a,b,c,d,det_a
            ! a b
            ! c d
            ! [a1 a2]
            ! [dual_a1
            !  dual_a2]

            a = a1(1)
            b = a2(1)
            c = a1(2)
            d = a2(2)

            det_a = a*d-b*c

            dual_a1(1) = d
            dual_a2(1) = -c
            dual_a1(2) = -b
            dual_a2(2) = a

            dual_a1 = dual_a1/det_a
            dual_a2 = dual_a2/det_a
            return

        end subroutine dual_2d
