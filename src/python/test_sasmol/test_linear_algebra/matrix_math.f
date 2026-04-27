      subroutine matrix_multiply(a, b, c, dim_a1, dim_a2, dim_b2)
      integer dim_a1, dim_a2, dim_b2
      real*4 a(dim_a1, dim_a2)
      real*4 b(dim_a2, dim_b2)
      real*4 c(dim_a1, dim_b2)
      real*4 cij

cf2py  intent(in) :: a, b, dim_a1, dim_a2, dim_b2
cf2py  intent(out):: c
cf2py  intent(hide):: i, j, k

      do 98, i = 1, dim_a1
          do 99, j = 1, dim_b2
              c(i, j) = 0.0
  99      continue
  98  continue

      do 300, i = 1, dim_a1
          do 200, j = 1, dim_b2
              cij = 0.0
              do 100, k = 1, dim_a2
c                  print *, 'a(', i, ',', k, ') = ', a(i, k), ', b(', k, ',', j, ') = '
c     &                b(k, j)
                  cij = cij + a(i, k) * b(k, j)
  100         continue
              c(i, j) = cij
c              print *, 'c(', i, ',', j, ') = ', c(i, j)
  200     continue
  300 continue

      end subroutine matrix_multiply