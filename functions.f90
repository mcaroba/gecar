module functions 

contains 





subroutine cross_product(u,v,cross)

  implicit none

  real*8 :: u(1:3), v(1:3), cross(1:3)

  cross(1:3) = (/ u(2)*v(3)-u(3)*v(2), &
                 -u(1)*v(3)+u(3)*v(1), &
                  u(1)*v(2)-u(2)*v(1) /)

end subroutine







pure function matinv3(A) result(B)
! Performs a direct calculation of the inverse of a 3×3 matrix.
! Copy-pasted from Fortran Wiki
  real*8, intent(in) :: A(3,3)
  real*8 :: B(3,3)
  real*8 :: detinv

! Calculate the inverse determinant of the matrix
  detinv = 1.d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function







end module
