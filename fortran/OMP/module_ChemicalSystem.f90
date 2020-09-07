module ChemicalSystem
    implicit none

    integer(kind=4), parameter                  :: num_reactions = 9
    integer(kind=4), parameter                  :: num_species   = 5
    integer(kind=4), parameter, dimension(5,9)  :: nu = reshape((/ -2, 1, 0, 0, 0, 2, -1, 0, 0, 0, 0, -1, -1, 1, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1 /),(/5,9/))
       real(kind=8), dimension(5)               :: X0 = (/ 10, 0, 20, 0, 0 /)
       real(kind=8), dimension(9)               :: c = (/ 50, 1000, 50, 1000, 1, 10, 3, 1, 6 /)

    contains

    function propensities(X) result(a)
        real(kind=8), dimension(num_species), intent(in) :: X
        real(kind=8), dimension(num_reactions)           :: a

        a(1) = c(1) * 0.5d0 * X(1) * (X(1)-1)
        a(2) = c(2) * X(2)
        a(3) = c(3) * X(2) * X(3)
        a(4) = c(4) * X(4)
        a(5) = c(5) * X(3)
        a(6) = c(6) * X(4)
        a(7) = c(7) * X(5)
        a(8) = c(8) * X(1)
        a(9) = c(9) * X(5)
    end function

    function propJacobian(X) result(Ja)
        real(kind=8), dimension(num_species), intent(in)   :: X
        real(kind=8), dimension(num_reactions,num_species) :: Ja

        Ja(1,1) = c(1) * 0.5d0 * ( 2*X(1) - 1 )
        Ja(2,1) = 0
        Ja(3,1) = 0
        Ja(4,1) = 0
        Ja(5,1) = 0
        Ja(6,1) = 0
        Ja(7,1) = 0
        Ja(8,1) = c(8)
        Ja(9,1) = 0
        Ja(1,2) = 0
        Ja(2,2) = c(2)
        Ja(3,2) = c(3) * X(3)
        Ja(4,2) = 0
        Ja(5,2) = 0
        Ja(6,2) = 0
        Ja(7,2) = 0
        Ja(8,2) = 0
        Ja(9,2) = 0
        Ja(1,3) = 0
        Ja(2,3) = 0
        Ja(3,3) = c(3) * X(2)
        Ja(4,3) = 0
        Ja(5,3) = c(5)
        Ja(6,3) = 0
        Ja(7,3) = 0
        Ja(8,3) = 0
        Ja(9,3) = 0
        Ja(1,4) = 0
        Ja(2,4) = 0
        Ja(3,4) = 0
        Ja(4,4) = c(4)
        Ja(5,4) = 0
        Ja(6,4) = c(6)
        Ja(7,4) = 0
        Ja(8,4) = 0
        Ja(9,4) = 0
        Ja(1,5) = 0
        Ja(2,5) = 0
        Ja(3,5) = 0
        Ja(4,5) = 0
        Ja(5,5) = 0
        Ja(6,5) = 0
        Ja(7,5) = c(7)
        Ja(8,5) = 0
        Ja(9,5) = c(9)
    end function

    function relaxRates(X) result(lambda)
        real(kind=8), dimension(num_species), intent(in) :: X
        real(kind=8), dimension(num_reactions)           :: lambda

        lambda(1) =  + c(1) * (2*X(1)-1) + c(2)
        lambda(2) =  + c(1) * (2*X(1)-1) + c(2)
        lambda(3) =  + c(3) * X(3) + c(3) * X(2) + c(4)
        lambda(4) =  + c(3) * X(3) + c(3) * X(2) + c(4)
        lambda(5) = c(5)
        lambda(6) = c(6)
        lambda(7) = c(7)
        lambda(8) =  + c(8)
        lambda(9) =  + c(9)
    end function
end module ChemicalSystem
