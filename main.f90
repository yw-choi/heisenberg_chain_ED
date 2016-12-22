program main
    implicit none

    integer, parameter :: L = 2  ! chain length
    integer, parameter :: NS = 2 ! 2*S+1 

    integer :: dimH  ! size of the Hilbert space = NS^L
    double precision, allocatable :: Hmat(:,:)
    
    integer :: i,j,is,S_i
    integer, allocatable :: basis(:)
    double precision :: S, Jex, val

    dimH = NS**L
    S = (NS-1.d0)/2.d0
    Jex = 1.0d0 ! antiferro

    print *, "Chain Length  = ", L
    print *, "S  = ", S
    print *, "2*S+1  = ", NS
    print *, "Hamiltonian Size = ", dimH
    allocate(Hmat(dimH,dimH))
    allocate(basis(L+1)) ! one basis state. basis(L+1)==basis(1)

    Hmat = 0.0d0

    ! calculate the matrix element
    do j=1,dimH
        ! find j-th basis state
        do is=1,L
            S_i = mod(int((j-1)/(NS**(is-1))), NS)
            basis(is) = S_i 
        enddo
        basis(L+1) = basis(1)

        ! diagonal terms
        do is=1,L
            Hmat(j,j) = Hmat(j,j) + Jex*(basis(is)-S)*(basis(is+1)-S)
        enddo

        ! off-diagonal terms
        do is=1,L
            ! S^(+)_(i) * S^(-)_(i+1)
            if (basis(is+1) .ne. 0 .and. basis(is) .ne. (NS-1)) then
                ! 1/2*normalization factor 
                val = Jex*0.5d0*sqrt((S-(basis(is)-S))*(S+(basis(is)-S)+1)* &
                      (S+(basis(is+1)-S))*(S-(basis(is+1)-S)+1))

                ! find target basis index
                i = find_target_basis_index(basis, is, is+1)
                Hmat(i,j) = Hmat(i,j) + val
            endif 

            ! S^(-)_(i) * S^(+)_(i+1)
            if (basis(is+1) .ne. (NS-1) .and. basis(is) .ne. 0) then
                ! 1/2*normalization factor 
                val = Jex*0.5d0*sqrt((S+(basis(is)-S))*(S-(basis(is)-S)+1)* &
                      (S-(basis(is+1)-S))*(S+(basis(is+1)-S)+1))

                ! find target basis index
                i = find_target_basis_index(basis, is+1, is)
                Hmat(i,j) = Hmat(i,j) + val
            endif 
        enddo
    enddo

    do i=1,dimH
        do j=1,dimH
            write(*,"(F8.2)", advance="no") Hmat(i,j)
        enddo
        write(*,*)
    enddo
    
contains
    integer function find_target_basis_index(basis, ip, im) result(i)
        integer :: basis(:), ip, im
        integer :: is
        if (ip==L+1) then
            ip = 1
        endif
        if (im==L+1) then
            im = 1
        endif
        i = 1
        do is=1,L
            if (is.eq.ip) then
                i = i + NS**(is-1)*(basis(is)+1)
            else if (is.eq.im) then
                i = i + NS**(is-1)*(basis(is)-1)
            else
                i = i + NS**(is-1)*basis(is)
            endif
        enddo
        
    end function find_target_basis_index
end program main