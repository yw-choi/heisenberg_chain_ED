program main
    implicit none

    !!! input parameters  !!!
    integer :: L ! chain length
    integer :: NS ! 2*S+1 

    integer :: dimH  ! size of the Hilbert space = NS^L
    double precision, allocatable :: Hmat(:,:), energies(:)
    
    integer :: i,j,is,S_i
    integer, allocatable :: basis(:)
    double precision :: S, Jex, val, t1, t2
    character(len=32) :: arg

    ! read command line argument
    if (iargc().gt.1) then
        call getarg(1,arg)
        read(arg,"(I)") L
        call getarg(2,arg)
        read(arg,"(I)") NS
    else
        ! default values
        L = 4 ! 4 sites
        NS = 2 ! s=1/2
    endif

    dimH = NS**L
    S = (NS-1.d0)/2.d0
    Jex = 1.0d0 ! antiferro

    print *, "Chain Length  = ", L
    print *, "S  = ", S
    print *, "2*S+1  = ", NS
    print *, "Hamiltonian Size = ", dimH
    allocate(Hmat(dimH,dimH))
    allocate(energies(dimH))
    allocate(basis(L+1)) ! one basis state. basis(L+1)==basis(1)

    Hmat = 0.0d0

    ! calculate the matrix element
    call cpu_time(t1)
    print *, ">> Bulding Hamiltonian..."
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
    call cpu_time(t2)
    print "(1x,a,f5.2,a)", "Elapsed time = ", (t2-t1)/60.d0, " min."
    print *, ">> End of Bulding Hamiltonian."


    ! Dump Hamiltonian
    ! open(unit=11, file="H.dat", form="formatted")
    ! do i=1,dimH
    !     do j=1,dimH
    !         write(11, "(F8.2)", advance="no") Hmat(i,j) 
    !     enddo
    !     write(11,*)
    ! enddo
    ! close(11)

    ! Diagonalization
    call diag_lapack(Hmat, energies)

    ! Dump eigenvalues & eigenvectors
    open(unit=11,file="energies.dat",form="formatted")
    do i=1,dimH
        write(11,"(F20.16)") energies(i)
    enddo
    close(11)
    ! ground state only
    open(unit=11,file="ground_state.dat",form="formatted")
    do i=1,dimH
        write(11,*) Hmat(i,1)
    enddo
    close(11)

    print *, "End of run"
contains
    subroutine diag_lapack(Hmat, energies)
        double precision :: Hmat(dimH,dimH), energies(dimH)

        integer :: info, lwork
        double precision, allocatable :: work(:)
        double precision :: t1, t2
        integer :: i
        lwork = 3*dimH-1
        allocate(work(lwork))

        call cpu_time(t1)
        print *, ">> Start of diagonalization"
        call DSYEV( 'V', 'U', dimH, Hmat, dimH, energies, work, lwork, info)
        call cpu_time(t2)
        print "(1x,a,f5.2,a)", "Elapsed time = ", (t2-t1)/60.d0, " min."
        print *, "<< End of diagonalization"

        print *, "DSYEV info = ", info

    end subroutine diag_lapack

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
