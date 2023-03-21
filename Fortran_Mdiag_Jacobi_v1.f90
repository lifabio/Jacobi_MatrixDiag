Module data_mod
	implicit none
	save
	!Double/Single precision definition
    integer, parameter :: single = selected_real_kind(6,37)
	integer, parameter :: double = selected_real_kind(15,307)

    !Global const.
    real(kind=double), parameter    :: pi=3.1415926535897932384626433832795d0       !Def. pi
    complex(kind=double), parameter :: iu=(0.d0, 1.d0)                              !Def. Imaginary Unit
	integer, parameter              :: L_MAX=1000000
    integer(kind=double), parameter :: Nmax=10000
    real(kind=double)               :: toll=1.d-16
    
	!Global Variables
	complex(kind=double)			 :: c0
    integer						     :: i0
    real(kind=double),allocatable    :: V(:,:),psi_Im(:,:),a(:,:),a2(:,:)
    complex(kind=double),allocatable :: H(:,:),Id(:,:),Rd(:,:),Rd2(:,:),Rd3(:,:),D(:,:),U(:,:),D2(:,:),Udag(:,:)
    
    
end module data_mod

module subs_mod
    use data_mod
    implicit none
    contains
    
    subroutine diagonalM(N,diagval,M)
        complex(kind=double), intent(in) :: diagval
        integer,intent(in)               :: N
        complex(kind=double),intent(out) :: M(N,N)
        integer :: i,j
        
        do i=1,N
            do j=1,N
                if (i==j) then
                    M(i,j)=diagval
                else
                    M(i,j)=0.d0
                endif
            enddo
        enddo
    end subroutine diagonalM
    
    subroutine jacobiRot_theta(M,N,k,l,theta)
        integer,intent(in)               :: N,k,l
        real(kind=double),intent(in)     :: theta
        complex(kind=double),intent(out) :: M(N,N)
        integer :: i,j
        real(kind=double) :: c,s
        
        call diagonalM(N,(1.d0,0d0),M)
        c=cos(theta)
        s=sin(theta)
        M(k,k)=c
        M(k,l)=s
        M(l,k)=-s
        M(l,l)=c
    end subroutine jacobiRot_theta
    
    subroutine jacobiRot(M,N,k,l,c,s)
        integer,intent(in)               :: N,k,l
        real(kind=double),intent(in)     :: c,s
        complex(kind=double),intent(out) :: M(N,N)
        integer :: i,j
        
        call diagonalM(N,(1.d0,0d0),M)
        M(k,k)=c
        M(k,l)=s
        M(l,k)=-s
        M(l,l)=c
    end subroutine jacobiRot
    
    subroutine jacobiDiag(A,N,D,U)    !For Tridiagonal Matrix
        integer,intent(in)               :: N
        complex(kind=double),intent(in) :: A(N,N)           !Input Matrix
        complex(kind=double),intent(out) :: D(N,N)          !Output Diagonal Matrix D=U*A*Ut
        complex(kind=double),intent(out) :: U(N,N)          !Output Diagonalization Matrix U (Eigenvectors on columns)
        integer :: i,j,k,l,r_piv,c_piv,mm
        real(kind=double) :: c,s,t,Akk,All,Akl,w,pivot,tau
        complex(kind=double) :: Rj(N,N),Dk(N),Dl(N),Ucol_k(N),Ucol_l(N)
                    
        call diagonalM(N,(1.d0,0d0),U)
        D=A     !Initialization of D
        call ABSmax(D,N,pivot,r_piv,c_piv)
        do while ( ( abs(pivot) ) >toll )
            k=r_piv
            l=c_piv
            Akl=D(k,l)    !element to be zeroed
            Akk=D(k,k)
            All=D(l,l)
            
            !Compute t=tan(theta), c=cos(theta) and s=sin(theta)
            w=(All - Akk)/(2.d0*Akl)
            if (w>=0) then
                t= 1/( dabs(w) + sqrt(w**2.d0 + 1.d0) ) !t=-w + sqrt(w**2.d0 + 1.d0)
            else
                t=-1/( dabs(w) + sqrt(w**2.d0 + 1.d0) ) !t=-w - sqrt(w**2.d0 + 1.d0)
            endif
            c=1/sqrt(1 + t**2.d0)
            s=t/sqrt(1 + t**2.d0)
            tau=s/(1 + c)
            !Update D (kth and lth rows and columns)
            do mm=1,N
                Dk(mm) = D(k,mm) - s*(D(l,mm) + tau*D(k,mm)) !c*D(k,mm) - s*D(l,mm)
                Dl(mm) = D(l,mm) + s*(D(k,mm) - tau*D(l,mm)) !c*D(l,mm) - s*D(k,mm)
            enddo
            do mm=1,N
                D(k,mm)=Dk(mm)
                D(l,mm)=Dl(mm)
                !Column update
                D(mm,k) = D(k,mm)
                D(mm,l) = D(l,mm)
            enddo
            D(k,l)= 0!(c**2.d0 - s**2.d0)*Akl + c*s*(Akk - All)
            D(l,k)=D(k,l)
            D(k,k) = Akk - t*Akl 
            D(l,l) = All + t*Akl
            
            !Update U
            !call jacobiRot(Rj,N,k,l,c,s)               !Compute/Save Jacobi rotation matrix
            !U = MATMUL(U, Rj)
            !Save Colums k,l of old U
            do j=1,N
                Ucol_k(j)=U(j,k)
                Ucol_l(j)=U(j,l)
            enddo
            !Update 
            do j=1,N
                U(j,k)=c*Ucol_k(j) - s*Ucol_l(j)
                U(j,l)=s*Ucol_k(j) + c*Ucol_l(j)
            enddo
            
            !Debug test
            !print*,'piv',pivot,r_piv,c_piv
            !print*,'Akl',Akk,Akl,All
            !print*,'w',w
            !print*,'t',t,' s',s,' c',c
            !write(*,*) 'D='
            !call printRealMatrix(D,N)
            !pause
            
            !Compute the pivot for the next iteration
            call ABSmax(D,N,pivot,r_piv,c_piv)
        enddo
    end subroutine jacobiDiag
    
    subroutine ABSmax(A,N,max,r_max,c_max)
        integer,intent(in)               :: N
        complex(kind=double),intent(in)  :: A(N,N)
        real(kind=double),intent(out)    :: max
        integer,intent(out)              :: r_max,c_max 
        integer :: i,j
        
        !Look for ||max|| in the upper triangular part of the matrix A 
        max=0.d0
        do i=1,N-1
            do j=i+1,N
                if ( abs(A(i,j))>=max ) then
                    max=abs(A(i,j))
                    r_max=i
                    c_max=j
                endif
            enddo
        enddo 
    end subroutine ABSmax
    
    subroutine printRealMatrix(M,N)
        integer,intent(in)               :: N
        complex(kind=double),intent(in) :: M(N,N)
        integer :: i,j
        character(len=20) :: exFmt
        
        write(exFmt,'("(",I0,"(F8.3))")') N
        do i=1,N
            Write( *, exFmt) ( real(M(i,j)), j = 1, N )
        enddo
        write(*,*)
    end subroutine printRealMatrix
    
end module subs_mod
    
program example
    !use omp_lib            !Load OpenMP library
    use data_mod
    use subs_mod
    implicit none

    integer :: n,i,j,k,Np
    real(kind=double) :: dt,dt2,tstart,tend

    
    Np=3
    allocate(H(Np,Np))
    allocate(Id(Np,Np))
    allocate(Rd(Np,Np))
    allocate(Rd2(Np,Np))
    allocate(Rd3(Np,Np))
    allocate(D(Np,Np))
    allocate(U(Np,Np))
    allocate(D2(Np,Np))
    allocate(Udag(Np,Np))
    
    !Creating a tridiagonal matrix
    do i=1,Np
        do j=1,Np
            if (i==j) then
                H(i,j)=-2.d0
            else
                if ( (i==j+1).OR.(i==j-1) ) then
                    H(i,j)=1.d0
                else
                    H(i,j)=0.d0
                endif
            endif
        enddo
    enddo
    
    
    !print H to verify 
    write(*,*) 'H='
    if (Np<200) call printRealMatrix(H,Np)
    
    !call diagonalM(Np,(1.d0,0d0),Id)
    !!print Id to verify 
    !if (Np<200) call printRealMatrix(Id,Np)
    !
    !call jacobiRot_theta(Rd,Np,1,3,pi/3)
    !!print Rd to verify 
    !if (Np<200) call printRealMatrix(Rd,Np)
    !
    !call jacobiRot(Rd2,Np,1,3,cos(pi/3),sin(pi/3))
    !!print Rd2 to verify 
    !if (Np<200) call printRealMatrix(Rd2,Np)
    !
    !Rd3=matmul(2*Id,3*Id)
    !!print Rd2 to verify 
    !if (Np<200) call printRealMatrix(Rd3,Np)
    
    write(*,*) 'Matrix diagonalization:'
    call jacobiDiag(H,Np,D,U)
    write(*,*) 'D='
    if (Np<200) call printRealMatrix(D,Np)
    write(*,*) 'U='
    if (Np<200) call printRealMatrix(U,Np)
    !write(*,*) D(1,5),D(5,1)
    
    Udag=transpose(U)
    D2=matmul(Udag,H)
    D2=matmul(D2,U)
    write(*,*) 'D2='
    if (Np<200) call printRealMatrix(D2,Np)
    pause
    
    
    

end program

