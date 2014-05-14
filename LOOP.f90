!INCLUDE 'MATHIO.F90'
!INCLUDE 'TENSOR.F90'
! ############### CONST ###################
MODULE CONST
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
END MODULE CONST
! ############### MODEL ###################
MODULE MODEL
	USE CONST
	INTEGER :: L = 5 
END MODULE MODEL
! ############## PHYSICS ###################
MODULE PHYSICS
	USE TENSORIAL
	INTEGER, ALLOCATABLE :: PTS(:,:) ! coordinate of lattice point
	INTEGER, ALLOCATABLE :: IBK(:), IBD(:) ! bulk and boundary point inds
	INTEGER, ALLOCATABLE :: LKS(:,:) ! links
	TYPE(TENSOR) :: H
CONTAINS
! set up lattice
SUBROUTINE SET_LATTICE()
	USE MODEL
	! local variables
	INTEGER :: NPTS, NBD, NBK, NLK
	INTEGER :: X, Y, Z, I, I0, I1, X0(3), X1(3), A(3)
	
	NPTS = (L+1)**3 ! number of points
	NBD = 2+6*L**2 ! num of boundary pts
	NBK = NPTS - NBD ! cal num of bulk pts
	NLK = 3*L*(L+1)**2 ! num of links
	ALLOCATE(PTS(NPTS,3), IBK(NBK), IBD(NBD), LKS(NLK,2))
	I = 1
	I0 = 1
	I1 = 1
	DO X = 0, L ! X runs slowest
	DO Y = 0, L
	DO Z = 0, L ! Z runs fastest
		PTS(I,:) = [X,Y,Z] ! keep coordinate
		IF (X==0.OR.Y==0.OR.Z==0.OR.X==L.OR.Y==L.OR.Z==L) THEN
			! at boundary
			IBD(I0) = I
			I0 = I0 + 1
		ELSE
			! in the bulk
			IBK(I1) = I
			I1 = I1 + 1
		END IF
		I = I + 1 ! inc
	END DO
	END DO
	END DO
	! set links
	A = [(L+1)**2,L+1,1]
	I1 = 1
	DO I = 1, NPTS
		X0 = PTS(I,:)
		X1 = X0 + [1,0,0]
		IF (X1(1) <= L) THEN
			LKS(I1,1) = I
			LKS(I1,2) = DOT_PRODUCT(X1,A)+1
			I1 = I1 + 1
		END IF
		X1 = X0 + [0,1,0]
		IF (X1(2) <= L) THEN
			LKS(I1,1) = I
			LKS(I1,2) = DOT_PRODUCT(X1,A)+1
			I1 = I1 + 1
		END IF
		X1 = X0 + [0,0,1]
		IF (X1(3) <= L) THEN
			LKS(I1,1) = I
			LKS(I1,2) = DOT_PRODUCT(X1,A)+1
			I1 = I1 + 1
		END IF
	END DO
END SUBROUTINE SET_LATTICE
! make Hamiltonian
SUBROUTINE MAKE_H()
	USE CONST
	! local variables
	INTEGER :: NLK, ILK, REC, IP1, IP2, PT1(3), PT2(3)
	INTEGER, ALLOCATABLE :: DIMS(:), INDS(:), VALS(:)
	
	NPTS = (L+1)**3 ! cal num of sites
	DIMS = [NPTS, NPTS, 5, 5]
	NLK = SIZE(LKS,1) ! get num of links
	ALLOCATE(INDS(10*NLK), VALS(10*NLK))
	DO ILK = 1, NLK ! for each link
		! get end points indices
		IP1 = LKS(ILK,1)
		IP2 = LKS(ILK,2)
		REC = 10*(ILK-1) ! record starting pointer
		! cal sparse matrix inds in the Hamiltonian
		INDS(REC+[1:5]) = (IP1-1) + (IP2-1)*NPTS + [0:4]*NPTS**2*6
		INDS(REC+[6:10])= (IP2-1) + (IP1-1)*NPTS + [0:4]*NPTS**2*6
		! initialize vals by 1
		VALS(REC+[1:10])= 1.
		! get end points coordinates
		PT1 = PTS(IP1,:)
		PT2 = PTS(IP2,:)
		! 
	END DO
END SUBROUTINE MAKE_H
! end of module PHYSICS
END MODULE PHYSICS
! ################ TASK ####################
MODULE TASK
	USE PHYSICS
CONTAINS
! ------------ Data --------------
! ------------ Tests -------------
! test routine
SUBROUTINE TEST()
	USE MATHIO
	
	CALL SET_LATTICE()
	CALL MAKE_H()
END SUBROUTINE TEST
! test SET_LATTICE
SUBROUTINE TEST_LATT()
	USE MODEL
	USE MATHIO
	
	L = 5
	CALL SET_LATTICE()
	CALL EXPORT('PTS', PTS)
	CALL EXPORT('IBK', IBK)
	CALL EXPORT('IBD', IBD)
	CALL EXPORT('LKS', LKS)
END SUBROUTINE TEST_LATT
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	PRINT *, '------------ NLSM -------------'

	CALL TEST()
!	CALL TEST_LATT()
END PROGRAM MAIN