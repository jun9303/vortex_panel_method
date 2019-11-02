!====================================================================
    MODULE MOD_VARS
!====================================================================
    REAL*8, PARAMETER :: PI = 4*ATAN(1.D0)

    INTEGER*8         :: NPTS  ! NUMBER OF NODE POINTS
    INTEGER*8         :: ISSYM ! FLOW SYMMETRY (Y:1, N:0 - KUTTA AT THETA = 150DEG.)
    REAL*8            :: ALPHA ! FREESTREAM FLOW ANGLE

    REAL*8, DIMENSION(:),  ALLOCATABLE :: NODEX, NODEY, CTRLX, CTRLY ! COORDINATES
    REAL*8, DIMENSION(:),  ALLOCATABLE :: THETA, SGMNT               ! GEOMETRY
    REAL*8, DIMENSION(:),  ALLOCATABLE :: GAMMA, FRVEL               ! VRTX STRENGTH, VEL. BY FREESTREAM
    REAL*8, DIMENSION(:),  ALLOCATABLE :: VLCTY, PCOEF, THWTS        ! VELOCITY, C_P, BL SEP. CRITERION
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: CN1, CN2, AN              ! MATRICES FOR BC (NO NORMAL PENETRATION)

    CONTAINS
!====================================================================
! SUBROUTINES========================================================
!====================================================================
    SUBROUTINE READ
!====================================================================
    IMPLICIT NONE
    CHARACTER*10 :: DUMMY
    INTEGER*8 :: N

    OPEN(10,FILE='input.in')
    
    READ(10,*) DUMMY
    READ(10,*) NPTS, ALPHA, ISSYM
    READ(10,*) DUMMY

    ALLOCATE(NODEX(NPTS+1))
    ALLOCATE(NODEY(NPTS+1))

    DO N = 1, NPTS
        READ(10,*) NODEX(N), NODEY(N)
    ENDDO

    CLOSE(10)

    RETURN
    END SUBROUTINE READ
!====================================================================
    SUBROUTINE ALLOCATE
!====================================================================
    IMPLICIT NONE

    ALLOCATE(CTRLX(NPTS))
    ALLOCATE(CTRLY(NPTS))
    ALLOCATE(THETA(NPTS))
    ALLOCATE(SGMNT(NPTS))  
    ALLOCATE(GAMMA(NPTS+1))
    ALLOCATE(FRVEL(NPTS+1))
    ALLOCATE(VLCTY(NPTS))  
    ALLOCATE(PCOEF(NPTS))
    ALLOCATE(THWTS(NPTS))
    ALLOCATE(CN1(NPTS,NPTS))
    ALLOCATE(CN2(NPTS,NPTS))
    ALLOCATE(AN(NPTS+1,NPTS+1))

    RETURN
    END SUBROUTINE ALLOCATE
!====================================================================
    SUBROUTINE KUTTACON
!====================================================================
    IMPLICIT NONE
    REAL*8 :: KT 
    REAL*8 :: ANG1, ANG2, GAM
    REAL*8 :: NX(NPTS), NY(NPTS)
    INTEGER*8 :: M, N, NEXT

    NX = NODEX(1:NPTS)
    NY = NODEY(1:NPTS)

    IF (ISSYM .EQ. 0) THEN
        KT = -5./6. * PI
    ELSE
        KT = PI
    ENDIF

    DO N = 1, NPTS
        ANG1 = ATAN2(NODEY(N),NODEX(N)) - KT
        ANG2 = ATAN2(NODEY(N+1),NODEX(N+1)) - KT
        IF (ANG1 * ANG2 .LT. 0.) THEN
            EXIT
        ENDIF
    ENDDO

    NX = CSHIFT(NX, N-1)
    NY = CSHIFT(NY, N-1)

    NODEX(1:NPTS) = NX
    NODEY(1:NPTS) = NY
    NODEX(NPTS+1) = NX(1)
    NODEY(NPTS+1) = NY(1)

    RETURN
    END SUBROUTINE KUTTACON
!====================================================================
    SUBROUTINE CTRLPNTS
!====================================================================
    IMPLICIT NONE
    INTEGER*8 :: N
    REAL*8, DIMENSION(NPTS) :: ANGLE

    DO N = 1, NPTS
        CTRLX(N) = (NODEX(N)+NODEX(N+1))/2.
        CTRLY(N) = (NODEY(N)+NODEY(N+1))/2.
        ANGLE(N) = ATAN2(CTRLY(N),CTRLX(N))*180./PI
    ENDDO

    WRITE(*,*) '***************ANGLE***************'
    WRITE(*,*) ANGLE

    RETURN
    END SUBROUTINE CTRLPNTS
!====================================================================
    SUBROUTINE GEOMETRY
!====================================================================
    IMPLICIT NONE
    INTEGER*8 :: N
    REAL*8    :: DX, DY

    DO N = 1, NPTS
        DX = NODEX(N+1)-NODEX(N)
        DY = NODEY(N+1)-NODEY(N)
        SGMNT(N) = (DX**2.+DY**2.)**.5
        THETA(N) = ATAN2(DY,DX)
    ENDDO

    RETURN
    END SUBROUTINE GEOMETRY
!====================================================================
    SUBROUTINE BCMATRIX
!====================================================================
    IMPLICIT NONE
    INTEGER*8 :: M, N

    CN1 = 0.
    CN2 = 0.

    DO M = 1, NPTS
        DO N = 1, NPTS
            IF (M .EQ. N) THEN
                CN1(M,N) = -1.
                CN2(M,N) = 1.
            ELSE
                CN2(M,N) = D(M,N) + .5*Q(M,N)*F(M,N)/SGMNT(N) &
                           - (A(M,N)*C(M,N) + D(M,N)*E(M,N))*G(M,N)/SGMNT(N)
                CN2(M,N) = NINT(CN2(M,N)*1E6)/1.E6 ! REMOVING FLOATING POINT ERROR
                CN1(M,N) = .5*D(M,N)*F(M,N) + C(M,N)*G(M,N) - CN2(M,N)
                CN1(M,N) = NINT(CN1(M,N)*1E6)/1.E6 ! REMOVING FLOATING POINT ERROR
            ENDIF
        ENDDO
    ENDDO

    DO M = 1, NPTS
        AN(M,1) = CN1(M,1)
        AN(M,NPTS+1) = CN2(M,NPTS)
        DO N = 2, NPTS
            AN(M,N) = CN1(M,N) + CN2(M,N-1)
        ENDDO
        FRVEL(M) = SIN(THETA(M)-ALPHA)
    ENDDO
    AN(NPTS+1,1) = 1.
    AN(NPTS+1,NPTS+1) = 1.
    DO N = 2, NPTS
        AN(NPTS+1, N) = 0.
    ENDDO
    FRVEL(NPTS+1) = 0. 

    AN = INV(AN)

    RETURN
    END SUBROUTINE BCMATRIX
!==================================================================== 
    SUBROUTINE VELOCITY
!==================================================================== 
    IMPLICIT NONE
    INTEGER*8, PARAMETER :: DIVM = 100
    INTEGER*8 :: M, N, DIV
    REAL*8 :: X, Y, G, VX, VY

    DO M = 1, NPTS
        VX = -COS(ALPHA)
        VY = -SIN(ALPHA)
        DO N = 1, NPTS
            DO DIV = 1, DIVM
                X = (NODEX(N+1)-NODEX(N))/DIVM*(DIV-1) + NODEX(N)
                Y = (NODEY(N+1)-NODEY(N))/DIVM*(DIV-1) + NODEY(N)
                G = (GAMMA(N+1)-GAMMA(N))/DIVM*(DIV-1) + GAMMA(N)

                IF (((CTRLX(M)-X)**2.+(CTRLY(M)-Y)**2.).GT.1.E-32) THEN 
                    VX = VX + (G*SGMNT(N)/DIVM)/((CTRLX(M)-X)**2. &
                         +(CTRLY(M)-Y)**2.)*(-(CTRLY(M)-Y))
                    VY = VY + (G*SGMNT(N)/DIVM)/((CTRLX(M)-X)**2. &
                         +(CTRLY(M)-Y)**2.)*(CTRLX(M)-X)
                ENDIF
            ENDDO
        ENDDO
        VX = VX
        VY = VY
        VLCTY(M) = 2.*(VX**2.+VY**2.)**.5
    ENDDO

    WRITE(*,*) '***************VELOCITY***************'
    WRITE(*,*) VLCTY

    RETURN
    END SUBROUTINE VELOCITY
!==================================================================== 
    SUBROUTINE PRESCOEF
!==================================================================== 
    IMPLICIT NONE
    INTEGER*8 :: N

    DO N = 1, NPTS
        PCOEF(N) = 1 - VLCTY(N)**2.
    ENDDO

    WRITE(*,*) '***************PCOEF****************'
    WRITE(*,*) PCOEF

    RETURN
    END SUBROUTINE PRESCOEF
!==================================================================== 
    SUBROUTINE THWAITES
!==================================================================== 
    IMPLICIT NONE
    REAL*8, DIMENSION(NPTS/2+1) :: S, VU, VL, KU, KL
    REAL*8 :: DVDX, INTD5, SS, SEPANGU, SEPANGL
    INTEGER*8 :: M, N

    SS = ATAN2(NODEY(1),NODEX(1))
    SS = -180.-SS*180./PI

    DO N = 1, NPTS/2+1
        S(N) = .5*(2.*PI/NPTS)*(N-1)
        VU(N) = VLCTY(N)
        VL(N) = VLCTY(NPTS-N+1)
    ENDDO

    KU = 0.
    KL = 0.

    DO N = 2, NPTS/2+1
        DVDX = ((VU(N+1)-VU(N))/(S(N+1)-S(N))+(VU(N)-VU(N-1))/(S(N)-S(N-1)))/2.
        INTD5 = 0.
        DO M = 1, N-1
            INTD5 = INTD5 + (VU(M)**5.+VU(M+1)**5.)/2.*(S(M+1)-S(M))
        ENDDO
        KU(N) = .45/VU(N)**6.*DVDX*INTD5 + .09
        DVDX = ((VL(N+1)-VL(N))/(S(N+1)-S(N))+(VL(N)-VL(N-1))/(S(N)-S(N-1)))/2.
        INTD5 = 0.
        DO M = 1, N-1
            INTD5 = INTD5 + (VL(M)**5.+VL(M+1)**5.)/2.*(S(M+1)-S(M))
        ENDDO
        KL(N) = .45/VL(N)**6.*DVDX*INTD5 + .09
    ENDDO

    WRITE(*,*) '***************KU****************'
    WRITE(*,*) KU

    DO N = 1, NPTS/2 + 1
        IF (KU(N)*KU(N+1).LT.0.) THEN
            SEPANGU = 2.*PI/NPTS*(N-1) + 2.*PI/NPTS*(ABS(KU(N))/(ABS(KU(N)+ABS(KU(N+1)))))
            SEPANGU = SEPANGU*180./PI
            EXIT
        ENDIF
    ENDDO

    DO N = 1, NPTS/2 + 1
        IF (KL(N)*KL(N+1).LT.0.) THEN
            SEPANGL = 2.*PI/NPTS*(N-1) + 2.*PI/NPTS*(ABS(KL(N))/(ABS(KL(N)+ABS(KL(N+1)))))
            SEPANGL = SEPANGL*180./PI
            EXIT
        ENDIF
    ENDDO

    WRITE(*,*) '***************SEPANG (UPPER, LOWER)****************'
    WRITE(*,*) 'STARTING FROM FRONT STAGNATION POINT'
    WRITE(*,*) SEPANGU, SEPANGL
    WRITE(*,*) SS


    RETURN
    END SUBROUTINE THWAITES
!====================================================================
    SUBROUTINE DEALLOCATE
!====================================================================
    IMPLICIT NONE

    DEALLOCATE(CTRLX)
    DEALLOCATE(CTRLY)
    DEALLOCATE(THETA)
    DEALLOCATE(SGMNT)  
    DEALLOCATE(GAMMA)
    DEALLOCATE(FRVEL)
    DEALLOCATE(VLCTY)  
    DEALLOCATE(PCOEF)
    DEALLOCATE(THWTS)
    DEALLOCATE(CN1)
    DEALLOCATE(CN2)
    DEALLOCATE(AN)

    RETURN
    END SUBROUTINE DEALLOCATE
!====================================================================
      SUBROUTINE PRINT_CURRENT_TIME
!====================================================================
    IMPLICIT NONE
    CHARACTER*8 ::  DATE
    CHARACTER*10::  NOW
    CHARACTER*5 ::  ZONE
    INTEGER*8   ::  VALS(8)

    CALL DATE_AND_TIME(DATE,NOW,ZONE,VALS)

    WRITE(*,101) VALS(1),VALS(2),VALS(3),VALS(5),VALS(6),VALS(7)
    WRITE(*,*) ''
101 FORMAT(' @ 'I0.4,'-',I0.2,'-',I0.2,' ',I0.2,':',I0.2,':',I0.2)

    RETURN
    END SUBROUTINE PRINT_CURRENT_TIME
!=====================================================================
! FUNCTIONS ==========================================================
!=====================================================================
    FUNCTION A(I,J) RESULT(AA)
!=====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: AA 

    AA = -(CTRLX(I)-NODEX(J))*COS(THETA(J))-(CTRLY(I)-NODEY(J))*SIN(THETA(J))

    END FUNCTION A
!====================================================================
    FUNCTION B(I,J) RESULT(BB)
!====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: BB

    BB = (CTRLX(I)-NODEX(J))**2. + (CTRLY(I)-NODEY(J))**2.

    END FUNCTION B
!====================================================================
    FUNCTION C(I,J) RESULT(CC)
!====================================================================
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: CC

    CC = SIN(THETA(I)-THETA(J))

    END FUNCTION C
!====================================================================
!====================================================================
    FUNCTION D(I,J) RESULT(DD)
!====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: DD

    DD = COS(THETA(I)-THETA(J))

    END FUNCTION D
!====================================================================
    FUNCTION E(I,J) RESULT(EE)
!====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: EE

    EE = (CTRLX(I)-NODEX(J))*SIN(THETA(J))-(CTRLY(I)-NODEY(J))*COS(THETA(J))

    END FUNCTION E
!====================================================================
    FUNCTION F(I,J) RESULT(FF)
!====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: FF

    FF = LOG(1.+(SGMNT(J)**2.+2*A(I,J)*SGMNT(J))/B(I,J))

    END FUNCTION F
!====================================================================
    FUNCTION G(I,J) RESULT(GG)
!====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: GG

    GG = ATAN2(E(I,J)*SGMNT(J),B(I,J)+A(I,J)*SGMNT(J))

    END FUNCTION G
!====================================================================
    FUNCTION P(I,J) RESULT(PP)
!====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: PP

    PP = (CTRLX(I)-NODEX(J))*SIN(THETA(I)-2.*THETA(J)) &
         + (CTRLY(I)-NODEY(J))*COS(THETA(I)-2.*THETA(J))

    END FUNCTION P
!====================================================================
    FUNCTION Q(I,J) RESULT(QQ)
!====================================================================
    IMPLICIT NONE
    INTEGER*8, INTENT(IN) :: I, J
    REAL*8 :: QQ

    QQ = (CTRLX(I)-NODEX(J))*COS(THETA(I)-2.*THETA(J)) &
         - (CTRLY(I)-NODEY(J))*SIN(THETA(I)-2.*THETA(J))

    END FUNCTION Q
!====================================================================
    FUNCTION INV(M) RESULT(MINV)
!====================================================================
    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN) :: M
    REAL*8, DIMENSION(SIZE(M,1), SIZE(M,2)) :: MINV

    REAL*8, DIMENSION(SIZE(M,1))    :: WK
    INTEGER*8, DIMENSION(SIZE(M,1)) :: IP    
    INTEGER*8 :: N, INFO

    EXTERNAL DGETRF ! LAPACK LIB.
    EXTERNAL DGETRI ! LAPACK LIB.

    MINV = M
    N = SIZE(M,1)

    CALL DGETRF(N,N,MINV,N,IP,INFO)

    IF (INFO .NE. 0) THEN
        STOP 'SINGULAR'
    ENDIF

    CALL DGETRI(N,MINV,N,IP,WK,N,INFO)

    IF (INFO .NE. 0) THEN
        STOP 'FAILURE'
    ENDIF

    END FUNCTION INV
!====================================================================
    END MODULE MOD_VARS
!====================================================================

!====================================================================
    PROGRAM PROJECT_TWO
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE
    INTEGER *8 :: I, M, N

    !=== READ NODE COORDINATES & BEGIN ALLOCATION
    CALL READ()
    CALL ALLOCATE()

    !=== REARRANGE POINT ORDER FOR KUTTA CONDITION
    CALL KUTTACON()

    !=== SET VALUES AND MATRICES BASED ON GIVEN NODES
    CALL CTRLPNTS()
    CALL GEOMETRY()
    CALL BCMATRIX()

    !=== ACQUIRE GAMMA USING INVERSE MATRIX MULTIPLICATION
    DO M = 1, NPTS + 1
            GAMMA(M) = 0.
        DO N = 1, NPTS + 1
            GAMMA(M) = GAMMA(M) + AN(M,N) * FRVEL(N)
        ENDDO
    ENDDO

    !=== CALCULATE TANGENTIAL VELOCITY AT EACH CONTROL POINT
    CALL VELOCITY()

    !=== CALCULATE PRESSURE COEFFICIENT ALONG SURFACE
    CALL PRESCOEF()

    !=== CALCULATE SEPARATION POINTS USING THWAITES' METHOD
    CALL THWAITES()

    !=== PRINT ENDTIME OF EXECUTION
    CALL PRINT_CURRENT_TIME()
    CALL DEALLOCATE()

    STOP
    END PROGRAM PROJECT_TWO
!====================================================================