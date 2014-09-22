      PROGRAM FCC

      IMPLICIT NONE
      INTEGER, PARAMETER :: NPART = 256
      DOUBLE PRECISION, PARAMETER :: DNS = 0.8
      INTEGER          :: G, I, K, NUC, IX, IY, IZ, M, IREF
      DOUBLE PRECISION :: BOXLH, UCL, UCLH, RROOT3, BOXL, VOLUME
      DOUBLE PRECISION :: RANF, DUMMY
      DOUBLE PRECISION :: Q(NPART,4), R(NPART,3) 
      CHARACTER (LEN=20) :: FILENAME
      VOLUME = DFLOAT(NPART)*DNS
      BOXL = VOLUME**(1.D0/3.D0)  
      BOXLH = 0.5D0*BOXL

!     NUMBER OF UNIT CELLS

      I = 1

      DO WHILE (4 * I ** 3 < NPART)
         I = I + 1
      ENDDO

      NUC   = I 
        
!     UNIT CELL LENGTH

      UCL     = 2.D0 * BOXLH / DFLOAT(NUC) 
      UCLH    = 0.5D0 * UCL
      RROOT3  = 1.D0 / DSQRT(3.D0)

!     BUILD THE UNIT CELL 

!     SUBLATTICE A 

      DO K = 1, 3

         R(1,K) = 0.D0

      ENDDO

      DO K = 1, 4

         Q(1,K) = RROOT3

      END DO
 
!     SUBLATTICE B 

      R(2,1) =  UCLH
      R(2,2) =  UCLH
      R(2,3) =  0.D0
      Q(2,1) =  1.D0 
      Q(2,2) =  RROOT3 
      Q(2,3) = -RROOT3 
      Q(2,4) = -RROOT3
!     SUBLATTICE C 

      R(3,1) =  0.D0
      R(3,2) =  UCLH
      R(3,3) =  UCLH
      Q(3,1) =  1.D0 
      Q(3,2) = -RROOT3 
      Q(3,3) =  RROOT3
      Q(3,4) = -RROOT3 

!     SUBLATTICE D 

      R(4,1) =  UCLH
      R(4,2) =  0.D0
      R(4,3) =  UCLH
      Q(4,1) =  1.D0 
      Q(4,2) = -RROOT3 
      Q(4,3) = -RROOT3
      Q(4,4) =  RROOT3 

!     CONSTRUCT THE LATTICE FROM THE UNIT CELL 

      M = 0

      DO IZ = 1, NUC

         DO IY = 1, NUC

            DO IX = 1, NUC

               DO IREF = 1, 4

                  R(IREF+M,1) = R(IREF,1) + UCL * DFLOAT(IX - 1)
                  R(IREF+M,2) = R(IREF,2) + UCL * DFLOAT(IY - 1)
                  R(IREF+M,3) = R(IREF,3) + UCL * DFLOAT(IZ - 1)
                  Q(IREF+M,1) = Q(IREF,1) 
                  Q(IREF+M,2) = Q(IREF,2)
                  Q(IREF+M,3) = Q(IREF,3)
                  Q(IREF+M,4) = Q(IREF,4)
               ENDDO

               M = M + 4

            ENDDO 

         ENDDO

      ENDDO

!     SHIFT THE CENTER OF THE BOX TO THE ORIGIN

      DO I = 1, NPART
         R(I,:) = R(I,:) - BOXLH
      ENDDO

      FILENAME = 'COORDS'
      OPEN(UNIT=1, FILE=FILENAME, STATUS='UNKNOWN')
      DO I = 1, NPART
        WRITE(1,*) R(I,1), R(I,2), R(I,3)
      END DO
      DO I = 1, NPART
        WRITE(1,*) Q(I,1), Q(I,2), Q(I,3), Q(I,4)
      END DO

      END PROGRAM FCC
!     ----------------------------------------------------------------------------------------------

!     DRAWING A UNIFORM RANDOM VARIATE BETWEEN 0 AND 1 

      FUNCTION RANF(DUMMY)

      INTEGER, PARAMETER :: L = 1029, C = 221591, X = 1048576
      INTEGER ::          SEED
      DOUBLE PRECISION :: DUMMY
      SAVE             SEED
      DATA             SEED / 0 /

      SEED = MOD(SEED * L + C, X)
      RANF = DFLOAT(SEED) / DFLOAT(X)

      
      END FUNCTION RANF

