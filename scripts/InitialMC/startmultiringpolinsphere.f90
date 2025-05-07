PROGRAM MULTIRINGPOL2SPHEREINITIALTXT
    IMPLICIT NONE
    INTEGER,PARAMETER:: nring=6, npart= 300, nsmall=300, mcniter = 20000
    REAL*8,PARAMETER:: lb = 1.0, pi = 3.1417, radsphere=8.250d0, sigma=0.80d0, T = 1.0d0, eps = 1.0d0, sprconst= 200.0d0
    REAL*8,PARAMETER:: radius= 5.0, theta = (2*pi)/npart, dela = lb/10, rs = (1.122 + 2.0d0)
    INTEGER:: i, j, k, p, q, s, mm, time, maxneigh, b, c, d, nn, size, count
    REAL*8:: lx, ly, lz, dx, dy, dz, dr, x1, x2, y1, y2, z1, z2, r1, r2, E, Ei, Ef, dE, Energy, LJ, rpos, E1, E2
    REAL*8:: n, p1, p2, p3, s1, s2, s3, h, u, x3, y3, z3, r3
    REAL*8,DIMENSION(:,:),ALLOCATABLE:: polarpos,pos
    INTEGER,DIMENSION(:),ALLOCATABLE:: no_neigh
    INTEGER,DIMENSION(:,:),ALLOCATABLE:: ne_list
    INTEGER,ALLOCATABLE:: seed(:)
    CHARACTER(len=23)::c1='ringpol3spherestart300_'
    CHARACTER(len=32)::c2
    CHARACTER(len=5)::c3
    CHARACTER(len=4)::c4= '.xyz'

    maxneigh = npart*nring -1
    ALLOCATE(pos((npart*nring),3))
    ALLOCATE(polarpos((npart*nring),2))
    ALLOCATE(no_neigh(npart*nring))
    ALLOCATE(ne_list((npart*nring),maxneigh))

    CALL SYSTEM_CLOCK(count)
    CALL RANDOM_SEED(size=size)
    ALLOCATE(seed(size), source= count+37*[(i,i=0,size-1)])
    CALL RANDOM_SEED(put=seed)
!----------------------------------------------------------------------------------------------------------

    !OPEN(unit=41,file='ringpolarc3spherestart400n6_    0.xyz',status='unknown')
    !WRITE(41,*) npart*nring
    !WRITE(41,*) ""
!-----------------------------------------------------------------------------------------------

    OPEN(unit=40,file='ringpol_300_spherestart_n6_vf0.2.txt',status='unknown')
    WRITE(40,*) "LAMMPS data file for ring polymer of length 6*600"
    WRITE(40,*) ""
    WRITE(40,*) nring*npart, "atoms"
    WRITE(40,*) nring*npart,"bonds"
    WRITE(40,*) ""
    WRITE(40,*) nring, "atom types"
    WRITE(40,*) 1, "bond types"
    WRITE(40,*) ""
    WRITE(40,*) "Masses"
    WRITE(40,*) ""
    DO i = 1, nring
        WRITE(40,*) i, 1
    ENDDO
    WRITE(40,*) ""
    WRITE(40,*) "Atoms"
    WRITE(40,*) ""
!-----------------------------------------------------------------------------------------

    DO i = 1,nring
        DO j= 1,npart
            k = (i-1)*npart + j
            polarpos(k,1) = radius
            polarpos(k,2) = theta*(k-1)
        ENDDO
    ENDDO

    DO i= 1,nring
        DO j = 1, npart
            k = (i-1)*npart + j
            pos(k,3) = -6.0d0 + 2*(i-1)
            pos(k,1) = (polarpos(k,1))*COS(polarpos(k,2))
            pos(k,2) = (polarpos(k,1))*SIN(polarpos(k,2))
        ENDDO
    ENDDO


    DO i= 1,nring
        DO j = 1, npart
            k = (i-1)*npart + j
            !WRITE(41,*) pos(k,1), pos(k,2), pos(k,3)
        ENDDO
    ENDDO
!-------------------------------------------------------------------------------------------------------
    CALL neighbour_list(pos,no_neigh,ne_list,npart,maxneigh,rs,nring)

    E = 0.0d0
    DO i=1, npart*nring
        IF (mod(i,npart).ne.0) THEN
          dx = pos(i+1,1) - pos(i,1)
          dy = pos(i+1,2) - pos(i,2)
          dz = pos(i+1,3) - pos(i,3)
          dr = dsqrt(dx**2 + dy**2 + dz**2)
        ELSE
          dx = pos(i,1) - pos((i-npart+1),1)
          dy = pos(i,2) - pos((i-npart+1),2)
          dz = pos(i,3) - pos((i-npart+1),3)
          dr = dsqrt(dx**2 + dy**2 + dz**2)
        ENDIF
        E = E + Energy(dr,sprconst,lb) + LJ(i,npart,pos,eps,sigma,no_neigh,ne_list,maxneigh,nring)
    ENDDO

    !#MC steps#!
      DO time= 1, mcniter  
        IF (mod(time,10)==0) THEN
          CALL neighbour_list(pos,no_neigh,ne_list,npart,maxneigh,rs,nring) 
        END IF
        E1= 0.0d0; E2 = 0.0d0

        !#energycalc#
        !OPEN(2,file="EnergyLJP60N200000T1.0sig0.8.csv",status="unknown")
        !OPEN(4, file="LJbondlength5&6T3.0.csv")
        DO mm= 1,(npart*nring)
          CALL RANDOM_NUMBER(n); i = int(n*dfloat(npart*nring)) + 1
          if ((mod(i,npart))==1) then
            x1 = pos(i,1)-pos((i+npart-1),1)
            x2 = pos((i+1),1)-pos(i,1)

            y1 = pos(i,2)-pos((i+npart-1),2)
            y2 = pos((i+1),2)-pos(i,2)

            z1 = pos(i,3)-pos((i+npart-1),3)
            z2 = pos((i+1),3)-pos(i,3)

          else if (mod(i,npart)==0) then
            x1 = pos(i,1)-pos((i-1),1)
            x2 = pos((i+1-npart),1)-pos(i,1)

            y1 = pos(i,2)-pos((i-1),2)
            y2 = pos((i+1-npart),2)-pos(i,2)

            z1 = pos(i,3)-pos((i-1),3)
            z2 = pos((i+1-npart),3)-pos(i,3)

          else
            x1 = pos(i,1)-pos((i-1),1)
            x2 = pos((i+1),1)-pos(i,1)

            y1 = pos(i,2)-pos((i-1),2)          
            y2 = pos((i+1),2)-pos(i,2)

            z1 = pos(i,3)-pos((i-1),3)
            z2 = pos((i+1),3)-pos(i,3)
          end if  
          r1= dsqrt(x1*x1 + y1*y1 + z1*z1)
          r2= dsqrt(x2*x2 + y2*y2 + z2*z2)

          Ei = Energy(r1,sprconst,lb)+Energy(r2,sprconst,lb)+LJ(i,npart,pos,eps,sigma,no_neigh,ne_list,maxneigh,nring)        

          !Displacement
          CALL RANDOM_NUMBER(p1); s1 = (p1-0.50d0)*(dela/0.50d0)
          CALL RANDOM_NUMBER(p2); s2 = (p2-0.50d0)*(dela/0.5d0)
          CALL RANDOM_NUMBER(p3); s3 = (p3-0.50d0)*(dela/0.50d0)

          pos(i,1) = pos(i,1)+s1
          pos(i,2) = pos(i,2)+s2
          pos(i,3) = pos(i,3)+s3

          rpos= dsqrt(pos(i,1)**2 + pos(i,2)**2 + pos(i,3)**2)

          if(rpos.ge.radsphere) then
            pos(i,1) = pos(i,1)-s1
            pos(i,2) = pos(i,2)-s2
            pos(i,3) = pos(i,3)-s3
          else

            if ((mod(i,npart)==1)) then
                x1 = pos(i,1)-pos((i+npart-1),1)
                x2 = pos((i+1),1)-pos(i,1)

                y1 = pos(i,2)-pos((i+npart-1),2)
                y2 = pos((i+1),2)-pos(i,2)

                z1 = pos(i,3)-pos((i+npart-1),3)
                z2 = pos((i+1),3)-pos(i,3)
  
            else if (mod(i,npart)==0) then
                x1 = pos(i,1)-pos((i-1),1)
                x2 = pos((i-npart+1),1)-pos(i,1)

                y1 = pos(i,2)-pos((i-1),2)
                y2 = pos((i-npart+1),2)-pos(i,2)

                z1 = pos(i,3)-pos((i-1),3)
                z2 = pos((i-npart+1),3)-pos(i,3)
            else    
                x1 = pos(i,1)-pos((i-1),1)
                x2 = pos((i+1),1)-pos(i,1)

                y1 = pos(i,2)-pos((i-1),2)          
                y2 = pos((i+1),2)-pos(i,2)

                z1 = pos(i,3)-pos((i-1),3)
                z2 = pos((i+1),3)-pos(i,3)
            endif

            r1= dsqrt(x1*x1 + y1*y1 + z1*z1)
            r2= dsqrt(x2*x2 + y2*y2 + z2*z2)

            Ef = Energy(r1,sprconst,lb) + Energy(r2,sprconst,lb) + LJ(i,npart,pos,eps,sigma,no_neigh,ne_list,maxneigh,nring)

            dE = Ef - Ei
            !print*, time, mm, dE

            if (dE.le.0.0d0) then
                E = E + dE
            else
                u = exp(-dE/T)
                CALL RANDOM_NUMBER(h)
                if (h.lt.u) then
                    E = E + dE
                else
                    pos(i,1) = pos(i,1)-s1
                    pos(i,2) = pos(i,2)-s2
                    pos(i,3) = pos(i,3)-s3
                end if        
            end if
          endif
        ENDDO

        DO i=1, npart*nring
            IF ((mod(i,npart)).ne.0) THEN
                dx = pos(i+1,1) - pos(i,1)
                dy = pos(i+1,2) - pos(i,2)
                dz = pos(i+1,3) - pos(i,3)
                dr = dsqrt(dx**2 + dy**2 + dz**2)
            ELSE
                dx = pos(i,1) - pos((i-npart+1),1)
                dy = pos(i,2) - pos((i-npart+1),2)
                dz = pos(i,3) - pos((i-npart+1),3)
                dr = dsqrt(dx**2 + dy**2 + dz**2)
            ENDIF
            E1 = E1 + Energy(dr,sprconst,lb)
            E2 = E2 + LJ(i,npart,pos,eps,sigma,no_neigh,ne_list,maxneigh,nring)
        ENDDO
        PRINT*, "Iteration",time, E1, E2

       !PRINT*, "Iteration:",",", time, E,"Energy per monomer:", ",",E/dfloat(npart)
       !WRITE(2,*) time, "," , E/dfloat(npart)

        if (mod(time,1000)==0) then
            !write(c3,'(i5)') time
            !c2 = c1//c3//c4
            !OPEN(file= c2, unit=3, status="unknown")
            !WRITE(3,*) npart
            !WRITE(3,*) ""
            !Do i = 1, npart
                !WRITE(3,*) pos(i,1), pos(i,2), pos(i,3)
            !End do
            !CLOSE(3)
        end if
    ENDDO
    !OPEN(50, file='test.dat', status='unknown')
    DO i= 1, nring
        DO j = 1, npart
            k = npart*(i-1)+j
            !WRITE(50,*) pos(i,1), pos(i,2), pos(i,3)
            WRITE(40,*) k, 0, i, pos(k,1), pos(k,2), pos(k,3), 0, 0, 0
        ENDDO
    ENDDO

    WRITE(40,*) ""
    WRITE(40,*) "Velocities"
    WRITE(40,*) ""

    DO i = 1, npart*nring
        WRITE(40,*) i, 0, 0, 0
    ENDDO

    WRITE(40,*) ""
    WRITE(40,*) "Bonds"
    WRITE(40,*) ""

    DO i = 1, npart*nring
            IF ((mod(i,npart)).ne.0) THEN
                WRITE(40,*) i, 1, i, i+1
            ELSE
                WRITE(40,*) i, 1, i, (i-npart+1)
            ENDIF
    ENDDO

END PROGRAM MULTIRINGPOL2SPHEREINITIALTXT


REAL*8 FUNCTION Energy(p,sprconst,a) 
     IMPLICIT NONE
     REAL*8,INTENT(IN):: p
     REAL*8,INTENT(IN):: sprconst, a     
      Energy = 0.50d0*sprconst*((p-a)**2)
      RETURN
  END FUNCTION Energy

 
  REAL*8 FUNCTION LJ(c,npart,pos,eps,sigma,no_neigh,ne_list,maxneigh,nring) RESULT(V)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::c, npart, maxneigh, nring
    REAL*8,INTENT(IN):: eps,sigma
    REAL*8,DIMENSION(npart*nring,3),INTENT(IN):: pos
    INTEGER,DIMENSION(npart*nring),INTENT(IN):: no_neigh
    INTEGER,DIMENSION(npart*nring,maxneigh),INTENT(IN):: ne_list
    INTEGER:: j, p
    REAL*8:: x,y,z,r
    V=0.0d0
      DO j= 1, no_neigh(c)
        p= ne_list(c,j)
        x= pos(p,1)-pos(c,1)
        y= pos(p,2)-pos(c,2)
        z= pos(p,3)-pos(c,3)

        r= dsqrt(x*x + y*y + z*z)
        IF (r.lt.((2**(1/6))*sigma)) THEN 
          V = V + 4*eps*(((sigma/r)**12) - ((sigma/r)**6)) + eps
        ELSE
          V= V+ 0.0d0
        END IF
      END DO
  END FUNCTION LJ

  SUBROUTINE neighbour_list(pos,no_neigh,ne_list,npart,maxneigh,rs,nring)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::npart,nring,maxneigh
    REAL*8,INTENT(IN):: rs
    REAL*8,DIMENSION(npart*nring,3),INTENT(IN):: pos
    INTEGER,DIMENSION(npart*nring),INTENT(INOUT):: no_neigh
    INTEGER,DIMENSION(npart*nring,maxneigh),INTENT(INOUT):: ne_list
    REAL*8:: x,y,z,r
    INTEGER:: ij, jk
    no_neigh=0; ne_list=0
    DO ij = 1, npart*nring
      DO jk = 1, npart*nring
        IF(ij.ne.jk) THEN
        x= pos(ij,1) - pos(jk,1)
        y= pos(ij,2) - pos(jk,2)
        z= pos(ij,3) - pos(jk,3)

        r= dsqrt(x*x + y*y + z*z)

        IF (r.lt.rs) THEN
          no_neigh(ij)= no_neigh(ij) + 1
          ne_list(ij,no_neigh(ij)) = jk
        END IF
        ENDIF
      END DO
    END DO 
  END SUBROUTINE