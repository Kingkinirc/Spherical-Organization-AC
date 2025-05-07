PROGRAM ANALYSISARC1INSPHERE
    IMPLICIT NONE
    INTEGER,PARAMETER:: npart= 200, niter= 10001, nring = 1, nsmall = 50
    REAL*8,PARAMETER:: mass = 1.0d0, sphererad= 4.0d0, pi = 3.1417
    REAL*8,PARAMETER:: shellsize=0.250d0, shellsizemon= 0.250d0
    REAL*8:: pos((niter*npart*nring),3), posCMloop1(niter,3), posCMloop2(niter,3), posCMloop3(niter,3)
    REAL*8:: posCMbig(niter,3)
    INTEGER,ALLOCATABLE:: CMdistsmall(:), CMdistbig(:)
    REAL*8,ALLOCATABLE:: mondistbig(:), mondistsmall(:), timeavgmondistsmall(:), timeavgmondistbig(:)
    REAL*8,ALLOCATABLE:: mondistsmall1(:), mondistsmall2(:), mondistsmall3(:)
    REAL*8,ALLOCATABLE:: timeavgmondistsmall1(:), timeavgmondistsmall2(:), timeavgmondistsmall3(:)
    REAL*8,ALLOCATABLE:: cummondistbig(:), cummondistsmall(:), tavgrmonbig(:), tavgrmonsmall(:)
    REAL*8,ALLOCATABLE:: tavgcummondistbig(:), tavgcummondistsmall(:)
    REAL*8:: rCMbig, r2CMbig, rCMsmall, r2CMsmall, rmon, densitybig, densitysmall, expno, cumexpno
    REAL*8:: rCMloop1, rCMloop2, rCMloop3, r2CMloop1, r2CMloop2, r2CMloop3, denbig, densmall, expno2
    INTEGER:: i, j, k, time, m, n1, n2, p, q, w, ab

    
    !denbig = densitybig/nring
    !densmall = densitysmall/nring

    m = (int(sphererad/shellsize)) + 1
    ALLOCATE(CMdistbig(m))
    ALLOCATE(CMdistsmall(m))

    p = (int(sphererad/shellsizemon)) + 1
    print*, p
    densitybig = ((npart-(nsmall))*nring*3)/(4*pi*((sphererad)**3))
    densitysmall =(nsmall*nring*3)/(4*pi*((sphererad)**3))

    ALLOCATE(mondistsmall(p))
    ALLOCATE(mondistsmall1(p))
    ALLOCATE(mondistsmall2(p))
    ALLOCATE(mondistsmall3(p))
    ALLOCATE(mondistbig(p))
    ALLOCATE(timeavgmondistsmall(p))
    ALLOCATE(timeavgmondistbig(p))
    ALLOCATE(timeavgmondistsmall1(p))
    ALLOCATE(timeavgmondistsmall2(p))
    ALLOCATE(timeavgmondistsmall3(p))
    ALLOCATE(tavgrmonsmall(p))
    ALLOCATE(tavgrmonbig(p))

    ALLOCATE(cummondistsmall(p))
    ALLOCATE(cummondistbig(p))
    ALLOCATE(tavgcummondistsmall(p))
    ALLOCATE(tavgcummondistbig(p))

    OPEN(unit=30, file="arc1_200_posdata_n1_vf0.2.dat", status="old")
    !OPEN(unit=31, file="rinsph_xyzCM_n100.csv",status="unknown")
    !OPEN(unit=32, file="arc1_200_n1_rCM_vf0.2.csv", status="unknown")
    OPEN(unit=33, file="arc1_200_n1_CoMdist_vf0.2.csv", status="unknown")
    OPEN(unit=34, file="arc1_200_n1_mondist_vf0.2.csv", status="unknown")
    OPEN(unit=36, file="arc1_200_n1_rmonbignsmall_vf0.2.csv", status='unknown')
    !OPEN(unit=37, file='arc1_200_n1_rmonsmall1_vf0.2.csv',status='unknown')
    !OPEN(unit=38, file='arc1_200_n1_rmonsmall2_vf0.2.csv',status='unknown')
    !OPEN(unit=39, file='arc1_200_n1_rmonsmall3_binsize0.05.csv',status='unknown')
    OPEN(unit=51, file="arc1_200_n1_cummondist_vf0.2.csv",status='unknown')
    OPEN(unit=52, file='arc1_200_n1_extracalc_vf0.2.csv',status='unknown')

    DO time = 1, niter
        ab = (time-1)*npart*nring
        DO i= 1,9
            READ(30,*)
        ENDDO

        DO i = 1, nring*npart
            READ(30,*) k, pos((ab+k),1), pos((ab+k),2), pos((ab+k),3)
        ENDDO
    ENDDO

    DO time = 1, niter
        DO i = 1, nring*npart
            k = (time-1)*npart*nring + i
            !PRINT*, pos(k,1), pos(k,2), pos(k,3)
        ENDDO
    ENDDO

    posCMbig= 0.0d0; CMdistbig = 0; CMdistsmall = 0
    posCMloop1= 0.0d0; !posCMloop2= 0.0d0; posCMloop3= 0.0d0

    DO time = 1, niter
        rCMbig = 0.0d0; r2CMbig= 0.0d0; rCMsmall=0.0d0; r2CMsmall= 0.0d0; rCMloop1=0.0d0
        r2CMloop1=0.0d0; !rCMloop2=0.0d0; r2CMloop2=0.0d0; rCMloop3=0.0d0; r2CMloop3=0.0d0

        DO j = 1, 3
            DO i = 1, nring
                DO k = 1, npart
                ab = (time-1)*nring*npart + (i-1)*npart + k
                !if (i.eq.1) then
                if ((k.ge.51).and.(k.le.100)) then
                    posCMloop1(time,j) = posCMloop1(time,j) + mass*pos(ab,j)
                !elseif ((k.ge.101).and.(k.le.150)) then
                    !posCMloop2(time,j) = posCMloop2(time,j) + mass*pos(ab,j)
                !elseif ((k.ge.151).and.(k.le.200)) then
                    !posCMloop3(time,j) = posCMloop3(time,j) + mass*pos(ab,j)
                else
                    posCMbig(time,j) = posCMbig(time,j) + mass*pos(ab,j)
                endif
                !endif
                ENDDO
            ENDDO
            posCMloop1(time,j) = posCMloop1(time,j)/(dfloat(nring*nsmall)*mass)
            !posCMloop2(time,j) = posCMloop2(time,j)/(dfloat(nring*nsmall)*mass)
            !posCMloop3(time,j) = posCMloop3(time,j)/(dfloat(nring*nsmall)*mass)
            posCMbig(time,j) = posCMbig(time,j)/(dfloat(nring*(npart-(nsmall)))*mass)
            r2CMbig = r2CMbig + (posCMbig(time,j))**2
            r2CMloop1 = r2CMloop1 + (posCMloop1(time,j))**2
            !r2CMloop2 = r2CMloop2 + (posCMloop2(time,j))**2
            !r2CMloop3 = r2CMloop3 + (posCMloop3(time,j))**2
            r2CMsmall = r2CMsmall+(posCMloop1(time,j))**2
        ENDDO

        rCMloop1 = dsqrt(r2CMloop1)
        !rCMloop2 = dsqrt(r2CMloop2)
        !rCMloop3 = dsqrt(r2CMloop3)
        rCMsmall = dsqrt(r2CMsmall)
        rCMbig = dsqrt(r2CMbig)
        n1 = int(rCMsmall/shellsize) + 1
        CMdistsmall(n1) = CMdistsmall(n1) + 1
        n2 = int(rCMbig/shellsize) + 1
        CMdistbig(n2) = CMdistbig(n2) + 1
        !WRITE(32,*) "timestep",",","rCMbig",",","rCMsmall",",","rCMloop1",",","rCMloop2",",","rCMloop3"
        !WRITE(32,*) (((time-1)*10000)),",",rCMbig,",",rCMsmall,",",rCMloop1,",",rCMloop2,",",rCMloop3
    ENDDO 

    !WRITE(31,*) "Co-ordinate:", ",", "x", ",", "y", ",", "z"
    !DO time = 1, niter            
        !WRITE(31,*) (((time-1)*1000)), ",", posCM(time,1), ",", posCM(time,2), ",", posCM(time,3)
    !ENDDO

    DO i = 1, m
        WRITE(33,*) (2*i-1)*shellsize/2, ",", dfloat(CMdistbig(i))/dfloat(niter), ",", dfloat(CMdistsmall(i))/dfloat(niter)
    ENDDO

    timeavgmondistbig = 0; timeavgmondistsmall = 0; timeavgmondistsmall1 = 0; timeavgmondistsmall2 = 0
    timeavgmondistsmall3 = 0
    DO time= 1, niter
        mondistbig = 0; mondistsmall = 0; mondistsmall1=0; mondistsmall2=0; mondistsmall3=0
        cummondistbig= 0; cummondistsmall=0
        DO i = 1, nring
            DO j = 1, npart
                k =(time-1)*npart*nring + (i-1)*npart + j
                !if (i.eq.6) then
                if ((j.ge.51).and.(j.le.100)) then
                    rmon = sqrt(pos(k,1)**2 + pos(k,2)**2 + pos(k,3)**2)
                    q = int(rmon/shellsizemon) + 1
                    mondistsmall1(q) = mondistsmall1(q) + 1
                    do w = q, p
                        cummondistsmall(w) = cummondistsmall(w) + 1.0
                    enddo
                else
                    rmon = sqrt(pos(k,1)**2 + pos(k,2)**2 + pos(k,3)**2)
                    q = int(rmon/shellsizemon) + 1
                    mondistbig(q) = mondistbig(q) + 1
                    do w = q, p
                        cummondistbig(w) = cummondistbig(w) + 1.0
                    enddo
                endif
                !endif
            ENDDO
        ENDDO

        DO i = 1, p
            tavgrmonbig(i) = tavgrmonbig(i) + mondistbig(i)
            tavgrmonsmall(i)= tavgrmonsmall(i) + mondistsmall1(i)
        ENDDO

        DO i = 1, p
            expno=densitybig*(4.0d0/3.0d0)*pi*(3.0*shellsizemon*((i*shellsizemon)**2)-(3.0*i*shellsizemon**3)+(shellsizemon**3))
        expno2=dfloat(150*nring)*(4.0d0/3.0d0)*pi*(3.0*shellsizemon*((i*shellsizemon)**2)-(3.0*i*shellsizemon**3)+(shellsizemon**3))
            !expno = expno/nring
            cumexpno = dfloat((npart-nsmall)*nring)
            mondistsmall2(i) = mondistbig(i)/expno2
            mondistbig(i) = mondistbig(i)/expno
            cummondistbig(i) = cummondistbig(i)/cumexpno
        ENDDO

        DO i = 1, p                                     !monomerdistribution for small loop
            expno=densitysmall*(4.0d0/3.0d0)*pi*(3.0*shellsizemon*((i*shellsizemon)**2)-(3.0*i*shellsizemon**3)+(shellsizemon**3))
            expno2=dfloat(50*nring)
            expno2= expno2*(4.0d0/3.0d0)*pi*(3.0*shellsizemon*((i*shellsizemon)**2)-(3.0*i*shellsizemon**3)+(shellsizemon**3))
            !expno = expno/nring
            cumexpno = dfloat((nsmall)*nring)
            mondistsmall3(i) = mondistsmall1(i)/expno2
            mondistsmall1(i) = mondistsmall1(i)/expno
            cummondistsmall(i) = cummondistsmall(i)/cumexpno
        ENDDO

        !if (mod(time,500)==0) then
            !print*, time
        !do i = 1, p
            !print*, cummondistbig(i), cummondistsmall(i)
        !enddo
        !endif

        DO i = 1, p
            timeavgmondistbig(i) = timeavgmondistbig(i) + mondistbig(i)
            timeavgmondistsmall1(i)= timeavgmondistsmall1(i) + mondistsmall1(i)
            timeavgmondistsmall2(i) = timeavgmondistsmall2(i) + mondistsmall2(i)
            timeavgmondistsmall3(i)= timeavgmondistsmall3(i) + mondistsmall3(i)
            tavgcummondistbig(i) = tavgcummondistbig(i) + cummondistbig(i)
            tavgcummondistsmall(i) = tavgcummondistsmall(i) + cummondistsmall(i)
        ENDDO
    ENDDO
    timeavgmondistbig = timeavgmondistbig/(niter)
    timeavgmondistsmall1 = timeavgmondistsmall1/(niter)
    timeavgmondistsmall2 = timeavgmondistsmall2/(niter)
    timeavgmondistsmall3 = timeavgmondistsmall3/(niter)
    timeavgmondistsmall = (timeavgmondistsmall1)
    tavgrmonbig = tavgrmonbig/dfloat(150*nring*niter)
    tavgrmonsmall = tavgrmonsmall/dfloat(50*nring*niter)
    tavgcummondistbig = tavgcummondistbig/(niter)
    tavgcummondistsmall = tavgcummondistsmall/(niter)
    DO i = 1, p
    WRITE(36,*) (2*i-1)*shellsizemon/2,",",tavgrmonbig(i),",",tavgrmonsmall(i)
    WRITE(34,*) (2*i-1)*shellsizemon/2,",",timeavgmondistbig(i),",",timeavgmondistsmall(i)
    WRITE(52,*) (2*i-1)*shellsizemon/2,",",timeavgmondistsmall2(i),",",timeavgmondistsmall3(i)
    WRITE(51,*) (2*i-1)*shellsizemon/2,",",tavgcummondistbig(i),",",tavgcummondistsmall(i)
    ENDDO

    !print*, "timeavg:"
    DO i = 1,p
        !print*, tavgcummondistbig(i), tavgcummondistsmall(i)
    ENDDO

    CLOSE(30)
    !CLOSE(31)
    CLOSE(32)
    CLOSE(33)
    CLOSE(34)
    CLOSE(35)
    CLOSE(36)
    CLOSE(51)
    CLOSE(52)

END PROGRAM ANALYSISARC1INSPHERE