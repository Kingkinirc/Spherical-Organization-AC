PROGRAM ANALYSISMULTIRINGINSPHERE
    IMPLICIT NONE
    INTEGER,PARAMETER:: npart= 300, niter= 10001, nring = 6, nsmall = 50
    REAL*8,PARAMETER:: mass = 1.0d0, sphererad= 8.250d0, pi = 3.1417
    REAL*8,PARAMETER:: shellsize=0.250d0
    REAL*8:: pos((niter*npart*nring),3)
    REAL*8,ALLOCATABLE:: mondist50(:), mondist100(:), timeavgmondist50(:), timeavgmondist100(:)
    REAL*8,ALLOCATABLE:: mondist150(:), mondist200(:), mondist250(:), mondist125(:)
    REAL*8,ALLOCATABLE:: timeavgmondist150(:),timeavgmondist200(:),timeavgmondist125(:),timeavgmondist250(:)
    REAL*8:: densitybig, densitysmall, n2, a, b, c, d, e, f
    REAL*8:: rmon50, rmon100, rmon125, rmon150, rmon200, rmon250, denbig, densmall
    INTEGER:: i, k, time, m, n1, p, ab

    m = (int(sphererad/shellsize)) + 1

    p = (int(sphererad/shellsize)) + 1

    densitybig = dfloat((npart-(3*nsmall))*nring*3)/(4*pi*(sphererad**3))
    densitysmall = dfloat(nsmall*nring*3)/(4*pi*(sphererad**3))
    denbig = densitybig/nring
    densmall = densitysmall/nring

    ALLOCATE(mondist50(p))
    ALLOCATE(mondist100(p))
    ALLOCATE(mondist150(p))
    ALLOCATE(mondist125(p))
    ALLOCATE(mondist200(p))
    ALLOCATE(mondist250(p))
    ALLOCATE(timeavgmondist50(p))
    ALLOCATE(timeavgmondist100(p))
    ALLOCATE(timeavgmondist125(p))
    ALLOCATE(timeavgmondist150(p))
    ALLOCATE(timeavgmondist200(p))
    ALLOCATE(timeavgmondist250(p))

    OPEN(unit=30, file="arc3_300_posdata_n6_vf0.2.dat", status="old")
    !OPEN(unit=31, file="rinsph_xyzCM_n100.csv",status="unknown")
    !OPEN(unit=32, file="arc3_300_n4_rCM.csv", status="unknown")
    OPEN(unit=33, file="arc3_300_n6_cldist_vf0.2.csv", status="unknown")

    DO time = 1, niter
        ab = (time-1)*npart*nring
        DO i= 1,9
            READ(30,*)
        ENDDO

        DO i = 1, nring*npart
            READ(30,*) k, pos((ab+k),1), pos((ab+k),2), pos((ab+k),3)
        ENDDO
    ENDDO

    timeavgmondist50 = 0.0d0; timeavgmondist100= 0.0d0; timeavgmondist150=0.0d0
    timeavgmondist200= 0.0d0; timeavgmondist250=0.0d0; timeavgmondist125=0.0d0
    DO time = 1, niter
        mondist50 = 0.0d0; mondist100= 0.0d0; mondist150=0.0d0; mondist200= 0.0d0; mondist250=0.0d0;
        mondist125 = 0.0d0
        !r2CMloop1=0.0d0; rCMloop2=0.0d0; r2CMloop2=0.0d0; rCMloop3=0.0d0; r2CMloop3=0.0d0
            DO i = 1, nring
                DO k = 1, npart
                ab = (time-1)*nring*npart + (i-1)*npart + k
                !if (i.eq.1) then
                if (k.eq.50) then
                    rmon50 = dsqrt((pos(ab,1))**2 + (pos(ab,2))**2 + (pos(ab,3))**2)
                    n1 = int(rmon50/shellsize)
                    mondist50(n1) = mondist50(n1) + 1
                elseif (k.eq.100) then
                    rmon100 = dsqrt((pos(ab,1))**2 + (pos(ab,2))**2 + (pos(ab,3))**2)
                    n1 = int(rmon100/shellsize)
                    mondist100(n1) = mondist100(n1) + 1
                elseif (k.eq.150) then
                    rmon150 = dsqrt((pos(ab,1))**2 + (pos(ab,2))**2 + (pos(ab,3))**2)
                    n1 = int(rmon150/shellsize)
                    mondist150(n1) = mondist150(n1) + 1
                elseif (k.eq.200) then
                    rmon200 = dsqrt((pos(ab,1))**2 + (pos(ab,2))**2 + (pos(ab,3))**2)
                    n1 = int(rmon200/shellsize)
                    mondist200(n1) = mondist200(n1) + 1
                elseif (k.eq.250) then
                    rmon250 = dsqrt((pos(ab,1))**2 + (pos(ab,2))**2 + (pos(ab,3))**2)
                    n1 = int(rmon250/shellsize)
                    mondist250(n1) = mondist250(n1) + 1
                elseif (k.eq.125) then
                    rmon125 = dsqrt((pos(ab,1))**2 + (pos(ab,2))**2 + (pos(ab,3))**2)
                    n1 = int(rmon125/shellsize)
                    mondist125(n1) = mondist125(n1) + 1
                endif
                !endif
                ENDDO
            ENDDO
            timeavgmondist50 = timeavgmondist50 + mondist50
            timeavgmondist100 = timeavgmondist100 + mondist100
            timeavgmondist125 = timeavgmondist125 + mondist125
            timeavgmondist150 = timeavgmondist150 + mondist150
            timeavgmondist200 = timeavgmondist200 + mondist200
            timeavgmondist250 = timeavgmondist250 + mondist250
    ENDDO 

    !WRITE(31,*) "Co-ordinate:", ",", "x", ",", "y", ",", "z"
    !DO time = 1, niter            
        !WRITE(31,*) (((time-1)*1000)), ",", posCM(time,1), ",", posCM(time,2), ",", posCM(time,3)
    !ENDDO

    timeavgmondist50 = timeavgmondist50/dfloat(niter*nring)
    timeavgmondist100 = timeavgmondist100/dfloat(niter*nring)
    timeavgmondist125 = timeavgmondist125/dfloat(niter*nring)
    timeavgmondist150 = timeavgmondist150/dfloat(niter*nring)
    timeavgmondist200 = timeavgmondist200/dfloat(niter*nring)
    timeavgmondist250 = timeavgmondist250/dfloat(niter*nring)

DO i = 1, p
    n2 = (2*i-1)*shellsize/2
    a = timeavgmondist50(i)
    b = timeavgmondist100(i)
    c = timeavgmondist150(i)
    d = timeavgmondist200(i)
    e = timeavgmondist250(i)
    f = timeavgmondist125(i)
    WRITE(33,*) n2,",",a,",",b,",",c,",",d,",",e, ",", f
ENDDO

END PROGRAM ANALYSISMULTIRINGINSPHERE