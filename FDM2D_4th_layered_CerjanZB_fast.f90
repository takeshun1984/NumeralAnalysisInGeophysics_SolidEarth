!! ------------------------------------------------------------------------ !!
!! Fast version of FDM2D_4th_layered_CerjanZB.f90
!!
!! Time-invariant coefficients are computed once before the time loop
!! (same idea as the precomputation in the Python notebooks 08-10):
!!
!!  Memory variable (Crank-Nicolson), e.g. for r_xx:
!!    dr/dt = -(1/tau_s) [ r + F ],
!!    F_xx = (lam+2mu)(tau_eP/tau_s - 1) div v - 2mu (tau_eS/tau_s - 1) dvz/dz
!!    -> r^(n+1) = C1 * r^n + C2 * F
!!       C1 = (1 + TU*DT/2)/(1 - TU*DT/2),  C2 = TU*DT/(1 - TU*DT/2),  TU = -1/tau_s
!!    P_DIV, P_2MU, P_MU : coefficients of F multiplied by C2
!!  Stress: S_DIV, S_2MU, S_MU : coefficients multiplied by DT
!!  Velocity: BX_DT, BZ_DT = DT * (inverse of the mean density)
!!  rigxz (harmonic mean of rigidity at the SXZ position) is included in
!!  P_MU and S_MU.  FD coefficients (4th/2nd order switch) depend only on k
!!  for a flat free surface -> C40X(k), C41X(k), C40Z(k), C41Z(k).
!!
!! This removes the divisions and most multiplications from the time loop
!! (about 1.5x faster than the original with -O3). The results are the same
!! as the original code within float32 round-off.
!! ------------------------------------------------------------------------ !!
module params
    use iso_fortran_env, only: real32
    implicit none
    public
    integer , parameter :: SP  = real32 !! Single Precision 
    real(SP), parameter :: PI  = atan(1.0_SP)*4.0_SP
    real(SP), parameter :: EPS = epsilon(1.0)
end module params

program FDM2D_2nd_layered
    use params
    implicit none

    ! model region
    integer , parameter :: NX  = 800
    integer , parameter :: NZ  = 250
    integer , parameter :: KFS = 25
    real(SP), parameter :: DX  = 0.4_SP
    real(SP), parameter :: DZ  = 0.4_SP
    real(SP), parameter :: DT  = 0.02_SP
    real(SP), parameter :: DTXZ = DT/(DX*DZ)
    integer , parameter :: NTMAX = 8000

    ! absorbing boundary condition
    real(SP), parameter :: beta = 0.09_SP
    integer , parameter :: NA = 20

    ! source informaiton 
    !! location
    real(SP), parameter :: X0 = 40.0_SP
    real(SP), parameter :: Z0 = 15.0_SP
    integer , parameter :: I0 = int(X0/DX)
    integer , parameter :: K0 = int(Z0/DZ)+KFS !! truncated if Z0/DZ is not an integer (actual depth = (K0-KFS)*DZ)
    !! duration & rupture-starting time from simulation start
    real(SP), parameter :: T0 = 3.0_SP  
    real(SP), parameter :: TS = 0.0_SP 
    !! momet tensor (45-deg reverse faulting)
    real(SP), parameter :: mxx = -1.00_SP
    real(SP), parameter :: mzz =  1.00_SP
    real(SP), parameter :: mxz =  0.00_SP
    real(SP), parameter :: MO  =  1.00_SP

    ! FDM coefficients
    real(SP), parameter :: rc40x = 17.0_SP / 16.0_SP / DX
    real(SP), parameter :: rc40z = 17.0_SP / 16.0_SP / DZ
    real(SP), parameter :: rc41x =  1.0_SP / 48.0_SP / DX
    real(SP), parameter :: rc41z =  1.0_SP / 48.0_SP / DZ
    real(SP), parameter :: rd40x = -1.0_SP / 16.0_SP / DX
    real(SP), parameter :: rd40z = -1.0_SP / 16.0_SP / DZ
    real(SP), parameter :: rd41x = -1.0_SP / 48.0_SP / DX
    real(SP), parameter :: rd41z = -1.0_SP / 48.0_SP / DZ  

    ! output 
    integer , parameter :: NTDEC  = 100
    integer , parameter :: NXD = 2
    integer , parameter :: NZD = 2

    real(SP), parameter :: DST = 10.0_SP
    integer , parameter :: NST = int(NX*DX/DST)-1
    integer , parameter :: NSKIP = 2
    integer , parameter :: NWMAX = NTMAX/NSKIP

    ! attenuation
    real(SP), parameter :: F0 = 1.0_SP/T0
    real(SP), parameter :: Q0 = 300.0_SP

    ! basin
    real(SP), parameter :: XB0 = 150.0_SP
    real(SP), parameter :: XB1 = 180.0_SP
    real(SP), parameter :: Bth = 4.0_SP

    ! Stress & Velocity field
    real(SP) :: SXX(-1:NX+2,-1:NZ+2), SZZ(-1:NX+2,-1:NZ+2), SXZ(-1:NX+2,-1:NZ+2)
    real(SP) :: VX (-1:NX+2,-1:NZ+2), VZ (-1:NX+2,-1:NZ+2)
    !! Memory variables
    real(SP) :: RXX(-1:NX+2,-1:NZ+2), RZZ(-1:NX+2,-1:NZ+2), RXZ(-1:NX+2,-1:NZ+2)
    real(SP) :: RXXN, RZZN, RXZN

    ! Absrobing boundary of Cerjan (1985)
    real(SP) :: gx1(NX), gx2(NX)
    real(SP) :: gz1(NZ), gz2(NZ)

    ! Spatial derivatives
    real(SP) :: DXVX , DZVX , DXVZ , DZVZ
    real(SP) :: DXSXX, DZSZZ, DXSXZ, DZSXZ

    ! Physical properties
    real(SP) :: LAM(-1:NX+2,-1:NZ+2), RIG(-1:NX+2,-1:NZ+2), RHO(-1:NX+2,-1:NZ+2)
    real(SP) :: QP (-1:NX+2,-1:NZ+2), QS (-1:NX+2,-1:NZ+2), TU (-1:NX+2,-1:NZ+2)

    ! topography data
    integer  :: kfs_bot(-1:NX+2), kfs_top(-1:NX+2)

    ! output files
    integer  :: ISTX(NST), ISTZ(NST)
    real(SP) :: VWX(NWMAX, NST), VWZ(NWMAX, NST)

    character(256) :: WNAME , ONAME
    character(6)   :: WNAME0, ONAME0
    real(SP) :: div, rot, vxmax

    ! other variables
    integer  :: i , k , is, it, io, ier
    integer  :: ii, kk, it1
    real(SP) :: Z, ZT, VP, VS, RO
    real(SP) :: TAUP, TAUS, TAU, W0
    real(SP) :: T, SDROP, texp, kupper

    real(SP) :: rigxz, rig00, rig01, rig11, rig10

    integer  :: isign

    ! time-invariant coefficients (computed once before the time loop)
    real(SP) :: C1(NX,NZ), P_DIV(NX,NZ), P_2MU(NX,NZ), P_MU(NX,NZ)
    real(SP) :: S_DIV(NX,NZ), S_2MU(NX,NZ), S_MU(NX,NZ)
    real(SP) :: BX_DT(NX,NZ), BZ_DT(NX,NZ)
    real(SP) :: C40X(NZ), C41X(NZ), C40Z(NZ), C41Z(NZ)
    real(SP) :: C2, PIV, shear, HDT

    ! 6 characters  ======
    ONAME0       = "psv.l."
    WNAME0       = "wav.l."

    ! initilization
    VX (:,:) = 0.00_SP
    VZ (:,:) = 0.00_SP
    SXX(:,:) = 0.00_SP
    SZZ(:,:) = 0.00_SP
    SXZ(:,:) = 0.00_SP

    RXX(:,:) = 0.00_SP
    RZZ(:,:) = 0.00_SP
    RXZ(:,:) = 0.00_SP

    ! Layered strucutre (F-net 1D model by Kubo et al. 2002)
    ZT = 0.00_SP
    do i = 0, NX+1
        do k = 0, NZ+1
            Z  = real(k-KFS)*DZ ! km
            if( Z.le.ZT ) then
                VP = 0.00_SP  ! km/s
                VS = 0.00_SP  ! km/s 
                RO = 0.001_SP ! g/cm3
            else if( Z.lt.3.00_SP ) then
                VP = 5.50_SP
                VS = 3.14_SP
                RO = 2.30_SP
            else if( Z.lt.18.0_SP ) then
                VP = 6.00_SP
                VS = 3.55_SP
                RO = 2.40_SP
            else if( Z.lt.33.0_SP ) then
                VP = 6.70_SP
                VS = 3.83_SP
                RO = 2.80_SP
            else 
                VP = 7.80_SP
                VS = 4.46_SP
                RO = 3.20_SP
            end if
            RIG(i,k) = RO*VS*VS
            LAM(i,k) = RO*VP*VP-2.0_SP*RO*VS*VS
            RHO(i,k) = RO 
            QP (i,k) = Q0*2.0_SP
            QS (i,k) = Q0
        end do
    end do

    do i = 1, NX
        kfs_top(i) = max(kfs - 2, 1 )
        kfs_bot(i) = min(kfs + 2, NZ)
    end do

    do is = 1, NST
        ISTX(is) = int(real(is)*DST/DX)
        ISTZ(is) = KFS+1
    end do

    ! model checking
    open(newunit=io,file='model.dat',action='write',status='unknown')
    do k = 1, NZ
    do i = 1, NX
        Z = real(k-kfs)*DZ
        VS = sqrt(RIG(i,k)/RHO(i,k))
        VP = sqrt((LAM(i,k)+2.0_SP*RIG(i,k))/RHO(i,k))
        write(io,'(2F9.3,3ES12.3E3)') &
            real(i)*dx, Z, VP, VS, RHO(i,k)
    end do
    end do
    close(io)

    ! Zener body Blach er al. (1995)
    W0 = 2.0_SP*PI*F0
    do k = 1, NZ
    do i = 1, NX
        TAU  = 1.0_SP/W0*(sqrt(1.0_SP+QP(i,k)**(-2))-1.0_SP/QP(i,k))
        TAUP = 1.0_SP/W0**2/TAU
        TAUS = (1.0_SP+W0*TAU*QS(i,k))/(W0*QS(i,k)-W0**2*TAU)
        
        QP(i,k) = TAUP/TAU
        QS(i,k) = TAUS/TAU
        TU(i,k) = -1.0_SP/TAU
    end do
    end do

    ! absorbing boundary condition by Cerjan (1985)
    ! initiation
    gx1(:) = 1.00_SP
    gx2(:) = 1.00_SP
    gz1(:) = 1.00_SP
    gz2(:) = 1.00_SP
    do i = 1, NA
        gx1(i) = exp(-beta*(1.0_SP-real(i-0.5)/real(na))**2)
        gz1(i) = exp(-beta*(1.0_SP-real(i-0.5)/real(na))**2)
        gx2(i) = exp(-beta*(1.0_SP-real(i)    /real(na))**2)
        gz2(i) = exp(-beta*(1.0_SP-real(i)    /real(na))**2)

        gx1(nx-i+1) = gx2(i)
        gx2(nx-i+1) = gx1(i)
        gz1(nz-i+1) = gz2(i)
        gz2(nz-i+1) = gz1(i)
    end do

    i = NX/2
    open(newunit=io,file='2ndFD.k',action='write',status='unknown')
    do k = 1, NZ
        isign = sign(1, (k - kfs_top(i)) * (kfs_bot(i) - k))
        write(io,'(2I4)') k-KFS, isign
    end do
    close(io)

    ! precompute time-invariant coefficients
    !! FD coefficients: 2nd order near the free surface, 4th order elsewhere
    !! (flat free surface: kfs_top, kfs_bot are the same for all i)
    do k = 1, NZ
        isign = sign(1, (k - kfs_top(1)) * (kfs_bot(1) - k))
        C40X(k) = rc40x + isign * rd40x
        C41X(k) = rc41x + isign * rd41x
        C40Z(k) = rc40z + isign * rd40z
        C41Z(k) = rc41z + isign * rd41z
    end do

    HDT = 0.5_SP*DT
    do k = 1, NZ
    do i = 1, NX
        !! harmonic mean of rigidity at the SXZ position
        rig00 = RIG(i  ,k  )
        rig10 = RIG(i+1,k  )
        rig01 = RIG(i  ,k+1)
        rig11 = RIG(i+1,k+1)
        rigxz = 4.0_SP * rig00*rig01*rig10*rig11/ &
                (rig00*rig01*rig10 + rig00*rig01*rig11 + &
                 rig00*rig10*rig11 + rig01*rig10*rig11 + EPS)

        !! memory variables: r^(n+1) = C1*r^n + C2*F  (QP, QS hold tau_e/tau_s)
        PIV = LAM(i,k) + 2.0_SP*RIG(i,k)
        C1(i,k) = (1.0_SP + TU(i,k)*DT*0.5_SP) / (1.0_SP - TU(i,k)*DT*0.5_SP)
        C2      =  TU(i,k)*DT              / (1.0_SP - TU(i,k)*DT*0.5_SP)
        P_DIV(i,k) = C2 * PIV*(QP(i,k)-1.0_SP)
        P_2MU(i,k) = C2 * 2.0_SP*RIG(i,k)*(QS(i,k)-1.0_SP)
        P_MU (i,k) = C2 * rigxz*(QS(i,k)-1.0_SP)

        !! stress: coefficients multiplied by DT
        S_DIV(i,k) = PIV*QP(i,k)*DT
        S_2MU(i,k) = 2.0_SP*RIG(i,k)*QS(i,k)*DT
        S_MU (i,k) = rigxz*QS(i,k)*DT

        !! velocity: DT * inverse of the mean density at the VX, VZ positions
        BX_DT(i,k) = 2.0_SP/(RHO(i,k)+RHO(i+1,k))*DT
        BZ_DT(i,k) = 2.0_SP/(RHO(i,k)+RHO(i,k+1))*DT
    end do
    end do

    ! Time steps
    do it = 1, NTMAX
        T = real(it)*DT
        ! spatial derivative for velocity field
        ! stress update
        do k = 1, NZ
        do i = 1, NX
            ! spatial derivatives of velocity (4th order, 2nd order near the surface)
            DXVX = (VX(i  ,k  ) - VX(i-1,k  )) * C40X(k) &
                 - (VX(i+1,k  ) - VX(i-2,k  )) * C41X(k)
            DZVZ = (VZ(i  ,k  ) - VZ(i  ,k-1)) * C40Z(k) &
                 - (VZ(i  ,k+1) - VZ(i  ,k-2)) * C41Z(k)
            DXVZ = (VZ(i+1,k  ) - VZ(i  ,k  )) * C40X(k) &
                 - (VZ(i+2,k  ) - VZ(i-1,k  )) * C41X(k)
            DZVX = (VX(i  ,k+1) - VX(i  ,k  )) * C40Z(k) &
                 - (VX(i  ,k+2) - VX(i  ,k-1)) * C41Z(k)
            div   = DXVX + DZVZ
            shear = DXVZ + DZVX

            ! update memory variable: r^(n+1) = C1*r^n + C2*F
            RXXN = RXX(i,k)
            RZZN = RZZ(i,k)
            RXZN = RXZ(i,k)
            RXX(i,k) = C1(i,k)*RXXN + P_DIV(i,k)*div - P_2MU(i,k)*DZVZ
            RZZ(i,k) = C1(i,k)*RZZN + P_DIV(i,k)*div - P_2MU(i,k)*DXVX
            RXZ(i,k) = C1(i,k)*RXZN + P_MU(i,k)*shear

            ! update stress (with absorbing boundary)
            SXX(i,k) = (SXX(i,k) + S_DIV(i,k)*div - S_2MU(i,k)*DZVZ &
                     + HDT*(RXX(i,k)+RXXN)) * gx1(i)*gz1(k)
            SZZ(i,k) = (SZZ(i,k) + S_DIV(i,k)*div - S_2MU(i,k)*DXVX &
                     + HDT*(RZZ(i,k)+RZZN)) * gx1(i)*gz1(k)
            SXZ(i,k) = (SXZ(i,k) + S_MU(i,k)*shear &
                     + HDT*(RXZ(i,k)+RXZN)) * gx2(i)*gz2(k)
        end do
        end do

        SDROP = MO*kupper(T,TS,T0)*DTXZ
        SXX(I0,K0) = SXX(I0,K0)-MXX*SDROP
        SZZ(I0,K0) = SZZ(I0,K0)-MZZ*SDROP
        SXZ(I0  ,K0  ) = SXZ(I0  ,K0  )-MXZ*SDROP*0.25_SP
        SXZ(I0-1,K0  ) = SXZ(I0-1,K0  )-MXZ*SDROP*0.25_SP
        SXZ(I0  ,K0-1) = SXZ(I0  ,K0-1)-MXZ*SDROP*0.25_SP
        SXZ(I0-1,K0-1) = SXZ(I0-1,K0-1)-MXZ*SDROP*0.25_SP

        ! spatial derivative for stress field
        ! velocity update
        do k = 1, NZ
        do i = 1, NX
            ! spatial derivatives of stress
            DXSXX = (SXX(i+1,k  )-SXX(i  ,k  ))*C40X(k) - (SXX(i+2,k  )-SXX(i-1,k  ))*C41X(k)
            DXSXZ = (SXZ(i  ,k  )-SXZ(i-1,k  ))*C40X(k) - (SXZ(i+1,k  )-SXZ(i-2,k  ))*C41X(k)
            DZSZZ = (SZZ(i  ,k+1)-SZZ(i  ,k  ))*C40Z(k) - (SZZ(i  ,k+2)-SZZ(i  ,k-1))*C41Z(k)
            DZSXZ = (SXZ(i  ,k  )-SXZ(i  ,k-1))*C40Z(k) - (SXZ(i  ,k+1)-SXZ(i  ,k-2))*C41Z(k)

            ! update velocity (with absorbing boundary)
            VX(i,k) = (VX(i,k) + (DXSXX+DZSXZ)*BX_DT(i,k)) * gx2(i)*gz1(k)
            VZ(i,k) = (VZ(i,k) + (DXSXZ+DZSZZ)*BZ_DT(i,k)) * gx1(i)*gz2(k)
        end do
        end do

        if( mod(it,NSKIP) == 0 ) then
            do is = 1, NST
                ii = ISTX(is)
                kk = ISTZ(is)

                it1 = IT/NSKIP

                VWX(it1,is) =  VX(ii,kk)
                VWZ(it1,is) = -VZ(ii,kk) ! Convert “up” to “positive”
            end do
        end if

        if( mod(it,ntdec) == 0 ) then
            vxmax = maxval(VX)

            write(6,'(I5,A1,I5,A4,F6.2,A3,1X,ES12.3E3)') it,'/',NTMAX,': T=',real(it)*DT,'[s]', vxmax
            write(ONAME,'(A,I5.5,A4)') ONAME0,it,'.out'
            open(newunit=io, file=ONAME, action='write',status='unknown')
            do k = 1, NZ, NZD
            do i = 1, NX, NXD
                DXVX = (VX(i  ,k  )-VX(i-1,k  ))/DX
                DXVZ = (VZ(i+1,k  )-VZ(i  ,k  ))/DX

                DZVX = (VX(i  ,k+1)-VX(i  ,k  ))/DZ
                DZVZ = (VZ(i,  k  )-VZ(i  ,k-1))/DZ

                div = DXVX+DZVZ
                rot = DXVZ-DZVX
                write(io,'(2F9.3,4ES12.3E3)') &
                    real(i)*dx, real(k-KFS)*dz, VX(i,k), VZ(i,k), div, rot
            end do
            end do
            close(io)
        end if
    end do

    ! waveform output
    do is = 1, NST
        write(WNAME,'(A,I4.4,A4)') WNAME0,int(ISTX(is)*DX),".dat"
        open(newunit=io,file=WNAME,action='write',status='unknown')
        do it = 1, NWMAX
            write(io,'(F9.3,2ES14.5E3)')  real(it)*DT*NSKIP, VWX(it,is), VWZ(it,is)
        end do
        close(io)
    end do

end program FDM2D_2nd_layered

function texp(t, ts, tr)
    use params
    real(SP), intent(in) :: t  !! time
    real(SP), intent(in) :: ts !! rupture start time
    real(SP), intent(in) :: tr !! characteristic time (duration)
    real(SP) :: texp
    real(SP) :: tt

    if (ts <= t) then
        tt = t - ts
        texp = (2 * PI)**2 * tt / (tr * tr) * exp(-2 * PI * tt / tr)
    else
        texp = 0.0
    end if

end function texp

function kupper(t, ts, tr)
    use params
    real(SP), intent(in) :: t  !< time
    real(SP), intent(in) :: ts !< rupture start time
    real(SP), intent(in) :: tr !< rise time (duration)
    real(SP) :: kupper

    !! ----
    if ( ts <= t .and. t <= ts + tr ) then
      kupper = 3 * Pi * ( sin( Pi*(t-ts)/tr ) )**3 / ( 4 * tr )
    else
      kupper = 0.0
    end if

end function kupper
