!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine psadj(plev    ,plevp1  ,plat    ,plon    ,t       , &
                       press_m ,press_i ,phis_old,phis_new,ps_old  , &
                       ps_new  )
!
!-----------------------------------------------------------------------
!
! Adjust Ps based on difference between "analysis" phis and model phis.
! Also uses T and P arrays
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plev      ! vertical dimension
      integer plevp1    ! "plev+1"
      integer plat      ! latitude dimension
      integer plon      ! longitude dimension
!
      real*8 t       (plon,plat,plev)    ! analysis Temperatures
      real*8 press_m (plon,plat,plev)    ! analysis pressures
      real*8 press_i (plon,plat,plevp1)  ! analysis pressures (interfaces)
      real*8 phis_old(plon,     plat)    ! analysis phis
      real*8 phis_new(plon,     plat)    ! model phis
      real*8 ps_old  (plon,     plat)    ! analysis Ps (horizontally
!                                        ! interpolated to model grid)
!
!     OUTPUTS
!
      real*8 ps_new  (plon,     plat)    ! adjusted model Ps
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 tsurf
      real*8 t_ref1
      real*8 t_ref2
      real*8 t0
      real*8 tbot
      real*8 pbot
      real*8 dtdz
      real*8 lapse
      real*8 boltz
      real*8 avogad
      real*8 mwdair
      real*8 rgas
      real*8 rdair
      real*8 gravit
      real*8 del_phis
      real*8 x
      real*8 threshold
      real*8 tmp
      real*8 z
      real*8 z_min
      real*8 z_incr
      real*8 hkk
      integer i, j, k, kk         ! Indices
!
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2
      boltz   = 1.38065e-23       ! boltzmann's constant ~ J/k/molecule
      avogad  = 6.02214e26        ! avogadro's number ~ molecules/kmole
      mwdair  = 28.966            ! molecular weight dry air ~ kg/kmole

      rgas    = avogad*boltz      ! universal gas constant ~ J/k/kmole
      rdair   = rgas/mwdair       ! constant for dry air   ~ J/k/kg

      t_ref1    = 290.5
      t_ref2    = 255.0
      threshold = 0.001
      
      do j = 1,plat
        do i = 1,plon

          del_phis = phis_old(i,j) - phis_new(i,j)
!
! If difference between analysis and model phis is negligible,
! then set model Ps = analysis
!
          if(abs(del_phis) .le. threshold) then
            ps_new(i,j) = ps_old(i,j)
!
! Else, go nuts...
!
          else
!
! Tbot and Pbot are determined from the first model level that is at
! least 150m above the surface
!
            z_min = 150.
            z     = 0.

            do k = plev,1,-1
              kk     = k
              hkk    = 0.5*( press_i(i,j,k+1) - press_i(i,j,k) )/ &
                                                          press_m(i,j,k)
              z_incr = (rdair/gravit)*t(i,j,k)*hkk
              z      = z + z_incr
              if(z .gt. z_min) go to 10
              z      = z + z_incr
            end do

            write(6,*) 'Error:  could not find model level above ',z_min
            call abort

   10       continue
            lapse = -dtdz
            k     = kk
!
! Define Tbot & Pbot
!
            tbot  = t      (i,j,k)
            pbot  = press_m(i,j,k)
            tmp   = lapse*(rdair/gravit)*(ps_old(i,j)/pbot - 1.)
            tsurf = tbot*(1. + tmp)
            t0    = tsurf + lapse*phis_old(i,j)/gravit
!
! Based on heuristic equations:
!
            if     (t0 .gt. t_ref1 .and. tsurf .le. t_ref1) then
              lapse = (t_ref1 - tsurf)*gravit/phis_old(i,j)
            elseif (t0 .gt. t_ref1 .and. tsurf .gt. t_ref1) then
              lapse = 0.
              tsurf = (t_ref1 + tsurf)*0.5
            endif

            if(tsurf .lt. t_ref2) then
              lapse = -dtdz
              tsurf = (t_ref2 + tsurf)*0.5
            end if              

            x   = lapse*del_phis/(gravit*tsurf)
            tmp = 1. - x/2. + x**2./3.
            tmp = del_phis/(rdair*tsurf)*tmp
            ps_new(i,j) = ps_old(i,j)*exp(tmp)

          endif
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine psadj_se(plev    ,plevp1  ,pcol    ,t       , &
                       press_m ,press_i ,phis_old,phis_new,ps_old  , &
                       ps_new  )
!
!-----------------------------------------------------------------------
!
! Adjust Ps based on difference between "analysis" phis and model phis.
! Also uses T and P arrays
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plev      ! vertical dimension
      integer plevp1    ! "plev+1"
      integer pcol      
!
      real*8 t       (pcol,plev)    ! analysis Temperatures
      real*8 press_m (pcol,plev)    ! analysis pressures
      real*8 press_i (pcol,plevp1)  ! analysis pressures (interfaces)
      real*8 phis_old(pcol)         ! analysis phis
      real*8 phis_new(pcol)         ! model phis
      real*8 ps_old  (pcol)         ! analysis Ps (horizontally
!                                   ! interpolated to model grid)
!
!     OUTPUTS
!
      real*8 ps_new  (pcol)         ! adjusted model Ps
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 tsurf
      real*8 t_ref1
      real*8 t_ref2
      real*8 t0
      real*8 tbot
      real*8 pbot
      real*8 dtdz
      real*8 lapse
      real*8 boltz
      real*8 avogad
      real*8 mwdair
      real*8 rgas
      real*8 rdair
      real*8 gravit
      real*8 del_phis
      real*8 x
      real*8 threshold
      real*8 tmp
      real*8 z
      real*8 z_min
      real*8 z_incr
      real*8 hkk
      integer n, k, kk         ! Indices
!
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2
      boltz   = 1.38065e-23       ! boltzmann's constant ~ J/k/molecule
      avogad  = 6.02214e26        ! avogadro's number ~ molecules/kmole
      mwdair  = 28.966            ! molecular weight dry air ~ kg/kmole

      rgas    = avogad*boltz      ! universal gas constant ~ J/k/kmole
      rdair   = rgas/mwdair       ! constant for dry air   ~ J/k/kg

      t_ref1    = 290.5
      t_ref2    = 255.0
      threshold = 0.001
      
      do n = 1,pcol

          del_phis = phis_old(n) - phis_new(n)
!
! If difference between analysis and model phis is negligible,
! then set model Ps = analysis
!
          if(abs(del_phis) .le. threshold) then
            ps_new(n) = ps_old(n)
!
! Else, go nuts...
!
          else
!
! Tbot and Pbot are determined from the first model level that is at
! least 150m above the surface
!
            z_min = 150.
            z     = 0.

            do k = plev,1,-1
              kk     = k
              hkk    = 0.5*( press_i(n,k+1) - press_i(n,k))/press_m(n,k)
              z_incr = (rdair/gravit)*t(n,k)*hkk
              z      = z + z_incr
              if(z .gt. z_min) go to 10
              z      = z + z_incr
            end do

            write(6,*) 'Error:  could not find model level above ',z_min
            call abort

   10       continue
            lapse = -dtdz
            k     = kk
!
! Define Tbot & Pbot
!
            tbot  = t      (n,k)
            pbot  = press_m(n,k)
            tmp   = lapse*(rdair/gravit)*(ps_old(n)/pbot - 1.)
            tsurf = tbot*(1. + tmp)
            t0    = tsurf + lapse*phis_old(n)/gravit
!
! Based on heuristic equations:
!
            if     (t0 .gt. t_ref1 .and. tsurf .le. t_ref1) then
              lapse = (t_ref1 - tsurf)*gravit/phis_old(n)
            elseif (t0 .gt. t_ref1 .and. tsurf .gt. t_ref1) then
              lapse = 0.
              tsurf = (t_ref1 + tsurf)*0.5
            endif

            if(tsurf .lt. t_ref2) then
              lapse = -dtdz
              tsurf = (t_ref2 + tsurf)*0.5
            end if              

            x   = lapse*del_phis/(gravit*tsurf)
            tmp = 1. - x/2. + x**2./3.
            tmp = del_phis/(rdair*tsurf)*tmp
            ps_new(n) = ps_old(n)*exp(tmp)

          endif
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! NCLFORTSTART
      subroutine tsadj(plat    ,plon    ,phis_old,phis_new,ts      )
!
!-----------------------------------------------------------------------
!
! Adjust Ts based on difference between old and new phis.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plat !  latitude dimension
      integer plon !  longitude dimension
!
      real*8 phis_old(plon,plat) ! analysis phis (e.g., ECMWF)
      real*8 phis_new(plon,plat) ! model phis
!
!     INPUT/OUTPUT
!
      real*8 ts      (plon,plat) ! Surface Temp
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 dtdz
      real*8 gravit
      real*8 del_z
      integer i, j, k             ! Indices
!
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2

      do j = 1,plat
        do i = 1,plon

          del_z = ( phis_new(i,j) - phis_old(i,j) )/gravit
          ts(i,j) = ts(i,j) + dtdz*del_z

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine tsadj_se(pcol, phis_old, phis_new, ts)
!
!-----------------------------------------------------------------------
!
! Adjust Ts based on difference between old and new phis.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer pcol
!
      real*8 phis_old(pcol) ! analysis phis (e.g., ECMWF)
      real*8 phis_new(pcol) ! model phis
!
!     INPUT/OUTPUT
!
      real*8 ts      (pcol) ! Surface Temp
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 dtdz
      real*8 gravit
      real*8 del_z
      integer n, k             ! Indices
!
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2

      do n = 1,pcol

          del_z = ( phis_new(n) - phis_old(n) )/gravit
          ts(n) = ts(n) + dtdz*del_z

      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
