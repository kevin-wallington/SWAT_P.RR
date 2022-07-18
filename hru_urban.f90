      subroutine hru_urban
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine computes loadings from urban areas using the
!!    USGS regression equations or a build-up/wash-off algorithm

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dirtmx(:)    |kg/curb km     |maximum amount of solids allowed to
!!                                 |build up on impervious surfaces
!!    hru_km(:)    |km^2           |area of HRU in square kilometers
!!    ihru         |none           |HRU number
!!    nsweep(:)    |none           |sequence number of street sweeping operation
!!                                 |within the year
!!    peakr        |m^3/s          |peak runoff rate
!!    sedyld(:)    |metric tons    |daily soil loss caused by water erosio
!!    surfq(:)     |mm H2O         |surface runoff for the day in HRU
!!    tconc(:)     |hr             |time of concentration
!!    thalf(:)     |days           |time for the amount of solids on
!!                                 |impervious areas to build up to 1/2
!!                                 |the maximum level
!!    tnconc(:)    |mg N/kg sed    |concentration of total nitrogen in
!!                                 |suspended solid load from impervious
!!                                 |areas
!!    tno3conc(:)  |mg NO3-N/kg sed|concentration of NO3-N in suspended
!!                                 |solid load from impervious areas
!!    tpconc(:)    |mg P/kg sed    |concentration of total phosphorus in
!!                                 |suspended solid load from impervious
!!                                 |areas
!!    twash(:)     |days           |time that solids have built-up on streets
!!    urbcoef(:)   |1/mm           |wash-off coefficient for removal of
!!                                 |constituents from an impervious surface
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    twash(:)    |
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max, Log
!!    SWAT: Regres, sweep

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use hru_module, only : hru, ihru, sedyld, surfq, ulu, sanyld, silyld, clayld, sagyld, lagyld, sedorgn, sedorgp,  &
        surqno3, surqsolp, twash, tconc, peakr, precipday, sedminpa, sedminps ! sedminpa and sedminps added by KDW 08/10
      use urban_data_module
      use hydrograph_module
      use climate_module
      
      implicit none

      real :: cod          !kg            |carbonaceous biological oxygen demand of 
                           !              |surface runoff from urban area
      real :: sus_sol      !kg            |suspended solid loading in surface runoff
                           !              |from urban area
      real :: tn           !kg            |total nitrogen in surface runoff from
                           !              |urban area
      real :: tp           !kg            |total phosphorus in surface runoff from 
                           !              |urban area
      real :: urbk         !1/hr          |
      real :: turo         !              |
      real :: dirto        !kg/ha         |amount of solids built up on impervious
                           !              |surfaces at beginning of day
      real :: durf         !hr            |duration of rainfall
      real :: rp1          !none          |variable to hold intermediate calculation
                           !              |value 
      real :: dirt         !kg/curb km    |amount of solids built up on impervious
                           !              |surfaces 
      integer :: j         !none          |HRU number
      integer :: iob       !              |
      real :: regres       !              |
      real :: xx           !              |
      real :: exp          !              |
      real :: tno3         !              |

      j = ihru
      ulu = hru(j)%luse%urb_lu
      iob = hru(j)%obj_no
      iwst = ob(iob)%wst

      select case (hru(j)%luse%urb_ro)

        case ("usgs_reg")                         !! USGS regression equations
        if (precipday > .1 .and. surfq(j) > .1) then
          cod = 0.
          sus_sol = 0.
          tn = 0.
          tp = 0.
          cod = Regres(1)
          sus_sol = Regres(2)
          tn = Regres(3)
          tp = Regres(4)

      !     sedyld(j) = (.001 * sus_sol) * urbdb(ulu)%fimp + sedyld(j)    & 
      !                                        * (1. - urbdb(ulu)%fimp)

      ! ! !! The sediment loading from urban imprevious area is assumed 
	    ! ! !! to be all sitly particles
      !     silyld(j) = (.001 * sus_sol) * urbdb(ulu)%fimp                &
      !                             + silyld(j) * (1. - urbdb(ulu)%fimp)
      !     sanyld(j) = sanyld(j) * (1. - urbdb(ulu)%fimp)
      !     clayld(j) = clayld(j) * (1. - urbdb(ulu)%fimp)
      !     sagyld(j) = sagyld(j) * (1. - urbdb(ulu)%fimp)
      !     lagyld(j) = lagyld(j) * (1. - urbdb(ulu)%fimp)

      !     sedorgn(j)=(.7 * tn / (hru(j)%km *100.))*urbdb(ulu)%fimp         &
      !                          + sedorgn(j) * (1. - urbdb(ulu)%fimp)
      !     surqno3(j)=(.3 * tn / (hru(j)%km *100.))*urbdb(ulu)%fimp         &
      !                          + surqno3(j) * (1. - urbdb(ulu)%fimp)
      !     sedorgp(j)=(.25 * tp / (hru(j)%km*100.))*urbdb(ulu)%fimp         & ! changed from .75 to .25 so mass balance correct with addition below, KDW, 08/10
      !                         +  sedorgp(j) * (1. - urbdb(ulu)%fimp)
      !     surqsolp(j)=(.25 * tp / (hru(j)%km *100.))*urbdb(ulu)%fimp         &
      !                         + surqsolp(j) * (1. - urbdb(ulu)%fimp)
      !     sedminpa(j) =(.25 * tp / (hru(j)%km*100.))*urbdb(ulu)%fimp         & ! KDW, 08/10
      !                         +  sedminpa(j) * (1. - urbdb(ulu)%fimp)
      !     sedminps(j) =(.25 * tp / (hru(j)%km*100.))*urbdb(ulu)%fimp         & ! KDW, 08/10
      !                         +  sedminps(j) * (1. - urbdb(ulu)%fimp)
      !! All above yields replaced below by KDW, 11/22/21, 01/14/21
      sedyld(j) = (.001 * sus_sol) * urbdb(ulu)%fimp + sedyld(j)
      silyld(j) = (.001 * sus_sol) * urbdb(ulu)%fimp + silyld(j)

      sedorgn(j)=(.7 * tn / (hru(j)%km *100.))*urbdb(ulu)%fimp + sedorgn(j)
      surqno3(j)=(.3 * tn / (hru(j)%km *100.))*urbdb(ulu)%fimp + surqno3(j)
      sedorgp(j)=(.25 * tp / (hru(j)%km*100.))*urbdb(ulu)%fimp +  sedorgp(j)        ! changed from .75 to .25 so mass balance correct with addition below, KDW, 08/10
      surqsolp(j)=(.25 * tp / (hru(j)%km *100.))*urbdb(ulu)%fimp + surqsolp(j)
      sedminpa(j) =(.25 * tp / (hru(j)%km*100.))*urbdb(ulu)%fimp +  sedminpa(j)
      sedminps(j) =(.25 * tp / (hru(j)%km*100.))*urbdb(ulu)%fimp +  sedminps(j)
      !! finish new yields from KDW
        endif

        case ("buildup_washoff")                         !! build-up/wash-off algorithm

        if (surfq(j) > 0.1) then
          !! rainy day: no build-up, street cleaning allowed

          !! calculate amount of dirt on streets prior to wash-off
          dirt = 0.
          dirto = 0.
          dirt = urbdb(ulu)%dirtmx * twash(j) /                         &                       
                                   (urbdb(ulu)%thalf + twash(j))
          dirto = dirt

          !! calculate wash-off of solids
          urbk = 0.
          urbk = urbdb(ulu)%urbcoef * (peakr * 3.6 / hru(j)%km)
                                     !! expression in () peakr in mm/hr
          rp1 = 0.
          durf = 0.
          turo = 0.
          rp1 = -2. * Log(1.- wst(iwst)%weat%precip_half_hr)
          durf = 4.605 / rp1         
          turo = durf + tconc(j)
          if (turo > 24.) turo = 24.
          xx = urbk * turo
          if (xx > 24.) xx = 24.
          dirt = dirt * Exp (-xx)
          if (dirt < 1.e-6) dirt = 0.0

          !! set time to correspond to lower amount of dirt
          twash(j) = 0.
          twash(j) = urbdb(ulu)%thalf * dirt / (urbdb(ulu)%dirtmx -dirt)

          sus_sol = 0.
          tn = 0.
          tp = 0.
          tno3 = 0.
          !! amounts are kg/ha
          sus_sol = Max(0., (dirto - dirt) * urbdb(ulu)%curbden)
          tn = urbdb(ulu)%tnconc * sus_sol / 1.e6
          tp = urbdb(ulu)%tpconc * sus_sol / 1.e6
          tno3 = urbdb(ulu)%tno3conc * sus_sol / 1.e6
      !! All yields replaced below by KDW, 11/22/21
      !     sedyld(j) = (.001 * sus_sol*hru(j)%area_ha)*urbdb(ulu)%fimp +            &
      !                              sedyld(j) * (1. - urbdb(ulu)%fimp)

      !     !! The sediment loading from urban imprevious area is assumed 
	    ! !! to be all sitly particles
      !     silyld(j) = (.001 * sus_sol * hru(j)%area_ha) *                         &
      !                urbdb(ulu)%fimp + silyld(j) * (1. -                     &
      !                                                urbdb(ulu)%fimp)
      !     sanyld(j) = sanyld(j) * (1. - urbdb(ulu)%fimp)
      !     clayld(j) = clayld(j) * (1. - urbdb(ulu)%fimp)
      !     sagyld(j) = sagyld(j) * (1. - urbdb(ulu)%fimp)
      !     lagyld(j) = lagyld(j) * (1. - urbdb(ulu)%fimp)

      !     surqno3(j) = tno3 * urbdb(ulu)%fimp + surqno3(j) *           &
      !                                        (1. - urbdb(ulu)%fimp)
      !     sedorgn(j) = (tn - tno3) * urbdb(ulu)%fimp + sedorgn(j) *    &
      !                                        (1. - urbdb(ulu)%fimp)
      !     sedorgp(j) = .25 * tp * urbdb(ulu)%fimp + sedorgp(j)         & ! changed from .75 to .25 so mass balance correct with addition below, KDW, 08/10
      !                                      * (1. - urbdb(ulu)%fimp)
      !     surqsolp(j) = .25 * tp * urbdb(ulu)%fimp + surqsolp(j) *     &
      !                                        (1. - urbdb(ulu)%fimp)
      !     sedminpa(j) = .25 * tp * urbdb(ulu)%fimp + sedminpa(j) *     & ! KDW, 08/10
      !                                        (1. - urbdb(ulu)%fimp)
      !     sedminps(j) = .25 * tp * urbdb(ulu)%fimp + sedminps(j) *     & ! KDW, 08/10
      !                                        (1. - urbdb(ulu)%fimp)             
          sedyld(j) = (.001 * sus_sol * hru(j)%area_ha) * urbdb(ulu)%fimp + sedyld(j) ! KDW 01/14/21, added area, had forgetten
          silyld(j) = (.001 * sus_sol * hru(j)%area_ha) * urbdb(ulu)%fimp + silyld(j) ! KDW 01/14/21, added area, had forgetten
    
          sedorgn(j)= .7 * tn * urbdb(ulu)%fimp + sedorgn(j) ! KDW 01/14/21, deleted division by area, my earlier mistake, here and all below
          surqno3(j)= .3 * tn * urbdb(ulu)%fimp + surqno3(j)
          sedorgp(j)= .25 * tp * urbdb(ulu)%fimp +  sedorgp(j)        ! changed from .75 to .25 so mass balance correct with addition below, KDW, 08/10
          surqsolp(j)= .25 * tp * urbdb(ulu)%fimp + surqsolp(j)
          sedminpa(j) = .25 * tp * urbdb(ulu)%fimp +  sedminpa(j)
          sedminps(j) = .25 * tp * urbdb(ulu)%fimp +  sedminps(j)
          !! finish new yields from KDW                                
        else
          !! dry day
          twash(j) = twash(j) + 1.

        end if

      end select

      return
      end subroutine hru_urban