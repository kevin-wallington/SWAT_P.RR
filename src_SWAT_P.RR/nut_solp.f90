      subroutine nut_solp
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates the movement of soluble phosphorus in the soil
!!    profile and theamount of phosphorus lost from the soil in surface runoff,
!!    lateral flow, and tile flow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihru          |none         |HRU number
!!    pperco        |none         |phosphorus percolation coefficient (0-1)
!!                                |0:concentration of phosphorus in surface runoff
!!                                |  is zero
!!                                |1:surface runoff has same concentration of
!!                                |  phosphorus as percolate      
!!    surfq(:)      |mm H2O       |surface runoff generated on day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min, Max

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use basin_module
      use organic_mineral_mass_module
      use hru_module, only : hru, surqsolp, percp, tilesolp, latsolp, surfq, ihru, qtile, i_sep, surfq_pervious
      use soil_module
      use hydrograph_module, only : ht1 ! - KDW
      
      implicit none 

      integer :: j         !none          |HRU number
      integer :: jj        !none          |counter   
      real :: ssfplyr      !kg P/ha       |phos transported in lateral flow from layer
      real :: percplyr     !kg P/ha       |phos leached to next lower layer with
                           !              |percolation
      real :: cosurf       !kg N/mm       |concentration of phos in surface runoff 
      real :: colatperc    !kg N/mm       |concentration of phos in lateral flow and percolation

      real :: xx           !varies        |variable to hold intermediate calculation
      real :: sro          !mm H2O        |surface runoff 
      real :: tileqout     !mm H2O        |total percolation out of and through layer induced by tile flow, KDW 06/10 
                           !              |result

      ssfplyr = 0.; percplyr = 0. ; cosurf = 0. ; colatperc = 0. ; xx = 0.; sro = 0. ! KDW 04/07
      tileqout = 0. ! KDW 06/10
      j = ihru

      percplyr = ht1%solp ! changed from 0, to add solp from runon - KDW
      latsolp(j) = 0

  if (bsn_cc%sol_P_model == 0 .or. bsn_cc%sol_P_model == 1) then  
            soil1(j)%mp(1)%lab = soil1(j)%mp(1)%lab + percplyr  !! add soluble P from runon, KDW 08/13
    !! compute soluble P lost in surface runoff
            xx = soil(j)%phys(1)%bd * soil(j)%phys(1)%d * bsn_prm%phoskd
            if (hru(j)%luse%urb_lu > 0) then ! KDW 11/22/21, corrected accounting for urban HRU, here rather than in hru_urban subroutine
              surqsolp(j) = soil1(j)%mp(1)%lab  * surfq_pervious(j) / (xx + 1.) 
            else
              surqsolp(j) = soil1(j)%mp(1)%lab  * surfq(j) / (xx + 1.)   !dont merge; KDW replaced rsd1(j)%mp%lab with soil1(j)%mp(1)%lab here and below in this if section, 08/13
                !!units ==> surqsolp = [kg/ha * mm] / [t/m^3 * mm * m^3/t] = kg/ha
            endif
            surqsolp(j) = Min(surqsolp(j), soil1(j)%mp(1)%lab)
            surqsolp(j) = Max(surqsolp(j), 0.)
            soil1(j)%mp(1)%lab = soil1(j)%mp(1)%lab - surqsolp(j)


      !! compute soluble P leaching
            percplyr = soil1(j)%mp(1)%lab * soil(j)%ly(1)%prk /                   &
              ((soil(j)%phys(1)%conv_wt/ 1000.) * bsn_prm%pperco + .1)   !dont merge
              percplyr = Min(percplyr, .5 * soil1(j)%mp(1)%lab)
              soil1(j)%mp(1)%lab = soil1(j)%mp(1)%lab - percplyr
            if (soil(j)%nly >= 2) then
              soil1(j)%mp(2)%lab = soil1(j)%mp(2)%lab + percplyr
            end if
        
            do jj = 2, soil(j)%nly-1
              percplyr = 0.
              if (jj/=i_sep(j)) then
                percplyr = soil1(j)%mp(jj)%lab * soil(j)%ly(jj)%prk /                &
                  ((soil(j)%phys(jj)%conv_wt / 1000.) * bsn_prm%pperco + .1)
                  percplyr = Min(percplyr, .2 * soil1(j)%mp(jj)%lab)
                soil1(j)%mp(jj)%lab = soil1(j)%mp(jj)%lab - percplyr
                soil1(j)%mp(jj+1)%lab = soil1(j)%mp(jj+1)%lab + percplyr
                if (jj == hru(j)%lumv%ldrain) then
                  percplyr = soil1(j)%mp(jj)%lab * qtile /                        &
                    (soil(j)%phys(jj)%conv_wt / 1000. * bsn_prm%pperco + .1)  !dont merge
                    percplyr = Min(percplyr, soil1(j)%mp(jj)%lab)
                  soil1(j)%mp(jj)%lab = soil1(j)%mp(jj)%lab - percplyr
                endif
              endif
            end do
        percp(j) = percplyr
        tilesolp(j) = 0. ! KDW 08/13
        latsolp(j) = 0.  ! KDW 08/13

  else ! SWAT+P formulation by KDW, soluble P leaching throughout the soil profile and transport with lateral and tile flow.
      do jj = 1, soil(j)%nly
        !! add phos leached from layer above
        soil1(j)%mp(jj)%lab = soil1(j)%mp(jj)%lab + percplyr
	      if (soil1(j)%mp(jj)%lab < 1.e-6) soil1(j)%mp(jj)%lab = 0.0
        if (jj == 1) then ! this if/else statement copied from "nut_nlch" to include surface runoff when calculating concentration - KDW
          if (hru(j)%luse%urb_lu > 0) then ! KDW 11/22/21, corrected accounting for urban HRU, here rather than in hru_urban subroutine
            sro = surfq_pervious(j) ! units: mm
          else
            sro = surfq(j) ! units: mm
          endif
        else
          sro = 0.
        endif
        !! determine concentration of phos in mobile water
        xx = (soil(j)%phys(jj)%st + soil(j)%phys(jj)%wpmm + sro) ! sro added by KDW, changed wp to wpmm on 11/12/21
        cosurf =  (soil1(j)%mp(1)%lab / bsn_prm%phoskd)  / (xx + 0.0001)  !new formulation - KDW
        ! in the above, phoskd has been repurposed, original range and default value not applicable
        cosurf = Max(cosurf, 0.)     !kg/ha/mm (if * 100 = ppm)  
        colatperc = (soil1(j)%mp(jj)%lab * bsn_prm%pperco)  / (xx + 0.0001)
        colatperc = Max(colatperc, 0.)     !kg/ha/mm (if * 100 = ppm)

        !! calculate phos in surface runoff
        if (jj == 1) then
          if (hru(j)%luse%urb_lu > 0) then ! KDW 11/22/21, corrected accounting for urban HRU, here rather than in hru_urban subroutine
            surqsolp(j) = surfq_pervious(j) * cosurf ! units: mm
          else
            surqsolp(j) = surfq(j) * cosurf
          endif
          surqsolp(j) = Min(surqsolp(j), soil1(j)%mp(jj)%lab)
          soil1(j)%mp(jj)%lab = soil1(j)%mp(jj)%lab - surqsolp(j)
        endif
        !! calculate phos in tile flow
        if (hru(j)%lumv%ldrain == jj) then
          tilesolp(j) = colatperc * qtile
          tilesolp(j) = Min(tilesolp(j), soil1(j)%mp(jj)%lab)
          soil1(j)%mp(jj)%lab = soil1(j)%mp(jj)%lab - tilesolp(j)
          tileqout = 0 ! KDW 06/10
        else 
          if (hru(j)%lumv%ldrain > jj) then !that is, layer above tile drain - KDW 06/10
            tileqout = tileqout + soil(j)%ly(jj)%qtilelyr ! tile induced flow passing through (from layers above) or leaving layer - KDW 06/10
          else
            tileqout = 0 ! KDW 06/10
          end if
        end if               

        !! calculate phos in lateral flow
        ssfplyr = colatperc * soil(j)%ly(jj)%flat
        ssfplyr = Min(ssfplyr, soil1(j)%mp(jj)%lab)
        latsolp(j) = latsolp(j) + ssfplyr
        soil1(j)%mp(jj)%lab = soil1(j)%mp(jj)%lab - ssfplyr

        !! calculate phos in percolate
        percplyr = colatperc * (soil(j)%ly(jj)%prk + tileqout) ! added tileqout for downward movement of P induced by tile flow KDW, 06/10
        percplyr = Min(percplyr, soil1(j)%mp(jj)%lab)
        soil1(j)%mp(jj)%lab = soil1(j)%mp(jj)%lab - percplyr

        if (soil1(j)%mp(jj)%lab < 1.e-6) soil1(j)%mp(jj)%lab = 0.0 ! KDW, 03/24
      end do
      if (surqsolp(j) < 1.e-6) surqsolp(j) = 0.0 ! KDW, 03/26
      if (tilesolp(j) < 1.e-6) tilesolp(j) = 0.0 ! KDW, 03/26
      if (latsolp(j) < 1.e-6) latsolp(j) = 0.0 ! KDW, 03/26

      !! calculate phos leaching from soil profile
      percp(j) = percplyr
  endif

      return
      end subroutine nut_solp