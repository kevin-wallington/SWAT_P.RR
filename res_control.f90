      subroutine res_control (jres)
      
      use basin_module
      use reservoir_data_module 
      use time_module
      use reservoir_module
      use climate_module
      use hydrograph_module
      use conditional_module
      use water_body_module
      
      implicit none

      integer :: ii                   !none          |counter 
      integer, intent(in) :: jres                 !none          
      integer :: idat                 !              |
      integer :: ihyd                 !none          |counter
      integer :: ised                 !none          |counter
      integer :: irel                 !              |
      integer :: inut                 !none          |counter
      integer :: iob                  !none          |counter
      integer :: irch                  !none          |counter
      real :: pvol_m3
      real :: evol_m3
      integer :: ch_par = 0
      real :: mxsa                    !m^2              ! maximum surface area, i.e. esa - KDW
      integer :: iswet                !binary         |is this a wetland (1) or a reservoir (0), KDW 03/30
      
      iswet = 0 ! KDW 03/30

      iob = res_ob(jres)%ob
      mxsa = res_ob(jres)%esa * 10000 ! - KDW
      
      !! set water body pointer to res
      wbody => res(jres)
      wbody_wb => res_wat_d(jres)
      wbody_bed => res_bed(jres)
      
      ht1 = ob(icmd)%hin    !! set incoming flow
      ht2 = resz            !! zero outgoing flow

      if (time%yrs > pco%nyskip) then ! moved up to precede trapping, KDW 04/08/22
        res_in_d(jres) = ht1 
      end if    

      !! add incoming flow to reservoir
      if (bsn_cc%wq == 1 .or. bsn_cc%wq == 0) then
          res(jres) = res(jres) + ht1
      else
          res(jres)%flo = res(jres)%flo + ht1%flo !made to add "flo" only instead of all of ht1 because sediment and nutrients are added within sub-routines - KDW
      endif

      if (time%yrc > res_hyd(jres)%iyres .or. (time%mo >= res_hyd(jres)%mores   &
                                   .and. time%yrc == res_hyd(jres)%iyres)) then
        !! perform reservoir water/sediment balance
        idat = res_ob(jres)%props
        ihyd = res_dat(idat)%hyd
        ised = res_dat(idat)%sed
        if(time%step == 0) then
          !! determine reservoir outflow
          irel = res_dat(idat)%release
          d_tbl => dtbl_res(irel)
          pvol_m3 = res_ob(jres)%pvol
          evol_m3 = res_ob(jres)%evol
          call conditions (jres, irel)
          call res_hydro (jres, irel, ihyd, pvol_m3, evol_m3)
          if (bsn_cc%wq == 1 .or. bsn_cc%wq == 0) then
              call res_sediment (jres, ihyd, ised)
              !! below lines moved up from lines ~125 - KDW
              inut = res_dat(jres)%nut 
              ch_par = res_nut(inut)%ch_params !this indicates which *channel* nutcha to use - KDW
              call res_nutrient (jres, inut, iob) !KDW 08/18
              call res_pest (jres)
              !! above lines moved up from lines ~125 - KDW
          endif
          if (bsn_cc%wq == 2) then
            call res_trap (jres, ihyd, ised, mxsa)
            !! below lines moved up from lines ~125 - KDW
            inut = res_dat(jres)%nut 
            ch_par = res_nut(inut)%ch_params !this indicates which *channel* nutcha to use - KDW
            call res_nutrient2 (jres, ch_par, iob, mxsa, iswet) !KDW 03/30
            call res_pest (jres)
            !! above lines moved up from lines ~125 - KDW
            call res_sediment2 (jres, ihyd, ised, mxsa)
          endif
	      else
	      !call res_hourly
        endif

        
        !! calculate water balance for day
        iwst = ob(iob)%wst
        res_wat_d(jres)%evap = 10. * res_hyd(ihyd)%evrsv * wst(iwst)%weat%pet * res_wat_d(jres)%area_ha
        res_wat_d(jres)%seep = 240. * res_hyd(ihyd)%k * res_wat_d(jres)%area_ha
        res_wat_d(jres)%precip = 10. * wst(iwst)%weat%precip * res_wat_d(jres)%area_ha

        !! add precip to reservoir storage
        res(jres)%flo = res(jres)%flo + res_wat_d(jres)%precip

        !! subtract outflow from reservoir storage
        !res(jres)%flo = res(jres)%flo - ht2%flo
        !subtract all nutrients, etc. not just flo, replaced above, KDW 10/15/21
        res(jres)%flo = res(jres)%flo - ht2%flo
        res(jres)%sed = res(jres)%sed - ht2%sed        
        res(jres)%orgn = res(jres)%orgn - ht2%orgn        
        res(jres)%orgp = res(jres)%orgp - ht2%orgp
        res(jres)%minpa = res(jres)%minpa - ht2%minpa
        res(jres)%minps = res(jres)%minps - ht2%minps
        res(jres)%no3 = res(jres)%no3 - ht2%no3
        res(jres)%solp = res(jres)%solp - ht2%solp
        res(jres)%chla = res(jres)%chla - ht2%chla
        res(jres)%nh3 = res(jres)%nh3 - ht2%nh3
        res(jres)%no2 = res(jres)%no2 - ht2%no2
        res(jres)%cbod = res(jres)%cbod - ht2%cbod
        res(jres)%dox = res(jres)%dox - ht2%dox
        res(jres)%san = res(jres)%san - ht2%san
        res(jres)%sil = res(jres)%sil - ht2%sil
        res(jres)%cla = res(jres)%cla - ht2%cla
        res(jres)%sag = res(jres)%sag - ht2%sag
        res(jres)%lag = res(jres)%lag - ht2%lag
        res(jres)%grv = res(jres)%grv - ht2%grv
        
        if (res(jres)%flo < 0.) then
          ht2%flo = ht2%flo + res(jres)%flo
          res(jres)%flo = 0.
        end if

        !! Below is reservoir withdrawal (m^3), specific to Lake Decatur - KDW
        !! Need to remove nutrients/sediment with this withdrawal, based on percent removed? - KDW
          if (time%mo == 1) then
            res(jres)%flo = res(jres)%flo - 126971
          elseif (time%mo == 2) then
            res(jres)%flo = res(jres)%flo - 133923
          elseif (time%mo == 3) then
            res(jres)%flo = res(jres)%flo - 141113
          elseif (time%mo == 4) then
            res(jres)%flo = res(jres)%flo - 151718
          elseif (time%mo == 5) then
            res(jres)%flo = res(jres)%flo - 153763
          elseif (time%mo == 6) then
            res(jres)%flo = res(jres)%flo - 140319
          elseif (time%mo == 7) then
            res(jres)%flo = res(jres)%flo - 120973
          elseif (time%mo == 8) then
            res(jres)%flo = res(jres)%flo - 114286
          elseif (time%mo == 9) then
            res(jres)%flo = res(jres)%flo - 111097
          elseif (time%mo == 10) then
            res(jres)%flo = res(jres)%flo - 111802
          elseif (time%mo == 11) then
            res(jres)%flo = res(jres)%flo - 120630
          elseif (time%mo == 12) then
            res(jres)%flo = res(jres)%flo - 128882
          end if
        !! Above is reservoir withdrawal (m^3), specific to Lake Decatur - KDW

          if (res(jres)%flo < 0.) then ! KDW 04/25/22
            res(jres)%flo = 0.
          end if

        !! subtract evaporation from reservoir storage
        res(jres)%flo = res(jres)%flo - res_wat_d(jres)%evap
        if (res(jres)%flo < 0.) then
          res_wat_d(jres)%evap = res_wat_d(jres)%evap + res(jres)%flo
          res(jres)%flo = 0.
        end if
      
        !! subtract seepage from reservoir storage
        res(jres)%flo = res(jres)%flo - res_wat_d(jres)%seep
        if (res(jres)%flo < 0.) then
          res_wat_d(jres)%seep = res_wat_d(jres)%seep + res(jres)%flo
          res(jres)%flo = 0.
        end if

        !! update surface area
        if (res(jres)%flo > 0.) then
          res_wat_d(jres)%area_ha = res_ob(jres)%br1 * res(jres)%flo ** res_ob(jres)%br2
        else
          res_wat_d(jres)%area_ha = 0.
        end if

        !! subtract sediment leaving from reservoir
        !res(jres)%sed = res(jres)%sed - ht2%sed
        !res(jres)%sil = res(jres)%sil - ht2%sil
        !res(jres)%cla = res(jres)%cla - ht2%cla

        !! set values for outflow variables
        ob(icmd)%hd(1) = ht2

        if (time%step > 0) then
          do ii = 1, time%step
            ob(icmd)%ts(1,ii) = ht2 / real(time%step)
          end do
        end if

        !! set inflow and outflow variables for reservoir_output
        if (time%yrs > pco%nyskip) then
          !res_in_d(jres) = ht1  !inflow variables moved up before trapping, KDW 04/08/22
          res_out_d(jres) = ht2
          !res_in_d(jres)%flo = res_in_d(jres)%flo / 10000.          !m^3 -> ha-m
          !res_out_d(jres)%flo = res_out_d(jres)%flo / 10000.        !m^3 -> ha-m
          !res_wat_d(jres)%evap = res_wat_d(jres)%evap / 10000.      !m^3 -> ha-m
          !res_wat_d(jres)%seep = res_wat_d(jres)%seep / 10000.      !m^3 -> ha-m
          !res_wat_d(jres)%precip = res_wat_d(jres)%precip / 10000.  !m^3 -> ha-m
        end if             
        
      else
        !! reservoir has not been constructed yet
        ob(icmd)%hd(1) = ob(icmd)%hin
      end if

      return
      end subroutine res_control