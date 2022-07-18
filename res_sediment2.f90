      subroutine res_sediment2 (jres, ihyd, ised, mxsa)

      use reservoir_data_module
      use reservoir_module
      use conditional_module
      use climate_module
      use time_module
      use hydrograph_module
      use water_body_module
      
      implicit none

      integer, intent (in) :: jres          !none          |reservoir number
      integer, intent (in) :: ihyd          !none          |res hydrologic data pointer
      integer, intent (in) :: ised          !none          |res sediment data pointer
      real :: trapres                       !              |
      real :: susp                          !              |
      real :: velofl                        !              |  
      real :: xx, sed_ppm, sil_ppm, cla_ppm, sag_ppm, san_ppm, lag_ppm, grv_ppm
      real :: sed_settle, sil_settle, cla_settle, fines_settle, settle_remain            ! tons - KDW
      real :: san_settle, sag_settle, lag_settle, grv_settle    ! tons - KDW
      real :: fines_inflow, srfarea_inflow, srfarea_settle, srfarea_wtclm
      real :: conv_mass2vol, vol_h2o_aer, vol_h2o_anaer, vol_h2o_newdep                !m^3          |volume of pore water in each layer - KDW
      real :: aer_mass, anaer_mass     
      real :: minpa_aer, minps_aer, orgp_aer, disp_aer                !kg          |kg P by form and location - KDW
      real :: minpa_anaer, minps_anaer, orgp_anaer, disp_anaer                !kg          |kg P by form and location - KDW
      real :: minpa_wtclm, minps_wtclm, orgp_wtclm, disp_wtclm                !kg          |kg P by form and location - KDW
      real :: percsfar_settle, percmass_settle, percvol_settle, perc_settle
      real :: orgp_settle, disp_settle, minpa_settle, minps_settle      !kg          |kg P settling - KDW
      real :: new_aer_minpa, new_aer_minps, new_aer_orgp, new_aer_disp                !kg          |kg P by form and location, end of time step - KDW
      real :: new_anaer_minpa, new_anaer_minps, new_anaer_orgp, new_anaer_disp                !kg          |kg P by form and location, end of time step - KDW
      real :: new_aer_srfarea, new_anaer_srfarea
      real :: mass_wtclm, fines_wtclm
      real, intent (in) :: mxsa ! - KDW

      sil_settle = 0. ; cla_settle = 0. ; sed_settle = 0. ; fines_settle = 0.
      san_settle = 0. ; sag_settle = 0. ; lag_settle = 0. ; grv_settle = 0.
      fines_inflow = 0. ; srfarea_inflow = 0. ; srfarea_settle = 0. ; srfarea_wtclm = 0.
      conv_mass2vol = 0. ; vol_h2o_aer = 0. ; vol_h2o_anaer = 0. ; vol_h2o_newdep = 0.
      aer_mass = 0. ; anaer_mass = 0.     
      minpa_aer = 0. ; minps_aer = 0. ; orgp_aer = 0. ; disp_aer = 0.
      minpa_anaer = 0. ; minps_anaer = 0. ; orgp_anaer = 0. ; disp_anaer = 0.
      minpa_wtclm = 0. ; minps_wtclm = 0. ; orgp_wtclm = 0. ; disp_wtclm = 0.
      percsfar_settle = 0. ; percmass_settle = 0. ; percvol_settle = 0.
      orgp_settle = 0. ; disp_settle = 0. ; minpa_settle = 0. ; minps_settle = 0.
      new_aer_minpa = 0. ; new_aer_minps = 0. ; new_aer_orgp = 0. ; new_aer_disp = 0.
      new_anaer_minpa = 0. ; new_anaer_minps = 0. ; new_anaer_orgp = 0. ; new_anaer_disp = 0.
      new_aer_srfarea = 0. ; new_anaer_srfarea = 0.
      mass_wtclm = 0. ; fines_wtclm = 0. ; srfarea_wtclm = 0.
      
      ! aer_mass = .001 * 1 * (mxsa * wbody_bed%bed_bd) ! depth of aer layer is 1mm
      ! anaer_mass = .1 * (mxsa * wbody_bed%bed_bd)
      aer_mass = .002 * 1 * (mxsa * wbody_bed%bed_bd) ! depth of aer layer is 1mm! KDW debug 04/28/22, test thicker layers, replacing above
      anaer_mass = .2 * (mxsa * wbody_bed%bed_bd)! KDW debug 04/28/22, test thicker layers, replacing above
      if (wbody%flo < 1.e-6) then
        ! reservoir is empty
        san_settle = ht1%san + wbody%san
        sag_settle = ht1%sag + wbody%sag
        lag_settle = ht1%lag + wbody%lag
        grv_settle = ht1%grv + wbody%grv
        sil_settle = ht1%sil + wbody%cla
        cla_settle = ht1%cla + wbody%san
        wbody%sil = 0
        wbody%cla = 0
        wbody%sag = 0
        wbody%san = 0
        wbody%lag = 0
        wbody%grv = 0
        wbody%sed = 0
        sed_ppm = 0
        mass_wtclm = wbody%sed
        fines_wtclm = wbody%cla + wbody%sil + wbody%sag + wbody%san
        !srfarea_wtclm = wbody%cla + wbody%sil * .2 + wbody%sag * (1./15.) + wbody%san * .01
        srfarea_wtclm = 3. * ((wbody%cla/0.000002) + (wbody%sil/0.00001) + (wbody%sag/0.00003) + (wbody%san/0.0002)) ! KDW 04/14/22

      else
        !! compute reservoir water column concentrations after adding unsettled inflow 
          wbody%sil = ht1%sil + wbody%sil
          wbody%cla = ht1%cla + wbody%cla
          wbody%sag = ht1%sag + wbody%sag
          wbody%san = ht1%san + wbody%san
          wbody%lag = ht1%lag + wbody%lag
          wbody%grv = ht1%grv + wbody%grv
          if (wbody%sil < 1.e-9) wbody%sil = 0.0 ! KDW, 03/26, 03/29
          if (wbody%cla < 1.e-9) wbody%cla = 0.0 ! KDW, 03/26, 03/29
          if (wbody%sag < 1.e-9) wbody%sag = 0.0 ! KDW, 03/26, 03/29
          if (wbody%san < 1.e-9) wbody%san = 0.0 ! KDW, 03/26, 03/29
          if (wbody%lag < 1.e-9) wbody%lag = 0.0 ! KDW, 03/26, 03/29
          if (wbody%grv < 1.e-9) wbody%grv = 0.0 ! KDW, 03/26 , 03/29
          wbody%sed = wbody%cla + wbody%sil + wbody%sag + wbody%san + wbody%lag + wbody%grv
          sed_ppm = 1000000. * (wbody%sed) / wbody%flo
          mass_wtclm = wbody%sed
          fines_wtclm = wbody%cla + wbody%sil + wbody%sag + wbody%san
          !srfarea_wtclm = wbody%cla + wbody%sil * .2 + wbody%sag * (1./15.) + wbody%san * .01
          srfarea_wtclm = 3. * ((wbody%cla/0.000002) + (wbody%sil/0.00001) + (wbody%sag/0.00003) + (wbody%san/0.0002)) ! KDW 04/14/22

      !! compute settling within reservoir toward equilibrium concentration
      if (sed_ppm > res_sed(ised)%nsed) then
        xx = (sed_ppm - res_sed(ised)%nsed) * res_sed(ised)%sed_stlr + res_sed(ised)%nsed
        sed_settle = (sed_ppm - xx) * wbody%flo / 1000000
        !sed_ppm = xx !! KDW
        !wbody%sed = sed_ppm * wbody%flo / 1000000
        if (sed_settle <= wbody%grv) then
          grv_settle = sed_settle
          settle_remain = 0
        else
          grv_settle = wbody%grv
          settle_remain = sed_settle - grv_settle
          if (settle_remain <= wbody%lag) then
            lag_settle = settle_remain
            settle_remain = 0
          else
            lag_settle = wbody%lag
            settle_remain = settle_remain - lag_settle
            if (settle_remain <= wbody%san) then
              san_settle = settle_remain
              settle_remain = 0
            else
              san_settle = wbody%san
              settle_remain = settle_remain - san_settle
              if (settle_remain <= wbody%sag) then
                sag_settle = settle_remain
                settle_remain = 0
              else
                sag_settle = wbody%sag
                settle_remain = settle_remain - sag_settle
                if (settle_remain <= wbody%sil) then ! KDW 08/02, allow silt to settle before clay
                  sil_settle = settle_remain
                  settle_remain = 0
                else
                  sil_settle = wbody%sil
                  settle_remain = settle_remain - sil_settle
                  if (0.1*settle_remain <= wbody%cla) then ! KDW 08/02, allow silt to settle before clay
                    cla_settle = 0.1*settle_remain ! KDW debug/temp 04/08/22, multiple by 0.1, here and above
                    settle_remain = 0
                  else
                    cla_settle = wbody%cla
                    settle_remain = settle_remain - cla_settle
                  endif
                endif
              endif
            endif
          endif
        endif
        ! if (wbody%sil + wbody%cla < 1.e-9) then ! KDW, 03/29
        !   sil_settle = 0! KDW, 03/29
        !   cla_settle = 0! KDW, 03/29
        ! else! KDW, 03/29
        !   perc_settle = settle_remain / (wbody%sil + wbody%cla)
        !   perc_settle = Min(1.0, perc_settle)
        !   sil_settle = perc_settle * wbody%sil
        !   cla_settle = perc_settle * wbody%cla
        ! endif
      else
        !no settling
        san_settle = 0
        sag_settle = 0
        lag_settle = 0
        grv_settle = 0
        sil_settle = 0
        cla_settle = 0
      endif

      wbody%sil = wbody%sil - sil_settle
      sil_ppm = 1000000. * (wbody%sil) / wbody%flo
      wbody%cla = wbody%cla - cla_settle
      cla_ppm = 1000000. * (wbody%cla) / wbody%flo
      wbody%sag = wbody%sag - sag_settle
      sag_ppm = 1000000. * (wbody%sag) / wbody%flo
      wbody%san = wbody%san - san_settle
      san_ppm = 1000000. * (wbody%san) / wbody%flo
      wbody%lag = wbody%lag - lag_settle
      lag_ppm = 1000000. * (wbody%lag) / wbody%flo
      wbody%grv = wbody%grv - grv_settle
      grv_ppm = 1000000. * (wbody%grv) / wbody%flo
      wbody%sed = wbody%sil + wbody%cla + wbody%sag + wbody%san + wbody%lag + wbody%grv ! KDW 04/13/22
      sed_ppm = 1000000. * (wbody%sed) / wbody%flo! KDW 04/13/22
      endif
      !sed_settle = sil_settle + cla_settle + san_settle + sag_settle + lag_settle + grv_settle
      !wbody%sed = wbody%cla + wbody%sil + wbody%sag + wbody%san + wbody%lag + wbody%grv
      !srfarea_settle = cla_settle + sil_settle * .2 + sag_settle * (1./15.) + san_settle * .01
      srfarea_settle = 3. * ((cla_settle/0.000002) + (sil_settle/0.00001) + (sag_settle/0.00003) + (san_settle/0.0002)) ! KDW 04/14/22

      fines_settle = sil_settle + cla_settle + san_settle + sag_settle

!! Phosphorus movement associated with sediment deposition - KDW
  ! What are water volumes in layers (and in newly deposited sediment)
      conv_mass2vol = wbody_bed%bed_por / wbody_bed%bed_bd ! convert from mass of sediment to volume of h2o in pore space
      vol_h2o_aer = wbody_bed%bed_por * (mxsa * .001)
      vol_h2o_anaer = wbody_bed%bed_por * (mxsa * .1)
      vol_h2o_newdep = wbody_bed%bed_por * (fines_settle / wbody_bed%bed_bd)
! How much P mass of each form is in each layer and the inflow wtclm, initially
      if (wbody_bed%aer_disp < 1.e-6) wbody_bed%aer_disp = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_disp < 1.e-6) wbody_bed%anaer_disp = 0.0 ! KDW, 03/26
      if (wbody_bed%aer_orgp < 1.e-6) wbody_bed%aer_orgp = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_orgp < 1.e-6) wbody_bed%anaer_orgp = 0.0 ! KDW, 03/26
      if (wbody_bed%aer_minpa < 1.e-6) wbody_bed%aer_minpa = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_minpa < 1.e-6) wbody_bed%anaer_minpa = 0.0 ! KDW, 03/26
      if (wbody_bed%aer_minps < 1.e-6) wbody_bed%aer_minps = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_minps < 1.e-6) wbody_bed%anaer_minps = 0.0 ! KDW, 03/26
      if (wbody%orgp < 1.e-6) wbody%orgp = 0.0 ! KDW, 03/26
      if (wbody%minpa < 1.e-6) wbody%minpa = 0.0 ! KDW, 03/26
      if (wbody%minps < 1.e-6) wbody%minps = 0.0 ! KDW, 03/26
      if (wbody%solp < 1.e-6) wbody%solp = 0.0 ! KDW, 03/26
      minpa_aer = wbody_bed%aer_minpa * aer_mass / 1000 ! kg
      minps_aer = wbody_bed%aer_minps * aer_mass / 1000 ! kg
      minpa_anaer = wbody_bed%anaer_minpa * anaer_mass / 1000 ! kg
      minps_anaer = wbody_bed%anaer_minps * anaer_mass / 1000 ! kg
      minpa_wtclm = wbody%minpa ! kg
      minps_wtclm = wbody%minps ! kg
      orgp_aer = wbody_bed%aer_orgp * aer_mass / 1000 ! kg
      disp_aer =wbody_bed%aer_disp * vol_h2o_aer / 1000 ! kg 
      orgp_anaer = wbody_bed%anaer_orgp * anaer_mass / 1000 ! kg
      disp_anaer = wbody_bed%anaer_disp * vol_h2o_anaer / 1000 ! kg 
      orgp_wtclm = wbody%orgp ! kg
      disp_wtclm = wbody%solp ! kg 

! How much phosphorus is associated with the different pools of depositing sediment
      if (wbody%flo < 1.e-6) then
        percvol_settle = 1
      else
        percvol_settle = vol_h2o_newdep / wbody%flo
      endif
      if (srfarea_wtclm < 1.e-6) then      
        percmass_settle = 1
      else
        percsfar_settle = srfarea_settle / srfarea_wtclm
      endif
      if (fines_wtclm < 1.e-6) then   
        percsfar_settle = 1
      else
        percmass_settle = fines_settle / fines_wtclm
      endif
      orgp_settle = percsfar_settle * wbody%orgp ! KDW, changed from "percmass *" to "percsfar *"  02/17/22
      wbody%orgp = (1- percsfar_settle)* wbody%orgp
      disp_settle = percvol_settle * wbody%solp
      wbody%solp = (1- percvol_settle)* wbody%solp
      minpa_settle = percsfar_settle * wbody%minpa
      wbody%minpa = (1- percsfar_settle)* wbody%minpa
      minps_settle = percsfar_settle * wbody%minps
      wbody%minps = (1- percsfar_settle)* wbody%minps

! How much phosphorus is in each layer after erosion and deposition 

      if (fines_settle <= aer_mass) then! = added to inequality on 02/19/22, KDW
          new_aer_srfarea = srfarea_settle + ((aer_mass - fines_settle) / aer_mass) * wbody_bed%aer_srfarea
          new_anaer_srfarea = (fines_settle / aer_mass) * wbody_bed%aer_srfarea + &
            ((anaer_mass - fines_settle) / anaer_mass) * wbody_bed%anaer_srfarea
          new_aer_minpa = minpa_settle + ((aer_mass - fines_settle) / aer_mass) * minpa_aer
          new_aer_minps = minps_settle + ((aer_mass - fines_settle) / aer_mass) * minps_aer
          new_aer_orgp = orgp_settle + ((aer_mass - fines_settle) / aer_mass) * orgp_aer
          new_aer_disp = disp_settle + ((aer_mass - fines_settle) / aer_mass) * disp_aer
          new_anaer_minpa = (fines_settle / aer_mass) * minpa_aer + &
            ((anaer_mass - fines_settle) / anaer_mass) * minpa_anaer
          new_anaer_minps = (fines_settle / aer_mass) * minps_aer + &
            ((anaer_mass - fines_settle) / anaer_mass) * minps_anaer
          new_anaer_orgp = (fines_settle / aer_mass) * orgp_aer + &
            ((anaer_mass - fines_settle) / anaer_mass) * orgp_anaer
          new_anaer_disp = (fines_settle / aer_mass) * disp_aer + &
            ((anaer_mass - fines_settle) / anaer_mass) * disp_anaer
      else
          if (fines_settle <= anaer_mass) then! = added to inequality on 02/19/22, KDW
              new_aer_srfarea = srfarea_settle * (aer_mass / fines_settle)
              new_anaer_srfarea = ((fines_settle - aer_mass)/ fines_settle) * srfarea_settle + &
              wbody_bed%aer_srfarea + ((anaer_mass - fines_settle) / anaer_mass) * wbody_bed%anaer_srfarea
              new_aer_minpa = minpa_settle * (aer_mass / fines_settle)
              new_aer_minps = minps_settle * (aer_mass / fines_settle)
              new_aer_orgp = orgp_settle * (aer_mass / fines_settle)
              new_aer_disp = disp_settle * (aer_mass / fines_settle)
              new_anaer_minpa = ((fines_settle - aer_mass)/ fines_settle) * minpa_settle + &
                minpa_aer + ((anaer_mass - fines_settle) / anaer_mass) * minpa_anaer
              new_anaer_minps = ((fines_settle - aer_mass)/ fines_settle) * minps_settle + &
                minps_aer + ((anaer_mass - fines_settle) / anaer_mass) * minps_anaer
              new_anaer_orgp = ((fines_settle - aer_mass)/ fines_settle) * orgp_settle + &
                orgp_aer + ((anaer_mass - fines_settle) / anaer_mass) * orgp_anaer
              new_anaer_disp = ((fines_settle - aer_mass)/ fines_settle) * disp_settle + &
                disp_aer + ((anaer_mass - fines_settle) / anaer_mass) * disp_anaer
          else
              if (fines_settle <= aer_mass + anaer_mass) then! = added to inequality on 02/19/22, KDW
                  new_aer_srfarea = srfarea_settle * (aer_mass / fines_settle)
                  new_anaer_srfarea = ((fines_settle - aer_mass)/ fines_settle) * srfarea_settle + &
                    ((anaer_mass + aer_mass - fines_settle)/ aer_mass) * wbody_bed%aer_srfarea
                  new_aer_minpa = minpa_settle * (aer_mass / fines_settle)
                  new_aer_minps = minps_settle * (aer_mass / fines_settle)
                  new_aer_orgp = orgp_settle * (aer_mass / fines_settle)
                  new_aer_disp = disp_settle * (aer_mass / fines_settle)
                  new_anaer_minpa = ((fines_settle - aer_mass)/ fines_settle) * minpa_settle + &
                    ((anaer_mass + aer_mass - fines_settle)/ aer_mass) * minpa_aer
                  new_anaer_minps = ((fines_settle - aer_mass)/ fines_settle) * minps_settle + &
                    ((anaer_mass + aer_mass - fines_settle)/ aer_mass) * minps_aer
                  new_anaer_orgp = ((fines_settle - aer_mass)/ fines_settle) * orgp_settle + &
                    ((anaer_mass + aer_mass - fines_settle)/ aer_mass) * orgp_aer
                  new_anaer_disp = ((fines_settle - aer_mass)/ fines_settle) * disp_settle + &
                    ((anaer_mass + aer_mass - fines_settle)/ aer_mass) * disp_aer
              else
                  new_aer_srfarea = srfarea_settle * (aer_mass / fines_settle)
                  new_anaer_srfarea = srfarea_settle * (anaer_mass / fines_settle)
                  new_aer_minpa = minpa_settle * (aer_mass / fines_settle)
                  new_aer_minps = minps_settle * (aer_mass / fines_settle)
                  new_aer_orgp = orgp_settle * (aer_mass / fines_settle)
                  new_aer_disp = disp_settle * (aer_mass / fines_settle)
                  new_anaer_minpa = minpa_settle * (anaer_mass / fines_settle)
                  new_anaer_minps = minps_settle * (anaer_mass / fines_settle)
                  new_anaer_orgp = orgp_settle * (anaer_mass / fines_settle)
                  new_anaer_disp = disp_settle * (anaer_mass / fines_settle)
              end if
          end if
      end if

      wbody_bed%aer_orgp = 1000 * new_aer_orgp / aer_mass
      wbody_bed%aer_disp = 1000 * new_aer_disp / vol_h2o_aer ! corrected from "/ aer_mass" to "/ vol_h2o_aer" on 02/19/22, KDW
      wbody_bed%aer_minpa = 1000 * new_aer_minpa / aer_mass
      wbody_bed%aer_minps = 1000 * new_aer_minps / aer_mass
      wbody_bed%aer_srfarea = new_aer_srfarea

      wbody_bed%anaer_orgp = 1000 * new_anaer_orgp / anaer_mass
      wbody_bed%anaer_disp = 1000 * new_anaer_disp / vol_h2o_anaer ! corrected from "/ anaer_mass" to "/ vol_h2o_anaer" on 02/19/22, KDW
      wbody_bed%anaer_minpa = 1000 * new_anaer_minpa / anaer_mass
      wbody_bed%anaer_minps = 1000 * new_anaer_minps / anaer_mass
      wbody_bed%anaer_srfarea = new_anaer_srfarea

      wbody_wb%phos_dep = wbody_wb%phos_dep+ orgp_settle + disp_settle + minpa_settle + minps_settle ! KDW, 06/07

      if (wbody_bed%aer_disp < 1.e-6) wbody_bed%aer_disp = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_disp < 1.e-6) wbody_bed%anaer_disp = 0.0 ! KDW, 03/26
      if (wbody_bed%aer_orgp < 1.e-6) wbody_bed%aer_orgp = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_orgp < 1.e-6) wbody_bed%anaer_orgp = 0.0 ! KDW, 03/26
      if (wbody_bed%aer_minpa < 1.e-6) wbody_bed%aer_minpa = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_minpa < 1.e-6) wbody_bed%anaer_minpa = 0.0 ! KDW, 03/26
      if (wbody_bed%aer_minps < 1.e-6) wbody_bed%aer_minps = 0.0 ! KDW, 03/26
      if (wbody_bed%anaer_minps < 1.e-6) wbody_bed%anaer_minps = 0.0 ! KDW, 03/26
      if (wbody%orgp < 1.e-6) wbody%orgp = 0.0 ! KDW, 03/26
      if (wbody%minpa < 1.e-6) wbody%minpa = 0.0 ! KDW, 03/26
      if (wbody%minps < 1.e-6) wbody%minps = 0.0 ! KDW, 03/26
      if (wbody%solp < 1.e-6) wbody%solp = 0.0 ! KDW, 03/26

      !! compute sediment and nutrients leaving reservoir - ppm -> t
      
      if (wbody%flo < 1.e-6) then
        ht2%sed = 0
        ht2%sil = 0
        ht2%cla = 0
        ht2%san = 0
        ht2%sag = 0
        ht2%lag = 0
        ht2%grv = 0
        ht2%no3 = 0
        ht2%orgn = 0
        ht2%orgp = 0
        ht2%minpa = 0
        ht2%minps = 0
        ht2%solp = 0
        ht2%chla = 0
        ht2%cbod = 0
        ht2%dox = 0
        ht2%nh3 = 0
        ht2%no2 = 0
      else
        ht2%sed = sed_ppm * ht2%flo / 1000000. ! KDW 04/08/22, added 0.7 factor due to stratification? if so, do for all sediment associated WQ
        ht2%sil = sil_ppm * ht2%flo / 1000000. ! KDW 04/08/22
        ht2%cla = cla_ppm * ht2%flo / 1000000. ! KDW 04/08/22
        ht2%san = san_ppm * ht2%flo / 1000000. ! KDW 04/08/22
        ht2%sag = sag_ppm * ht2%flo / 1000000. ! KDW 04/08/22
        ht2%lag = lag_ppm * ht2%flo / 1000000. ! KDW 04/08/22
        ht2%grv = grv_ppm * ht2%flo / 1000000. ! KDW 04/08/22

        ht2%no3 = wbody%no3 * ht2%flo / wbody%flo !N not part of sediment routine, but just calculated here with P
        ht2%orgn = wbody%orgn * ht2%flo / wbody%flo ! KDW 04/08/22
        ht2%orgp = wbody%orgp * ht2%flo / wbody%flo ! KDW 04/08/22
        ht2%minpa = wbody%minpa * ht2%flo / wbody%flo ! KDW 04/08/22
        ht2%minps = wbody%minps * ht2%flo / wbody%flo ! KDW 04/08/22
        ht2%solp = wbody%solp * ht2%flo / wbody%flo ! KDW 04/08/22, add 0.7 factor due to stratification?
        ht2%chla = wbody%chla * ht2%flo / wbody%flo ! KDW 04/08/22
        ht2%cbod = wbody%cbod * ht2%flo / wbody%flo ! KDW 04/08/22
        ht2%dox = wbody%dox * ht2%flo / wbody%flo ! KDW 04/08/22
        ht2%nh3 = wbody%nh3 * ht2%flo / wbody%flo
        ht2%no2 = wbody%no2 * ht2%flo / wbody%flo
      endif
      
      return
      end subroutine res_sediment2