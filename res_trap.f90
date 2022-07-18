      subroutine res_trap (jres, ihyd, ised, mxsa)

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
      real :: sil_settle, cla_settle, sed_settle       ! tons - KDW
      real :: san_settle, sag_settle, lag_settle, grv_settle    ! tons - KDW
      real :: fines_inflow, srfarea_inflow, srfarea_settle, srfarea_wtclm
      real :: conv_mass2vol, vol_h2o_aer, vol_h2o_anaer, vol_h2o_newdep                !m^3          |volume of pore water in each layer - KDW
      real :: aer_mass, anaer_mass     
      real :: minpa_aer, minps_aer, orgp_aer, disp_aer                !kg          |kg P by form and location - KDW
      real :: minpa_anaer, minps_anaer, orgp_anaer, disp_anaer                !kg          |kg P by form and location - KDW
      real :: minpa_wtclm, minps_wtclm, orgp_wtclm, disp_wtclm                !kg          |kg P by form and location - KDW
      real :: percsfar_settle, percmass_settle, percvol_settle
      real :: orgp_settle, disp_settle, minpa_settle, minps_settle      !kg          |kg P settling - KDW
      real :: new_aer_minpa, new_aer_minps, new_aer_orgp, new_aer_disp               !kg          |kg P by form and location, end of time step - KDW
      real :: new_anaer_minpa, new_anaer_minps, new_anaer_orgp, new_anaer_disp                !kg          |kg P by form and location, end of time step - KDW
      real :: new_aer_srfarea, new_anaer_srfarea, settle_remain
      real :: trap_gra, trap_san, trap_sil, trap_cla, trap_sag, trap_lag ! KDW 04/08/22
      real, intent (in) :: mxsa ! - KDW

      sil_settle = 0. ; cla_settle = 0. ; sed_settle = 0.
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
      new_aer_srfarea = 0. ; new_anaer_srfarea = 0.; settle_remain = 0.
      trap_gra = 0. ; trap_san = 0. ; trap_sil = 0. ; trap_cla = 0. ; trap_sag = 0. ; trap_lag = 0.
 


      ! aer_mass = .001 * 1 * (mxsa * wbody_bed%bed_bd) ! depth of aer layer is 1mm
      ! anaer_mass = .1 * (mxsa * wbody_bed%bed_bd)
      aer_mass = .002 * 1 * (mxsa * wbody_bed%bed_bd) ! depth of aer layer is 1mm! KDW debug 04/28/22, test thicker layers, replacing above
      anaer_mass = .2 * (mxsa * wbody_bed%bed_bd)! KDW debug 04/28/22, test thicker layers, replacing above
      fines_inflow = ht1%cla + ht1%sil + ht1%sag + ht1%san
      !srfarea_inflow = ht1%cla + ht1%sil * .2 + ht1%sag * (1./15.) + ht1%san * .01
      srfarea_inflow = 3. * ((ht1%cla/0.000002) + (ht1%sil/0.00001) + (ht1%sag/0.00003) + (ht1%san/0.0002)) ! KDW 04/14/22

      if (wbody%flo < 1.e-6) then
        ! reservoir is empty
        san_settle = ht1%san
        sag_settle = ht1%sag
        lag_settle = ht1%lag
        grv_settle = ht1%grv
        sil_settle = ht1%sil
        cla_settle = ht1%cla
        ht1%sil = 0
        ht1%cla = 0
        ht1%sag = 0
        ht1%san = 0
        ht1%lag = 0
        ht1%grv = 0
        ht1%sed = 0
      else

        !! compute new sediment concentration in reservoir
	     if (ht1%sed < 1.e-6) ht1%sed = 0.0      
        !! velsetl = 1.35 for clay particle m/d
	     if (wbody_wb%area_ha > 1.e-6 .and. ht2%flo > 1.e-6) then
          velofl = (ht2%flo / wbody_wb%area_ha) / 10000.  ! m3/d / ha * 10000. = m/d - changed from wbody%flo to ht2%flo to match old formulation - KDW
	        !trapres = res_sed(ised)%velsetlr / velofl ! KDW debug/temp 04/08/22
          trap_gra = 10000 * res_sed(ised)%velsetlr * (2.00**2.) / velofl ! KDW debug/temp 04/08/22
          trap_san = 10000 * res_sed(ised)%velsetlr *(0.20**2.) / velofl! KDW debug/temp 04/08/22
          trap_sil = 10000 * res_sed(ised)%velsetlr *(0.01**2.) / velofl! KDW debug/temp 04/08/22
          trap_cla = 10000 * res_sed(ised)%velsetlr *(0.002**2.) / velofl! KDW debug/temp 04/08/22
          trap_sag = 10000 * res_sed(ised)%velsetlr *(0.03**2.) / velofl! KDW debug/temp 04/08/22
          trap_lag = 10000 * res_sed(ised)%velsetlr *(0.50**2.) / velofl! KDW debug/temp 04/08/22
          !trapres = setlr / velofl! KDW debug/temp 04/08/22
       else
          !trapres = 1 ! - KDW 05/12
          trap_gra = 1.! KDW debug/temp 04/08/22
          trap_san = 1.! KDW debug/temp 04/08/22
          trap_sil = 1.! KDW debug/temp 04/08/22
          trap_cla = 1.! KDW debug/temp 04/08/22
          trap_sag = 1.! KDW debug/temp 04/08/22
          trap_lag = 1.! KDW debug/temp 04/08/22
       end if
       !if (trapres > 1.) trapres = 1.
       if (trap_gra > 1.) trap_gra = 1.! KDW debug/temp 04/08/22
       if (trap_san > 1.) trap_san = 1.! KDW debug/temp 04/08/22
       if (trap_sil > 1.) trap_sil = 1.! KDW debug/temp 04/08/22
       if (trap_cla > 1.) trap_cla = 1.! KDW debug/temp 04/08/22
       if (trap_sag > 1.) trap_sag = 1.! KDW debug/temp 04/08/22
       if (trap_lag > 1.) trap_lag = 1.! KDW debug/temp 04/08/22
       !susp =  1. - res_sed(ised)%trapeff * trapres !! editted so not overly trapping,(was just "1. - trapres") - KDW

       grv_settle = ht1%grv * res_sed(ised)%trapeff * trap_gra! KDW debug/temp 04/08/22
       san_settle = ht1%san * res_sed(ised)%trapeff * trap_san! KDW debug/temp 04/08/22
       sil_settle = ht1%sil * res_sed(ised)%trapeff * trap_sil! KDW debug/temp 04/08/22
       cla_settle = ht1%cla * res_sed(ised)%trapeff * trap_cla! KDW debug/temp 04/08/22
       sag_settle = ht1%sag * res_sed(ised)%trapeff * trap_sag! KDW debug/temp 04/08/22
       lag_settle = ht1%lag * res_sed(ised)%trapeff * trap_lag! KDW debug/temp 04/08/22

       !!! Allow larger sediment to settle first - KDW, 08/02
      !  settle_remain = ht1%sed * (1-susp)
      !  if (settle_remain <= ht1%grv) then
      !     grv_settle = settle_remain
      !     settle_remain = 0
      !  else
      !     grv_settle = ht1%grv
      !     settle_remain = settle_remain - grv_settle
      !     if (settle_remain <= ht1%lag) then
      !       lag_settle = settle_remain
      !       settle_remain = 0
      !     else
      !       lag_settle = ht1%lag
      !       settle_remain = settle_remain - lag_settle
      !       if (settle_remain <= ht1%san) then
      !         san_settle = settle_remain
      !         settle_remain = 0
      !       else
      !         san_settle = ht1%san
      !         settle_remain = settle_remain - san_settle
      !         if (settle_remain <= ht1%sag) then
      !           sag_settle = settle_remain
      !           settle_remain = 0
      !         else
      !           sag_settle = ht1%sag
      !           settle_remain = settle_remain - sag_settle
      !           if (settle_remain <= ht1%sil) then
      !             sil_settle = settle_remain
      !             settle_remain = 0
      !           else
      !             sil_settle = ht1%sil
      !             settle_remain = settle_remain - sil_settle
      !             if (settle_remain <= ht1%cla) then
      !               cla_settle = settle_remain
      !               settle_remain = 0
      !             else
      !               cla_settle = ht1%cla
      !               settle_remain = settle_remain - cla_settle
      !             end if
      !           end if
      !         end if
      !       end if
      !     end if
      !  end if
       ht1%grv = ht1%grv - grv_settle
       ht1%lag = ht1%lag - lag_settle
       ht1%san = ht1%san - san_settle
       ht1%sag = ht1%sag - sag_settle
       ht1%sil = ht1%sil - sil_settle
       ht1%cla = ht1%cla - cla_settle
       ht1%sed = ht1%sed - (grv_settle + lag_settle + san_settle + sag_settle + &
            sil_settle + cla_settle) ! KDW 08/02


        !! old formulation, all sizes settled in equal proportion - KDW
        ! san_settle = ht1%san * (1-susp)
        ! sag_settle = ht1%sag * (1-susp)
        ! lag_settle = ht1%lag * (1-susp)
        ! grv_settle = ht1%grv * (1-susp)
        ! sil_settle = ht1%sil * (1-susp)
        ! cla_settle = ht1%cla * (1-susp)
        ! ht1%sil = ht1%sil * susp
        ! ht1%cla = ht1%cla * susp
        ! ht1%sag = ht1%sag * susp
        ! ht1%san = ht1%san * susp
        ! ht1%lag = ht1%lag * susp
        ! ht1%grv = ht1%grv * susp
        ! ht1%sed = ht1%sed * susp
      end if

      sed_settle = sil_settle + cla_settle + san_settle + sag_settle !actually just the "fines" doesn't include lag_settle, grv_settle
      !srfarea_settle = cla_settle + sil_settle * .2 + sag_settle * (1./15.) + san_settle * .01
      srfarea_settle = 3. * ((cla_settle/0.000002) + (sil_settle/0.00001) + (sag_settle/0.00003) + (san_settle/0.0002)) ! KDW 04/14/22

!! Phosphorus movement associated with sediment deposition - KDW
        
  ! What are water volumes in layers (and in newly deposited sediment)
        conv_mass2vol = wbody_bed%bed_por / wbody_bed%bed_bd ! convert from mass of sediment to volume of h2o in pore space
        vol_h2o_aer = wbody_bed%bed_por * (mxsa * .001)
        vol_h2o_anaer = wbody_bed%bed_por * (mxsa * .1)
        vol_h2o_newdep = wbody_bed%bed_por * (sed_settle / wbody_bed%bed_bd)
  ! How much P mass of each form is in each layer and the inflow wtclm, initially
        minpa_aer = wbody_bed%aer_minpa * aer_mass / 1000 ! kg
        minps_aer = wbody_bed%aer_minps * aer_mass / 1000 ! kg
        minpa_anaer = wbody_bed%anaer_minpa * anaer_mass / 1000 ! kg
        minps_anaer = wbody_bed%anaer_minps * anaer_mass / 1000 ! kg
        minpa_wtclm = ht1%minpa ! kg
        minps_wtclm = ht1%minps ! kg
        orgp_aer = wbody_bed%aer_orgp * aer_mass / 1000 ! kg
        disp_aer = wbody_bed%aer_disp * vol_h2o_aer / 1000 ! kg 
        orgp_anaer = wbody_bed%anaer_orgp * anaer_mass / 1000 ! kg
        disp_anaer = wbody_bed%anaer_disp * vol_h2o_anaer / 1000 ! kg 
        orgp_wtclm = ht1%orgp ! kg
        disp_wtclm = ht1%solp ! kg 

  ! How much phosphorus is associated with the different pools of depositing sediment
        if (srfarea_inflow < 1.e-6) then
          percvol_settle = 1
        else
          percsfar_settle = srfarea_settle / srfarea_inflow
        endif
        if (fines_inflow < 1.e-6) then
          percmass_settle = 1
        else
          percmass_settle = sed_settle / fines_inflow !recall, sed_settle is actually fines settle
        endif
        if (ht1%flo < 1.e-6) then
          percsfar_settle = 1
        else
          percvol_settle = vol_h2o_newdep / ht1%flo
        endif
        orgp_settle = percsfar_settle * ht1%orgp ! KDW, changed from "percmass *" to "percsfar *"  02/17/22
        ht1%orgp = (1- percsfar_settle)* ht1%orgp
        disp_settle = percvol_settle * ht1%solp
        ht1%solp = (1- percvol_settle)* ht1%solp
        minpa_settle = percsfar_settle * ht1%minpa
        ht1%minpa = (1- percsfar_settle)* ht1%minpa
        minps_settle = percsfar_settle * ht1%minps
        ht1%minps = (1- percsfar_settle)* ht1%minps

  ! How much phosphorus is in each layer after erosion and deposition 

        if (sed_settle <= aer_mass) then ! = added to inequality on 02/19/22, KDW
            new_aer_srfarea = srfarea_settle + ((aer_mass - sed_settle) / aer_mass) * wbody_bed%aer_srfarea
            new_anaer_srfarea = (sed_settle / aer_mass) * wbody_bed%aer_srfarea + &
              ((anaer_mass - sed_settle) / anaer_mass) * wbody_bed%anaer_srfarea
            new_aer_minpa = minpa_settle + ((aer_mass - sed_settle) / aer_mass) * minpa_aer
            new_aer_minps = minps_settle + ((aer_mass - sed_settle) / aer_mass) * minps_aer
            new_aer_orgp = orgp_settle + ((aer_mass - sed_settle) / aer_mass) * orgp_aer
            new_aer_disp = disp_settle + ((aer_mass - sed_settle) / aer_mass) * disp_aer
            new_anaer_minpa = (sed_settle / aer_mass) * minpa_aer + &
              ((anaer_mass - sed_settle) / anaer_mass) * minpa_anaer
            new_anaer_minps = (sed_settle / aer_mass) * minps_aer + &
              ((anaer_mass - sed_settle) / anaer_mass) * minps_anaer
            new_anaer_orgp = (sed_settle / aer_mass) * orgp_aer + &
              ((anaer_mass - sed_settle) / anaer_mass) * orgp_anaer
            new_anaer_disp = (sed_settle / aer_mass) * disp_aer + &
              ((anaer_mass - sed_settle) / anaer_mass) * disp_anaer
        else
            if (sed_settle <= anaer_mass) then! = added to inequality on 02/19/22, KDW
                new_aer_srfarea = srfarea_settle * (aer_mass / sed_settle)
                new_anaer_srfarea = ((sed_settle - aer_mass)/ sed_settle) * srfarea_settle + &
                  wbody_bed%aer_srfarea + ((anaer_mass - sed_settle) / anaer_mass) * wbody_bed%anaer_srfarea
                new_aer_minpa = minpa_settle * (aer_mass / sed_settle)
                new_aer_minps = minps_settle * (aer_mass / sed_settle)
                new_aer_orgp = orgp_settle * (aer_mass / sed_settle)
                new_aer_disp = disp_settle * (aer_mass / sed_settle)
                new_anaer_minpa = ((sed_settle - aer_mass)/ sed_settle) * minpa_settle + &
                  minpa_aer + ((anaer_mass - sed_settle) / anaer_mass) * minpa_anaer
                new_anaer_minps = ((sed_settle - aer_mass)/ sed_settle) * minps_settle + &
                  minps_aer + ((anaer_mass - sed_settle) / anaer_mass) * minps_anaer
                new_anaer_orgp = ((sed_settle - aer_mass)/ sed_settle) * orgp_settle + &
                  orgp_aer + ((anaer_mass - sed_settle) / anaer_mass) * orgp_anaer
                new_anaer_disp = ((sed_settle - aer_mass)/ sed_settle) * disp_settle + &
                  disp_aer + ((anaer_mass - sed_settle) / anaer_mass) * disp_anaer
            else
                if (sed_settle <= aer_mass + anaer_mass) then! = added to inequality on 02/19/22, KDW
                    new_aer_srfarea = srfarea_settle * (aer_mass / sed_settle)
                    new_anaer_srfarea = ((sed_settle - aer_mass)/ sed_settle) * srfarea_settle + &
                      ((anaer_mass + aer_mass - sed_settle)/ aer_mass) * wbody_bed%aer_srfarea
                    new_aer_minpa = minpa_settle * (aer_mass / sed_settle)
                    new_aer_minps = minps_settle * (aer_mass / sed_settle)
                    new_aer_orgp = orgp_settle * (aer_mass / sed_settle)
                    new_aer_disp = disp_settle * (aer_mass / sed_settle)
                    new_anaer_minpa = ((sed_settle - aer_mass)/ sed_settle) * minpa_settle + &
                      ((anaer_mass + aer_mass - sed_settle)/ aer_mass) * minpa_aer
                    new_anaer_minps = ((sed_settle - aer_mass)/ sed_settle) * minps_settle + &
                      ((anaer_mass + aer_mass - sed_settle)/ aer_mass) * minps_aer
                    new_anaer_orgp = ((sed_settle - aer_mass)/ sed_settle) * orgp_settle + &
                      ((anaer_mass + aer_mass - sed_settle)/ aer_mass) * orgp_aer
                    new_anaer_disp = ((sed_settle - aer_mass)/ sed_settle) * disp_settle + &
                      ((anaer_mass + aer_mass - sed_settle)/ aer_mass) * disp_aer
                else
                    new_aer_srfarea = srfarea_settle * (aer_mass / sed_settle)
                    new_anaer_srfarea = srfarea_settle * (anaer_mass / sed_settle)
                    new_aer_minpa = minpa_settle * (aer_mass / sed_settle)
                    new_aer_minps = minps_settle * (aer_mass / sed_settle)
                    new_aer_orgp = orgp_settle * (aer_mass / sed_settle)
                    new_aer_disp = disp_settle * (aer_mass / sed_settle)
                    new_anaer_minpa = minpa_settle * (anaer_mass / sed_settle)
                    new_anaer_minps = minps_settle * (anaer_mass / sed_settle)
                    new_anaer_orgp = orgp_settle * (anaer_mass / sed_settle)
                    new_anaer_disp = disp_settle * (anaer_mass / sed_settle)
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

        wbody_wb%phos_dep = orgp_settle + disp_settle + minpa_settle + minps_settle ! KDW, 06/07

        if (wbody_bed%aer_disp < 1.e-6) wbody_bed%aer_disp = 0.0 ! KDW, 03/26
        if (wbody_bed%anaer_disp < 1.e-6) wbody_bed%anaer_disp = 0.0 ! KDW, 03/26
        if (wbody_bed%aer_orgp < 1.e-6) wbody_bed%aer_orgp = 0.0 ! KDW, 03/26
        if (wbody_bed%anaer_orgp < 1.e-6) wbody_bed%anaer_orgp = 0.0 ! KDW, 03/26
        if (wbody_bed%aer_minpa < 1.e-6) wbody_bed%aer_minpa = 0.0 ! KDW, 03/26
        if (wbody_bed%anaer_minpa < 1.e-6) wbody_bed%anaer_minpa = 0.0 ! KDW, 03/26
        if (wbody_bed%aer_minps < 1.e-6) wbody_bed%aer_minps = 0.0 ! KDW, 03/26
        if (wbody_bed%anaer_minps < 1.e-6) wbody_bed%anaer_minps = 0.0 ! KDW, 03/26
        if (ht1%orgp < 1.e-6) ht1%orgp = 0.0 ! KDW, 03/26
        if (ht1%minpa < 1.e-6) ht1%minpa = 0.0 ! KDW, 03/26
        if (ht1%minps < 1.e-6) ht1%minps = 0.0 ! KDW, 03/26
        if (ht1%solp < 1.e-6) ht1%solp = 0.0 ! KDW, 03/26

      return
      end subroutine res_trap