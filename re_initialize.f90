      subroutine re_initialize
    
      use hru_module, only : hru, hru_init, bss
      use soil_module
      use plant_module
      use organic_mineral_mass_module
      use mgt_operations_module
      use hydrograph_module
      use hru_lte_module
      use sd_channel_module
      use channel_module
      use aquifer_module
      use water_body_module ! KDW
      use reservoir_data_module, only : wet_dat, wet_hyd
      use reservoir_module, only : wet_ob, res_ob
      
      implicit none

      integer :: iihru, iires, iprop, ihyd! KDW
      real :: x1, wet_h, wet_h1, wet_fr ! KDW 04/08

      !! re-initialize all hru parameters
      if (sp_ob%hru > 0) then
        hru = hru_init
        soil = soil_init
        soil1 = soil1_init
        rsd1 = rsd1_init
        pcom = pcom_init
        pl_mass = pl_mass_init
        wet = wet_om_init
        bss = 0.
        do iihru = 1, sp_ob%hru
          iprop = hru(iihru)%dbs%surf_stor! KDW 04/08
          if (iprop > 0) then ! KDW 04/08, and all below except as noted
            wet_bed(iihru) = wbody_bed_init  ! KDW 03/29
            wet_wat_d(iihru)%area_ha = 0.
            ihyd = wet_dat(iprop)%hyd
            if (wet(iihru)%flo > 0.) then
              x1 = wet_hyd(ihyd)%bcoef ** 2 + 4. * wet_hyd(ihyd)%ccoef * (1. - wet(iihru)%flo / wet_ob(iihru)%pvol)
              if (x1 < 1.e-6) then
                wet_h = 0.
              else
                wet_h1 = (-wet_hyd(ihyd)%bcoef - sqrt(x1)) / (2. * wet_hyd(ihyd)%ccoef)
                wet_h = wet_h1 + wet_hyd(ihyd)%bcoef
              end if
              wet_fr = (1. + wet_hyd(ihyd)%acoef * wet_h)
              wet_fr = min(wet_fr,1.)
              wet_wat_d(iihru)%area_ha = hru(iihru)%area_ha * wet_fr
            endif
          end if ! end KDW 04/08
        end do
      end if
      
      !! re-initialize hru_lte parameters
      if (sp_ob%hru_lte > 0) then
        hlt = hlt_init
      end if
      
      !! re-initialize channel lte storage and dimensions
      if (sp_ob%chandeg > 0) then
        sdch_init = sd_ch
        ch_stor = ch_om_water_init
      end if

      !! Added by KDW
      !! re-initialize channel lte storage and dimensions
      if (sp_ob%chan > 0) then
        ch = channel_init ! KDW
      end if
      !! End addition by KDW

      !! re-initialize reservoir storage
      if (sp_ob%res > 0) then
        res = res_om_init
        do iires = 1, sp_ob%res
          res_wat_d(iires)%area_ha = res_ob(iires)%br1 * res(iires)%flo ** res_ob(iires)%br2 ! KDW 04/08
          res_bed(iires) = wbody_bed_init ! KDW 03/29
        end do
      end if
      
      !! re-initialize aquifer storage
      if (sp_ob%aqu > 0) then
        aqu_om_init = aqu_om_init
      end if

      return
      end subroutine re_initialize