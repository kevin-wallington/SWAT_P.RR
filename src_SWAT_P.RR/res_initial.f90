      subroutine res_initial
      
      use reservoir_module
      use maximum_data_module
      use reservoir_data_module
      use hydrograph_module
      use constituent_mass_module
      use pesticide_data_module
      use water_body_module
      use channel_module
      
      implicit none
      
      integer :: ires        !none          |counter
      integer :: iprop       !              |     
      integer :: ihyd        !none          |counter 
      integer :: ised        !none          |counter 
      integer :: lnvol       !              |
      real :: resdif         !              |
      integer :: i           !none          |counter
      integer :: idat        !none          |counter
      integer :: init        !              | 
      integer :: ipest       !none          |counter
      integer :: ipath       !              |
      integer :: isalt       !              |
      integer :: ipest_db    !none      |counter
      !integer :: iob ! - KDW
      integer :: ch_par ! - KDW
      real :: aer_mass ! - KDW
      real :: anaer_mass ! - KDW

      do ires = 1, sp_ob%res
        !! set initial volumes for res and hru types
        !! convert units
        iprop = res_ob(ires)%props
        ihyd = res_dat(iprop)%hyd
        res_ob(ires)%evol = res_hyd(ihyd)%evol * 10000.       !! ha-m => m**3
        res_ob(ires)%pvol = res_hyd(ihyd)%pvol * 10000.       !! ha-m => m**3
        res_ob(ires)%esa = res_hyd(ihyd)%esa
        res_ob(ires)%psa = res_hyd(ihyd)%psa
        
        !! calculate shape parameters for surface area equation
        resdif = res_hyd(ihyd)%evol - res_hyd(ihyd)%pvol
        if ((res_hyd(ihyd)%esa - res_hyd(ihyd)%psa) > 0. .and. resdif > 0.) then
          lnvol = Log10(res_ob(ires)%evol) - Log10(res_ob(ires)%pvol)
          if (lnvol > 1.e-4) then
            res_ob(ires)%br2 = (Log10(res_ob(ires)%esa) - Log10(res_ob(ires)%psa)) / lnvol
          else  
            res_ob(ires)%br2 = (Log10(res_ob(ires)%esa) - Log10(res_ob(ires)%psa)) / 0.001
          end if
          if (res_ob(ires)%br2 > 0.9) then
            res_ob(ires)%br2 = 0.9
            res_ob(ires)%br1 = (res_ob(ires)%psa / res_ob(ires)%pvol) ** 0.9
          else
            res_ob(ires)%br1 = (res_ob(ires)%esa / res_ob(ires)%evol) ** res_ob(ires)%br2
          end if  
        else
          res_ob(ires)%br2 = 0.9
          if (res_ob(ires)%pvol > 1.e-6) then
            res_ob(ires)%br1 = (res_ob(ires)%psa / res_ob(ires)%pvol) ** 0.9
          else
            res_ob(ires)%br1 = .1
          end if
        end if
        
      end do
      
      do ires = 1, sp_ob%res
        idat = res_ob(ires)%props
        i = res_dat(idat)%init
        
        !! initialize surface area of aerobic and anaerobic pools - KDW
        !iob = res_ob(ires)%ob
        !ch_par = ob_out(iob)%objno !this indicates which *channel* to use for bed sediment concentrations below - KDW
        ch_par = 1 ! - temporary, b/c code is currently structured to run res_intital and wet_initial before hyd_connect_out - KDW
        ! aer_mass = .001 * 1 * (10000 * res_ob(ires)%esa * res_bed(ires)%bed_bd)
        aer_mass = .002 * 1 * (10000 * res_ob(ires)%esa * res_bed(ires)%bed_bd)! KDW debug 04/28/22, test thicker layers, replacing above
        !res_bed(ires)%aer_srfarea = aer_mass * (ch(ch_par)%bed_cla + ch(ch_par)%bed_sil * .2 + ch(ch_par)%bed_san * .01)
        res_bed(ires)%aer_srfarea = 3 * aer_mass * (ch(ch_par)%bed_cla/0.000002 + ch(ch_par)%bed_sil/0.00001 + ch(ch_par)%bed_san/0.0002) ! KDW 04/19/22
        !anaer_mass = .1 * 1 * (10000 * res_ob(ires)%esa * res_bed(ires)%bed_bd)
        anaer_mass = .2 * 1 * (10000 * res_ob(ires)%esa * res_bed(ires)%bed_bd)! KDW debug 04/28/22, test thicker layers, replacing above
        !res_bed(ires)%anaer_srfarea = anaer_mass * (ch(ch_par)%bed_cla + ch(ch_par)%bed_sil * .2 + ch(ch_par)%bed_san * .01) 
        res_bed(ires)%anaer_srfarea = 3 * anaer_mass * (ch(ch_par)%bed_cla/0.000002 + ch(ch_par)%bed_sil/0.00001 + ch(ch_par)%bed_san/0.0002) ! KDW 04/19/22

        !! initialize org-min in reservoir
        init = res_init(i)%org_min
        res(ires) = om_init_water(init)
        call res_convert_mass (res(ires), res_ob(ires)%pvol)

        !! initialize reservoir bed concentrations, KDW 11/04/21
        init = res_init(i)%bed
        res_bed(ires)%aer_orgp = streambed_init(init)%aer_orgp
        res_bed(ires)%aer_minpa = streambed_init(init)%aer_minpa
        res_bed(ires)%aer_minps = streambed_init(init)%aer_minps
        res_bed(ires)%aer_disp = streambed_init(init)%aer_disp
        res_bed(ires)%anaer_orgp = streambed_init(init)%anaer_orgp
        res_bed(ires)%anaer_minpa = streambed_init(init)%anaer_minpa
        res_bed(ires)%anaer_minps = streambed_init(init)%anaer_minps
        res_bed(ires)%anaer_disp = streambed_init(init)%anaer_disp
        
        !! set initial reservoir org-min to reset for soft calibration
        res_om_init(ires) = res(ires)

        !! initialize pesticides in reservoir water and benthic from input data
        init = res_init(i)%pest
        do ipest = 1, cs_db%num_pests
          ipest_db = cs_db%pest_num(ipest)
          res_water(ires)%pest(ipest) = pest_water_ini(init)%water(ipest)
          res_benthic(ires)%pest(ipest) = pest_water_ini(init)%benthic(ipest)
          !! calculate mixing velocity using molecular weight and porosity
          ised = res_dat(idat)%sed
          res_ob(ires)%aq_mix(ipest) = pestdb(ipest_db)%mol_wt * (1. - res_sed(ised)%bd / 2.65)
        end do
                  
        !! initialize pathogens in reservoir water and benthic from input data
        init = res_init(i)%path
        do ipath = 1, cs_db%num_paths
          res_water(ires)%path(ipath) = path_water_ini(init)%water(ipath)
          res_benthic(ires)%path(ipath) = path_water_ini(init)%benthic(ipath)
        end do
                        
        !! initialize salts in reservoir water and benthic from input data
        init = res_init(i)%salt
        do isalt = 1, cs_db%num_salts
          res_water(ires)%salt(isalt) = salt_water_ini(init)%water(isalt)
          res_benthic(ires)%salt(isalt) = salt_water_ini(init)%benthic(isalt)
        end do
        
        !! calculate initial surface area       
        res_wat_d(ires)%area_ha = res_ob(ires)%br1 * res(ires)%flo ** res_ob(ires)%br2

      end do
      close(105)

      return
      end subroutine res_initial