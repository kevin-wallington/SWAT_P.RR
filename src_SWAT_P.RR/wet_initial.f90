      subroutine wet_initial
      
      use reservoir_module
      use reservoir_data_module
      !use reservoir_data_module
      use hydrograph_module
      use hru_module, only : hru, ihru
      use maximum_data_module
      use water_body_module
      use channel_module
      
      implicit none
      
      integer :: iprop          !              | 
      integer :: ihyd           !none          |counter 
      integer :: init           !              |
      real :: cnv               !none          |conversion factor (mm => m^3)
      real :: x1                !              |
      real :: wet_h             !              |
      real :: wet_h1            !              | 
      real :: wet_fr            !              | 
      integer :: init_om
      integer :: init_str ! KDW 11/04/21
      !integer :: iob ! - KDW
      integer :: ch_par ! - KDW
      real :: aer_mass ! - KDW
      real :: anaer_mass ! - KDW
      
  
      do ihru = 1, sp_ob%hru
        !! set initial volumes and convert units
        iprop = hru(ihru)%dbs%surf_stor
        if (iprop > 0) then
          ihyd = wet_dat(iprop)%hyd
          !! ha*mm*10. => m**3
          wet_ob(ihru)%evol = wet_hyd(ihyd)%esa * hru(ihru)%area_ha * wet_hyd(ihyd)%edep * 10.
          wet_ob(ihru)%pvol = wet_hyd(ihyd)%psa * hru(ihru)%area_ha * wet_hyd(ihyd)%pdep * 10.
          wet_ob(ihru)%psa = wet_hyd(ihyd)%psa * hru(ihru)%area_ha 
          wet_ob(ihru)%esa = wet_hyd(ihyd)%esa * hru(ihru)%area_ha 

          !!set initial n and p concentrations --> (ppm) * (m^3) / 1000 = kg    !! ppm = t/m^3 * 10^6
          init = wet_dat(iprop)%init
          init_om = wet_init(init)%org_min
          cnv = om_init_water(init_om)%flo * wet_ob(ihru)%pvol / 1000.
          wet(ihru) = cnv * om_init_water(init_om)
          wet(ihru)%flo = om_init_water(init_om)%flo * wet_ob(ihru)%pvol
          wet_om_init(ihru) = wet(ihru)

          !! initialize wetland "bed" concentrations, KDW 11/04/21
          init_str = wet_init(init)%bed
          wet_bed(ihru)%aer_orgp = streambed_init(init_str)%aer_orgp
          wet_bed(ihru)%aer_minpa = streambed_init(init_str)%aer_minpa
          wet_bed(ihru)%aer_minps = streambed_init(init_str)%aer_minps
          wet_bed(ihru)%aer_disp = streambed_init(init_str)%aer_disp
          wet_bed(ihru)%anaer_orgp = streambed_init(init_str)%anaer_orgp
          wet_bed(ihru)%anaer_minpa = streambed_init(init_str)%anaer_minpa
          wet_bed(ihru)%anaer_minps = streambed_init(init_str)%anaer_minps
          wet_bed(ihru)%anaer_disp = streambed_init(init_str)%anaer_disp

          !! update surface area
          !! wetland on hru - solve quadratic to find new depth
          wet_wat_d(ihru)%area_ha = 0.
          if (wet(ihru)%flo > 0.) then
            x1 = wet_hyd(ihyd)%bcoef ** 2 + 4. * wet_hyd(ihyd)%ccoef * (1. - wet(ihru)%flo / wet_ob(ihru)%pvol)
            if (x1 < 1.e-6) then
              wet_h = 0.
            else
              wet_h1 = (-wet_hyd(ihyd)%bcoef - sqrt(x1)) / (2. * wet_hyd(ihyd)%ccoef)
              wet_h = wet_h1 + wet_hyd(ihyd)%bcoef
            end if
            wet_fr = (1. + wet_hyd(ihyd)%acoef * wet_h)
            wet_fr = min(wet_fr,1.)
            wet_wat_d(ihru)%area_ha = hru(ihru)%area_ha * wet_hyd(ihyd)%psa * wet_fr ! should the psa term be removed to reconcile with line 144 wetland_control? KDW
          end if 

          !! initialize surface area of aerobic and anaerobic pools - KDW
          !iob = hru(ihru)%obj_no
          !ch_par = ob_out(iob)%objno !this indicates which *channel* to use for bed sediment concentrations below - KDW
          ch_par = 1 ! - temporary, b/c code is currently structured to run res_intital and wet_initial before hyd_connect_out - KDW
          !aer_mass = .001 * 1 * (10000 * wet_ob(ihru)%esa * wet_bed(ihru)%bed_bd)
          aer_mass = .002 * 1 * (10000 * wet_ob(ihru)%esa * wet_bed(ihru)%bed_bd)! KDW debug 04/28/22, test thicker layers, replacing above
          !wet_bed(ihru)%aer_srfarea = aer_mass * (ch(ch_par)%bed_cla + ch(ch_par)%bed_sil * .2 + ch(ch_par)%bed_san * .01)
          wet_bed(ihru)%aer_srfarea = 3 * aer_mass * (ch(ch_par)%bed_cla/0.000002 + ch(ch_par)%bed_sil/0.00001 + ch(ch_par)%bed_san/0.0002) ! KDW 04/19/22
          !anaer_mass = .1 * 1 * (10000 * wet_ob(ihru)%esa * wet_bed(ihru)%bed_bd)
          anaer_mass = .2 * 1 * (10000 * wet_ob(ihru)%esa * wet_bed(ihru)%bed_bd)! KDW debug 04/28/22, test thicker layers, replacing above
          !wet_bed(ihru)%anaer_srfarea = anaer_mass * (ch(ch_par)%bed_cla + ch(ch_par)%bed_sil * .2 + ch(ch_par)%bed_san * .01)
          wet_bed(ihru)%anaer_srfarea = 3 * anaer_mass * (ch(ch_par)%bed_cla/0.000002 + ch(ch_par)%bed_sil/0.00001 + ch(ch_par)%bed_san/0.0002) ! KDW 04/19/22

        end if
      
      end do
      close(105)

      return
      end subroutine wet_initial