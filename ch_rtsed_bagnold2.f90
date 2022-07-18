subroutine ch_rtsed_bagnold2
      
      !!    ~ ~ ~ PURPOSE ~ ~ ~
      !!    this subroutine routes sediment (and sediment attached phosphorus - KDW)
      !!    from subbasin to basin outlets
      !!    deposition is based on fall velocity and degradation on stream power
      
      !!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
      !!    name        |units         |definition
      !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      !!    ch_d(:)     |m             |average depth of main channel
      !!    ch_li(:)    |km            |initial length of main channel
      !!    ch_n(2,:)   |none          |Manning"s "n" value for the main channel
      !!    ch_s(2,:)   |m/m           |average slope of main channel
      !!    ch_si(:)    |m/m           |initial slope of main channel
      !!    ch_wdr(:)   |m/m           |channel width to depth ratio
      !!    rchdep      |m             |depth of flow on day
      !!    sdti        |m^3/s         |average flow on day in reach
      !!    sedst(:)    |metric tons   |amount of sediment stored in reach
      !!                               |reentrained in channel sediment routing
      !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      !!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
      !!    name        |units         |definition
      !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      !!    ch_d(:)     |m             |average depth of main channel
      !!    ch_s(2,:)   |m/m           |average slope of main channel
      !!    peakr       |m^3/s         |peak runoff rate in channel
      !!    sedst(:)    |metric tons   |amount of sediment stored in reach
      !!    sedrch      |metric tons   |sediment transported out of channel
      !!                               |during time step
      !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      
      !!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
      !!    name        |units         |definition
      !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      !!    dat2        |m             |change in channel depth during time step
      !!    deg         |metric tons   |sediment reentrained in water by channel
      !!                               |degradation
      !!    dep         |metric tons   |sediment deposited on river bottom
      !!    depdeg      |m             |depth of degradation/deposition from original
      !!    depnet      |metric tons   |
      !!    dot         |
      !!    jrch        |none          |reach number
      !!    qdin        |m^3 H2O       |water in reach during time step
      !!    vc          |m/s           |flow velocity in reach
      !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      
      !!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
      !!    Intrinsic: Max
      !!    SWAT: ttcoef
      
      !!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
      !!	Modification to the original SWAT sediment routine
      !!	By Balaji Narasimhan and Peter Allen
      !!    Bagnolds strempower approach combined with Einsteins deposition equation
      !!    Plus particle size tracking.
      
            use basin_module
            use channel_module
            use hydrograph_module
            use channel_data_module
            use channel_velocity_module
            
            implicit none
      
            integer :: ch_d50type          !units         |description
            real :: qdin                   !m^3 H2O       |water in reach during time step
            real :: sedin                  !units         |description
            real :: vc                     !m/s           |flow velocity in reach         
            real :: cyin                   !units         |description
            real :: cych                   !units         |description
            real :: depnet                 !metric tons   |
            real :: deg                    !metric tons   |sediment reentrained in water by channel
                                           !              |degradation 
            real :: dep                    !metric tons   |sediment deposited on river bottom
            real :: tbase                  !units         |description
            real :: depdeg                 !m             |depth of degradation/deposition from original
            real :: dot                    !mm            |actual depth from impermeable layer to water level
                                           !              |above drain during subsurface irrigation
            real :: vs                     !units         |description
            real :: x                      !units         |description
            real :: SC                     !units         |description
            real :: Tcbnk                  !units         |description
            real :: Tcbed                  !units         |description
            real :: Tbank                  !units         |description
            real :: Tbed                   !units         |description
            real :: asinea                 !units         |description
            real :: Tou                    !units         |description
            real :: sanin                  !units         |description
            real :: silin                  !units         |description
            real :: clain                  !units         |description
            real :: sagin                  !units         |description
            real :: lagin                  !units         |description
            real :: grain                  !units         |description
            real :: outfract               !units         |description
            real :: depsan                 !units         |description
            real :: depsil                 !units         |description
            real :: depcla                 !units         |description
            real :: depsag                 !units         |description
            real :: deplag                 !units         |description 
            real :: depgra                 !units         |description
            real :: degsan                 !units         |description
            real :: degsil                 !units         |description
            real :: degcla                 !units         |description
            real :: deggra                 !units         |description
            real :: bnksan                 !units         |description
            real :: bnksil                 !units         |description
            real :: bnkcla                 !units         |description
            real :: bnkgra                 !units         |description
            real :: pdep                   !units         |description
            real :: pdepbed                !units         |description
            real :: bedsize                !units         |description
            real :: USpower                !units         |description
            real :: adddep                 !units         |description
            real :: fpratio                !units         |description
            real :: watdep                 !units         |description
            real :: bnkrt                  !units         |description
            real :: bedrt                  !units         |description
            real :: effbnkbed              !units         |description
            real :: sedinorg               !units         |description
            real :: deg1                   !units         |description
            real :: deg1san                !units         |description
            real :: deg1sil                !units         |description
            real :: deg1cla                !units         |description
            real :: deg1sag                !units         |description
            real :: deg1lag                !units         |description
            real :: deg1gra                !units         |description
            real :: degremain              !units         |description
            real :: c                      !units         |description 
            real :: pbed                   !units         |description
            real :: pbank                  !units         |description
            real :: rh                     !m             |hydraulic radius
            real :: topw                   !m             |top width of main channel
            real :: sfbank                 !units         |description
            real :: vgra                   !units         |description
            real :: vsan                   !units         |description
            real :: vsil                   !units         |description
            real :: vcla                   !units         |description
            real :: vsag                   !units         |description
            real :: vlag                   !units         |description
            real :: dat2                   !m             |change in channel depth during time step 
            real :: fines, fines_wtclm, aer_mass, anaer_mass                  !tons (Mg)          |mass of previously deposited fines - KDW
            real :: aer_mass_remain, anaer_mass_remain                  !tons (Mg)          |mass of sedment from layer which ends still in layer - KDW
            real :: topmass, topcla, topsil, topsag, topsan                !tons          |previously deposited mass in top 1mm (aerobic layer) - KDW
            real :: andep_cla, andep_sil, andep_sag, andep_san, andep_tot                !tons          |previously deposited mass below top 1mm (anaerobic layer) - KDW
            real :: cla_aer, sil_aer, sag_aer, san_aer                !tons          |mass in top 1mm (anaerobic layer) - KDW 
            real :: cla_anaer, sil_anaer, sag_anaer, san_anaer                !tons          |mass below top 1mm (anaerobic layer) - KDW 
            real :: sfarea_anaer, sfarea_aer, sfarea_fines, sfarea_wtclm                !tons          |surface area of sediment, registered as equivilent in mass of clay - KDW
            real :: minpa_aer, minps_aer, orgp_aer, disp_aer                !kg          |kg P by form and location - KDW
            real :: minpa_anaer, minps_anaer, orgp_anaer, disp_anaer                !kg          |kg P by form and location - KDW
            real :: minpa_wtclm, minps_wtclm, orgp_wtclm, disp_wtclm                !kg          |kg P by form and location - KDW
            real :: vol_h2o_aer, vol_h2o_anaer, vol_h2o_newdep                !m^3          |volume of pore water in each layer - KDW
            real :: cla_ero, sil_ero, sag_ero, san_ero, bed_ero                !tons          |mass of eroded sediment by  - KDW
            real :: fines_newdep, bed_netero, ext_mass                !tons          |mass of deposited fines (this time step), net erosion, mass in from external - KDW
            real :: sfarea_ero, sfarea_dep                !tons          |surface area of sediment, registered as equivilent in mass of clay - KDW
            real :: rmn_aer_mass, rmn_anaer_mass, rmn_finmass_aer, rmn_finmass_anaer                !tons          |remaining mass of sediment or previously deposited "fines" in layer - KDW
            real :: rmn_aer_sfar, rmn_anaer_sfar, rmn_finsfar_aer, rmn_finsfar_anaer                !tons          |remaining surface area (equivilent clay mass) of sediment or previously deposited "fines" in layer - KDW
            real :: rmn_old_aer, rmn_old_anaer                !tons          |remaining mass of non-previously deposited sediment in layer - KDW
            real :: percmass_ero_aer, percmass_ero_anaer, percsfar_ero_aer, percsfar_ero_anaer                !%/100          |percent of layer which is eroded, by mass or surface area - KDW
            real :: orgp_ero_aer, disp_ero_aer, minpa_ero_aer, minps_ero_aer                !kg          |kg P eroded by form and location - KDW
            real :: orgp_ero_anaer, disp_ero_anaer, minpa_ero_anaer, minps_ero_anaer                !kg          |kg P eroded by form and location - KDW
            real :: rmn_orgp_aer, rmn_disp_aer, rmn_minpa_aer, rmn_minps_aer                !kg          |kg P which was initially in layer and remains, by form - KDW
            real :: rmn_orgp_anaer, rmn_disp_anaer, rmn_minpa_anaer, rmn_minps_anaer                !kg          |kg P which was initially in layer and remains, by form - KDW
            real :: percmass_dep, percsfar_dep, percvol_dep                !%/100          |percent of water column which is deposited, by mass, volume, or surface area - KDW
            real :: orgp_dep, disp_dep, minpa_dep, minps_dep                !kg          |kg P deposited by form - KDW
            real :: new_minpa_aer, new_minps_aer, new_orgp_aer, new_disp_aer                !kg          |kg P by form and location, end of time step - KDW
            real :: new_minpa_anaer, new_minps_anaer, new_orgp_anaer, new_disp_anaer                !kg          |kg P by form and location, end of time step - KDW
            real :: fines_still_aer, fines_from_aer, fines_still_anaer, fines_from_anaer                !tons          |mass of prev. deposited "fines' which end in layer, according to initial layer - KDW
            real :: old_still_aer, old_from_aer, old_still_anaer, old_from_anaer                !tons          |mass of non-prev.-dep-"fines' which end in layer, according to initial layer - KDW
            real :: minpa_from_aer, minps_from_aer, orgp_from_aer, disp_from_aer                !kg          |kg P from aer to anaer by form - KDW
            real :: minpa_still_aer, minps_still_aer, orgp_still_aer, disp_still_aer                !kg          |kg P from aer and still in aer by form - KDW
            real :: minpa_still_anaer, minps_still_anaer, orgp_still_anaer, disp_still_anaer                !kg          |kg P from anaer and still in anaer by form - KDW
            real :: minpa_from_anaer, minps_from_anaer, orgp_from_anaer, disp_from_anaer                !kg          |kg P from anaer to aer by form - KDW
            real :: conv_mass2vol, conv_conc              !units   |conversions necessary between pools as needed - KDW
            real :: mass_from_aer, mass_still_anaer, mass_still_aer, mass_from_anaer
      
            fines = 0. ; fines_wtclm = 0. ; aer_mass = 0. ; anaer_mass = 0. ! KDW
            aer_mass_remain = 0. ; anaer_mass_remain = 0. ! KDW
            cla_aer = 0. ; sil_aer = 0. ; sag_aer = 0. ; san_aer = 0.  ! KDW
            cla_anaer = 0. ; sil_anaer = 0. ; sag_anaer = 0. ; san_anaer = 0.  ! KDW
            sfarea_anaer = 0. ; sfarea_aer = 0. ; sfarea_fines = 0. ; sfarea_wtclm = 0. ! KDW
            minpa_aer = 0. ; minps_aer = 0. ; orgp_aer = 0. ; disp_aer = 0.                 ! KDW
            minpa_anaer = 0. ; minps_anaer = 0. ; orgp_anaer = 0. ; disp_anaer = 0.                 ! KDW
            minpa_wtclm = 0. ; minps_wtclm = 0. ; orgp_wtclm = 0. ; disp_wtclm = 0.                 ! KDW
            vol_h2o_aer = 0. ; vol_h2o_anaer = 0. ; vol_h2o_newdep = 0. ! KDW
            cla_ero = 0. ; sil_ero = 0. ; sag_ero = 0. ; san_ero = 0. ; bed_ero = 0.  ! KDW
            rmn_aer_mass = 0. ; rmn_anaer_mass = 0. ; rmn_finmass_aer = 0. ; rmn_finmass_anaer = 0.  ! KDW
            rmn_aer_sfar = 0. ; rmn_anaer_sfar = 0. ; rmn_finmass_aer = 0. ; rmn_finsfar_anaer = 0.  ! KDW
            rmn_old_aer = 0. ; rmn_old_anaer = 0. ! KDW
            percmass_ero_aer = 0. ; percmass_ero_anaer = 0. ; percsfar_ero_aer = 0. ; percsfar_ero_anaer = 0. ! KDW
            orgp_ero_aer = 0. ; disp_ero_aer = 0. ; minpa_ero_aer = 0. ; minps_ero_aer = 0. ! KDW
            orgp_ero_anaer = 0. ; disp_ero_anaer = 0. ; minpa_ero_anaer = 0. ; minps_ero_anaer = 0. ! KDW
            rmn_orgp_aer = 0. ; rmn_disp_aer = 0. ; rmn_minpa_aer = 0. ; rmn_minps_aer = 0. ! KDW
            rmn_orgp_anaer = 0. ; rmn_disp_anaer = 0. ; rmn_minpa_anaer = 0. ; rmn_minps_anaer = 0. ! KDW
            percmass_dep = 0. ; percsfar_dep = 0. ; percvol_dep = 0. ! KDW
            orgp_dep = 0. ; disp_dep = 0. ; minpa_dep = 0. ; minps_dep = 0. ! KDW
            new_minpa_aer = 0. ; new_minps_aer = 0. ; new_orgp_aer = 0. ; new_disp_aer = 0. ! KDW
            new_minpa_anaer = 0. ; new_minps_anaer = 0. ; new_orgp_anaer = 0. ; new_disp_anaer = 0. ! KDW
            fines_still_aer = 0. ; fines_from_aer = 0. ; fines_still_anaer = 0. ; fines_from_anaer = 0. ! KDW 
            old_still_aer = 0. ; old_from_aer = 0. ; old_still_anaer = 0. ; old_from_anaer = 0. ! KDW 
            minpa_from_aer = 0. ; minps_from_aer = 0. ; orgp_from_aer = 0. ; disp_from_aer = 0.  ! KDW
            minpa_from_anaer = 0. ; minps_from_anaer = 0. ; orgp_from_anaer = 0. ; disp_from_anaer = 0.  ! KDW
            minpa_still_aer = 0. ; minps_still_aer = 0. ; orgp_still_aer = 0. ; disp_still_aer = 0.  ! KDW
            minpa_still_anaer = 0. ; minps_still_anaer = 0. ; orgp_still_anaer = 0. ; disp_still_anaer = 0.  ! KDW
            conv_mass2vol = 0. ; conv_conc = 0. ! KDW
            mass_from_aer = 0. ; mass_still_anaer = 0. ; mass_still_aer = 0. ; mass_from_anaer = 0.
            topmass = 0. ; topcla = 0 ; topsil = 0. ; topsag = 0. ; topsan = 0. ! KDW
            andep_cla = 0. ; andep_sil = 0. ; andep_sag = 0. ; andep_san = 0. ; andep_tot = 0. ! KDW
            fines_newdep = 0. ; bed_netero = 0. ;  ext_mass = 0.    ! KDW
            sfarea_ero = 0. ; sfarea_dep = 0. ! KDW
      
      !! initialize water in reach during time step
            qdin = 0.
            qdin = rtwtr + ch(jrch)%rchstor
            if (qdin < 1.e-6) qdin = 0.0 ! KDW, 03/26
      
      !! initialize sediment in reach during time step
            sedin = 0.
              sanin = 0.
              silin = 0.
              clain = 0.
              sagin = 0.
              lagin = 0.
            sedin = ob(icmd)%hin%sed  + ch(jrch)%sedst
            sanin = ob(icmd)%hin%san  + ch(jrch)%sanst
            silin = ob(icmd)%hin%sil  + ch(jrch)%silst
            clain = ob(icmd)%hin%cla  + ch(jrch)%clast
            sagin = ob(icmd)%hin%sag  + ch(jrch)%sagst
            lagin = ob(icmd)%hin%lag  + ch(jrch)%lagst
            grain = ob(icmd)%hin%grv  + ch(jrch)%grast
            if (sedin < 1.e-6) sedin = 0.0 ! KDW, 03/26
            if (sanin < 1.e-6) sanin = 0.0 ! KDW, 03/26
            if (silin < 1.e-6) silin = 0.0 ! KDW, 03/26
            if (clain < 1.e-6) clain = 0.0 ! KDW, 03/26
            if (sagin < 1.e-6) sagin = 0.0 ! KDW, 03/26
            if (lagin < 1.e-6) lagin = 0.0 ! KDW, 03/26
            if (grain < 1.e-6) grain = 0.0 ! KDW, 03/26 
            sedinorg = sedin
            
      !!    Embeddment of stored sediment into "permanent" streambed - KDW
            ch(jrch)%depch = ch(jrch)%depch * exp(-bsn_prm%sbdk)
            ch(jrch)%depsanch = ch(jrch)%depsanch * exp(-bsn_prm%sbdk)
            ch(jrch)%depsilch = ch(jrch)%depsilch * exp(-bsn_prm%sbdk)
            ch(jrch)%depclach = ch(jrch)%depclach * exp(-bsn_prm%sbdk)
            ch(jrch)%depsagch = ch(jrch)%depsagch * exp(-bsn_prm%sbdk)
            ch(jrch)%deplagch = ch(jrch)%deplagch * exp(-bsn_prm%sbdk)
            ch(jrch)%depgrach = ch(jrch)%depgrach * exp(-bsn_prm%sbdk)
            if (ch(jrch)%depch < 1.e-6) ch(jrch)%depch = 0.0 ! KDW, 03/26
            if (ch(jrch)%depsanch < 1.e-6) ch(jrch)%depsanch = 0.0 ! KDW, 03/26
            if (ch(jrch)%depsilch < 1.e-6) ch(jrch)%depsilch = 0.0 ! KDW, 03/26
            if (ch(jrch)%depclach < 1.e-6) ch(jrch)%depclach = 0.0 ! KDW, 03/26
            if (ch(jrch)%depsagch < 1.e-6) ch(jrch)%depsagch = 0.0 ! KDW, 03/26
            if (ch(jrch)%deplagch < 1.e-6) ch(jrch)%deplagch = 0.0 ! KDW, 03/26
            if (ch(jrch)%depgrach < 1.e-6) ch(jrch)%depgrach = 0.0 ! KDW, 03/26
      !!    Embeddment of stored sediment into "permanent" streambed - KDW
      
      !!    Set apart the top 1mm of previously deposited sediment to erode first - KDW
            fines = ch(jrch)%depclach + ch(jrch)%depsilch + ch(jrch)%depsagch + ch(jrch)%depsanch
            aer_mass = .001 * 1 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd) ! depth of aer layer is 1
            anaer_mass = .1 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd)
          if (fines >= aer_mass) then
              aer_mass_remain = 0
              topcla  = aer_mass * (ch(jrch)%depclach / fines)
              topsil  = aer_mass * (ch(jrch)%depsilch / fines)
              topsag  = aer_mass * (ch(jrch)%depsagch / fines)
              topsan  = aer_mass * (ch(jrch)%depsanch / fines)
          else
                  aer_mass_remain = aer_mass - fines
                  topcla  = ch(jrch)%depclach
                  topsil  = ch(jrch)%depsilch
                  topsag  = ch(jrch)%depsagch
                  topsan  = ch(jrch)%depsanch
          end if
          cla_aer = topcla +  aer_mass_remain * ch(jrch)%bnk_cla/(1-ch(jrch)%bnk_gra) ! Note, final division term adjusts for not including gravel as part of aerboic/anaerobic mass
          sil_aer = topsil +  aer_mass_remain * ch(jrch)%bnk_sil/(1-ch(jrch)%bnk_gra)
          sag_aer = topsag
          san_aer = topsan +  aer_mass_remain * ch(jrch)%bnk_san/(1-ch(jrch)%bnk_gra)
      
          if (fines <= aer_mass + anaer_mass) then
                  andep_cla = ch(jrch)%depclach - topcla
                  andep_sil = ch(jrch)%depsilch - topsil
                  andep_sag = ch(jrch)%depsagch - topsag
                  andep_san = ch(jrch)%depsanch - topsan
                  andep_tot = andep_cla + andep_sil + andep_sag + andep_san
          else
                  andep_cla  = anaer_mass * (ch(jrch)%depclach / fines)
                  andep_sil  = anaer_mass * (ch(jrch)%depsilch / fines)
                  andep_sag  = anaer_mass * (ch(jrch)%depsagch / fines)
                  andep_san  = anaer_mass * (ch(jrch)%depsanch / fines)
                  andep_tot  = anaer_mass
          end if
          anaer_mass_remain = anaer_mass - andep_tot
          cla_anaer = andep_cla + anaer_mass_remain * ch(jrch)%bnk_cla/(1-ch(jrch)%bnk_gra)
          sil_anaer = andep_sil + anaer_mass_remain * ch(jrch)%bnk_sil/(1-ch(jrch)%bnk_gra)
          sag_anaer = andep_sag
          san_anaer = andep_san + anaer_mass_remain * ch(jrch)%bnk_san/(1-ch(jrch)%bnk_gra)
      
          sfarea_anaer = cla_anaer + sil_anaer * .2 + sag_anaer * (1/15) + san_anaer * .01
          sfarea_aer = cla_aer + sil_aer * .2 + sag_aer * (1/15) + san_aer * .01
          sfarea_fines = ch(jrch)%depclach + ch(jrch)%depsilch * .2 + & 
                  ch(jrch)%depsagch * (1/15) + ch(jrch)%depsanch * .01
          fines_wtclm = clain + silin + sagin + sanin + 0.0001
          sfarea_wtclm = clain + silin * .2 + sagin * (1/15) + sanin * .01 + 0.0001
      !!    Set apart the top 1mm of previously deposited sediment to erode first - KDW     
      
      
      !      if (rtwtr > 0. .and. rchdep > 0.) then !! combined with qdin if statement in line 202      
      
      !! do not perform sediment routing if no water in reach
            if (rtwtr > 0. .and. rchdep > 0. .and. qdin > 0.01) then ! NEW - changed from below, KDW
            !if (qdin > 0.01) then
      
      !! initialize reach peak runoff rate
            peakr = bsn_prm%prf * sdti
      
      !! calculate peak flow velocity
            vc = 0.
            if (rcharea < .010) then
              vc = 0.01
            else
              vc = peakr / rcharea
            end if
            
            if (vc > 5.) vc = 5.
      
            tbase = 0.
            tbase = ch_hyd(jhyd)%l * 1000. / (3600. * 24. * vc)
            if (tbase > 1.) tbase = 1.
      
      !! JIMMY"S NEW IMPROVED METHOD for sediment transport
            cyin = 0.
            cych = 0.
            depnet = 0.
            deg = 0.
      
            deg1 = 0.
            deg1san = 0.
            deg1sil = 0.
            deg1cla = 0.
            deg1sag = 0.
            deg1lag = 0.
            deg1gra = 0.
      
            degrte = 0.
            degremain = 0.
            deggra = 0.
            degsan = 0.
            degsil = 0.
            degcla = 0.
            bnksan = 0.
            bnksil = 0.
            bnkcla = 0.
            bnkgra = 0.
            bnkrte = 0.
            dep = 0.
            depsan = 0.
            depsil = 0.
            depcla = 0.
            depsag = 0.
            deplag = 0.
            depgra = 0.
            watdep = 0.
            bnkrt = 0.
            bedrt = 0.
            effbnkbed = 0.
      
            c = ch_hyd(jhyd)%side
              pbed = ch_vel(jrch)%wid_btm
            pbank = 2. * rchdep * Sqrt(1. + c * c)
            rh = rcharea / (pbed + pbank)
      
            topw = 0.
            if (rchdep <= ch_hyd(jhyd)%d) then
              topw = ch_vel(jrch)%wid_btm + 2. * rchdep * c
              fpratio = 0.
              watdep = rchdep
            else
              topw = 5 * ch_hyd(jhyd)%w + 2. * (rchdep - ch_hyd(jhyd)%d) * 4.
              adddep = rchdep - ch_hyd(jhyd)%d
              !! Area Ratio of water in flood plain to total cross sectional area
              fpratio = (rcharea - ch_vel(jrch)%area -                          &                 
                                          ch_hyd(jhyd)%w*adddep)/rcharea
              fpratio = max(0.,fpratio)
              watdep = ch_hyd(jhyd)%d
            end if
      
      !!	Applied Bank Shear Stress
      !!    Equations from Eaton and Millar (2004)
            SFbank = 10**(-1.4026 * log10((pbed/pbank) + 1.5) + 2.247)
      
            Tou = 9800. * rchdep * ch_hyd(jhyd)%s
      
            asinea = 1. / sqrt((1.**2) + (c**2))
      
            Tbank = Tou * (SFbank/100.) * (topw + pbed) * asinea/ (4.*rchdep)
      
            Tbed  = Tou * (1. - (SFbank/100.)) * (topw/(2.*pbed) + 0.5)
      
      !!    Potential Bank Erosion rate in metric tons per day
      !!    Assumed on an average Only one bank eroding due to meandering of channel
            bnkrte = ch_sed(jsed)%bnk_kd * (Tbank - ch_sed(jsed)%tc_bnk)*1e-06
            if (bnkrte < 0.) bnkrte = 0.
            bnkrte = bnkrte * ch_hyd(jhyd)%l * 1000.* (watdep *               &        
                           Sqrt(1. + c * c)) * ch_sed(jsed)%bnk_bd * 86400.
      
      !!    Potential Bed degradation rate in metric tons per day
            degrte = ch_sed(jsed)%bed_kd * (Tbed - ch_sed(jsed)%tc_bed)*1e-06
            if (degrte < 0.) degrte = 0.
            degrte = degrte * ch_hyd(jhyd)%l * 1000.* ch_vel(jrch)%wid_btm    &  
                                               * ch_sed(jsed)%bed_bd * 86400.
      
      !!    Relative potential for bank/bed erosion
            if (bnkrte + degrte > 1.e-6) then
              bnkrt = bnkrte / (bnkrte + degrte)
            else
              bnkrt = 1.0
            end if
            bnkrt = Min(1.0, bnkrt)
      !!    Relative potential for bed erosion
            bedrt = 1. - bnkrt
      
      !!    Incoming sediment concentration
            cyin = sedin/qdin
      
      !!    Streampower for sediment calculated based on Bagnold (1977) concept
            cych = bsn_prm%spcon * vc ** bsn_prm%spexp
      
      !!    Potential sediment Transport capacity
            depnet = qdin * (cych - cyin) 
      
            if (depnet .LE. 1.e-6) then
              depnet = 0.
              bnkrte = 0.
              degrte = 0.
            else
              !! First the deposited material will be degraded before channel bed or bank erosion
              if (depnet >= ch(jrch)%depch) then
                !! Effective erosion
                effbnkbed = depnet - ch(jrch)%depch
                !! Effective bank erosion
                if (effbnkbed*bnkrt <= bnkrte) bnkrte = effbnkbed*bnkrt
                bnksan = bnkrte * ch(jrch)%bnk_san
                bnksil = bnkrte * ch(jrch)%bnk_sil
                bnkcla = bnkrte * ch(jrch)%bnk_cla
                bnkgra = bnkrte * ch(jrch)%bnk_gra
      
                !! Effective bed erosion
                if (effbnkbed*bedrt <= degrte) degrte = effbnkbed*bedrt
                degsan = degrte * ch(jrch)%bed_san
                degsil = degrte * ch(jrch)%bed_sil
                degcla = degrte * ch(jrch)%bed_cla
                deggra = degrte * ch(jrch)%bed_gra
      
                deg1 = ch(jrch)%depch
                deg1san = ch(jrch)%depsanch
                deg1sil = ch(jrch)%depsilch
                deg1cla = ch(jrch)%depclach
                deg1sag = ch(jrch)%depsagch
                deg1lag = ch(jrch)%deplagch
                deg1gra = ch(jrch)%depgrach
      
                ch(jrch)%depch = 0.
                ch(jrch)%depsanch = 0.
                ch(jrch)%depsilch = 0.
                ch(jrch)%depclach = 0.
                ch(jrch)%depsagch = 0.
                ch(jrch)%deplagch = 0.
                ch(jrch)%depgrach = 0.
      
              else
      
                bnkrte = 0.
                degrte = 0.
                degsan = 0.
                degsil = 0.
                degcla = 0.
                deggra = 0.
                bnksan = 0.
                bnksil = 0.
                bnkcla = 0.
                bnkgra = 0.
      
                ch(jrch)%depch = ch(jrch)%depch - depnet
                deg1 = depnet
      
      !! Top 1mm of previously deposited sediment is eroded first - KDW
                if (topcla >= depnet) then
                  ch(jrch)%depclach = ch(jrch)%depclach - depnet
                  deg1cla = depnet
                  degremain = 0.
                else
                  degremain = depnet - topcla
                  deg1cla = topcla
                  ch(jrch)%depclach = ch(jrch)%depclach - topcla
                  if (topsil >= degremain) then
                    ch(jrch)%depsilch = ch(jrch)%depsilch - degremain
                    deg1sil = degremain
                    degremain = 0.
                  else
                    degremain = degremain - topsil
                    deg1sil = topsil
                    ch(jrch)%depsilch = ch(jrch)%depsilch - topsil
                    if (topsag >= degremain) then
                      ch(jrch)%depsagch = ch(jrch)%depsagch - degremain
                      deg1sag = degremain
                      degremain = 0.
                    else
                      degremain = degremain - topsag
                      deg1sag = topsag
                      ch(jrch)%depsagch = ch(jrch)%depsagch - topsag
                      if (topsan >= degremain) then
                        ch(jrch)%depsanch = ch(jrch)%depsanch - degremain
                        deg1san = degremain
                        degremain = 0.
                      else
                        degremain = degremain - topsan
                        deg1san = topsan
                        ch(jrch)%depsanch = ch(jrch)%depsanch - topsan
                      endif
                    endif
                   endif
                endif
      !! Top 1mm of previously deposited sediment is eroded first - KDW          
      
                  if (ch(jrch)%depclach >= degremain) then !changed to "degremain" from "depnet" - KDW 
                  ch(jrch)%depclach = ch(jrch)%depclach - degremain ! - NEW, changed depnet to degremain
                  deg1cla = deg1cla + degremain !added prior deg1, from top 1mm - KDW ! - NEW, changed depnet to degremain
                  degremain = 0.
                else
                  degremain = degremain - ch(jrch)%depclach !changed to "degremain" from "depnet" - KDW
                  deg1cla = deg1cla + ch(jrch)%depclach !added prior deg1, from top 1mm - KDW
                  ch(jrch)%depclach = 0.
                  if (ch(jrch)%depsilch >= degremain) then
                    ch(jrch)%depsilch = ch(jrch)%depsilch - degremain
                    deg1sil = deg1sil + degremain !added prior deg1, from top 1mm - KDW
                    degremain = 0.
                  else
                    degremain = degremain - ch(jrch)%depsilch
                    deg1sil = deg1sil + ch(jrch)%depsilch !added prior deg1, from top 1mm - KDW
                    ch(jrch)%depsilch = 0.
                    if (ch(jrch)%depsagch >= degremain) then
                      ch(jrch)%depsagch = ch(jrch)%depsagch - degremain
                      deg1sag = deg1sag + degremain !added prior deg1, from top 1mm - KDW
                      degremain = 0.
                    else
                      degremain = degremain - ch(jrch)%depsagch
                      deg1sag = deg1sag + ch(jrch)%depsagch !added prior deg1, from top 1mm - KDW
                      ch(jrch)%depsagch = 0.
                      if (ch(jrch)%depsanch >= degremain) then
                        ch(jrch)%depsanch = ch(jrch)%depsanch - degremain
                        deg1san = deg1san + degremain !added prior deg1, from top 1mm - KDW
                        degremain = 0.
                      else
                        degremain = degremain - ch(jrch)%depsanch
                        deg1san = deg1san + ch(jrch)%depsanch !added prior deg1, from top 1mm - KDW
                        ch(jrch)%depsanch = 0.
                        if (ch(jrch)%deplagch >= degremain) then
                          ch(jrch)%deplagch = ch(jrch)%deplagch - degremain
                          deg1lag = degremain
                          degremain = 0.
                        else
                          degremain = degremain - ch(jrch)%deplagch
                          deg1lag = ch(jrch)%deplagch
                          ch(jrch)%deplagch = 0.
                          if (ch(jrch)%depgrach >= degremain) then
                            ch(jrch)%depgrach = ch(jrch)%depgrach - degremain
                            deg1gra = degremain
                            degremain = 0.
                          else
                            degremain = degremain - ch(jrch)%depgrach
                            deg1gra = ch(jrch)%depgrach
                            ch(jrch)%depgrach = 0.
                          endif
                        endif
                      endif
                    endif
                   endif
                endif
      
              endif !depnet > ch(jrch)%depch loop - KDW
       
            end if !depnet .LE. 0 loop - KDW
      
            if (ch(jrch)%depch < 1.e-6) then
              ch(jrch)%depch = 0.
              ch(jrch)%depsanch = 0.
              ch(jrch)%depsilch = 0.
              ch(jrch)%depclach = 0.
              ch(jrch)%depsagch = 0.
              ch(jrch)%deplagch = 0.
              ch(jrch)%depgrach = 0.
            end if
      
      !!	Fall velocity Based on equation 1.36 from SWRRB manual
              vgra = 411.0 * ((2.00)**2.) / (3600.)
              vsan = 411.0 * ((0.20)**2.) / (3600.)
              vsil = 411.0 * ((0.01)**2.) / (3600.)
              vcla = 411.0 * ((0.002)**2.) / (3600.)
              vsag = 411.0 * ((0.03)**2.) / (3600.)
              vlag = 411.0 * ((0.50)**2.) / (3600.)
      
      !!	Deposition calculated based on Einstein Equation
              x = 0.
      
      !!	Gravel deposition
              x = 1.055 * 1000. * ch_hyd(jhyd)%l * vgra / (vc * rchdep)
              if (x > 20.) x = 20.
              pdep = min((1. - exp(-x)), 1.)
              depgra = grain * pdep
      
      !!	sand deposition
              x = 1.055 * 1000. * ch_hyd(jhyd)%l * vsan / (vc * rchdep)
              if (x > 20.) x = 20.
              pdep = min((1. - exp(-x)), 1.)
              depsan = sanin * pdep
      
      !!	Silt deposition
              x = 1.055 * 1000. * ch_hyd(jhyd)%l * vsil / (vc * rchdep)
              if (x > 20.) x = 20.
              pdep = min((1. - exp(-x)), 1.)
      
              depsil = silin * pdep
      
      !!	Clay deposition
              x = 1.055 * 1000. * ch_hyd(jhyd)%l * vcla / (vc * rchdep)
              if (x > 20.) x = 20.
              pdep = min((1. - exp(-x)), 1.)
      
              depcla = clain * pdep
      
      !!	Small aggregates deposition
              x = 1.055 * 1000. * ch_hyd(jhyd)%l * vsag / (vc * rchdep)
              if (x > 20.) x = 20.
              pdep = min((1. - exp(-x)), 1.)
      
              depsag = sagin * pdep
      
      !!	Large aggregates deposition
              x = 1.055 * 1000. * ch_hyd(jhyd)%l * vlag / (vc * rchdep)
              if (x > 20.) x = 20.
              pdep = min((1. - exp(-x)), 1.)
      
              deplag = lagin * pdep
      
              dep = depsan + depsil + depcla + depsag + deplag + depgra
      
      !!    Particles deposited on Floodplain (only silt and clay type particles)
              ch(jrch)%depfp = ch(jrch)%depfp + (depsil + depcla) * fpratio
              ch(jrch)%depsilfp = ch(jrch)%depsilfp + depsil * fpratio
              ch(jrch)%depclafp = ch(jrch)%depclafp + depcla * fpratio
      
      !!    Remaining is deposited in the channel
              ch(jrch)%depch =ch(jrch)%depch + dep - (depsil + depcla)*fpratio
              ch(jrch)%depsilch = ch(jrch)%depsilch + depsil * (1. - fpratio)
              ch(jrch)%depclach = ch(jrch)%depclach + depcla * (1. - fpratio)
              ch(jrch)%depsanch = ch(jrch)%depsanch + depsan
              ch(jrch)%depsagch = ch(jrch)%depsagch + depsag
              ch(jrch)%deplagch = ch(jrch)%deplagch + deplag
              ch(jrch)%depgrach = ch(jrch)%depgrach + depgra
              if (ch(jrch)%depch < 1.e-6) ch(jrch)%depch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depsanch < 1.e-6) ch(jrch)%depsanch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depsilch < 1.e-6) ch(jrch)%depsilch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depclach < 1.e-6) ch(jrch)%depclach = 0.0 ! KDW, 03/26
              if (ch(jrch)%depsagch < 1.e-6) ch(jrch)%depsagch = 0.0 ! KDW, 03/26
              if (ch(jrch)%deplagch < 1.e-6) ch(jrch)%deplagch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depgrach < 1.e-6) ch(jrch)%depgrach = 0.0 ! KDW, 03/26
      
            sedin  = sedin + degrte + bnkrte + deg1    - dep
            grain  = grain + deggra + bnkgra + deg1gra - depgra
            sanin  = sanin + degsan + bnksan + deg1san - depsan
            silin  = silin + degsil + bnksil + deg1sil - depsil
            clain  = clain + degcla + bnkcla + deg1cla - depcla
            sagin  = sagin + deg1sag - depsag
            lagin  = lagin + deg1lag - deplag
      
          if (sedin  < 1.e-6) then
              sedin = 0.
              sanin = 0.
              silin = 0.
              clain = 0.
              sagin = 0.
              lagin = 0.
              grain = 0.
            end if
      
            outfract = rtwtr / qdin
            if (outfract > 1.) outfract = 1.
      
            sedrch =  sedin * outfract
            rch_san = sanin * outfract
            rch_sil = silin * outfract
            rch_cla = clain * outfract
            rch_sag = sagin * outfract
            rch_lag = lagin * outfract
            rch_gra = grain * outfract
      
          if (sedrch  < 1.e-6) then
              sedrch = 0.
              rch_san = 0.
              rch_sil = 0.
              rch_cla = 0.
              rch_sag = 0.
              rch_lag = 0.
              rch_gra = 0.
            end if
      
            ch(jrch)%sedst = sedin - sedrch
            ch(jrch)%sanst = sanin - rch_san
            ch(jrch)%silst = silin - rch_sil
            ch(jrch)%clast = clain - rch_cla
            ch(jrch)%sagst = sagin - rch_sag
            ch(jrch)%lagst = lagin - rch_lag
            ch(jrch)%grast = grain - rch_gra
      
          if (ch(jrch)%sedst < 1.e-6) then
              ch(jrch)%sedst = 0.
              ch(jrch)%sanst = 0.
              ch(jrch)%silst = 0.
              ch(jrch)%clast = 0.
              ch(jrch)%sagst = 0.
              ch(jrch)%lagst = 0.
              ch(jrch)%grast = 0.
            end if
      
      !!    Mass balance tests
      !!      ambalsed = sedinorg + degrte + bnkrte + deg1 - dep - sedrch       
      !!     &            - sedst(jrch))
      !!      ambalsed = depch(jrch) - depsanch(jrch)-depsilch(jrch)            
      !!     &-ch(jrch)%depclach-depsagch(jrch)-deplagch(jrch)-depgrach(jrch)
      !!      ambalsed = sedst(jrch) - sanst(jrch)-silst(jrch)-clast(jrch)      
      !!     &-sagst(jrch)-lagst(jrch)-grast(jrch)
      !!      ambalsed = (sedin-sanin-silin-clain-sagin-lagin-grain)/sedin
      !!      ambalsed = sedrch-rch_san-rch_sil-rch_cla-rch_sag-rch_lag-rch_gra
      !!      if (abs(ambalsed) .gt. 1e-3) write (*,*) time%day,jrch,ambalsed,sedrch
            
      !!    Deposition during the previous time step
            ch(jrch)%depprch = ch(jrch)%depch  !! Channel
            ch(jrch)%depprfp = ch(jrch)%depfp  !! Flood plain
            
      
      
      !! compute changes in channel dimensions
            if (bsn_cc%deg == 1) then
              depdeg = 0.
              depdeg = ch_hyd(jhyd)%d - ch(jrch)%di
              if (depdeg < ch(jrch)%si * ch(jrch)%li * 1000.) then
                if (qdin > 1400000.) then
                  dot = 0.
                  dot = 358.6 * rchdep * ch_hyd(jhyd)%s * ch_sed(jsed)%cov1 
                  dat2 = 1.
                  dat2 =  dat2 * dot
                  ch_hyd(jhyd)%d = ch_hyd(jhyd)%d + dat2
                  ch_hyd(jhyd)%w = ch_hyd(jhyd)%wdr * ch_hyd(jhyd)%d
                  ch_hyd(jhyd)%s = ch_hyd(jhyd)%s - dat2 / (ch_hyd(jhyd)%l *      &
                                                     1000.)
                  ch_hyd(jhyd)%s = Max(.0001, ch_hyd(jhyd)%s)
                  call ch_ttcoef(jrch)
                endif
              endif
            endif
      
            else ! for if qdin > .01 loop
      
              dep = sedin ; depsan = sanin ; depsil = silin ; depcla = clain
              depsag = sagin ; deplag = lagin ; depgra = grain
              degcla = 0. ; deg1cla = 0. ; degsil = 0. ; deg1sil = 0. 
              deg1sag = 0. ; degsan = 0. ; deg1san = 0. 
              ch(jrch)%sedst = 0. ; ch(jrch)%sanst = 0. ; ch(jrch)%silst = 0. 
              ch(jrch)%clast = 0. ; ch(jrch)%sagst = 0. ; ch(jrch)%lagst = 0. 
              ch(jrch)%grast = 0. ;
              ch(jrch)%depch = ch(jrch)%depch + dep
              ch(jrch)%depsilch = ch(jrch)%depsilch + depsil
              ch(jrch)%depclach = ch(jrch)%depclach + depcla
              ch(jrch)%depsanch = ch(jrch)%depsanch + depsan
              ch(jrch)%depsagch = ch(jrch)%depsagch + depsag
              ch(jrch)%deplagch = ch(jrch)%deplagch + deplag
              ch(jrch)%depgrach = ch(jrch)%depgrach + depgra
              if (ch(jrch)%depch < 1.e-6) ch(jrch)%depch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depsanch < 1.e-6) ch(jrch)%depsanch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depsilch < 1.e-6) ch(jrch)%depsilch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depclach < 1.e-6) ch(jrch)%depclach = 0.0 ! KDW, 03/26
              if (ch(jrch)%depsagch < 1.e-6) ch(jrch)%depsagch = 0.0 ! KDW, 03/26
              if (ch(jrch)%deplagch < 1.e-6) ch(jrch)%deplagch = 0.0 ! KDW, 03/26
              if (ch(jrch)%depgrach < 1.e-6) ch(jrch)%depgrach = 0.0 ! KDW, 03/26
              sedrch = 0. ; rch_san = 0. ; rch_sil = 0. ; rch_cla = 0. 
              rch_sag = 0. ; rch_lag = 0. ; rch_gra = 0. 
      
      
            end if !! end of qdin > 0.01 loop
      
          !end if  !! end of rtwtr and rchdep > 0 loop !!commented out by KDW, enough to have qdin > .01 loop
      
      !! Phosphorus movement associated with sediment erosion - KDW
            
      ! What are water volumes in layers (and in newly deposited sediment)
          conv_mass2vol = ch(jrch)%bed_por / ch_sed(jsed)%bed_bd ! convert from mass of sediment to volume of h2o in pore space
          vol_h2o_aer = ch(jrch)%bed_por * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * .001) ! m^3
          vol_h2o_anaer = ch(jrch)%bed_por * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * .1)! m^3
          vol_h2o_newdep = ch(jrch)%bed_por * (fines_newdep / ch_sed(jsed)%bed_bd)
      ! How much P mass of each form is in each layer and the wtclm, initially
          if (ch(jrch)%aer_orgp < 1.e-6) ch(jrch)%aer_orgp = 0.0 ! KDW, 03/26
          if (ch(jrch)%aer_disp < 1.e-6) ch(jrch)%aer_disp = 0.0 ! KDW, 03/26
          if (ch(jrch)%aer_minpa < 1.e-6) ch(jrch)%aer_minpa = 0.0 ! KDW, 03/26
          if (ch(jrch)%aer_minps < 1.e-6) ch(jrch)%aer_minps = 0.0 ! KDW, 03/26
          if (ch(jrch)%anaer_orgp < 1.e-6) ch(jrch)%anaer_orgp = 0.0 ! KDW, 03/26
          if (ch(jrch)%anaer_disp < 1.e-6) ch(jrch)%anaer_disp = 0.0 ! KDW, 03/26
          if (ch(jrch)%anaer_minpa < 1.e-6) ch(jrch)%anaer_minpa = 0.0 ! KDW, 03/26
          if (ch(jrch)%anaer_minps < 1.e-6) ch(jrch)%anaer_minps = 0.0 ! KDW, 03/26
          if (ch(jrch)%organicp < 1.e-6) ch(jrch)%organicp = 0.0 ! KDW, 03/26
          if (ch(jrch)%disolvp < 1.e-6) ch(jrch)%disolvp = 0.0 ! KDW, 03/26
          if (ch(jrch)%minpa < 1.e-6) ch(jrch)%minpa = 0.0 ! KDW, 03/26
          if (ch(jrch)%minps < 1.e-6) ch(jrch)%minps = 0.0 ! KDW, 03/26
          minpa_aer = ch(jrch)%aer_minpa * aer_mass / 1000 ! kg
          minps_aer = ch(jrch)%aer_minps * aer_mass / 1000 ! kg
          minpa_anaer = ch(jrch)%anaer_minpa * anaer_mass / 1000 ! kg
          minps_anaer = ch(jrch)%anaer_minps * anaer_mass / 1000 ! kg
          minpa_wtclm = ch(jrch)%minpa * qdin / 1000 ! kg
          minps_wtclm = ch(jrch)%minps * qdin / 1000 ! kg
          orgp_aer = ch(jrch)%aer_orgp * aer_mass / 1000 ! kg
          disp_aer = ch(jrch)%aer_disp * vol_h2o_aer / 1000 ! kg 
          orgp_anaer = ch(jrch)%anaer_orgp * anaer_mass / 1000 ! kg
          disp_anaer = ch(jrch)%anaer_disp * vol_h2o_anaer / 1000 ! kg 
          orgp_wtclm = ch(jrch)%organicp * qdin / 1000 ! kg
          disp_wtclm = ch(jrch)%disolvp * qdin / 1000 ! kg 
      ! how much "fines" erode and deposit, and what is the surface area
          cla_ero = degcla + deg1cla
          sil_ero = degsil + deg1sil
          sag_ero = deg1sag
          san_ero = degsan + deg1san
          bed_ero = cla_ero + sil_ero + sag_ero + san_ero
          fines_newdep = depcla + depsil + depsag + depsan + 0.0001
          bed_netero = bed_ero - fines_newdep
          sfarea_ero  = cla_ero + sil_ero * .2 + sag_ero * (1/15) + san_ero * .01
          sfarea_dep  = depcla + depsil * .2 + depsag * (1/15) + depsan * .01    
          
      ! How much of each layer (aerobic and anaerobic) is left after erosion, how much is previously deposited "fines", and surface areas
          if (bed_ero <= aer_mass) then
                rmn_aer_mass = aer_mass - bed_ero
                rmn_anaer_mass = anaer_mass
                rmn_aer_sfar = sfarea_aer - sfarea_ero
                rmn_anaer_sfar = sfarea_anaer
          else
                rmn_aer_mass = 0
                rmn_aer_sfar = 0
                if (bed_ero <= aer_mass + anaer_mass) then
                      rmn_anaer_mass = anaer_mass + aer_mass - bed_ero
                      rmn_anaer_sfar = sfarea_anaer + sfarea_aer - sfarea_ero
                else
                      rmn_anaer_mass = 0
                      rmn_anaer_sfar = 0
                end if
          end if
      ! How much phosphorus is associated with the different pools of eroding or depositing sediment
          if (bed_ero < aer_mass) then
                percmass_ero_aer = bed_ero / aer_mass ! percentage (by mass) of aerobic layer that is eroded
                percmass_ero_anaer = 0.
                percsfar_ero_aer = sfarea_ero/sfarea_aer
                percsfar_ero_anaer = 0
                else
                      if (bed_ero < aer_mass + anaer_mass) then
                      percmass_ero_aer = 1 ! percentage (by mass) of aerobic layer that is eroded
                      percmass_ero_anaer = (bed_ero - aer_mass) / anaer_mass
                      percsfar_ero_aer = 1
                      percsfar_ero_anaer = (sfarea_ero - sfarea_aer) /sfarea_anaer
                      else
                            percmass_ero_aer = 1 ! percentage (by mass) of aerobic layer that is eroded
                            percmass_ero_anaer = 1
                            percsfar_ero_aer = 1
                            percsfar_ero_anaer = 1
                      endif
          end if
      
          orgp_ero_aer = orgp_aer * percmass_ero_aer
          disp_ero_aer = disp_aer * percmass_ero_aer
          minpa_ero_aer = percsfar_ero_aer * minpa_aer
          minps_ero_aer = percsfar_ero_aer * minps_aer
          rmn_orgp_aer = orgp_aer - orgp_ero_aer
          rmn_disp_aer = disp_aer - disp_ero_aer
          rmn_minpa_aer = minpa_aer - minpa_ero_aer
          rmn_minps_aer = minps_aer - minps_ero_aer
      
          orgp_ero_anaer = orgp_anaer * percmass_ero_anaer
          disp_ero_anaer = disp_anaer * percmass_ero_anaer
          minpa_ero_anaer = percsfar_ero_anaer * minpa_anaer
          minps_ero_anaer = percsfar_ero_anaer * minps_anaer
          rmn_orgp_anaer = orgp_anaer - orgp_ero_anaer
          rmn_disp_anaer = disp_anaer - disp_ero_anaer
          rmn_minpa_anaer = minpa_anaer - minpa_ero_anaer
          rmn_minps_anaer = minps_anaer - minps_ero_anaer
      
          if (fines_wtclm > 0) then
              percmass_dep = fines_newdep / fines_wtclm ! percentage (by mass) of aerobic layer that is eroded
          else
              percmass_dep = 1
          endif
          if (sfarea_wtclm > 0) then
              percsfar_dep = sfarea_dep / sfarea_wtclm
          else
              percsfar_dep = 1
          endif
          if (qdin > 0) then
              percvol_dep = vol_h2o_newdep / qdin
          else
              percvol_dep = 1
          endif
      
          orgp_dep = percmass_dep * orgp_wtclm
          disp_dep = percvol_dep * disp_wtclm
          minpa_dep = percsfar_dep * minpa_wtclm
          minps_dep = percsfar_dep * minps_wtclm
          
      ! How much phosphorus is in each layer after erosion and deposition
          if (fines_newdep >= aer_mass + anaer_mass) then
                new_minpa_aer = minpa_dep * (aer_mass / fines_newdep)
                new_minps_aer = minps_dep * (aer_mass / fines_newdep)
                new_orgp_aer = orgp_dep * (aer_mass / fines_newdep)
                new_disp_aer = disp_dep * (aer_mass / fines_newdep)
      
                new_minpa_anaer = minpa_dep * (anaer_mass / fines_newdep)
                new_minps_anaer = minps_dep * (anaer_mass / fines_newdep)
                new_orgp_anaer = orgp_dep * (anaer_mass / fines_newdep)
                new_disp_anaer = disp_dep * (anaer_mass / fines_newdep)
          else
                if (fines_newdep >= aer_mass) then ! new aerobic layer is entirely new deposition
                      new_minpa_aer = minpa_dep * (aer_mass / fines_newdep)
                      new_minps_aer = minps_dep * (aer_mass / fines_newdep)
                      new_orgp_aer = orgp_dep * (aer_mass / fines_newdep)
                      new_disp_aer = disp_dep * (aer_mass / fines_newdep)
      
                      new_minpa_anaer = minpa_dep * (1 - (aer_mass / fines_newdep))
                      new_minps_anaer = minps_dep * (1 - (aer_mass / fines_newdep))
                      new_orgp_anaer = orgp_dep * (1 - (aer_mass / fines_newdep))
                      new_disp_anaer = disp_dep * (1 - (aer_mass / fines_newdep))
      
                    if (fines_newdep + rmn_aer_mass >= aer_mass + anaer_mass) then
                            mass_from_aer = (anaer_mass + aer_mass) - fines_newdep
                            new_minpa_anaer = new_minpa_anaer + rmn_minpa_aer * (mass_from_aer / rmn_aer_mass)
                            new_minps_anaer = new_minps_anaer + rmn_minps_aer * (mass_from_aer / rmn_aer_mass)
                            new_orgp_anaer = new_orgp_anaer + rmn_orgp_aer * (mass_from_aer / rmn_aer_mass)
                            new_disp_anaer = new_disp_anaer + rmn_disp_aer * (mass_from_aer / rmn_aer_mass) 
                    else
                            new_minpa_anaer = new_minpa_anaer + rmn_minpa_aer
                            new_minps_anaer = new_minps_anaer + rmn_minps_aer
                            new_orgp_anaer = new_orgp_anaer + rmn_orgp_aer
                            new_disp_anaer = new_disp_anaer + rmn_disp_aer
                            if (fines_newdep + rmn_aer_mass + rmn_anaer_mass >= aer_mass + anaer_mass) then 
                                mass_still_anaer = (anaer_mass + aer_mass) - (fines_newdep + rmn_aer_mass)
                                new_minpa_anaer = new_minpa_anaer + rmn_minpa_anaer * (mass_still_anaer / rmn_anaer_mass)
                                new_minps_anaer = new_minps_anaer + rmn_minps_anaer * (mass_still_anaer / rmn_anaer_mass)
                                new_orgp_anaer = new_orgp_anaer + rmn_orgp_aer * (mass_still_anaer / rmn_anaer_mass)
                                new_disp_anaer = new_disp_anaer + rmn_disp_aer * (mass_still_anaer / rmn_anaer_mass)
                            else ! this final stanza is net erosion
                                    ext_mass = anaer_mass + aer_mass - (fines_newdep + rmn_aer_mass + rmn_anaer_mass)
                                    new_minpa_anaer = new_minpa_anaer + rmn_minpa_anaer + ext_mass * ch(jrch)%ext_minpa / 1000 ! NEW - divide by 1000!
                                    new_minps_anaer = new_minps_anaer + rmn_minps_anaer + ext_mass * ch(jrch)%ext_minps / 1000
                                    new_orgp_anaer = new_orgp_anaer + rmn_orgp_anaer + ext_mass * ch(jrch)%ext_orgp / 1000
                                    new_disp_anaer = new_disp_anaer + rmn_disp_anaer + ext_mass * conv_mass2vol * &
                                        ch(jrch)%ext_disp / 1000
                            end if
                    end if
                else
                      new_minpa_aer = minpa_dep
                      new_minps_aer = minps_dep
                      new_orgp_aer = orgp_dep
                      new_disp_aer = disp_dep
                    if (fines_newdep + rmn_aer_mass >= aer_mass) then !!!! THIS IS WHERE PREVIOUS CODE WAS WRONG - NEW !!!!
                            mass_still_aer = aer_mass - fines_newdep
                            mass_from_aer = rmn_aer_mass - mass_still_aer
                            
                            new_minpa_aer = new_minpa_aer + rmn_minpa_aer * (mass_still_aer / rmn_aer_mass)
                            new_minps_aer = new_minps_aer + rmn_minps_aer * (mass_still_aer / rmn_aer_mass)
                            new_orgp_aer = new_orgp_aer + rmn_orgp_aer * (mass_still_aer / rmn_aer_mass)
                            new_disp_aer = new_disp_aer + rmn_disp_aer * (mass_still_aer / rmn_aer_mass)
      
                            new_minpa_anaer = rmn_minpa_aer * (mass_from_aer / rmn_aer_mass)
                            new_minps_anaer = rmn_minps_aer * (mass_from_aer / rmn_aer_mass)
                            new_orgp_anaer = rmn_orgp_aer * (mass_from_aer / rmn_aer_mass)
                            new_disp_anaer = rmn_disp_aer * (mass_from_aer / rmn_aer_mass)
      
                            if (fines_newdep + rmn_aer_mass + rmn_anaer_mass >= aer_mass + anaer_mass) then 
                                mass_still_anaer = (anaer_mass + aer_mass) - (fines_newdep + rmn_aer_mass)
                                new_minpa_anaer = new_minpa_anaer + rmn_minpa_anaer * (mass_still_anaer / rmn_anaer_mass)
                                new_minps_anaer = new_minps_anaer + rmn_minps_anaer * (mass_still_anaer / rmn_anaer_mass)
                                new_orgp_anaer = new_orgp_anaer + rmn_orgp_anaer * (mass_still_anaer / rmn_anaer_mass)
                                new_disp_anaer = new_disp_anaer + rmn_disp_anaer * (mass_still_anaer / rmn_anaer_mass)
                            else
                                ext_mass = anaer_mass + aer_mass - (fines_newdep + rmn_aer_mass + rmn_anaer_mass)
                                new_minpa_anaer = new_minpa_anaer + rmn_minpa_anaer + ext_mass * ch(jrch)%ext_minpa / 1000
                                new_minps_anaer = new_minps_anaer + rmn_minps_anaer + ext_mass * ch(jrch)%ext_minps / 1000
                                new_orgp_anaer = new_orgp_anaer + rmn_orgp_anaer + ext_mass * ch(jrch)%ext_orgp / 1000
                                new_disp_anaer = new_disp_anaer + rmn_disp_anaer + ext_mass * conv_mass2vol * &
                                    ch(jrch)%ext_disp / 1000
      
                            endif                      
                    else  ! below is net erosion, rather than net deposition (above)
                            new_minpa_aer = new_minpa_aer + rmn_minpa_aer
                            new_minps_aer = new_minps_aer + rmn_minps_aer
                            new_orgp_aer = new_orgp_aer + rmn_orgp_aer
                            new_disp_aer = new_disp_aer + rmn_disp_aer
      
                            if (fines_newdep + rmn_aer_mass + rmn_anaer_mass >= aer_mass) then ! sed moving from to aerobic is from anaerbic
                                    mass_from_anaer = aer_mass - (fines_newdep + rmn_aer_mass)
                                    mass_still_anaer = rmn_anaer_mass - mass_from_anaer
                                    ext_mass = anaer_mass + aer_mass - (fines_newdep + rmn_aer_mass + rmn_anaer_mass)
      
                                    new_minpa_aer = new_minpa_aer + rmn_minpa_anaer * (mass_from_anaer / rmn_anaer_mass)
                                    new_minps_aer = new_minps_aer + rmn_minps_anaer * (mass_from_anaer / rmn_anaer_mass)
                                    new_orgp_aer = new_orgp_aer + rmn_orgp_anaer * (mass_from_anaer / rmn_anaer_mass)
                                    new_disp_aer = new_disp_aer + rmn_disp_anaer * (mass_from_anaer / rmn_anaer_mass)
      
                                    new_minpa_anaer = rmn_minpa_anaer * (mass_still_anaer / rmn_anaer_mass) + &
                                        ext_mass * ch(jrch)%ext_minpa / 1000
                                    new_minps_anaer = rmn_minps_anaer * (mass_still_anaer / rmn_anaer_mass) + &
                                        ext_mass * ch(jrch)%ext_minps / 1000
                                    new_orgp_anaer = rmn_orgp_anaer * (mass_still_anaer / rmn_anaer_mass) + &
                                        ext_mass * ch(jrch)%ext_orgp / 1000
                                    new_disp_anaer = rmn_disp_anaer * (mass_still_anaer / rmn_anaer_mass) + &
                                        ext_mass * ch(jrch)%ext_disp / 1000
                            else ! sediment entering aerobic layer from external of system
                                    ext_mass = aer_mass - (fines_newdep + rmn_aer_mass + rmn_anaer_mass)
                                    new_minpa_aer = new_minpa_aer + rmn_minpa_anaer + ext_mass * ch(jrch)%ext_minpa / 1000
                                    new_minps_aer = new_minps_aer + rmn_minps_anaer + ext_mass * ch(jrch)%ext_minps / 1000
                                    new_orgp_aer = new_orgp_aer + rmn_orgp_anaer + ext_mass * ch(jrch)%ext_orgp / 1000
                                    new_disp_aer = new_disp_aer + rmn_disp_anaer + ext_mass * conv_mass2vol * &
                                        ch(jrch)%ext_disp / 1000
      
                                    new_minpa_anaer = anaer_mass * ch(jrch)%ext_minpa / 1000
                                    new_minps_anaer = anaer_mass * ch(jrch)%ext_minps / 1000
                                    new_orgp_anaer = anaer_mass * ch(jrch)%ext_orgp / 1000
                                    new_disp_anaer = anaer_mass * conv_mass2vol * ch(jrch)%ext_disp / 1000
                            end if
                    end if
                end if
          end if
      
          ch(jrch)%aer_orgp = 1000 * new_orgp_aer / aer_mass
          ch(jrch)%aer_disp = 1000 * new_disp_aer / aer_mass
          ch(jrch)%aer_minpa = 1000 * new_minpa_aer / aer_mass
          ch(jrch)%aer_minps = 1000 * new_minps_aer / aer_mass
      
          ch(jrch)%anaer_orgp = 1000 * new_orgp_anaer / anaer_mass
          ch(jrch)%anaer_disp = 1000 * new_disp_anaer / anaer_mass
          ch(jrch)%anaer_minpa = 1000 * new_minpa_anaer / anaer_mass
          ch(jrch)%anaer_minps = 1000 * new_minps_anaer / anaer_mass
      
          if (qdin > 0.001) then
              ch(jrch)%organicp = ch(jrch)%organicp + 1000 * (orgp_ero_aer + orgp_ero_anaer - orgp_dep) / qdin
              ch(jrch)%disolvp = ch(jrch)%disolvp + 1000 * (disp_ero_aer + disp_ero_anaer - disp_dep) / qdin
              ch(jrch)%minpa = ch(jrch)%minpa + 1000 * (minpa_ero_aer + minpa_ero_anaer - minpa_dep) / qdin
              ch(jrch)%minps =ch(jrch)%minps + 1000 * (minps_ero_aer + minps_ero_anaer - minps_dep) / qdin
          else
              ch(jrch)%organicp = 0.
              ch(jrch)%disolvp = 0.
              ch(jrch)%minpa = 0.
              ch(jrch)%minps = 0.
          endif
          
          if (ch(jrch)%aer_orgp < 1.e-6) ch(jrch)%aer_orgp = 0.0 ! KDW 03/24
          if (ch(jrch)%aer_disp < 1.e-6) ch(jrch)%aer_disp = 0.0 ! KDW 03/24
          if (ch(jrch)%aer_minpa < 1.e-6) ch(jrch)%aer_minpa = 0.0 ! KDW 03/24
          if (ch(jrch)%aer_minps < 1.e-6) ch(jrch)%aer_minps = 0.0 ! KDW 03/24
          if (ch(jrch)%anaer_orgp < 1.e-6) ch(jrch)%anaer_orgp = 0.0 ! KDW 03/24
          if (ch(jrch)%anaer_disp < 1.e-6) ch(jrch)%anaer_disp = 0.0 ! KDW 03/24
          if (ch(jrch)%anaer_minpa < 1.e-6) ch(jrch)%anaer_minpa = 0.0 ! KDW 03/24
          if (ch(jrch)%anaer_minps < 1.e-6) ch(jrch)%anaer_minps = 0.0 ! KDW 03/24
          if (ch(jrch)%organicp < 1.e-6) ch(jrch)%organicp = 0.0 ! KDW 03/24
          if (ch(jrch)%disolvp < 1.e-6) ch(jrch)%disolvp = 0.0 ! KDW 03/24
          if (ch(jrch)%minpa < 1.e-6) ch(jrch)%minpa = 0.0 ! KDW 03/24
          if (ch(jrch)%minps < 1.e-6) ch(jrch)%minps = 0.0 ! KDW 03/24
      
          ch_d(jrch)%phos_dep = orgp_dep + disp_dep + minpa_dep + minps_dep ! KDW, 06/04
          ch_d(jrch)%phos_ero = orgp_ero_aer + disp_ero_aer + minpa_ero_aer + minps_ero_aer + &
            orgp_ero_anaer + disp_ero_anaer + minpa_ero_anaer + minps_ero_anaer! KDW, 06/04
      
      !! Phosphorus movement associated with sediment erosion - KDW
      
      
      !!    Organic nitrogen and Organic Phosphorus contribution from channel erosion
      !!    Only bank erosion is assumed to contribute to channel erosion
         !!   ch_orgn(jrch) = bnkrte * ch_nut(jnut)%onco * 1000.
         !!   ch_orgp(jrch) = bnkrte * ch_nut(jnut)%opco * 1000.
            ch(jrch)%orgn = bnkrte * ch_nut(jnut)%onco / 1000.
            ch(jrch)%orgp = bnkrte * ch_nut(jnut)%opco / 1000. ! this isn't added to the reach store... KDW
      
            return
            end subroutine ch_rtsed_bagnold2