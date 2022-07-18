      subroutine ch_rtsed2
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine routes sediment from subbasin to basin outlets
!!    deposition is based on fall velocity and degradation on stream

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!                               |1 no vegetative cover on channel
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

      use basin_module
      use channel_data_module
      use channel_module
      use hydrograph_module, only : ob, jrch, icmd
      use time_module
      use channel_velocity_module
             
      implicit none

      real :: qdin             !m^3 H2O       |water in reach during time step
      real :: sedin            !units         |description
      real :: vc               !m/s           |flow velocity in reach
      real :: cyin             !units         |description
      real :: cych             !units         |description
      real :: depnet           !metric tons   |
      real :: deg1             !units         |description
      real :: deg2             !units         |description
      real :: dep              !metric tons   |sediment deposited on river bottom
      real :: depdeg           !m             |depth of degradation/deposition from original
      real :: dot              !mm            |actual depth from impermeable layer to water level
                               !              |above drain during subsurface irrigation
      real :: outfract         !units         |description
      real :: deg              !metric tons   |sediment reentrained in water by channel
                               !              |degradation
      real :: sedinorg         !units         |description
      real :: tbase            !none          |flow duration (fraction of 1 hr)
      real :: dat2             !m             |change in channel depth during time step
      real :: aer_mass, anaer_mass              !tons (Mg)   |sediment mass in layer - KDW
      real :: conv_mass2vol, conv_conc              !units   |conversions necessary between pools as needed - KDW
      real :: minpa_wtclm, minps_wtclm, orgp_wtclm, disp_wtclm              !kg   |kg P by form in layer - KDW
      real :: minpa_aer, minps_aer, orgp_aer, disp_aer              !kg   |kg P by form in layer - KDW
      real :: minpa_anaer, minps_anaer, orgp_anaer, disp_anaer              !kg   |kg P by form in layer - KDW
      real :: mass_ero_aer, mass_ero_anaer, mass_dep_aer, mass_dep_anaer              !kg   |kg P eroded or deposited by layer - KDW
      real :: mass_anaer2aer, mass_aer2anaer              !kg   |kg P moving with layer boundary layer - KDW
      real :: xx, yy, zz, xyz

      sedin = 0.0

      if (rtwtr > 0. .and. rchdep > 0.) then

!! initialize water in reach during time step
      qdin = 0.
      qdin = rtwtr + ch(jrch)%rchstor

!! do not perform sediment routing if no water in reach
      if (qdin > 0.01) then

!! initialize sediment in reach during time step
      sedin = 0.
      sedin = ob(icmd)%hin%sed  + ch(jrch)%sedst
      sedinorg = sedin
!! initialize reach peak runoff rate
      peakr = bsn_prm%prf * sdti

!! calculate flow velocity
      vc = 0.
      if (rchdep < .010) then
        vc = 0.01
      else
        vc = peakr / (rcharea + .1)  !dont merge
      end if
      if (vc > 5.) vc = 5.

      tbase = 0.
      tbase = ch_hyd(jhyd)%l * 1000. / (3600. * 24. * vc + .1)  !dont merge

      if (tbase > 1.) tbase = 1.


!! JIMMY"S NEW IMPROVED METHOD for sediment transport
      cyin = 0.
      cych = 0.
      depnet = 0.
	  deg = 0.
      deg1 = 0.
	  deg2 = 0.
      dep = 0.
      aer_mass = 0. ; anaer_mass = 0.
      conv_mass2vol = 0. ; conv_conc = 0. 
      minpa_wtclm = 0. ; minps_wtclm = 0. ; orgp_wtclm = 0. ; disp_wtclm = 0. 
      minpa_aer = 0. ; minps_aer = 0. ; orgp_aer = 0. ; disp_aer = 0.
      minpa_anaer = 0. ; minps_anaer = 0. ; orgp_anaer = 0. ; disp_anaer = 0.
      mass_ero_aer = 0. ; mass_ero_anaer = 0. ; mass_dep_aer = 0. ; mass_dep_anaer = 0. 
      mass_anaer2aer = 0. ; mass_aer2anaer = 0. 

      cyin = sedin / qdin
      cych = bsn_prm%spcon * vc ** bsn_prm%spexp
      depnet = qdin * (cych - cyin)
	if(abs(depnet) < 1.e-6) depnet = 0.
      !if (vc < bsn_prm%vcrit) depnet = 0.

!!  tbase is multiplied so that erosion is proportional to the traveltime, 
!!  which is directly related to the length of the channel
!!  Otherwise for the same discharge rate and sediment deficit
!!  the model will erode more sediment per unit length of channel 
!!  from a small channel than a larger channel. Modification made by Balaji Narasimhan

!!    Embeddment of stored sediment into "permanent" streambed - KDW
      ch(jrch)%depch = ch(jrch)%depch * exp(-bsn_prm%sbdk)
!!    Embeddment of stored sediment into "permanent" streambed - KDW 

      if (depnet > 1.e-6) then
        deg = depnet * tbase
	  !! First the deposited material will be degraded before channel bed
	  if (deg >= ch(jrch)%depch) then
	    deg1 = ch(jrch)%depch
        deg2 = (deg - deg1) * ch_sed(jsed)%erod(time%mo) * ch_sed(jsed)%cov2
	  else
	    deg1 = deg
	    deg2 = 0.
	  endif
        dep = 0.
      else
        dep = -depnet * tbase
        deg = 0.
	  deg1 = 0.
	  deg2 = 0.
      endif

!!    Phosphorus movement - KDW
      xx = ch_vel(jrch)%wid_btm
      yy = ch_hyd(jhyd)%l
      zz = ch_sed(jsed)%bed_bd
      aer_mass = .001 * 1 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd) ! depth of aer layer is 1
      anaer_mass = .1 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd)
      conv_mass2vol = 1000 * ch(jrch)%bed_por / ch_sed(jsed)%bed_bd ! convert from mass of sediment (Mg) to volume of h2o (L) in pore space
      conv_conc = (qdin / sedin) ! convert concentration from mg P/L to g P/Mg sediment

      minpa_wtclm = ch(jrch)%minpa * qdin ! g
      minps_wtclm = ch(jrch)%minps * qdin ! g
      orgp_wtclm = ch(jrch)%organicp * qdin ! g
      disp_wtclm = ch(jrch)%disolvp * qdin * 1000 ! mg 

      if (depnet > aer_mass + anaer_mass) then
            mass_ero_aer = aer_mass
            mass_ero_anaer = anaer_mass
            minpa_aer = aer_mass * ch(jrch)%ext_minpa ! g for all
            minps_aer = aer_mass * ch(jrch)%ext_minps
            orgp_aer = aer_mass * ch(jrch)%ext_orgp
            disp_aer = aer_mass * conv_mass2vol * ch(jrch)%ext_disp ! mg
            minpa_anaer = anaer_mass * ch(jrch)%ext_minpa
            minps_anaer = anaer_mass * ch(jrch)%ext_minps
            orgp_anaer = anaer_mass * ch(jrch)%ext_orgp
            disp_anaer = anaer_mass * conv_mass2vol * ch(jrch)%ext_disp
      else
            if (depnet > anaer_mass) then
                  mass_ero_aer = aer_mass
                  mass_ero_anaer = depnet - mass_ero_aer
                  mass_anaer2aer = anaer_mass - mass_ero_anaer
                  minpa_aer = mass_anaer2aer * ch(jrch)%anaer_minpa + (aer_mass - mass_anaer2aer) * ch(jrch)%ext_minpa
                  minps_aer = mass_anaer2aer * ch(jrch)%anaer_minps + (aer_mass - mass_anaer2aer) * ch(jrch)%ext_minps
                  orgp_aer = mass_anaer2aer * ch(jrch)%anaer_orgp + (aer_mass - mass_anaer2aer) * ch(jrch)%ext_orgp
                  disp_aer = mass_anaer2aer * conv_mass2vol * ch(jrch)%anaer_disp + & 
                        (aer_mass - mass_anaer2aer) * conv_mass2vol * ch(jrch)%ext_disp 
                  minpa_anaer = anaer_mass * ch(jrch)%ext_minpa
                  minps_anaer = anaer_mass * ch(jrch)%ext_minps
                  orgp_anaer = anaer_mass * ch(jrch)%ext_orgp
                  disp_anaer = anaer_mass * conv_mass2vol * ch(jrch)%ext_disp
            else
                  if (depnet > aer_mass) then
                        mass_ero_aer = aer_mass
                        mass_ero_anaer = depnet - mass_ero_aer
                        mass_anaer2aer = aer_mass
                        minpa_aer = mass_anaer2aer * ch(jrch)%anaer_minpa
                        minps_aer = mass_anaer2aer * ch(jrch)%anaer_minps
                        orgp_aer = mass_anaer2aer * ch(jrch)%anaer_orgp
                        disp_aer = mass_anaer2aer * conv_mass2vol * ch(jrch)%anaer_disp
                        minpa_anaer = (anaer_mass - depnet) * ch(jrch)%anaer_minpa + depnet * ch(jrch)%ext_minpa
                        minps_anaer = (anaer_mass - depnet) * ch(jrch)%anaer_minps + depnet * ch(jrch)%ext_minps
                        orgp_anaer = (anaer_mass - depnet) * ch(jrch)%anaer_orgp + depnet * ch(jrch)%ext_orgp
                        disp_anaer = (anaer_mass - depnet) * conv_mass2vol * ch(jrch)%anaer_disp + & 
                              depnet * conv_mass2vol * ch(jrch)%ext_disp
                  else
                        if (depnet > 1.e-6) then
                              mass_ero_aer = depnet
                              mass_ero_anaer = 0
                              mass_anaer2aer = depnet
                              minpa_aer = mass_anaer2aer * ch(jrch)%anaer_minpa + (aer_mass - mass_ero_aer) * ch(jrch)%aer_minpa
                              minps_aer = mass_anaer2aer * ch(jrch)%anaer_minps + (aer_mass - mass_ero_aer) * ch(jrch)%aer_minps
                              orgp_aer = mass_anaer2aer * ch(jrch)%anaer_orgp + (aer_mass - mass_ero_aer) * ch(jrch)%aer_orgp
                              disp_aer = mass_anaer2aer * conv_mass2vol * ch(jrch)%anaer_disp + & 
                                    (aer_mass - mass_ero_aer) * conv_mass2vol * ch(jrch)%aer_disp
                              minpa_anaer = (anaer_mass - depnet) * ch(jrch)%anaer_minpa + depnet * ch(jrch)%ext_minpa
                              minps_anaer = (anaer_mass - depnet) * ch(jrch)%anaer_minps + depnet * ch(jrch)%ext_minps
                              orgp_anaer = (anaer_mass - depnet) * ch(jrch)%anaer_orgp + depnet * ch(jrch)%ext_orgp
                              disp_anaer = (anaer_mass - depnet) * conv_mass2vol * ch(jrch)%anaer_disp + & 
                                    depnet * conv_mass2vol * ch(jrch)%ext_disp
                        else
                              mass_ero_aer = 0
                              mass_ero_anaer = 0 
                              if (depnet > - aer_mass) then
                                    mass_dep_aer = -depnet
                                    mass_dep_anaer = 0
                                    mass_aer2anaer = -depnet
                                    minpa_aer = mass_dep_aer * ch(jrch)%minpa * conv_conc + (aer_mass - mass_dep_aer) * ch(jrch)%aer_minpa
                                    minps_aer = mass_dep_aer * ch(jrch)%minps * conv_conc + (aer_mass - mass_ero_aer) * ch(jrch)%aer_minps
                                    orgp_aer = mass_dep_aer * ch(jrch)%organicp * conv_conc + (aer_mass - mass_ero_aer) * ch(jrch)%aer_orgp
                                    disp_aer = mass_dep_aer * conv_mass2vol * ch(jrch)%disolvp + & 
                                          (aer_mass - mass_ero_aer) * conv_mass2vol * ch(jrch)%aer_disp
                                    minpa_anaer = mass_aer2anaer * ch(jrch)%aer_minpa + (anaer_mass - mass_aer2anaer) * ch(jrch)%anaer_minpa
                                    minps_anaer = mass_aer2anaer * ch(jrch)%aer_minps + (anaer_mass - mass_aer2anaer) * ch(jrch)%anaer_minps
                                    orgp_anaer = mass_aer2anaer * ch(jrch)%aer_orgp + (anaer_mass - mass_aer2anaer) * ch(jrch)%anaer_orgp
                                    disp_anaer = mass_aer2anaer * conv_mass2vol * ch(jrch)%aer_disp + & 
                                          (anaer_mass - mass_aer2anaer) * conv_mass2vol * ch(jrch)%anaer_disp
                              else
                                    if (depnet > - anaer_mass) then
                                          mass_dep_aer = aer_mass
                                          mass_dep_anaer = -depnet - aer_mass
                                          mass_aer2anaer = aer_mass
                                          minpa_aer = mass_dep_aer * ch(jrch)%minpa * conv_conc
                                          minps_aer = mass_dep_aer * ch(jrch)%minps * conv_conc
                                          orgp_aer = mass_dep_aer * ch(jrch)%organicp * conv_conc
                                          disp_aer = mass_dep_aer * conv_mass2vol * ch(jrch)%disolvp
                                          minpa_anaer = mass_dep_anaer * ch(jrch)%minpa * conv_conc + mass_aer2anaer * ch(jrch)%aer_minpa + & 
                                                (anaer_mass + depnet) * ch(jrch)%anaer_minpa
                                          minps_anaer = mass_dep_anaer * ch(jrch)%minps * conv_conc + mass_aer2anaer * ch(jrch)%aer_minps + & 
                                                (anaer_mass + depnet) * ch(jrch)%anaer_minps
                                          orgp_anaer = mass_dep_anaer * ch(jrch)%organicp * conv_conc + mass_aer2anaer * ch(jrch)%aer_orgp + & 
                                                (anaer_mass + depnet) * ch(jrch)%anaer_orgp
                                          disp_anaer = mass_dep_anaer * conv_mass2vol * ch(jrch)%disolvp + & 
                                                mass_aer2anaer * conv_mass2vol * ch(jrch)%aer_disp + & 
                                                (anaer_mass + depnet) * conv_mass2vol * ch(jrch)%anaer_disp
                                    else
                                          if (depnet > - anaer_mass - aer_mass) then
                                                mass_dep_aer = aer_mass
                                                mass_dep_anaer = -depnet - aer_mass
                                                mass_aer2anaer = anaer_mass - mass_dep_anaer
                                                minpa_aer = mass_dep_aer * ch(jrch)%minpa * conv_conc
                                                minps_aer = mass_dep_aer * ch(jrch)%minps * conv_conc
                                                orgp_aer = mass_dep_aer * ch(jrch)%organicp * conv_conc
                                                disp_aer = mass_dep_aer * conv_mass2vol * ch(jrch)%disolvp
                                                minpa_anaer = mass_dep_anaer * ch(jrch)%minpa * conv_conc + mass_aer2anaer * ch(jrch)%aer_minpa
                                                minps_anaer = mass_dep_anaer * ch(jrch)%minps * conv_conc + mass_aer2anaer * ch(jrch)%aer_minps
                                                orgp_anaer = mass_dep_anaer * ch(jrch)%organicp * conv_conc + mass_aer2anaer * ch(jrch)%aer_orgp
                                                disp_anaer = mass_dep_anaer * conv_mass2vol * ch(jrch)%disolvp + & 
                                                      mass_aer2anaer * conv_mass2vol * ch(jrch)%aer_disp
                                          else
                                                mass_dep_aer = aer_mass
                                                mass_dep_anaer = anaer_mass
                                                mass_aer2anaer = 0
                                                minpa_aer = mass_dep_aer * ch(jrch)%minpa * conv_conc
                                                minps_aer = mass_dep_aer * ch(jrch)%minps * conv_conc
                                                orgp_aer = mass_dep_aer * ch(jrch)%organicp * conv_conc
                                                disp_aer = mass_dep_aer * conv_mass2vol * ch(jrch)%disolvp
                                                minpa_anaer = mass_dep_anaer * ch(jrch)%minpa * conv_conc
                                                minps_anaer = mass_dep_anaer * ch(jrch)%minps * conv_conc
                                                orgp_anaer = mass_dep_anaer * ch(jrch)%organicp * conv_conc
                                                disp_anaer = mass_dep_anaer * conv_mass2vol * ch(jrch)%disolvp
                                          end if
                                    end if
                              end if
                        end if 
                  end if
            end if
      end if

      ch(jrch)%aer_orgp = orgp_aer / aer_mass ! g P / Mg
      ch(jrch)%aer_disp = disp_aer / (aer_mass * conv_mass2vol) ! mg P / L
      ch(jrch)%aer_minpa = minpa_aer / aer_mass
      ch(jrch)%aer_minps = minps_aer / aer_mass

      ch(jrch)%anaer_orgp = orgp_anaer / anaer_mass
      ch(jrch)%anaer_disp = disp_anaer / (anaer_mass * conv_mass2vol)
      ch(jrch)%anaer_minpa = minpa_anaer / anaer_mass
      ch(jrch)%anaer_minps = minps_anaer / anaer_mass

      minpa_wtclm = minpa_wtclm + mass_ero_aer * ch(jrch)%aer_minpa + mass_ero_anaer * ch(jrch)%anaer_minpa - & 
            dep * ch(jrch)%minpa * conv_conc
      minps_wtclm = minps_wtclm + mass_ero_aer * ch(jrch)%aer_minps + mass_ero_anaer * ch(jrch)%anaer_minps - & 
            dep * ch(jrch)%minps * conv_conc
      orgp_wtclm = orgp_wtclm + mass_ero_aer * ch(jrch)%aer_orgp + mass_ero_anaer * ch(jrch)%anaer_orgp - & 
            dep * ch(jrch)%organicp * conv_conc
      disp_wtclm = disp_wtclm + mass_ero_aer * conv_mass2vol * ch(jrch)%aer_disp + &
            mass_ero_anaer * conv_mass2vol * ch(jrch)%anaer_disp - dep * conv_mass2vol * ch(jrch)%disolvp
      
      ch(jrch)%organicp = orgp_wtclm / qdin
      ch(jrch)%disolvp = disp_wtclm / (qdin * 1000)
      ch(jrch)%minpa = minpa_wtclm / qdin
      ch(jrch)%minps = minps_wtclm / qdin
      
!!    Phosphorus movement - KDW

	ch(jrch)%depch = ch(jrch)%depch + dep - deg1
      if (ch(jrch)%depch < 1.e-6) ch(jrch)%depch = 0.

	sedin = sedin + deg1 + deg2 - dep
	if (sedin < 1.e-6) sedin = 0.

	outfract = rtwtr / qdin
	if (outfract > 1.) outfract = 1.

      sedrch = sedin * outfract
      if (sedrch < 1.e-6) sedrch = 0.

      ch(jrch)%sedst = sedin - sedrch
      if (ch(jrch)%sedst < 1.e-6) ch(jrch)%sedst = 0.

!!    Mass balance tests
!!	ambalsed = sedinorg + deg1 + deg2 - dep - sedrch - sedst(jrch)
!!	if (ambalsed .gt. 1e-3) write (*,*) time%day, jrch, ambalsed

!!  In this default sediment routing sediment is not tracked by particle size
      rch_san = 0.
      rch_sil = sedrch  !! As particles are not tracked by size, the sediments 
      rch_cla = 0.      !! in reach is assumed to be silt for mass conservation
      rch_sag = 0.
      rch_lag = 0.
      rch_gra = 0.

!!    Organic nitrogen and Organic Phosphorus contribution from channel erosion
!!    ch_orgn(jrch) = deg2 * ch_nut(jnut)%onco * 1000.
!!    ch_orgp(jrch) = deg2 * ch_nut(jnut)%opco * 1000.

      ch(jrch)%orgn = deg2 * ch_nut(jnut)%onco / 1000.
!      ch(jrch)%orgp = deg2 * ch_nut(jnut)%opco / 1000.

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
            ch_hyd(jhyd)%s = ch_hyd(jhyd)%s - dat2 / (ch_hyd(jhyd)%l  *     &
                                                           1000.)
            ch_hyd(jhyd)%s = Max(.0001, ch_hyd(jhyd)%s)
            call ch_ttcoef(jrch)
          endif
        endif
      endif

	else
	  sedrch = 0.
	  rch_san = 0.
	  rch_sil = 0.
	  rch_cla = 0.
	  rch_sag = 0.
	  rch_lag = 0.
	  rch_gra = 0.
        ch(jrch)%sedst = sedin

	endif !! end of qdin > 0.01 loop

      endif  !! end of rtwtr and rchdep > 0 loop

      return
      end subroutine ch_rtsed2