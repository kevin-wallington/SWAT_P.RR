      subroutine ch_watqual5

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine performs in-stream nutrient transformations and water
!!    quality calculations

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ai0          |ug chla/mg alg|ratio of chlorophyll-a to algal biomass
!!    ai1          |mg N/mg alg   |fraction of algal biomass that is nitrogen
!!    ai2          |mg P/mg alg   |fraction of algal biomass that is phosphorus
!!    ai3          |mg O2/mg alg  |the rate of oxygen production per unit of
!!                                |algal photosynthesis
!!    ai4          |mg O2/mg alg  |the rate of oxygen uptake per unit of algae
!!                                |respiration
!!    ai5          |mg O2/mg N    |the rate of oxygen uptake per unit of NH3
!!                                |nitrogen oxidation
!!    ai6          |mg O2/mg N    |the rate of oxygen uptake per unit of NO2
!!                                |nitrogen oxidation
!!    algae(:)     |mg alg/L      |algal biomass concentration in reach
!!    ammonian(:)  |mg N/L        |ammonia concentration in reach
!!    bc1(:)       |1/day         |rate constant for biological oxidation of NH3
!!                                |to NO2 in reach at 20 deg C
!!    bc2(:)       |1/day         |rate constant for biological oxidation of NO2
!!                                |to NO3 in reach at 20 deg C
!!    bc3(:)       |1/day         |rate constant for hydrolysis of organic N to
!!                                |ammonia in reach at 20 deg C
!!    bc4(:)       |1/day         |rate constant for the decay of organic P to
!!                                |dissolved P in reach at 20 deg C
!!    chlora(:)    |mg chl-a/L    |chlorophyll-a concentration in reach
!!    dayl(:)      |hours         |day length for current day
!!    disolvp(:)   |mg P/L        |dissolved phosphorus concentration in reach
!!    hru_ra(:)    |MJ/m^2        |solar radiation for the day in HRU
!!    inum1        |none          |reach number
!!    inum2        |none          |inflow hydrograph storage location number
!!    k_l          |MJ/(m2*hr)    |half saturation coefficient for light
!!    k_n          |mg N/L        |michaelis-menton half-saturation constant
!!                                |for nitrogen
!!    k_p          |mg P/L        |michaelis-menton half saturation constant
!!                                |for phosphorus
!!    lambda0      |1/m           |non-algal portion of the light extinction
!!                                |coefficient
!!    lambda1      |1/(m*ug chla/L)|linear algal self-shading coefficient
!!    lambda2      |(1/m)(ug chla/L)**(-2/3)
!!                                |nonlinear algal self-shading coefficient
!!    mumax        |1/day         |maximum specific algal growth rate at 20 deg 
!!                                |C
!!    nitraten(:)  |mg N/L        |nitrate concentration in reach
!!    nitriten(:)  |mg N/L        |nitrite concentration in reach
!!    organicn(:)  |mg N/L        |organic nitrogen concentration in reach
!!    organicp(:)  |mg P/L        |organic phosphorus concentration in reach
!!    rch_cbod(:)  |mg O2/L       |carbonaceous biochemical oxygen demand in
!!                                |reach 
!!    rch_dox(:)   |mg O2/L       |dissolved oxygen concentration in reach
!!    rchdep       |m             |depth of flow on day
!!    rchwtr       |m^3 H2O       |water stored in reach at beginning of day
!!    rhoq         |1/day         |algal respiration rate at 20 deg C
!!    rttime       |hr            |reach travel time
!!    rtwtr        |m^3 H2O       |flow out of reach
!!    tfact        |none          |fraction of solar radiation computed in the
!!                                |temperature heat balance that is
!!                                |photosynthetically active
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    algae(:)    |mg alg/L      |algal biomass concentration in reach
!!    ammonian(:) |mg N/L        |ammonia concentration in reach
!!    chlora(:)   |mg chl-a/L    |chlorophyll-a concentration in reach
!!    disolvp(:)  |mg P/L        |dissolved phosphorus concentration in reach
!!    nitraten(:) |mg N/L        |nitrate concentration in reach
!!    nitriten(:) |mg N/L        |nitrite concentration in reach
!!    organicn(:) |mg N/L        |organic nitrogen concentration in reach
!!    organicp(:) |mg P/L        |organic phosphorus concentration in reach
!!    rch_cbod(:) |mg O2/L       |carbonaceous biochemical oxygen demand in
!!                               |reach
!!    rch_dox(:)  |mg O2/L       |dissolved oxygen concentration in reach
!!    soxy        |mg O2/L       |saturation concetration of dissolved oxygen
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    algcon      |mg alg/L      |initial algal biomass concentration in reach
!!    algi        |MJ/(m2*hr)    |daylight average, photosynthetically active,
!!                               |light intensity
!!    algin       |mg alg/L      |algal biomass concentration in inflow
!!    ammoin      |mg N/L        |ammonium N concentration in inflow
!!    bc1mod      |1/day         |rate constant for biological oxidation of NH3
!!                               |to NO2 modified to reflect impact of low 
!!                               |oxygen concentration
!!    bc2mod      |1/day         |rate constant for biological oxidation of NO2
!!                               |to NO3 modified to reflect impact of low
!!                               |oxygen concentration
!!    cbodcon     |mg/L          |initial carbonaceous biological oxygen demand
!!                               |concentration in reach
!!    cbodin      |mg/L          |carbonaceous biological oxygen demand 
!!                               |concentration in inflow
!!    chlin       |mg chl-a/L    |chlorophyll-a concentration in inflow
!!    cinn        |mg N/L        |effective available nitrogen concentration
!!    cordo       |none          |nitrification rate correction factor
!!    disoxin     |mg O2/L       |dissolved oxygen concentration in inflow
!!    dispin      |mg P/L        |soluble P concentration in inflow
!!    f1          |none          |fraction of algal nitrogen uptake from
!!                               |ammonia pool
!!    fl_1        |none          |growth attenuation factor for light, based on
!!                               |daylight-average light intensity
!!    fll         |none          |growth attenuation factor for light averaged
!!                               |over the diurnal cycle
!!    fnn         |none          |algal growth limitation factor for nitrogen
!!    fpp         |none          |algal growth limitation factor for phosphorus
!!    gra         |1/day         |local algal growth rate at 20 deg C
!!    jrch        |none          |reach number
!!    lambda      |1/m           |light extinction coefficient
!!    nh3con      |mg N/L        |initial ammonia concentration in reach
!!    nitratin    |mg N/L        |nitrate concentration in inflow
!!    nitritin    |mg N/L        |nitrite concentration in inflow
!!    no2con      |mg N/L        |initial nitrite concentration in reach
!!    no3con      |mg N/L        |initial nitrate concentration in reach
!!    o2con       |mg O2/L       |initial dissolved oxygen concentration in 
!!                               |reach
!!    orgncon     |mg N/L        |initial organic N concentration in reach
!!    orgnin      |mg N/L        |organic N concentration in inflow
!!    orgpcon     |mg P/L        |initial organic P concentration in reach
!!    orgpin      |mg P/L        |organic P concentration in inflow
!!    solpcon     |mg P/L        |initial soluble P concentration in reach
!!    tday        |none          |flow duration (fraction of 24 hr)
!!    thbc1       |none          |temperature adjustment factor for local
!!                               |biological oxidation of NH3 to NO2
!!    thbc2       |none          |temperature adjustment factor for local
!!                               |biological oxidation of NO2 to NO3
!!    thbc3       |none          |temperature adjustment factor for local
!!                               |hydrolysis of organic N to ammonia N
!!    thbc4       |none          |temperature adjustment factor for local
!!                               |decay of organic P to dissolved P
!!    thgra       |none          |temperature adjustment factor for local algal
!!                               |growth rate
!!    thrho       |none          |temperature adjustment factor for local algal
!!                               |respiration rate
!!    thrk1       |none          |temperature adjustment factor for local CBOD
!!                               |deoxygenation
!!    thrk2       |none          |temperature adjustment factor for local oxygen
!!                               |reaeration rate
!!    thrk3       |none          |temperature adjustment factor for loss of
!!                               |CBOD due to settling
!!    thrk4       |none          |temperature adjustment factor for local
!!                               |sediment oxygen demand
!!    thrs1       |none          |temperature adjustment factor for local algal
!!                               |settling rate
!!    thrs2       |none          |temperature adjustment factor for local
!!                               |benthos source rate for dissolved phosphorus
!!    thrs3       |none          |temperature adjustment factor for local
!!                               |benthos source rate for ammonia nitrogen
!!    thrs4       |none          |temperature adjustment factor for local
!!                               |organic N settling rate
!!    thrs5       |none          |temperature adjustment factor for local
!!                               |organic P settling rate
!!    wtmp        |deg C         |temperature of water in reach
!!    wtrin       |m^3 H2O       |water flowing into reach on day
!!    uu          |varies        |variable to hold intermediate calculation
!!                               |result
!!    vv          |varies        |variable to hold intermediate calculation
!!                               |result
!!    wtrtot      |m^3 H2O       |inflow + storage water
!!    ww          |varies        |variable to hold intermediate calculation
!!                               |result
!!    xx          |varies        |variable to hold intermediate calculation
!!                               |result
!!    yy          |varies        |variable to hold intermediate calculation
!!                               |result
!!    zz          |varies        |variable to hold intermediate calculation
!!                               |result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log, Exp, Min
!!    SWAT: Theta

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
      !use jrw_datalib_module            
      use channel_module
      use hydrograph_module      
      use climate_module
      use channel_data_module
      use channel_velocity_module
      use basin_module
      use time_module !KDW 10/05

      real :: tday, wtmp, fll, gra
      real :: lambda, fnn, fpp, algi, fl_1, xx, yy, zz, ww, qq, xyz, zyx
      real :: uu, vv, cordo, f1, algcon, orgncon, nh3con, no2con, no3con
      real :: orgpcon, minpacon, minpscon, solpcon, cbodcon, o2con, wtrtot, bc1mod, bc2mod, cbodocon
      real :: orgpcon2, minpacon2, minpscon2, solpcon2, cbodcon2, o2con2, cbodocon2
      real :: aer_disp2, anaer_disp2, aer_orgp2, anaer_orgp2, aer_minpa2, anaer_minpa2, aer_minps2, anaer_minps2
      real :: algcon2, orgncon2, nh3con2, no2con2, no3con2
      real :: thgra = 1.047, thrho = 1.047, thrs1 = 1.024, thw12 = 1.047
      real :: thrs2 = 1.074, thrs3 = 1.074, thrs4 = 1.024, thrs5 = 1.024
      real :: thbc1 = 1.083, thbc2 = 1.047, thbc3 = 1.047, thbc4 = 1.047
      real :: thrk1 = 1.047, thrk2 = 1.024, thrk3 = 1.024, thrk4 = 1.060
      real, parameter :: bk = .0006     !              | ! KDW debug 04/28/22 test faster solid diffusion, was 0.0006
      real :: rs3_s, rk4_s, rk1_k, rk1_m, rk3_k, rk2_m, rk2_k 
      real :: bc1_k, bc2_m, bc3_k, bc3_m, rs2_s, bc4_k, bc4_m
      real :: chlin            !mg chl-a/L    |chlorophyll-a concentration in inflow
      real :: algin            !mg alg/L      |algal biomass concentration in inflow
      real :: orgnin           !mg N/L        |organic N concentration in inflow
      real :: ammoin           !mg N/L        |ammonium N concentration in inflow
      real :: nitratin         !mg N/L        |nitrate concentration in inflow
      real :: nitritin         !mg N/L        |nitrite concentration in inflow
      real :: orgpin           !mg P/L        |organic P concentration in inflow 
      real :: minpain          !mg P/L        |active mineral P concentration in inflow - KDW
      real :: minpsin          !mg P/L        |stable mineral P concentration in inflow - KDW
      real :: dispin           !mg P/L        |soluble P concentration in inflow
      real :: cbodin           !mg/L          |carbonaceous biological oxygen demand
      real :: disoxin          !mg O2/L       |dissolved oxygen concentration in inflow
      real :: soxy             !mg O2/L       |saturation concetration of dissolved oxygen
      real :: factk, factm
      real :: alg_m, alg_md, alg_mg, alg_set, algcon_out
      real :: alg_m_o2, alg_m_disp, alg_m_orgp, alg_nh4_m, alg_no3_m
      real :: fines, aer_mass, anaer_mass, aer_mass_remain, anaer_mass_remain
      real :: sedin, sanin, silin, clain, sagin
      real :: topcla, topsil, topsag, topsan
      real :: cla_aer, san_aer, sil_aer, sag_aer
      real :: cla_anaer, san_anaer, sil_anaer, sag_anaer
      real :: andep_cla, andep_sil, andep_sag, andep_san, andep_tot
      real :: w12, diff_aer2anaer, diff_aer2wtclm, diff_anaer2wtclm
      real :: dox_factor, cbod_factor, btrb_orgp, btrb_minpa, btrb_minps
      real :: mnrlz_aer, mnrlz_anaer
      real :: conv_mass2vol, conv_conc, wtclm_por
      !real :: adsts_wtclm, adsts_aer, adsts_anaer
      !real :: frsts_wtclm, frsts_aer, frsts_anaer
      real :: pai_wtclm, pai_aer, pai_anaer
      real :: rto_wtclm, rto_aer, rto_anaer
      real :: ads_gradient, ssd_gradient
      real :: sfarea_anaer, sfarea_aer, sfarea_wtclm
      real :: rmp_aer, rmp_anaer, roc_aer, roc_anaer
      !real :: ads_m, ads_m1, ads_m2
      real :: ssd_m, ssd_m1, ssd_m2
      real :: qdin
      real :: aer_vol, anaer_vol
      real :: equil_conc ! KDW 06/02, 06/15
      real :: disox_ave
      real :: diff, max_diff ! KDW 02/23/22
      real :: wtclm_mixed_minpa, wtclm_mixed_minpa_conc, wtclm_mixed_solp, wtclm_mixed_solp_conc ! KDW 02/23/22
      real :: diffremain! KDW 02/23/22
      real :: P_adsorbed, P_min_new, P_tot, P_adsorbed_mass, P_desorbed_mass, ads_inflow! KDW 03/03/22
      real :: max_diff_minpa, max_diff_disp, max_disp_conc, max_Ptot_conc ! KDW 03/07/22, 03/21/22
      real :: min_disp_conc, min_Ptot_conc, max_diff1, max_diff2, max_diff3 ! KDw 03/22/22
      real :: rs2_ratio, diff_aer2wtclm2, Ptot_conc_balanced,aer_balanced ! KDw 03/29/22, 04/01/22
      real :: aer_disp_temp, aer_minpa_temp, anaer_disp_temp, anaer_minpa_temp ! KDw 04/01/22
      real :: wtclm_mixed_solp_conc_temp, wtclm_mixed_minpa_conc_temp ! KDw 04/01/22
      real :: SA_VOL_wtclm, SA_VOL_anaer, SA_VOL_aer ! KDW 04/15/22

      !! initialize new variables - KDW
      rs3_s = 0. ; rk4_s = 0. ; rk1_k = 0. ; rk1_m = 0. ; rk3_k = 0. ; rk2_m = 0. ; rk2_k = 0. 
      bc1_k = 0. ; bc2_m  = 0. ; bc3_k = 0. ; bc3_m = 0. ; rs2_s = 0. ; bc4_k = 0. ; bc4_m = 0.
      factk = 0. ; factm = 0.
      alg_m = 0. ; alg_md = 0. ; alg_mg = 0. ; alg_set = 0. ; algcon_out = 0.
      alg_m_o2 = 0. ; alg_m_disp = 0. ; alg_m_orgp = 0. ; alg_nh4_m = 0. ; alg_no3_m = 0.
      fines = 0. ; aer_mass = 0. ; anaer_mass = 0. ; aer_mass_remain = 0. ; anaer_mass_remain = 0.
      sedin = 0. ; sanin = 0. ; silin = 0. ; clain = 0. ; sagin = 0.
      topcla = 0. ; topsil = 0. ; topsag = 0. ; topsan = 0.
      cla_aer = 0. ; san_aer = 0. ; sil_aer = 0. ; sag_aer = 0.
      cla_anaer = 0. ; san_anaer = 0. ; sil_anaer = 0. ; sag_anaer = 0.
      andep_cla = 0. ; andep_sil = 0. ; andep_sag = 0. ; andep_san = 0. ; andep_tot = 0.
      w12 = 0. ; rs2_s = 0. ; diff_aer2anaer = 0. ; diff_aer2wtclm = 0. ; diff_anaer2wtclm = 0.
      dox_factor = 0. ; cbod_factor = 0. ; btrb_orgp = 0. ; btrb_minpa = 0. ; btrb_minps = 0.
      mnrlz_aer = 0. ; mnrlz_anaer = 0.
      conv_mass2vol = 0. ; conv_conc = 0. ; wtclm_por = 0.
      !adsts_wtclm = 0. ; adsts_aer = 0. ; adsts_anaer = 0.
      !frsts_wtclm = 0. ; frsts_aer = 0. ; frsts_anaer = 0.
      pai_wtclm = 0. ; pai_aer = 0. ; pai_anaer = 0.
      rto_wtclm = 0. ; rto_aer = 0. ; rto_anaer = 0.
      ads_gradient = 0. ; ssd_gradient = 0.
      sfarea_anaer = 0. ; sfarea_aer = 0. ; sfarea_wtclm = 0.
      rmp_aer = 0. ; rmp_anaer = 0. ; roc_aer = 0. ; roc_anaer = 0.
      ads_m = 0. ; ads_m1 = 0. ; ads_m2 = 0.
      ssd_m = 0. ; ssd_m1 = 0. ; ssd_m2 = 0.
      aer_vol = 0. ; anaer_vol = 0.
      equil_conc = 0. ! KDW 06/02, 06/15
      disox_ave = 0.
      xyz = 0. ; zyx = 0. ; diff = 0. ; max_diff = 0. ! KDW 02/23/22
      wtclm_mixed_minpa = 0. ; wtclm_mixed_minpa_conc = 0. ! KDW 02/23/22
      wtclm_mixed_solp = 0. ; wtclm_mixed_solp_conc = 0. ! KDW 02/23/22
      diffremain = 0.! KDW 02/23/22
      yyy = 0. ! 02/28/22
      P_adsorbed = 0.; P_min_new = 0.; P_tot = 0.; P_adsorbed_mass = 0.; P_desorbed_mass = 0.; ads_inflow = 0.  ! KDW 03/03/22
      max_diff2 = 0.; B_aer = 0.; B_anaer = 0.; B_wtclm = 0. ! KDW 03/05/22
      max_diff_minpa = 0.; max_diff_disp = 0.; max_disp_conc = 0.; max_Ptot_conc = 0. ! KDW 03/07/22, 03/21/22
      min_disp_conc = 0.; min_Ptot_conc = 0.; max_diff1 = 0.; max_diff2 = 0.; max_diff3 = 0. ! KDw 03/22/22
      rs2_ratio = 0.1 ! KDW ratio of rs2_s b/w wtclm and aerobic layer vs. b/w aerobic and anaerobic layer, 03/29/22
      diff_aer2wtclm2 = 0.; Ptot_conc_balanced = 0.; aer_balanced = 0.; ! KDw 03/29/22, 04/01/22
      aer_disp_temp = 0.; aer_minpa_temp = 0.; anaer_disp_temp = 0.; anaer_minpa_temp = 0. ! KDw 04/01/22
      wtclm_mixed_solp_conc_temp = 0.; wtclm_mixed_minpa_conc_temp = 0. ! KDw 04/01/22
      SA_VOL_wtclm = 0. ; SA_VOL_anaer = 0. ; SA_VOL_aer = 0. ! KDW 04/15/22

!!    rnum1        |none          |fraction of overland flow

         !! calculate flow duration
         !tday = rttime / 24.0
         wtrin = wtrin + qdbank ! added to add qdbank, KDW 04/04/22, moved up from lines 337-338 on 04/18/22
         qdin = wtrin + vol_0 ! changed from ch(jrch)%rchstor to vol_0, KDW 04/04/22, here and throughout routine
         if (wtrin > 0) then
            tday = vol_0/ wtrin ! KDW - this code really uses this ratio, not truly residence time ! changed from ch(jrch)%rchstor to vol_0, KDW 04/04/22
         else
            tday = 10000.
         endif
         tday = min(10000., tday) ! sets arbitratily large maximum residence time of water
         tday = max(.001, tday) ! sets arbitrarily small minimum residence time of water

        !! calculate temperature in stream Stefan and Preudhomme. 1993.  Stream temperature estimation 
        !! from air temperature.  Water Res. Bull. p. 27-45 SWAT manual equation 2.3.13
         wtmp = 5.0 + 0.75 * wst(iwst)%weat%tave
         if (wtmp <= 0.) wtmp = 0.1

!        benthic sources/losses in mg   
        
         rs3_s =  Theta(ch_nut(jnut)%rs3,thrs3,wtmp)*ch_hyd(jhyd)%l *ch_hyd(jhyd)%w * rt_delt ! [(g/day)/(m^2 * mg/L)] - KDW
         rk4_s =  Theta(ch_nut(jnut)%rk4,thrk4,wtmp)*ch_hyd(jhyd)%l *ch_hyd(jhyd)%w * rt_delt ! [(g/day)/(m^2 * mg/L)] - KDW

            !! initialize concentration of nutrient in reach
         wtrtot = 0.
         algcon = 0.
         orgncon = 0.
         nh3con = 0.
         no2con = 0.
         no3con = 0.
         orgpcon = 0.
         minpacon = 0. ! - KDW
         minpscon = 0. ! - KDW
         solpcon = 0.
         cbodcon = 0.
         o2con = 0.
         ch(jrch)%rch_cbod = amax1(1.e-6, ch(jrch)%rch_cbod)
         
        !wtrin = ob(icmd)%hin%flo
      !   wtrin = wtrin + qdbank ! added to add qdbank, KDW 04/04/22, moved up to lines 301-302 on 04/18/22
      !   qdin = wtrin + vol_0 ! changed from ch(jrch)%rchstor to vol_0, KDW 04/04/22, here and throughout routine
        if ( wtrin > 0.0001) then ! changed from ob(icmd)%hin%flo
         chlin = 1000. * ob(icmd)%hin%chla * rt_delt / wtrin 
         algin = 1000. * chlin / ch_nut(jnut)%ai0        !! QUAL2E equation III-1
         orgnin = 1000. * ob(icmd)%hin%orgn  *rt_delt/ wtrin  
         ammoin = 1000. * ob(icmd)%hin%nh3 * rt_delt / wtrin  
         nitritin = 1000. * ob(icmd)%hin%no2 * rt_delt / wtrin
         nitratin = 1000. * ob(icmd)%hin%no3 * rt_delt  / wtrin
         !orgpin = 1000. * ob(icmd)%hin%sedp *rt_delt / wtrin 
         orgpin = 1000. * ob(icmd)%hin%orgp *rt_delt / wtrin ! - KDW
         minpain = 1000. * ob(icmd)%hin%minpa *rt_delt / wtrin ! - KDW
         minpsin = 1000. * ob(icmd)%hin%minps *rt_delt / wtrin ! - KDW
         dispin = 1000. * ob(icmd)%hin%solp* rt_delt  / wtrin
         cbodin = 1000. * ob(icmd)%hin%cbod * rt_delt / wtrin
         disoxin = 1000. * ob(icmd)%hin%dox  * rt_delt / wtrin
         if (chlin < 1.e-6) chlin = 0.0 ! - KDW, 03/26
         if (algin < 1.e-6) algin = 0.0 ! - KDW, 03/26
         if (orgnin < 1.e-6) orgnin = 0.0 ! - KDW, 03/26
         if (ammoin < 1.e-6) ammoin = 0.0 ! - KDW, 03/26
         if (nitritin < 1.e-6) nitritin = 0.0 ! - KDW, 03/26
         if (nitratin < 1.e-6) nitratin = 0.0 ! - KDW, 03/26
         if (orgpin < 1.e-6) orgpin = 0.0 ! - KDW, 03/26
         if (minpain < 1.e-6) minpain = 0.0 ! - KDW, 03/26
         if (minpsin < 1.e-6) minpsin = 0.0 ! - KDW, 03/26
         if (dispin < 1.e-6) dispin = 0.0 ! - KDW, 03/26
         if (cbodin < 1.e-6) cbodin = 0.0 ! - KDW, 03/26
         if (disoxin < 1.e-6) disoxin = 0.0 ! - KDW, 03/26
        else
         chlin = 0
         algin = 0
         orgnin = 0  
         ammoin = 0 
         nitritin = 0
         nitratin = 0
         !orgpin = 0 
         orgpin = 0
         minpain = 0
         minpsin = 0
         dispin = 0
         cbodin = 0
         disoxin = 0
        end if    
                 
         !algcon =ch(jrch)%algae 
        algcon = 1000. * ch(jrch)%chlora / ch_nut(jnut)%ai0 ! changed to match watqual4
        orgncon = ch(jrch)%organicn
        nh3con = ch(jrch)%ammonian 
        no2con =  ch(jrch)%nitriten
        no3con = ch(jrch)%nitraten 
        orgpcon =ch(jrch)%organicp 
        minpacon =ch(jrch)%minpa ! - KDW
        minpscon =ch(jrch)%minps ! - KDW
        solpcon = ch(jrch)%disolvp 
        cbodcon= ch(jrch)%rch_cbod 
        o2con = ch(jrch)%rch_dox 
        if (orgncon < 1.e-6) orgncon = 0.0
	    if (nh3con < 1.e-6) nh3con = 0.0
	    if (no2con < 1.e-6) no2con = 0.0
	    if (no3con < 1.e-6) no3con = 0.0
        if (orgpcon < 1.e-6) orgpcon = 0.0
        if (minpacon < 1.e-6) minpacon = 0.0 ! - KDW
        if (minpscon < 1.e-6) minpscon = 0.0 ! - KDW
	    if (solpcon < 1.e-6) solpcon = 0.0
	    if (cbodcon < 1.e-6) cbodcon = 0.0
        if (o2con < 1.e-6) o2con = 0.0
         
        if (wtrin> 0.0001) then !changed from > 0 - KDW
          !! Update dispin here for diffusion between water column and aerobic layer
         disoxin = disoxin - rk4_s /wtrin
         if (disoxin < 0) then   ! subtract from o2con - KDW
            if (vol_0> 0) then
                o2con = ((o2con * vol_0) + (disoxin * wtrin)) / vol_0
            else
                o2con= 0
            endif
            disoxin = 0
         endif
                  
         if (o2con < 1.e-6) o2con = 0.0
         !disoxin = max(0., disoxin)
         !dispin = dispin + rs2_s / wtrin ! new formulation in P section below - KDW
         ammoin = ammoin + rs3_s / wtrin
         if (ammoin < 0) then   ! subtract from ammocon - KDW
             if (vol_0 > 0) then   
                nh3con = ((nh3con * vol_0) + (ammoin * wtrin)) / vol_0
             else
                 nh3con = 0
             endif
            ammoin = 0
         endif
         if (nh3con < 1.e-6) nh3con = 0.0

         !! calculate effective concentration of available nitrogen QUAL2E equation III-15
         !cinn = nh3con + no3con
         cinn = ((nh3con + no3con) * vol_0 + (ammoin +nitratin) * wtrin)/(wtrin + vol_0) ! KDW 12/13/21 in case no water in reach to start day, this formulation prevents instability in k2m
         cinp = (solpcon * vol_0 + dispin * wtrin)/(wtrin + vol_0)

         !! calculate saturation concentration for dissolved oxygen
         !! QUAL2E section 3.6.1 equation III-29
         ww = -139.34410 + (1.575701e05 / (wtmp + 273.15))
         xx = 6.642308e07 / ((wtmp + 273.15)**2)
         yy = 1.243800e10 / ((wtmp + 273.15)**3)
         zz = 8.621949e11 / ((wtmp + 273.15)**4)
         soxy = Exp(ww - xx + yy - zz)
         if (soxy < 1.e-6) soxy = 0. 
         !! end initialize concentrations

         !! O2 impact calculations

         !! calculate nitrification rate correction factor for low
         !! oxygen QUAL2E equation III-21
	     if (o2con.le.0.001) o2con=0.001
	     if (o2con.gt.30.) o2con=30.
         cordo = 1.0 - Exp(-0.6 * o2con)
         !! modify ammonia and nitrite oxidation rates to account for low oxygen
         bc1mod = ch_nut(jnut)%bc1 * cordo !!! not in watqual4
         bc2mod = ch_nut(jnut)%bc2 * cordo !!! not in watqual4
         !! end O2 impact calculations
       
         !! algal growth
         !! calculate light extinction coefficient 
         !! (algal self shading) QUAL2E equation III-12
         if (ch_nut(jnut)%ai0 * algcon > 1.e-6) then
           lambda = ch_nut(jnut)%lambda0 + (ch_nut(jnut)%lambda1 *      &
              ch_nut(jnut)%ai0 * algcon)                                &
              + ch_nut(jnut)%lambda2 * (ch_nut(jnut)%ai0 *              & 
              algcon) ** (.66667)
         else
           lambda = ch_nut(jnut)%lambda0
         endif

	     if (lambda > ch_nut(jnut)%lambda0) lambda = ch_nut(jnut)%lambda0
         !! calculate algal growth limitation factors for nitrogen
         !! and phosphorus QUAL2E equations III-13 & III-14
         fnn = cinn / (cinn + ch_nut(jnut)%k_n)
         !fpp = solpcon / (solpcon + ch_nut(jnut)%k_p) ! replaced with below by KDW in case no water in reach at beginning of time step, 12/13/21
         fpp = cinp / (cinp + ch_nut(jnut)%k_p)      

         !! calculate daylight average, photosynthetically active,
         !! light intensity QUAL2E equation III-8
         !! Light Averaging Option # 2
         iwgn = wst(iwst)%wco%wgn
         if (wst(iwst)%weat%daylength > 0.) then ! Corrected from wgn_pms(iwgn)%daylth (the dormancy threshold) to wst(iwst)%weat%daylength - KDW, 02/12/22
           algi = wst(iwst)%weat%solrad * ch_nut(jnut)%tfact /  wst(iwst)%weat%daylength ! Corrected from wgn_pms(iwgn)%daylth (the dormancy threshold) to wst(iwst)%weat%daylength - KDW, 02/12/22
         else
           algi = 0.00001
         end if

         !! calculate growth attenuation factor for light, based on
         !! daylight average light intensity QUAL2E equation III-7b
         if (rchdep > 0) then
            fl_1 = (1. / (lambda * rchdep)) *                               &                             
               Log((ch_nut(jnut)%k_l + algi) / (ch_nut(jnut)%k_l + algi *  &
               (Exp(-lambda * rchdep))))
         else
            fl_1 = 1
         endif
         fll = 0.92 * (wst(iwst)%weat%daylength / 24.) * fl_1 ! Corrected from wgn_pms(iwgn)%daylth (the dormancy threshold) to wst(iwst)%weat%daylength - KDW, 05/24

         !! calculcate local algal growth rate
         gra = 0. ! KDW 03/10/22
         if (algcon < 5000.) then
          select case (ch_nut(jnut)%igropt)
           case (1)
             !! multiplicative QUAL2E equation III-3a
             gra = ch_nut(jnut)%mumax * fll * fnn * fpp
           case (2)
             !! limiting nutrient QUAL2E equation III-3b
             gra = ch_nut(jnut)%mumax * fll * Min(fnn, fpp)
           case (3)
             !! harmonic mean QUAL2E equation III-3c
             if (fnn > 1.e-6 .and. fpp > 1.e-6) then
               gra = ch_nut(jnut)%mumax * fll * 2. / ((1. / fnn) + (1. / fpp))
             else
               gra = 0.
             endif
          end select
         end if

         !! calculate algal biomass concentration at end of day (phytoplanktonic algae)
         !! QUAL2E equation III-2
         !ch(jrch)%algae = 0. ! - eliinated to match watqual4 - KDW
         !factm = 0.! - eliinated to match watqual4 - KDW 
         factk = Theta(gra,thgra,wtmp)   ! added to correct growth rate calc, below - KDW     
         alg_mg = wq_k2m (tday, rt_delt, factk, algcon, algin) ! the growth rate ; changed from semianalyt to k2m - KDW
         if ((alg_mg * ch_nut(jnut)%ai2) > (cinp*(wtrin + vol_0))) then !KDW 02/12/22, so growth doesn't require more disp than available
            alg_mg  = (cinp*(wtrin + vol_0)) / ch_nut(jnut)%ai2
         endif
         if ((alg_mg * ch_nut(jnut)%ai1) > (cinn*(wtrin + vol_0))) then!KDW 02/12/22, so growth doesn't require more N than available
            alg_mg  = (cinn*(wtrin + vol_0)) / ch_nut(jnut)%ai1
         endif
         if (alg_mg  > 500) then  !KDW 02/12/22, so algal growth isn't unstable
            alg_mg  = 500
         endif
         if (alg_mg < 1.e-6) alg_mg = 0 ! KDW 12/13/21
         factk = Theta(gra,thgra,wtmp) - Theta(ch_nut(jnut)%rhoq, thrho, wtmp)      
         alg_m = wq_k2m(tday, rt_delt, factk, algcon, algin) ! net growth rate; changed from semianalyt to k2m - KDW
         alg_md = alg_mg-alg_m ! the death rate ; edited by negating RHS, was the negative of the death rate - KDW

         if (alg_md < 1.e-6) alg_md = 0 ! KDW 12/13/21
         if (ABS(alg_m) < 1.e-6) alg_m = 0 ! KDW 12/13/21
         alg_set = 0.
         if (rchdep > 0.001) alg_set = Theta (ch_nut(jnut)%rs1, thrs1, wtmp) / rchdep 
         
         algcon_out = wq_semianalyt (tday, rt_delt, alg_m, -alg_set, algcon, algin)           

         ch(jrch)%algae = algcon_out
         if (ch(jrch)%algae < 1.e-6) ch(jrch)%algae = 0.
	     if (ch(jrch)%algae > 5000.) ch(jrch)%algae = 5000.

         !! calculate chlorophyll-a concentration at end of day QUAL2E equation III-1
         ch(jrch)%chlora = ch(jrch)%algae * ch_nut(jnut)%ai0 / 1000.
         !! end algal growth 

         !! oxygen calculations
         !! calculate carbonaceous biological oxygen demand at end of day QUAL2E section 3.5 equation III-26
         !! adjust rk1 to m-term and BOD & O2 mass availability
         
         cbodocon = min (cbodcon,o2con)
         cbodoin = min (cbodin,disoxin)
         rk1_k = -Theta (ch_nut(jnut)%rk1, thrk1,wtmp)
         rk1_m = wq_k2m (tday, rt_delt, rk1_k, cbodocon, cbodoin)
         ! calculate corresponding m-term
         rk3_k=0.
         if (rchdep > 0.001)  rk3_k = -Theta (ch_nut(jnut)%rk3, thrk3, wtmp) / rchdep
         factm = rk1_m
         factk = rk3_k
         ch(jrch)%rch_cbod = 0.
         ch(jrch)%rch_cbod = wq_semianalyt (tday, rt_delt, factm, factk, cbodcon, cbodin)
         if (ch(jrch)%rch_cbod < 1.e-6) ch(jrch)%rch_cbod = 0.

         !! calculate dissolved oxygen concentration if reach at 
         !! end of day QUAL2E section 3.6 equation III-28

         rk2_m = Theta (ch_nut(jnut)%rk2, thrk2, wtmp) * soxy ! this is the piece multiplied by saturated oxy level - KDW
         rk2_k = -Theta (ch_nut(jnut)%rk2, thrk2, wtmp) 

         bc1_k = -Theta(ch_nut(jnut)%bc1,thbc1,wtmp)
         bc2_k = -Theta(ch_nut(jnut)%bc2,thbc2,wtmp)
         bc1_m = wq_k2m (tday, rt_delt, bc1_k, nh3con, ammoin) ! changed from rk2_k to bc1_k - KDW
         bc2_m = wq_k2m (tday, rt_delt, bc2_k, no2con, nitritin)

         alg_m_o2 = -ch_nut(jnut)%ai4 * alg_md + ch_nut(jnut)%ai3 * alg_mg
         disox_ave = (o2con * vol_0+ disoxin * wtrin) / (wtrin + vol_0) ! KDW 12/13/21 so o2 consumed in respiration doesn't exceed available o2
         if(alg_m_o2 + disox_ave < 0) then
             alg_m_o2 = -disox_ave  ! KDW 12/13/21 so o2 consumed in respiration doesn't exceed available o2 (including from photosynthesis in current step)
         endif

         factm = rk1_m + rk2_m + alg_m_o2 + bc1_m * ch_nut(jnut)%ai5 + bc2_m * ch_nut(jnut)%ai6 ! added alg_m_o2, removed rk4_m b/c not defined and already accounted for in doxin - KDW
         ch(jrch)%rch_dox = wq_semianalyt (tday, rt_delt, factm, rk2_k, o2con, disoxin)
         if (ch(jrch)%rch_dox.le.0.001) ch(jrch)%rch_dox=0.001 ! KDW 02/11/22, repeated from lines 465-66 so o2con doesn't grow unreasonably
         if (ch(jrch)%rch_dox.gt.30.) ch(jrch)%rch_dox=30. ! KDW 02/11/22, repeated from lines 465-66 so o2con doesn't grow unreasonably
         if (ch(jrch)%rch_dox < 1.e-6) ch(jrch)%rch_dox = 0. ! changed from < 0 by KDW 12/14/21
         
         !! end oxygen calculations 

         !! nitrogen calculations
         !! calculate fraction of algal nitrogen uptake from ammonia pool QUAL2E equation III-18
         f1 = ch_nut(jnut)%p_n * nh3con / (ch_nut(jnut)%p_n * nh3con +     &
         (1. - ch_nut(jnut)%p_n) * no3con + 1.e-6)

         !! calculate organic N concentration at end of day
         !! QUAL2E section 3.3.1 equation III-16
         bc3_k = Theta(ch_nut(jnut)%bc3,thbc3,wtmp) 
         rs4_k=0.
         if (rchdep > 0.001)  rs4_k = Theta (ch_nut(jnut)%rs4, thrs4, wtmp) / rchdep   

         bc3_m = wq_k2m (tday, rt_delt, -bc3_k, orgncon, orgnin)
         factk =-rs4_k
         alg_orgN_m = alg_md * ch_nut(jnut)%ai1 ! added algae respiration contribution - KDW
         factm = bc3_m + alg_orgN_m ! this formulation includes alg_orgN_m, organic N source from dying algea (only to nitrate/nitrite) - KDW
         ch(jrch)%organicn = wq_semianalyt (tday, rt_delt, factm, factk, orgncon, orgnin)
         if (ch(jrch)%organicn < 1.e-6) ch(jrch)%organicn = 0. 

         !! calculate ammonia nitrogen concentration at end of day QUAL2E section 3.3.2 equation III-17
         factk= 0. ! changed from -bc1_k to 0 for clarity since not used as k factor below - KDW
         alg_nh4_m = alg_mg * f1 * ch_nut(jnut)%ai1
         factm= bc1_m - bc3_m - alg_nh4_m
         ch(jrch)%ammonian= wq_semianalyt(tday,rt_delt,factm,0.,nh3con,ammoin) 
         if (ch(jrch)%ammonian < 1.e-6) ch(jrch)%ammonian = 0. 
  
         !! calculate concentration of nitrite at end of day QUAL2E section 3.3.3 equation III-19
         factm=-bc1_m + bc2_m
         ch(jrch)%nitriten = wq_semianalyt(tday,rt_delt,factm,0.,no2con,nitritin)
         if (ch(jrch)%nitriten < 1.e-6) ch(jrch)%nitriten = 0. 

         !! calculate nitrate concentration at end of day QUAL2E section 3.3.4 equation III-20
         factk = 0.
         alg_no3_m = alg_mg * (1. - f1) * ch_nut(jnut)%ai1
         factm = -bc2_m - alg_no3_m
        
         ch(jrch)%nitraten = wq_semianalyt (tday, rt_delt, factm,0., no3con, nitratin)
         if (ch(jrch)%nitraten < 1.e-6) ch(jrch)%nitraten = 0. 
         !! end nitrogen calculations

         !! phosphorus calculations

         !! Determine surface area of pools
         sedin = ob(icmd)%hin%sed  + ch(jrch)%sedst
         sanin = ob(icmd)%hin%san  + ch(jrch)%sanst
         silin = ob(icmd)%hin%sil  + ch(jrch)%silst
         clain = ob(icmd)%hin%cla  + ch(jrch)%clast
         sagin = ob(icmd)%hin%sag  + ch(jrch)%sagst

         fines = ch(jrch)%depclach + ch(jrch)%depsilch + ch(jrch)%depsagch + ch(jrch)%depsanch
         ! aer_mass = .001 * 1 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd) ! depth of aer layer is 1
         ! anaer_mass = .1 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd)
         ! aer_vol = .001 * ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * 1000 ! liters
         ! anaer_vol = .1 * ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * 1000 ! liters
         aer_mass = .002 * 1 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd) ! depth of aer layer is 1 ! KDW debug 04/28/22, test thicker layers, replacing above
         anaer_mass = .2 * (ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * ch_sed(jsed)%bed_bd) ! KDW debug 04/28/22, test thicker layers
         aer_vol = .002 * ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * 1000 ! liters ! KDW debug 04/28/22, test thicker layers
         anaer_vol = .2 * ch_vel(jrch)%wid_btm * 1000 * ch_hyd(jhyd)%l * 1000 ! liters ! KDW debug 04/28/22, test thicker layers
         vol_h2o_aer = ch(jrch)%bed_por * aer_vol ! KDW 03/21/22
         vol_h2o_anaer = ch(jrch)%bed_por * anaer_vol ! KDW 03/21/22

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
         cla_anaer = andep_cla + anaer_mass_remain * ch(jrch)%bnk_cla /(1-ch(jrch)%bnk_gra) ! Note, final division term adjusts for not including gravel as part of aerboic/anaerobic mass
         sil_anaer = andep_sil + anaer_mass_remain * ch(jrch)%bnk_sil /(1-ch(jrch)%bnk_gra)
         sag_anaer = andep_sag
         san_anaer = andep_san + anaer_mass_remain * ch(jrch)%bnk_san /(1-ch(jrch)%bnk_gra)

         ! sfarea_anaer = cla_anaer + sil_anaer * .2 + sag_anaer * (1./15.) + san_anaer * .01
         ! sfarea_aer = cla_aer + sil_aer * .2 + sag_aer * (1./15.) + san_aer * .01
         ! sfarea_wtclm = clain + silin * .2 + sagin * (1./15.) + sanin * .01
         SA_VOL_anaer = 3. * ((cla_anaer/0.000002) + (sil_anaer/0.00001) + (sag_anaer/0.00003) + (san_anaer/0.0002)) / anaer_mass ! KDW 04/14/22
         SA_VOL_aer = 3. * ((cla_aer/0.000002) + (sil_aer/0.00001) + (sag_aer/0.00003) + (san_aer/0.0002)) / aer_mass ! KDW 04/14/22
         SA_VOL_wtclm = 3. * ((clain/0.000002) + (silin/0.00001) + (sagin/0.00003) + (sanin/0.0002)) / sedin ! KDW 04/14/22
                 
         !! Update aerobic and anaerobic layers for mineralization
         bc4_k = Theta(ch_nut(jnut)%bc4,thbc4,wtmp)
         mnrlz_aer = bc4_k * ch(jrch)%aer_orgp * rt_delt
         mnrlz_anaer = bc4_k * ch(jrch)%anaer_orgp * rt_delt
         ch(jrch)%aer_orgp = ch(jrch)%aer_orgp - mnrlz_aer
         ch(jrch)%anaer_orgp = ch(jrch)%anaer_orgp - mnrlz_anaer
         ch(jrch)%aer_disp = ch(jrch)%aer_disp + mnrlz_aer * ch_sed(jsed)%bed_bd / ch(jrch)%bed_por !  "* ch_sed(jsed)%bed_bd" added 02/28/22, KDW, " / ch(jrch)%bed_por" added 03/21/22
         ch(jrch)%anaer_disp = ch(jrch)%anaer_disp + mnrlz_anaer * ch_sed(jsed)%bed_bd / ch(jrch)%bed_por !  "* ch_sed(jsed)%bed_bd" added 02/28/22, KDW, " / ch(jrch)%bed_por" added 03/21/22

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Add Runge-Kutta, hourly step version, for comparison!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
         !! Update aerobic and anaerobic layers for diffusion
         rs2_s =  Theta(ch_nut(jnut)%rs2,thrs2,wtmp)*ch_hyd(jhyd)%l *ch_hyd(jhyd)%w * 1000 ! repurposed rs2, units [(g_P/day)/(m^2 * mg_P/L)], rs2_s units [(g_p/day)/(mg_P/L)]
         ! KDW 2/14/22 added "* 1000" and changed unit of rs2 to [(g_P/day)/(m^2 * mg_P/L)] because issue with calibration max value
         !if (rs2_s > 1) rs2_s = 1 ! commented out KDW, 06/01
         !!!!KDW 02/21/22 !!!!
         ! dox_factor = 1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half)))
         ! adsts_wtclm = 1000000 * ch_nut(jnut)%spkgc * (sfarea_wtclm / (qdin * 1000)) * &
         !    (0.75 + 0.25 * dox_factor)! mg P/L KDW 07/23
         ! adsts_aer = 1000000 * ch_nut(jnut)%spkgc * (sfarea_aer / aer_vol) * & ! use entire aer_vol (not just pore volume) b/c porosity comes in via Keq equation, KDW
         !    (0.75 + 0.25 * dox_factor) ! mg P/L kg P/Mg sed  *  Mg clay equivilent / Liters = mg / L
         ! adsts_anaer = 1000000 * ch_nut(jnut)%spkgc * (sfarea_anaer / anaer_vol) * 0.75 
         !! Before diffusion, re-equilibrate dissolved v adsorbed because sediment erosion and deposition - KDW 03/03/22
         !Aerobic layer
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half)))) ! KDW debug 04/29/22 was 0.75 + 0.25*...
         P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
         xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
         yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
         P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
            4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
         if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
         if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
         if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
         ch(jrch)%aer_minpa = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
         ch(jrch)%aer_disp = (P_tot - (ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
         if (ch(jrch)%aer_disp < .0) ch(jrch)%aer_disp = 0.! kDW 04/28/22
         ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
         ! P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
         ! P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
         ! P_adsorbed = P_min_new - ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
         ! P_adsorbed = MIN(P_adsorbed, ch(jrch)%aer_disp * ch(jrch)%bed_por) !  "* ch(jrch)%bed_por" added 03/21/22
         ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd))
         ! ch(jrch)%aer_minpa =  ch(jrch)%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
         ! ch(jrch)%aer_disp = ch(jrch)%aer_disp - P_adsorbed / ch(jrch)%bed_por !  "/ ch(jrch)%bed_por" added 03/21/22
         !Anaerobic layer
         dox_factor = 0.5
         P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
         xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
         yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
         P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
            4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
         if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
         if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
         if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
         ch(jrch)%anaer_minpa = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
         ch(jrch)%anaer_disp = (P_tot - (ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
         if (ch(jrch)%anaer_disp < .0) ch(jrch)%anaer_disp = 0. ! kDW 04/28/22
         ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
         ! P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por
         ! P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
         ! P_adsorbed = P_min_new - ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd
         ! P_adsorbed = MIN(P_adsorbed, ch(jrch)%anaer_disp * ch(jrch)%bed_por)
         ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd))
         ! ch(jrch)%anaer_minpa =  ch(jrch)%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
         ! ch(jrch)%anaer_disp = ch(jrch)%anaer_disp - P_adsorbed / ch(jrch)%bed_por

         !! Before diffusion, re-equilibrate dissolved v adsorbed because inflow mixing
         wtclm_mixed_minpa = (minpacon * vol_0 + minpain * wtrin) ! KDW 02/24/22
         wtclm_mixed_minpa_conc = wtclm_mixed_minpa / qdin! KDW 02/24/22
         wtclm_mixed_minps = (minpscon * vol_0 + minpsin * wtrin) ! KDW 02/24/22
         wtclm_mixed_minps_conc = wtclm_mixed_minps / qdin! KDW 02/24/22
         wtclm_mixed_solp = (solpcon * vol_0 + dispin * wtrin)! KDW 02/24/22
         wtclm_mixed_solp_conc = wtclm_mixed_solp / qdin! KDW 02/24/22
         !wtclm_por = 1 - (sedin / qdin) / 1.2 ! KDW 04/18/22
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
         P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc ! wtclm_por approximately 1, so leave out here and below (in contrast to bed, for per TOTAL volume calc)
         xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2 ) ! KDW 04/15/22, KDW 04/18/22
         yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
         if(xx > 0.)then !KDW 04/18/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
         else
            P_min_new = 0.
         endif
         if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
         if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22 
         if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
         P_adsorbed_mass = (P_min_new * xx - wtclm_mixed_minpa_conc)* qdin ! KDW 04/15/22 
         wtclm_mixed_minpa_conc = P_min_new * xx ! KDW 04/15/22, back to g/ kg
         wtclm_mixed_solp_conc = P_tot - wtclm_mixed_minpa_conc ! KDW 04/15/22, back to g/m^3 pore water volume
         if (wtclm_mixed_solp_conc < .0) wtclm_mixed_solp_conc = 0.! kDW 04/28/22
         wtclm_mixed_minpa = wtclm_mixed_minpa_conc * qdin ! KDW 04/15/22
         wtclm_mixed_solp = wtclm_mixed_solp_conc * qdin ! KDW 04/15/22 
         ! P_min_new = (adsts_wtclm + P_tot + yyy - ((adsts_wtclm + P_tot + yyy)**2 - 4*P_tot*adsts_wtclm) **0.5) / 2
         ! P_adsorbed = P_min_new - wtclm_mixed_minpa_conc   
         ! P_adsorbed = MIN(P_adsorbed, wtclm_mixed_solp_conc )
         ! P_adsorbed = MAX(P_adsorbed, -wtclm_mixed_minpa_conc )
         ! P_adsorbed_mass = P_adsorbed * qdin  
         ! wtclm_mixed_minpa_conc = P_min_new ! KDW 03/03/22
         ! wtclm_mixed_minpa = wtclm_mixed_minpa_conc * qdin! KDW 03/03/22
         ! wtclm_mixed_solp_conc = wtclm_mixed_solp_conc - P_adsorbed ! KDW 03/03/22
         ! wtclm_mixed_solp = wtclm_mixed_solp_conc * qdin ! KDW 03/03/22
         if (P_adsorbed_mass >= 0.) then
            if (P_adsorbed_mass < solpcon * vol_0) then !all adsorption/desorption from stored dissolved and mineral P
               solpcon = (solpcon * vol_0 - P_adsorbed_mass) / (vol_0 + 0.0001)
               minpacon = (minpacon * vol_0 + P_adsorbed_mass) / (vol_0 + 0.0001)
            else
               ads_inflow = P_adsorbed_mass - solpcon * vol_0
               solpcon = 0.
               minpacon = (minpacon * vol_0 + (P_adsorbed_mass - ads_inflow)) / (vol_0 + 0.0001)
               dispin = (dispin * wtrin - ads_inflow) / wtrin
               minpain = (minpain * wtrin + ads_inflow) / wtrin
            endif
         else
            P_desorbed_mass = -P_adsorbed_mass
            if (P_desorbed_mass < minpacon * vol_0) then
               minpacon = (minpacon * vol_0 - P_desorbed_mass) / (vol_0 + 0.0001)
               solpcon = (solpcon * vol_0 + P_desorbed_mass) / (vol_0 + 0.0001)
            else
               ads_inflow = P_desorbed_mass - minpacon * vol_0
               minpacon = 0.
               solpcon = (solpcon * vol_0 + (P_desorbed_mass - ads_inflow)) / (vol_0 + 0.0001)
               minpain = (minpain * wtrin - ads_inflow) / wtrin
               dispin = (dispin * wtrin + ads_inflow) / wtrin
            endif
         endif
         
         !!Exchangable P via diffusion!!
         ! Initial PAI used to estimate how much diffusion sourced from dissolved versus labile/adsorbed P
         ! frsts_wtclm = adsts_wtclm - wtclm_mixed_minpa_conc
         ! frsts_aer = adsts_aer - ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd ! mg P/L - aer_minpa is mgP/kg and is multiplied by Mg/m^3 (i.e. mg/L)
         ! frsts_anaer = adsts_anaer - ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd
         ! if (frsts_wtclm > 0) then
         !    pai_wtclm = 1 - (1 / (1 + (wtclm_por/((ch_nut(jnut)%adskeq / 100) * frsts_wtclm)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
         ! else
         !    frsts_wtclm = 0
         !    pai_wtclm = 1
         ! endif
         ! if (frsts_aer > 0) then
         !    pai_aer = 1 - (1 / (1 + (ch(jrch)%bed_por/((ch_nut(jnut)%adskeq / 100) * frsts_aer)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
         ! else
         !    frsts_aer = 0
         !    pai_aer = 1
         ! endif
         ! if (frsts_anaer > 0) then
         !    pai_anaer = 1 - (1 / (1 + (ch(jrch)%bed_por/((ch_nut(jnut)%adskeq / 100) * frsts_anaer)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
         ! else
         !    frsts_anaer = 0
         !    pai_anaer = 1
         ! endif
         if((wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) > 0.)then ! KDW 04/18/22
            pai_wtclm = wtclm_mixed_solp_conc / (wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) !KDW 04/15/22
         else
            pai_wtclm = 1. !no sediment, so no minpa?
         endif
         if ((ch(jrch)%aer_disp * ch(jrch)%bed_por + ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd) > 0.) then! KDW 04/18/22
            pai_aer = ch(jrch)%aer_disp * ch(jrch)%bed_por / (ch(jrch)%aer_disp * ch(jrch)%bed_por + ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd) !KDW 04/15/22
         else
            pai_aer = 0.01
         endif
         if ((ch(jrch)%anaer_disp * ch(jrch)%bed_por + ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd) > 0.) then! KDW 04/18/22
            pai_anaer = ch(jrch)%anaer_disp * ch(jrch)%bed_por / (ch(jrch)%anaer_disp * ch(jrch)%bed_por + ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd) !KDW 04/15/22
         else
            pai_anaer = 0.01
         endif

         !!Diffusion between the aerobic layer and the water column
         disp_aer_init = ch(jrch)%aer_disp
         if (qdin > 0) then !change from ch(jrch)%rchstor to qdin, KDW 03/24/22
            diff_aer2wtclm = rs2_s * (ch(jrch)%aer_disp - wtclm_mixed_solp_conc) * rt_delt ! unit grams P, by convention positive is diffusion from aerobic layer to water column
            !Check for overshooting equilibrium
            if (diff_aer2wtclm <= 0) then  
               diff = abs(diff_aer2wtclm)
               max_diff = (wtclm_mixed_solp_conc - ch(jrch)%aer_disp) * & !labile/adsorbed P can compensate for strictly dissolved P
               (1/ ( (pai_aer/(vol_h2o_aer / 1000)) + (pai_wtclm/qdin) ) ) ! use vol_h2o_aer here b/c want true pore volume concentration, KDW, 03/21/22
               max_diff1 = (wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc) * qdin ! KDW 03/22/22

               max_diff_disp = pai_wtclm * max_diff ! KDW 03/22/22
               P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc ! KDW 03/22/22, g/m^3 bulk
               min_disp_conc = wtclm_mixed_solp_conc - (max_diff_disp/qdin) ! KDW 03/22/22, g/m^3 bulk
               dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
               xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
               yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
               min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
               max_diff_minpa = (P_tot - min_Ptot_conc)*qdin - max_diff_disp  ! KDW 04/15/22, grams
               max_diff_minpa = MIN(max_diff_minpa, ((1-pai_wtclm) * max_diff)) ! KDW 04/15/22
               max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
               ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
               ! min_Ptot_conc = min_disp_conc *(1+(adsts_wtclm/(yyy + min_disp_conc))) ! KDW 03/22/22
               ! max_diff_minpa = (P_tot - min_Ptot_conc)*qdin - max_diff_disp  ! KDW 03/22/22
               ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_wtclm) * max_diff)) ! KDW 03/22/22
               ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

               max_diff_disp = pai_aer * max_diff ! KDW 03/07/22
               P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! KDW 03/21/22, g/m^3 bulk
               max_disp_conc = ch(jrch)%aer_disp * ch(jrch)%bed_por + (max_diff_disp/(aer_vol/1000)) ! KDW 03/21/22, grams
               dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
               xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
               yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
               max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
               ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
               ! max_Ptot_conc = max_disp_conc *(1+(adsts_aer/(yyy + max_disp_conc))) ! KDW 03/21/22
               max_diff_minpa = (max_Ptot_conc - P_tot)*(aer_vol/1000) - max_diff_disp  ! KDW 03/21/22, grams
               max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 03/21/22
               max_diff3 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

               max_diff = MIN(max_diff, max_diff1, max_diff2, max_diff3) ! KDW 03/22/22
               max_diff = MAX(max_diff, 0.) ! KDW 03/22/22
               if(diff > max_diff) diff = max_diff
               if(diff < 0.) diff = 0. ! KDW 04/28/22
               diff_aer2wtclm = -diff
            endif
            if (diff_aer2wtclm >= 0) then  
               diff = abs(diff_aer2wtclm)
               max_diff = (ch(jrch)%aer_disp - wtclm_mixed_solp_conc) * & !labile/adsorbed P can compensate for strictly dissolved P
               (1/ ( (pai_wtclm/qdin) + (pai_aer/(vol_h2o_aer / 1000)) ) )
               max_diff1 = (ch(jrch)%aer_disp * ch(jrch)%bed_por + ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd) * aer_vol / 1000  ! KDW 03/21/22, added "* ch(jrch)%bed_por"
              
               max_diff_disp = pai_aer * max_diff ! KDW 03/22/22
               P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! KDW 03/22/22
               min_disp_conc = ch(jrch)%aer_disp * ch(jrch)%bed_por - (max_diff_disp/(aer_vol/1000)) ! KDW 03/22/22
               dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
               xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
               yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
               min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
               max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
               max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 04/15/22
               max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
               ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
               ! min_Ptot_conc = min_disp_conc *(1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/22/22
               ! max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 03/22/22
               ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 03/22/22
               ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22
              
               max_diff_disp = pai_wtclm * max_diff ! KDW 03/07/22
               P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc! KDW 03/21/22
               max_disp_conc = wtclm_mixed_solp_conc + (max_diff_disp/qdin) ! KDW 03/21/22
               dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
               xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
               yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22
               max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
               ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
               ! max_Ptot_conc = max_disp_conc *(1+(adsts_wtclm/(yyy + max_disp_conc))) ! KDW 03/21/22
               max_diff_minpa = (max_Ptot_conc - P_tot)*qdin - max_diff_disp  ! KDW 03/21/22
               max_diff_minpa = MIN(max_diff_minpa, ((1-pai_wtclm) * max_diff)) ! KDW 03/21/22

               max_diff3 =  max_diff_disp + max_diff_minpa ! KDW 03/07/22

               max_diff = MIN(max_diff, max_diff1, max_diff2, max_diff3) ! KDW 03/22/22
               max_diff = MAX(max_diff, 0.)
               if(diff > max_diff) diff = max_diff
               if(diff < 0.) diff = 0. ! KDW 04/28/22
               diff_aer2wtclm = diff
            endif
         else
            diff_aer2wtclm = 0
         endif

         !Update aerobic and water column concentrations
         if (diff_aer2wtclm < 0) then 
            diff = abs(diff_aer2wtclm)
            ch(jrch)%aer_disp = ch(jrch)%aer_disp + (diff*pai_aer/(vol_h2o_aer/1000))
            ch(jrch)%aer_minpa = ch(jrch)%aer_minpa + (diff*(1-pai_aer)/aer_mass)
            if((diff*pai_wtclm)<(dispin*wtrin))then
               dispin = dispin - (diff*pai_wtclm/wtrin)
               diffremain = diff - diff*pai_wtclm
            else
               diffremain = diff - (dispin*wtrin)
               dispin = 0
            endif
            if ((diff*(1-pai_wtclm))<(minpain*wtrin))then
               minpain = minpain - (diff*(1-pai_wtclm)/ wtrin)
               diffremain = diffremain - diff*(1-pai_wtclm)
            else
               diffremain = diffremain - (minpain*wtrin)
               minpain = 0
            endif
            diff = diffremain
            if(diff > 0.)then
               if((diff*pai_wtclm)<(solpcon*vol_0))then
                  solpcon = solpcon - (diff*pai_wtclm/(vol_0 + 0.0001))
                  diffremain = diff - diff*pai_wtclm
               else
                  diffremain = diff - (solpcon*vol_0)
                  solpcon = 0
               endif
               if ((diff*(1-pai_wtclm))<(minpacon*vol_0))then
                  minpacon = minpacon - (diff*(1-pai_wtclm)/ (vol_0 + 0.0001))
                  diffremain = diffremain - diff*(1-pai_wtclm)
               else
                  diffremain = diffremain - (minpacon*vol_0)
                  minpacon = 0
               endif
            endif
            if(diffremain > 0.) then !necessary because don't know fraction of P from inflow vs storage
               if(diffremain < (dispin*wtrin))then
                  dispin = dispin - (diffremain/wtrin)
                  diffremain = 0
               else
                  diffremain = diffremain - (dispin*wtrin)
                  dispin = 0
                  if(diffremain < (solpcon*vol_0))then
                     solpcon = solpcon - (diffremain/(vol_0+ 0.0001))
                     diffremain = 0
                  else
                     diffremain = diffremain - (solpcon*vol_0)
                     solpcon = 0
                     if(diffremain < (minpain*wtrin))then
                        minpain = minpain - (diffremain/wtrin)
                        diffremain = 0
                     else
                        diffremain = diffremain - (minpain*wtrin)
                        minpain = 0
                        if(diffremain < (minpacon*vol_0))then
                           minpacon = minpacon - (diffremain/(vol_0 + 0.0001))
                           diffremain = 0
                        else
                           diffremain = diffremain - (minpacon*vol_0)
                           minpacon = 0
                        endif
                     endif
                  endif
               endif
            endif 
         endif
         if (diff_aer2wtclm > 0) then
            diff = abs(diff_aer2wtclm)
            if ((diff*pai_wtclm/wtrin) < 10.) then ! in case wtrin is low, don't want excessive dispin causing instability
               dispin = dispin + (diff*pai_wtclm/wtrin)
               minpain = minpain + (diff*(1-pai_wtclm)/wtrin)
            else
               dispin = dispin + 10
               minpain = minpain + 10*(1-pai_wtclm)/pai_wtclm
               diffremain = diff - 10*(1/pai_wtclm) * wtrin
               solpcon = solpcon + (diffremain*pai_wtclm/(vol_0+ 0.0001))
               minpacon = minpacon + (diffremain*(1-pai_wtclm)/(vol_0 + 0.0001))
               diffremain = 0.
            endif
            ch(jrch)%aer_disp = ch(jrch)%aer_disp - (diff*pai_aer/(vol_h2o_aer/1000))
            ch(jrch)%aer_minpa = ch(jrch)%aer_minpa - (diff*(1-pai_aer)/aer_mass)
         endif
         
         !! Bring water column and aerobic layer to adsorption equilibrium before simulating anaerobic diffusion - KDW 03/22/22
         !Aerobic layer
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
            xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
            if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
            if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
            if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
            ch(jrch)%aer_minpa = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
            ch(jrch)%aer_disp = (P_tot - (ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
            if (ch(jrch)%aer_disp < .0) ch(jrch)%aer_disp = 0.! kDW 04/28/22
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
            ! P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
            ! P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
            ! P_adsorbed = P_min_new - ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
            ! P_adsorbed = MIN(P_adsorbed, ch(jrch)%aer_disp * ch(jrch)%bed_por) !  "* ch(jrch)%bed_por" added 03/21/22
            ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd))
            ! ch(jrch)%aer_minpa =  ch(jrch)%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
            ! ch(jrch)%aer_disp = ch(jrch)%aer_disp - P_adsorbed / ch(jrch)%bed_por !  "/ ch(jrch)%bed_por" added 03/21/22
         !Water column
            wtclm_mixed_minpa = (minpacon * vol_0 + minpain * wtrin) 
            wtclm_mixed_minpa_conc = wtclm_mixed_minpa / qdin
            wtclm_mixed_solp = (solpcon * vol_0 + dispin * wtrin)
            wtclm_mixed_solp_conc = wtclm_mixed_solp / qdin
            !yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) 
            P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc ! wtclm_por approximately 1, so leave out here and below (in contrast to bed, for per TOTAL volume calc)
            xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
            yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22
            if(xx > 0.)then ! KDW 04/18/22
               P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
                  4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
            else
               P_min_new = 0.
            endif
            if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
            if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22 
            if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
            P_adsorbed_mass = (P_min_new * xx - wtclm_mixed_minpa_conc)* qdin ! KDW 04/15/22 
            wtclm_mixed_minpa_conc = P_min_new * xx ! KDW 04/15/22, back to g/ kg
            wtclm_mixed_solp_conc = P_tot - wtclm_mixed_minpa_conc ! KDW 04/15/22, back to g/m^3 pore water volume
            if (wtclm_mixed_solp_conc < .0) wtclm_mixed_solp_conc = 0.! kDW 04/28/22
            wtclm_mixed_minpa = wtclm_mixed_minpa_conc * qdin ! KDW 04/15/22
            wtclm_mixed_solp = wtclm_mixed_solp_conc * qdin ! KDW 04/15/22             
            ! P_min_new = (adsts_wtclm + P_tot + yyy - ((adsts_wtclm + P_tot + yyy)**2 - 4*P_tot*adsts_wtclm) **0.5) / 2
            ! P_adsorbed = P_min_new - wtclm_mixed_minpa_conc   
            ! P_adsorbed = MIN(P_adsorbed, wtclm_mixed_solp_conc )
            ! P_adsorbed = MAX(P_adsorbed, -wtclm_mixed_minpa_conc )
            ! P_adsorbed_mass = P_adsorbed * qdin  
            ! wtclm_mixed_minpa_conc = P_min_new 
            ! wtclm_mixed_minpa = wtclm_mixed_minpa_conc * qdin
            ! wtclm_mixed_solp_conc = wtclm_mixed_solp_conc - P_adsorbed 
            ! wtclm_mixed_solp = wtclm_mixed_solp_conc * qdin 
            if (P_adsorbed_mass >= 0.) then
               if (P_adsorbed_mass < solpcon * vol_0) then !all adsorption/desorption from stored dissolved and mineral P
                  solpcon = (solpcon * vol_0 - P_adsorbed_mass) / (vol_0 + 0.0001)
                  minpacon = (minpacon * vol_0 + P_adsorbed_mass) / (vol_0 + 0.0001)
               else
                  ads_inflow = P_adsorbed_mass - solpcon * vol_0
                  solpcon = 0.
                  minpacon = (minpacon * vol_0 + (P_adsorbed_mass - ads_inflow)) / (vol_0 + 0.0001)
                  dispin = (dispin * wtrin - ads_inflow) / wtrin
                  minpain = (minpain * wtrin + ads_inflow) / wtrin
               endif
            else
               P_desorbed_mass = -P_adsorbed_mass
               if (P_desorbed_mass < minpacon * vol_0) then
                  minpacon = (minpacon * vol_0 - P_desorbed_mass) / (vol_0 + 0.0001)
                  solpcon = (solpcon * vol_0 + P_desorbed_mass) / (vol_0 + 0.0001)
               else
                  ads_inflow = P_desorbed_mass - minpacon * vol_0
                  minpacon = 0.
                  solpcon = (solpcon * vol_0 + (P_desorbed_mass - ads_inflow)) / (vol_0 + 0.0001)
                  minpain = (minpain * wtrin - ads_inflow) / wtrin
                  dispin = (dispin * wtrin + ads_inflow) / wtrin
               endif
            endif
                     
         !!Diffusion between the aerobic and anaerobic layers (or anaerobic through aerobic to water column)
         diff_aer2anaer = rs2_s * (ch(jrch)%aer_disp - ch(jrch)%anaer_disp) * rt_delt * rs2_ratio ! unit grams P, / 10 added 03/26/22, KDW
         !Update mixed wtclm concentrations before calculating diffusion with anaerobic
         ! wtclm_mixed_minpa = (minpacon * vol_0 + minpain * wtrin)
         ! wtclm_mixed_minpa_conc = wtclm_mixed_minpa / qdin
         ! wtclm_mixed_solp = (solpcon * vol_0 + dispin * wtrin)
         ! wtclm_mixed_solp_conc = wtclm_mixed_solp / qdin
         ! frsts_wtclm = adsts_wtclm - wtclm_mixed_minpa_conc
         ! frsts_aer = adsts_aer - ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd
         ! if (frsts_wtclm > 0) then
         !    pai_wtclm = 1 - (1 / (1 + (wtclm_por/((ch_nut(jnut)%adskeq / 100) * frsts_wtclm)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
         ! else
         !    frsts_wtclm = 0
         !    pai_wtclm = 1
         ! endif
         ! if (frsts_aer > 0) then
         !    pai_aer = 1 - (1 / (1 + (ch(jrch)%bed_por/((ch_nut(jnut)%adskeq / 100) * frsts_aer)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
         ! else
         !    frsts_aer = 0
         !    pai_aer = 1
         ! endif
         if ((wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) > 0.) then ! KDW 04/18/22
            pai_wtclm = wtclm_mixed_solp_conc / (wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) !KDW 04/15/22
         else
            pai_wtclm = 1.
         endif
         if ((ch(jrch)%aer_disp * ch(jrch)%bed_por + ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd) > 0.) then ! KDW 04/18/22
            pai_aer = ch(jrch)%aer_disp * ch(jrch)%bed_por / (ch(jrch)%aer_disp * ch(jrch)%bed_por + ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd) !KDW 04/15/22
         else
            pai_aer = 0.01
         endif
         
         !! Calculate aer conc where flows to/from wtclm and anaer are balance 
         aer_balanced = (wtclm_mixed_solp_conc + ch(jrch)%anaer_disp * rs2_ratio) / (1 +rs2_ratio) ! KDW 03/30/22
         !! if anaer conc. is between aer and wtclm, allow aer to rise/fall to level of anaer (if sufficient diffusion potential)
         if ((ch(jrch)%aer_disp < ch(jrch)%anaer_disp) .and. (ch(jrch)%anaer_disp <= wtclm_mixed_solp_conc)) then
            diff = abs(diff_aer2anaer)
            max_diff = (ch(jrch)%anaer_disp - ch(jrch)%aer_disp) * &
            (1/ ( (pai_aer/(vol_h2o_aer / 1000)) + (pai_aer/(vol_h2o_anaer / 1000)) ) )
            max_diff1 = (ch(jrch)%anaer_disp * ch(jrch)%bed_por + ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd) * anaer_vol / 1000 

            max_diff_disp = pai_anaer * max_diff ! KDW 03/22/22
            P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por ! KDW 03/22/22
            min_disp_conc = ch(jrch)%anaer_disp * ch(jrch)%bed_por - (max_diff_disp/(anaer_vol/1000)) ! KDW 03/22/22
            dox_factor = 0.5 
            xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            max_diff_minpa = (P_tot - min_Ptot_conc)*(anaer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 04/15/22
            max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
            ! min_Ptot_conc = min_disp_conc *(1+(adsts_anaer/(yyy + min_disp_conc))) ! KDW 03/22/22
            ! max_diff_minpa = (P_tot - min_Ptot_conc)*(anaer_vol/1000) - max_diff_disp  ! KDW 03/22/22
            ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 03/22/22
            ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

            max_diff_disp = pai_aer * max_diff ! KDW 03/07/22
            P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! KDW 03/21/22
            max_disp_conc = ch(jrch)%aer_disp * ch(jrch)%bed_por + (max_diff_disp/(aer_vol/1000)) ! KDW 03/21/22
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
            ! max_Ptot_conc = max_disp_conc *(1+(adsts_aer/(yyy + max_disp_conc))) ! KDW 03/21/22
            max_diff_minpa = (max_Ptot_conc - P_tot)*(aer_vol/1000) - max_diff_disp  ! KDW 03/21/22
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 03/21/22
            max_diff3 =  max_diff_disp + max_diff_minpa ! KDW 03/07/22

            max_diff = MIN(max_diff, max_diff1, max_diff2, max_diff3) ! KDW 03/22/22
            max_diff = MAX(max_diff, 0.)

            if(diff > max_diff) diff = max_diff 
            if(diff < 0.) diff = 0. ! KDW 04/28/22
            diff_anaer2wtclm = 0. ! KDW 03/30/22
            diff_aer2anaer = -diff             
         elseif ((ch(jrch)%aer_disp > ch(jrch)%anaer_disp) .and. (ch(jrch)%anaer_disp >= wtclm_mixed_solp_conc)) then
            diff = abs(diff_aer2anaer)
            max_diff = (ch(jrch)%aer_disp - ch(jrch)%anaer_disp) * &
            (1/ ( (pai_aer/(vol_h2o_aer / 1000)) + (pai_anaer/(vol_h2o_anaer / 1000)) ) )
            max_diff1 = (ch(jrch)%aer_disp * ch(jrch)%bed_por + ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd) * aer_vol / 1000

            max_diff_disp = pai_aer * max_diff ! KDW 03/22/22
            P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! KDW 03/22/22
            min_disp_conc = ch(jrch)%aer_disp * ch(jrch)%bed_por - (max_diff_disp/(aer_vol/1000)) ! KDW 03/22/22
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 04/15/22
            max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
            ! min_Ptot_conc = min_disp_conc *(1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/22/22
            ! max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 03/22/22
            ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 03/22/22
            ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

            max_diff_disp = pai_anaer * max_diff ! KDW 03/07/22
            P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por ! KDW 03/21/22
            max_disp_conc = ch(jrch)%anaer_disp * ch(jrch)%bed_por + (max_diff_disp/(anaer_vol/1000)) ! KDW 03/21/22
            dox_factor = 0.5
            xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
            ! max_Ptot_conc = max_disp_conc *(1+(adsts_anaer/(yyy + max_disp_conc))) ! KDW 03/21/22
            max_diff_minpa = (max_Ptot_conc - P_tot)*(anaer_vol/1000) - max_diff_disp  ! KDW 03/21/22
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 03/21/22
            max_diff3 =  max_diff_disp + max_diff_minpa ! KDW 03/07/22

            max_diff = MIN(max_diff, max_diff1, max_diff2, max_diff3) ! KDW 03/22/22
            max_diff = MAX(max_diff, 0.)
            if(diff > max_diff) diff = max_diff
            if(diff < 0.) diff = 0. ! KDW 04/28/22
            diff_anaer2wtclm = 0. ! KDW 03/30/22
            diff_aer2anaer = diff
         else
            !! bring aer conc toward balancing conc via exchange with anaer
            if (ch(jrch)%anaer_disp < wtclm_mixed_solp_conc) then
               if (ch(jrch)%aer_disp < aer_balanced) then
               !! disP will move through aer between anaer and wtclm
                  diff_anaer2wtclm = -diff_aer2anaer 
                  diff_aer2anaer  = 0.
               else
               !! allow aer_disp to be brought down as low as a aer_balanced
                  ! how much diffusion to bring aer to aer_balanced
                  P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! KDW 03/21/22
                  min_disp_conc = aer_balanced * ch(jrch)%bed_por ! KDW 03/21/22, 04/02/22
                  dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
                  xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
                  yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
                  min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
                  ! min_Ptot_conc = min_disp_conc *(1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/21/22
                  max_diff2 = (P_tot - min_Ptot_conc)*(aer_vol/1000)  ! KDW 03/21/22
                  ! make sure doesn't cause P in anaer to exceed aer_balance
                  P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por ! KDW 04/04/22
                  max_disp_conc = aer_balanced * ch(jrch)%bed_por ! KDW 04/04/22
                  dox_factor = 0.5
                  xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
                  yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
                  max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
                  ! max_Ptot_conc = max_disp_conc *(1+(adsts_anaer/(yyy + max_disp_conc))) !KDW 04/04/22
                  max_diff = (max_Ptot_conc - P_tot)*(anaer_vol/1000)  ! KDW 04/04/22
                  if (max_diff > max_diff2) max_diff = max_diff2 ! KDW 04/04/22
                  if(max_diff < 0.) max_diff = 0. ! KDW 04/28/22
                  if (diff_aer2anaer > max_diff) then
                     diff_anaer2wtclm = -(diff_aer2anaer - max_diff)
                     diff_aer2anaer = max_diff
                  else
                     diff_anaer2wtclm = 0.
                  endif
               endif
            else !i.e. ch(jrch)%anaer_disp >= wtclm_mixed_solp_conc
               if (ch(jrch)%aer_disp > aer_balanced) then
               !! disP will move through aer between anaer and wtclm
                  diff_anaer2wtclm = -diff_aer2anaer 
                  diff_aer2anaer  = 0.
               else
               !! allow aer_disp to be brought up as high as a aer_balanced
                  ! how much diffusion to bring aer to aer_balanced                 
                  P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! KDW 03/21/22
                  max_disp_conc = aer_balanced * ch(jrch)%bed_por ! KDW 03/21/22, 04/02/22
                  dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
                  xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
                  yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
                  max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
                  ! max_Ptot_conc = max_disp_conc *(1+(adsts_aer/(yyy + max_disp_conc))) ! KDW 03/21/22
                  max_diff2 = (max_Ptot_conc - P_tot)*(aer_vol/1000)  ! KDW 03/21/22
                  ! make sure enough P in anaer to supply diffusion
                  P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por ! KDW 04/04/22
                  min_disp_conc = aer_balanced * ch(jrch)%bed_por ! KDW 04/04/22
                  dox_factor = 0.5 
                  xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
                  yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
                  min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
                  ! min_Ptot_conc = min_disp_conc *(1+(adsts_anaer/(yyy + min_disp_conc))) !KDW 04/04/22
                  max_diff = (P_tot - min_Ptot_conc)*(anaer_vol/1000)  ! KDW 04/04/22
                  if (max_diff > max_diff2) max_diff = max_diff2 ! KDW 04/04/22
                  if(max_diff < 0.) max_diff = 0. ! KDW 04/28/22
                  if (-diff_aer2anaer > max_diff) then
                     diff_anaer2wtclm = -(diff_aer2anaer + max_diff)
                     diff_aer2anaer = -max_diff
                  else
                     diff_anaer2wtclm = 0
                  endif
               endif
            endif
            !! update aerobic and anaerobic conc (use temp variables)
            aer_disp_temp = ch(jrch)%aer_disp - diff_aer2anaer / (vol_h2o_aer/1000)
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + aer_disp_temp * ch(jrch)%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22, 04/27/22
            xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
            if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
            if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
            if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
            aer_minpa_temp = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
            aer_disp_temp = (P_tot - (aer_minpa_temp * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume  
            if (aer_disp_temp < .0) aer_disp_temp = 0.! kDW 04/28/22         
            ! aer_disp_temp = ch(jrch)%aer_disp - diff_aer2anaer / (vol_h2o_aer/1000)
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
            ! P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + aer_disp_temp * ch(jrch)%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
            ! P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
            ! P_adsorbed = P_min_new - ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
            ! P_adsorbed = MIN(P_adsorbed, aer_disp_temp * ch(jrch)%bed_por) !  "* ch(jrch)%bed_por" added 03/21/22
            ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd))
            ! aer_minpa_temp =  ch(jrch)%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
            ! aer_disp_temp = aer_disp_temp - P_adsorbed / ch(jrch)%bed_por !  "/ ch(jrch)%bed_por" added 03/21/22
            
            anaer_disp_temp = ch(jrch)%anaer_disp + diff_aer2anaer / (vol_h2o_anaer/1000)
            dox_factor = 0.5 
            P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + anaer_disp_temp * ch(jrch)%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22, 04/27/22
            xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
            if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
            if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
            if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
            anaer_minpa_temp = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
            anaer_disp_temp = (P_tot - (anaer_minpa_temp * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume  
            if (anaer_disp_temp < .0) anaer_disp_temp = 0.! kDW 04/28/22         
            ! anaer_disp_temp = ch(jrch)%anaer_disp + diff_aer2anaer / (vol_h2o_anaer/1000)
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
            ! P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + anaer_disp_temp * ch(jrch)%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
            ! P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
            ! P_adsorbed = P_min_new - ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
            ! P_adsorbed = MIN(P_adsorbed, anaer_disp_temp * ch(jrch)%bed_por) !  "* ch(jrch)%bed_por" added 03/21/22
            ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd))
            ! anaer_minpa_temp =  ch(jrch)%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
            ! anaer_disp_temp = anaer_disp_temp - P_adsorbed / ch(jrch)%bed_por !  "/ ch(jrch)%bed_por" added 03/21/22
         endif
         
         !! allow exchange between anaer and wtclm, if diffusion potential available
         if(diff_anaer2wtclm > 0) then
            !! check diffusion for overshooting equilibrium
            diff = abs(diff_anaer2wtclm)
            max_diff = (anaer_disp_temp - wtclm_mixed_solp_conc) * &
            (1/ ( (pai_wtclm/qdin) + (pai_anaer/(vol_h2o_anaer / 1000)) ) )
            max_diff1 = (anaer_disp_temp * ch(jrch)%bed_por + anaer_minpa_temp * ch_sed(jsed)%bed_bd) * anaer_vol / 1000

            max_diff_disp = pai_anaer * max_diff ! KDW 03/22/22
            P_tot = anaer_minpa_temp * ch_sed(jsed)%bed_bd + anaer_disp_temp * ch(jrch)%bed_por ! KDW 03/22/22
            min_disp_conc = anaer_disp_temp * ch(jrch)%bed_por - (max_diff_disp/(anaer_vol/1000)) ! KDW 03/22/22
            dox_factor = 0.5 
            xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            max_diff_minpa = (P_tot - min_Ptot_conc)*(anaer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 04/15/22
            max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
            ! min_Ptot_conc = min_disp_conc *(1+(adsts_anaer/(yyy + min_disp_conc))) ! KDW 03/22/22
            ! max_diff_minpa = (P_tot - min_Ptot_conc)*(anaer_vol/1000) - max_diff_disp  ! KDW 03/22/22
            ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 03/22/22
            ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

            max_diff_disp = pai_wtclm * max_diff ! KDW 03/07/22
            P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc! KDW 03/21/22
            max_disp_conc = wtclm_mixed_solp_conc + (max_diff_disp/qdin) ! KDW 03/21/22
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
            yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
            max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
            ! max_Ptot_conc = max_disp_conc *(1+(adsts_wtclm/(yyy + max_disp_conc))) ! KDW 03/21/22
            max_diff_minpa = (max_Ptot_conc - P_tot)*qdin - max_diff_disp  ! KDW 03/21/22
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_wtclm) * max_diff)) ! KDW 03/21/22
            max_diff3 =  max_diff_disp + max_diff_minpa ! KDW 03/07/22

            max_diff = MIN(max_diff, max_diff1, max_diff2, max_diff3) ! KDW 03/22/22
            max_diff = MAX(max_diff, 0.)
            if(diff > max_diff) diff = max_diff
            if(diff < 0.) diff = 0. ! KDW 04/28/22
            diff_anaer2wtclm = diff

         elseif (diff_anaer2wtclm < 0) then
            !! check diffusion for overshooting equilibrium
            diff = abs(diff_anaer2wtclm)
            max_diff = (wtclm_mixed_solp_conc - anaer_disp_temp) * &
            (1/ ( (pai_wtclm/qdin) + (pai_anaer/(vol_h2o_anaer / 1000)) ) )
            max_diff1 = (wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) * qdin

            max_diff_disp = pai_wtclm * max_diff ! KDW 03/22/22
            P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc ! KDW 03/22/22
            min_disp_conc = wtclm_mixed_solp_conc - (max_diff_disp/qdin) ! KDW 03/22/22
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
            yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
            min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
            ! min_Ptot_conc = min_disp_conc *(1+(adsts_wtclm/(yyy + min_disp_conc))) ! KDW 03/22/22
            max_diff_minpa = (P_tot - min_Ptot_conc)*qdin - max_diff_disp  ! KDW 03/22/22
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_wtclm) * max_diff)) ! KDW 03/22/22
            max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

            max_diff_disp = pai_anaer * max_diff ! KDW 03/07/22
            P_tot = anaer_minpa_temp * ch_sed(jsed)%bed_bd + anaer_disp_temp * ch(jrch)%bed_por ! KDW 03/21/22
            max_disp_conc = anaer_disp_temp * ch(jrch)%bed_por + (max_diff_disp/(anaer_vol/1000)) ! KDW 03/21/22
            dox_factor = 0.5
            xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
            ! max_Ptot_conc = max_disp_conc *(1+(adsts_anaer/(yyy + max_disp_conc))) ! KDW 03/21/22
            max_diff_minpa = (max_Ptot_conc - P_tot)*(anaer_vol/1000) - max_diff_disp  ! KDW 03/21/22
            max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 03/21/22
            max_diff3 =  max_diff_disp + max_diff_minpa ! KDW 03/07/22

            max_diff = MIN(max_diff, max_diff1, max_diff2, max_diff3) ! KDW 03/22/22
            max_diff = MAX(max_diff, 0.)
            if(diff > max_diff) diff = max_diff
            if(diff < 0.) diff = 0. ! KDW 04/28/22
            diff_anaer2wtclm = -diff
         endif

         !! update wtclm_mixed_solp_conc and ch(jrch)%anaer_disp (use temp variables)
         wtclm_mixed_solp_conc_temp = wtclm_mixed_solp_conc + diff_anaer2wtclm / qdin
         !yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) 
         P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc_temp ! wtclm_por approximately 1, so leave out here and below (in contrast to bed, for per TOTAL volume calc)         
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
         xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
         yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
         if(xx > 0.)then ! KDW 04/18/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
         else
            P_min_new = 0.
         endif
         if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
         if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22 
         if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
         P_adsorbed_mass = (P_min_new * xx - wtclm_mixed_minpa_conc)* qdin ! KDW 04/15/22 
         wtclm_mixed_minpa_conc_temp = P_min_new * xx ! KDW 04/15/22, back to g/ kg
         wtclm_mixed_solp_conc_temp = P_tot - wtclm_mixed_minpa_conc ! KDW 04/15/22, back to g/m^3 pore water volume     
         if (wtclm_mixed_solp_conc_temp < .0) wtclm_mixed_solp_conc_temp = 0.! kDW 04/28/22    
         ! P_min_new = (adsts_wtclm + P_tot + yyy - ((adsts_wtclm + P_tot + yyy)**2 - 4*P_tot*adsts_wtclm) **0.5) / 2
         ! P_adsorbed = P_min_new - wtclm_mixed_minpa_conc   
         ! P_adsorbed = MIN(P_adsorbed, wtclm_mixed_solp_conc_temp )
         ! P_adsorbed = MAX(P_adsorbed, -wtclm_mixed_minpa_conc )
         ! P_adsorbed_mass = P_adsorbed * qdin  
         ! wtclm_mixed_minpa_conc = P_min_new 
         ! wtclm_mixed_solp_conc_temp = wtclm_mixed_solp_conc_temp - P_adsorbed 

         anaer_disp_temp = anaer_disp_temp + diff_aer2anaer / (vol_h2o_anaer/1000)
         !yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
         P_tot = anaer_minpa_temp * ch_sed(jsed)%bed_bd + anaer_disp_temp * ch(jrch)%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
         dox_factor = 0.5
         xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
         yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
         P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
            4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
         if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
         if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
         if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
         anaer_minpa_temp = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
         anaer_disp_temp = (P_tot - (anaer_minpa_temp * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume     
         if (anaer_disp_temp < .0) anaer_disp_temp = 0.! kDW 04/28/22    
         ! P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
         ! P_adsorbed = P_min_new - ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
         ! P_adsorbed = MIN(P_adsorbed, anaer_disp_temp * ch(jrch)%bed_por) !  "* ch(jrch)%bed_por" added 03/21/22
         ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd))
         ! anaer_minpa_temp =  ch(jrch)%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
         ! anaer_disp_temp = anaer_disp_temp - P_adsorbed / ch(jrch)%bed_por !  "/ ch(jrch)%bed_por" added 03/21/22

         if(abs(diff_anaer2wtclm) > 0) then
            !! find new rs2_balanced aer conc
            aer_balanced = (wtclm_mixed_solp_conc_temp + anaer_disp_temp * rs2_ratio) / (1 +rs2_ratio) ! KDW 03/30/22

            !! allow aer_disp to reach aer_balanced again
            P_tot = aer_minpa_temp * ch_sed(jsed)%bed_bd + aer_disp_temp * ch(jrch)%bed_por ! KDW 03/21/22
            min_disp_conc = aer_balanced * ch(jrch)%bed_por ! KDW 03/21/22, 04/02/22
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            Ptot_conc_balanced = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
            ! Ptot_conc_balanced = min_disp_conc * (1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/21/22
            diff = (P_tot - Ptot_conc_balanced)*(aer_vol/1000)  ! KDW 03/21/22
            !! assume that change in aer_disp (to aer_balanced) is due to flow with whichever pool is in correct direction (anaerobic or watercolumn)
            !!!!!!!!! need to check if sufficient P available in source pool!!!!!!!!!!!!
            if (diff > 0) then ! i.e. P leaving aer
               if (wtclm_mixed_solp_conc_temp > anaer_disp_temp) then
                  ! make sure enough anaer doesn't overshoot aer_balanced
                  P_tot = anaer_minpa_temp * ch_sed(jsed)%bed_bd + anaer_disp_temp * ch(jrch)%bed_por ! KDW 04/04/22
                  max_disp_conc = aer_balanced * ch(jrch)%bed_por ! KDW 04/04/22
                  dox_factor = 0.5
                  xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
                  yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
                  max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
                  ! max_Ptot_conc = max_disp_conc *(1+(adsts_anaer/(yyy + max_disp_conc))) !KDW 04/04/22
                  max_diff = (max_Ptot_conc - P_tot)*(anaer_vol/1000)  ! KDW 04/04/22
                  if (diff > max_diff) diff = max_diff ! KDW 04/04/22
                  if(diff < 0.) diff = 0. ! KDW 04/28/22
                  diff_aer2anaer = diff_aer2anaer + diff
               else
                  ! make sure enough wtclm doesn't overshoot aer_balanced
                  P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc_temp ! KDW 04/04/22
                  max_disp_conc = aer_balanced ! KDW 04/04/22
                  dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
                  xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
                  yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
                  max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
                  ! max_Ptot_conc = max_disp_conc *(1+(adsts_wtclm/(yyy + max_disp_conc))) !KDW 04/04/22
                  max_diff = (max_Ptot_conc - P_tot)*qdin  ! KDW 04/04/22
                  if (diff > max_diff) diff = max_diff ! KDW 04/04/22
                  if(diff < 0.) diff = 0. ! KDW 04/28/22
                  diff_aer2wtclm2 = diff
               endif
            endif
            if (diff < 0) then ! i.e. P entering aer
               if (wtclm_mixed_solp_conc_temp < anaer_disp_temp) then
                  ! make sure enough P in anaer to supply diffusion
                  P_tot = anaer_minpa_temp * ch_sed(jsed)%bed_bd + anaer_disp_temp * ch(jrch)%bed_por ! KDW 04/04/22
                  min_disp_conc = aer_balanced * ch(jrch)%bed_por ! KDW 04/04/22
                  dox_factor = 0.5 
                  xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
                  yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
                  min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
                  ! min_Ptot_conc = min_disp_conc *(1+(adsts_anaer/(yyy + min_disp_conc))) !KDW 04/04/22
                  max_diff = (P_tot - min_Ptot_conc)*(anaer_vol/1000)  ! KDW 04/04/22
                  if (-diff > max_diff) diff = -max_diff ! KDW 04/04/22
                  if(diff > 0.) diff = 0. ! KDW 04/28/22
                  diff_aer2anaer = diff_aer2anaer + diff
               else
                  ! make sure enough P in wtclm to supply diffusion
                  P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc_temp ! KDW 04/04/22
                  min_disp_conc = aer_balanced ! KDW 04/04/22
                  dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
                  xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
                  yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
                  min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
                  ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
                  ! min_Ptot_conc = min_disp_conc *(1+(adsts_wtclm/(yyy + min_disp_conc))) !KDW 04/04/22
                  max_diff = (P_tot - min_Ptot_conc)*qdin  ! KDW 04/04/22
                  if (-diff > max_diff) diff = -max_diff ! KDW 04/04/22
                  if(diff > 0.) diff = 0. ! KDW 04/28/22
                  diff_aer2wtclm2 = diff
               endif
            endif
         endif
         
         
         !Update water column, aerobic, and anaerobic concentrations
         if (diff_anaer2wtclm + diff_aer2wtclm2 < 0) then ! KDW 04/01/22, added diff_aer2wtclm2
            diff = abs(diff_anaer2wtclm + diff_aer2wtclm2) ! KDW 04/01/22, added diff_aer2wtclm2
            ch(jrch)%anaer_disp = ch(jrch)%anaer_disp + (-diff_anaer2wtclm*pai_anaer/(vol_h2o_anaer/1000)) ! KDW 04/01/22, changed from "diff" to "-diff_anaer2wtclm"
            ch(jrch)%anaer_minpa = ch(jrch)%anaer_minpa + (-diff_anaer2wtclm*(1-pai_anaer)/anaer_mass)! KDW 04/01/22, changed from "diff" to "-diff_anaer2wtclm"
            ch(jrch)%aer_disp = ch(jrch)%aer_disp + (-diff_aer2wtclm2*pai_aer/(vol_h2o_aer/1000))! KDW 04/01/22
            ch(jrch)%aer_minpa = ch(jrch)%aer_minpa + (-diff_aer2wtclm2*(1-pai_aer)/aer_mass)! KDW 04/01/22

            if((diff*pai_wtclm)<(dispin*wtrin))then
               dispin = dispin - (diff*pai_wtclm/wtrin)
               diffremain = diff - diff*pai_wtclm
            else
               diffremain = diff - (dispin*wtrin)
               dispin = 0
            endif
            if ((diff*(1-pai_wtclm))<(minpain*wtrin))then
               minpain = minpain - (diff*(1-pai_wtclm)/ wtrin)
               diffremain = diffremain - diff*(1-pai_wtclm)
            else
               diffremain = diffremain - (minpain*wtrin)
               minpain = 0
            endif
            diff = diffremain
            if(diff > 0.)then
               if((diff*pai_wtclm)<(solpcon*vol_0))then
                  solpcon = solpcon - (diff*pai_wtclm/(vol_0 + 0.0001))
                  diffremain = diff - diff*pai_wtclm
               else
                  diffremain = diff - (solpcon*vol_0)
                  solpcon = 0
               endif
               if ((diff*(1-pai_wtclm))<(minpacon*vol_0))then
                  minpacon = minpacon - (diff*(1-pai_wtclm)/ (vol_0 + 0.0001))
                  diffremain = diffremain - diff*(1-pai_wtclm)
               else
                  diffremain = diffremain - (minpacon*vol_0)
                  minpacon = 0
               endif
            endif
            if(diffremain > 0.) then !necessary because don't know fraction of P from inflow vs storage
               if(diffremain < (dispin*wtrin))then
                  dispin = dispin - (diffremain/wtrin)
                  diffremain = 0
               else
                  diffremain = diffremain - (dispin*wtrin)
                  dispin = 0
                  if(diffremain < (solpcon*vol_0))then
                     solpcon = solpcon - (diffremain/(vol_0 + 0.0001))
                     diffremain = 0
                  else
                     diffremain = diffremain - (solpcon*vol_0)
                     solpcon = 0
                     if(diffremain < (minpain*wtrin))then
                        minpain = minpain - (diffremain/wtrin)
                        diffremain = 0
                     else
                        diffremain = diffremain - (minpain*wtrin)
                        minpain = 0
                        if(diffremain < (minpacon*vol_0))then
                           minpacon = minpacon - (diffremain/(vol_0 + 0.0001))
                           diffremain = 0
                        else
                           diffremain = diffremain - (minpacon*vol_0)
                           minpacon = 0
                        endif
                     endif
                  endif
               endif
            endif 
         endif
         if (diff_anaer2wtclm + diff_aer2wtclm2 > 0) then! KDW 04/01/22, added diff_aer2wtclm2
            diff = abs(diff_anaer2wtclm + diff_aer2wtclm2)! KDW 04/01/22, added diff_aer2wtclm2
            if ((diff*pai_wtclm/wtrin) < 10.) then ! in case wtrin is low, don't want excessive dispin causing instability
               dispin = dispin + (diff*pai_wtclm/wtrin)
               minpain = minpain + (diff*(1-pai_wtclm)/wtrin)
            else
               dispin = dispin + 10
               minpain = minpain + 10*(1-pai_wtclm)/pai_wtclm
               diffremain = diff - 10*(1/pai_wtclm) * wtrin
               solpcon = solpcon + (diffremain*pai_wtclm/(vol_0 + 0.0001))
               minpacon = minpacon + (diffremain*(1-pai_wtclm)/(vol_0 + 0.0001))
               diffremain = 0.
            endif
            ch(jrch)%anaer_disp = ch(jrch)%anaer_disp - (diff_anaer2wtclm*pai_anaer/(vol_h2o_anaer/1000))! KDW 04/01/22, changed from "diff" to "diff_anaer2wtclm"
            ch(jrch)%anaer_minpa = ch(jrch)%anaer_minpa - (diff_anaer2wtclm*(1-pai_anaer)/anaer_mass)! KDW 04/01/22, changed from "diff" to "diff_anaer2wtclm"
            ch(jrch)%aer_disp = ch(jrch)%aer_disp - (diff_aer2wtclm2*pai_aer/(vol_h2o_aer/1000))! KDW 04/01/22
            ch(jrch)%aer_minpa = ch(jrch)%aer_minpa - (diff_aer2wtclm2*(1-pai_aer)/aer_mass)! KDW 04/01/22
         endif
         if (diff_aer2anaer > 0) then
            diff = abs(diff_aer2anaer)
            ch(jrch)%anaer_disp = ch(jrch)%anaer_disp + (diff*pai_anaer/(vol_h2o_anaer/1000))
            ch(jrch)%anaer_minpa = ch(jrch)%anaer_minpa + (diff*(1-pai_anaer)/anaer_mass)
            ch(jrch)%aer_disp = ch(jrch)%aer_disp - (diff*pai_aer/(vol_h2o_aer/1000))
            ch(jrch)%aer_minpa = ch(jrch)%aer_minpa - (diff*(1-pai_aer)/aer_mass)
         endif
         if (diff_aer2anaer < 0) then
            diff = abs(diff_aer2anaer)
            ch(jrch)%anaer_disp = ch(jrch)%anaer_disp - (diff*pai_anaer/(vol_h2o_anaer/1000))
            ch(jrch)%anaer_minpa = ch(jrch)%anaer_minpa - (diff*(1-pai_anaer)/anaer_mass)
            ch(jrch)%aer_disp = ch(jrch)%aer_disp + (diff*pai_aer/(vol_h2o_aer/1000))
            ch(jrch)%aer_minpa = ch(jrch)%aer_minpa + (diff*(1-pai_aer)/aer_mass)
         endif

         !!Set all pools to minimum of zero
         dispin = MAX(0.,dispin)
         solpcon = MAX(0.,solpcon)
         minpain = MAX(0.,minpain)
         minpacon = MAX(0.,minpacon)
         ch(jrch)%aer_disp = MAX(0.,ch(jrch)%aer_disp)
         ch(jrch)%aer_minpa = MAX(0.,ch(jrch)%aer_minpa)
         ch(jrch)%anaer_disp = MAX(0.,ch(jrch)%anaer_disp)
         ch(jrch)%anaer_minpa = MAX(0.,ch(jrch)%anaer_minpa)

         !!!!KDW 02/21/22 - End Diffusion !!!!

         w12 = Theta(ch_nut(jnut)%w12,thw12,wtmp) * ch_hyd(jhyd)%l * ch_hyd(jhyd)%w ! [(g_P/day)/(m^2 * mg_P/kg)] - KDW
         !dox_factor = ch(jrch)%rch_dox / (ch(jrch)%rch_dox + ch_nut(jnut)%dox_half)
         dox_factor = 1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))) ! KDW 10/05/21, replace above to be more sensistive near DO threshold
         cbod_factor = ch(jrch)%rch_cbod / (ch(jrch)%rch_cbod + ch_nut(jnut)%cbod_half)
         w12 = w12 * dox_factor * cbod_factor
         !if (w12 > 1) w12 = 1 ! commented out KDW, 06/01
         btrb_orgp = w12 * (ch(jrch)%aer_orgp - ch(jrch)%anaer_orgp) * rt_delt
         equil_conc = ((ch(jrch)%aer_orgp * aer_mass) + (ch(jrch)%anaer_orgp * anaer_mass)) / & ! KDW 06/02
            (aer_mass + anaer_mass)
         !for section below, had "aer_mass / 1000" and "anaer_mass / 1000" when /1000 doesn't belong, fixed 02/24/22, KDW
         if (btrb_orgp > (ch(jrch)%aer_orgp - equil_conc) * aer_mass .and. btrb_orgp > 0) then !added ".and. btrb_orb > 0" here and in lines below, KDW 11/01/21
            btrb_orgp = (ch(jrch)%aer_orgp - equil_conc) * aer_mass 
         endif
         if (btrb_orgp < -(ch(jrch)%anaer_orgp - equil_conc) * anaer_mass .and. btrb_orgp < 0) then
            btrb_orgp = -(ch(jrch)%anaer_orgp - equil_conc) * anaer_mass 
         endif
         btrb_minpa = w12 * (ch(jrch)%aer_minpa - ch(jrch)%anaer_minpa) * rt_delt
         equil_conc = ((ch(jrch)%aer_minpa * aer_mass) + (ch(jrch)%anaer_minpa * anaer_mass)) / & ! KDW 06/02
            (aer_mass + anaer_mass)
         if (btrb_minpa > (ch(jrch)%aer_minpa - equil_conc) * aer_mass .and. btrb_minpa > 0) then
            btrb_minpa = (ch(jrch)%aer_minpa - equil_conc) * aer_mass 
         endif
         if (btrb_minpa < -(ch(jrch)%anaer_minpa - equil_conc) * anaer_mass .and. btrb_minpa < 0) then
            btrb_minpa = -(ch(jrch)%anaer_minpa - equil_conc) * anaer_mass 
         endif
         btrb_minps = w12 * (ch(jrch)%aer_minps - ch(jrch)%anaer_minps) * rt_delt
         equil_conc = ((ch(jrch)%aer_minps * aer_mass) + (ch(jrch)%anaer_minps * anaer_mass)) / & ! KDW 06/02
            (aer_mass + anaer_mass)
         if (btrb_minps > (ch(jrch)%aer_minps - equil_conc) * aer_mass .and. btrb_minps > 0) then
            btrb_minps = (ch(jrch)%aer_minps - equil_conc) * aer_mass 
         endif
         if (btrb_minps < -(ch(jrch)%anaer_minps - equil_conc) * anaer_mass .and. btrb_minps < 0) then
            btrb_minps = -(ch(jrch)%anaer_minps - equil_conc) * anaer_mass
         endif
         ! ch(jrch)%aer_disp = ch(jrch)%aer_disp - (diff_aer2anaer/(aer_vol/1000)) - (diff_aer2wtclm/(aer_vol/1000))
         ! ch(jrch)%anaer_disp = ch(jrch)%anaer_disp + (diff_aer2anaer/(anaer_vol/1000)) - (diff_anaer2wtclm/(anaer_vol/1000))
         !Commented out above with below because diffusion mass transfer accounted for above now, KDW 02/25/22
         ch(jrch)%aer_orgp = ch(jrch)%aer_orgp - btrb_orgp/aer_mass
         ch(jrch)%anaer_orgp = ch(jrch)%anaer_orgp + btrb_orgp/anaer_mass
         ch(jrch)%aer_minpa = ch(jrch)%aer_minpa - btrb_minpa/aer_mass
         ch(jrch)%anaer_minpa = ch(jrch)%anaer_minpa + btrb_minpa/anaer_mass
         ch(jrch)%aer_minps = ch(jrch)%aer_minps - btrb_minps/aer_mass
         ch(jrch)%anaer_minps = ch(jrch)%anaer_minps + btrb_minps/anaer_mass

         !! Update aerobic and anaerobic layers for adsorption
         !!Assume adsorbed and dissolved reach equilibrium - KDW 03/03/22
            !Aerobic layer
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
            xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
            if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
            if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
            if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
            ch(jrch)%aer_minpa = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
            ch(jrch)%aer_disp = (P_tot - (ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
            if (ch(jrch)%aer_disp < .0) ch(jrch)%aer_disp = 0.! kDW 04/28/22 
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
            ! P_tot = ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%aer_disp * ch(jrch)%bed_por
            ! P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
            ! P_adsorbed = P_min_new - ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd
            ! P_adsorbed = MIN(P_adsorbed, ch(jrch)%aer_disp * ch(jrch)%bed_por)
            ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%aer_minpa * ch_sed(jsed)%bed_bd))
            ! ch(jrch)%aer_minpa =  ch(jrch)%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
            ! ch(jrch)%aer_disp = ch(jrch)%aer_disp - P_adsorbed / ch(jrch)%bed_por

            !Anaerobic layer
            dox_factor = 0.5
            P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
            xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/15/22
            yy = ch(jrch)%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
            if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
            if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
            if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
            ch(jrch)%anaer_minpa = P_min_new * xx / ch_sed(jsed)%bed_bd ! KDW 04/15/22, back to g/ kg
            ch(jrch)%anaer_disp = (P_tot - (ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd)) / ch(jrch)%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
            if (ch(jrch)%anaer_disp < .0) ch(jrch)%anaer_disp = 0.! kDW 04/28/22 
            ! yyy = ch(jrch)%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
            ! P_tot = ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd + ch(jrch)%anaer_disp * ch(jrch)%bed_por
            ! P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
            ! P_adsorbed = P_min_new - ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd
            ! P_adsorbed = MIN(P_adsorbed, ch(jrch)%anaer_disp * ch(jrch)%bed_por)
            ! P_adsorbed = MAX(P_adsorbed, (-ch(jrch)%anaer_minpa * ch_sed(jsed)%bed_bd))
            ! ch(jrch)%anaer_minpa =  ch(jrch)%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
            ! ch(jrch)%anaer_disp = ch(jrch)%anaer_disp - P_adsorbed / ch(jrch)%bed_por
            
         roc_aer = bk * (7. * ch(jrch)%aer_minpa - ch(jrch)%aer_minps) ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
         if (roc_aer > 0.) then
            roc_aer = Min(roc_aer, ch(jrch)%aer_minpa)
         endif
         if (roc_aer < 0.) then
            roc_aer = roc_aer * .1
            roc_aer = Max(roc_aer, -ch(jrch)%aer_minps)
         endif
         roc_anaer = bk * (7. * ch(jrch)%anaer_minpa - ch(jrch)%anaer_minps) ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
         if (roc_anaer > 0.) then
            roc_anaer = Min(roc_anaer, ch(jrch)%anaer_minpa)
         endif
         if (roc_anaer < 0.) then
            roc_anaer = roc_anaer * .1
            roc_anaer = Max(roc_anaer, -ch(jrch)%anaer_minps)
         endif

         if (ch(jrch)%aer_minpa - roc_aer < 0)  roc_aer = ch(jrch)%aer_minpa ! KDW 03/03/22
         xx = SA_VOL_aer * (1-ch(jrch)%bed_por) ! KDW 04/18/22
         if ((ch(jrch)%aer_minpa - roc_aer) > (ch_nut(jnut)%srpm2 * xx / ch_sed(jsed)%bed_bd)) then ! KDW 04/18/22, adsts replaced by (ch_nut(jnut)%srpm2 * xx / ch_sed(jsed)%bed_bd), here through next 8 lines
            roc_aer = ch(jrch)%aer_minpa - (ch_nut(jnut)%srpm2 * xx / ch_sed(jsed)%bed_bd) ! KDW 03/03/22! KDW 04/18/22
         endif
         if (ch(jrch)%anaer_minpa - roc_anaer < 0)  roc_anaer = ch(jrch)%anaer_minpa ! KDW 03/03/22
         xx = SA_VOL_anaer * (1-ch(jrch)%bed_por) ! KDW 04/18/22
         if ((ch(jrch)%anaer_minpa - roc_anaer) > (ch_nut(jnut)%srpm2 * xx / ch_sed(jsed)%bed_bd)) then! KDW 04/18/22
            roc_anaer = ch(jrch)%anaer_minpa - (ch_nut(jnut)%srpm2 * xx / ch_sed(jsed)%bed_bd) ! KDW 03/03/22! KDW 04/18/22
         endif

         !ch(jrch)%minps = ch(jrch)%minps + roc_wtclm
         !if (ch(jrch)%minps < 0.) ch(jrch)%minps = 0.
         ch(jrch)%aer_minps = ch(jrch)%aer_minps + roc_aer
         if (ch(jrch)%aer_minps < 0.) ch(jrch)%aer_minps = 0.
         ch(jrch)%anaer_minps = ch(jrch)%anaer_minps + roc_anaer
         if (ch(jrch)%anaer_minps < 0.) ch(jrch)%anaer_minps = 0.

         !ch(jrch)%minpa = ch(jrch)%minpa - roc_wtclm + rmp_wtclm
         !if (ch(jrch)%minpa < 0.) ch(jrch)%minpa = 0.
         ch(jrch)%aer_minpa = ch(jrch)%aer_minpa - roc_aer  ! KDW 03/03/22
         if (ch(jrch)%aer_minpa < 0.) ch(jrch)%aer_minpa = 0.
         ch(jrch)%anaer_minpa = ch(jrch)%anaer_minpa - roc_anaer  ! KDW 03/03/22
         if (ch(jrch)%anaer_minpa < 0.) ch(jrch)%anaer_minpa = 0.

         !ch(jrch)%disolvp = ch(jrch)%disolvp - rmp_wtclm
         !if (ch(jrch)%disolvp < 0.) ch(jrch)%disolvp = 0.    
         
         !! Update water column pools, adding new m factors from above - KDW

         !! calculate organic phosphorus concentration at end of day QUAL2E section 3.3.6 equation III-24
         bc4_m = wq_k2m(tday,rt_delt,-bc4_k,orgpcon,orgpin) 
         !rs5_k=0. ! rs5 no longer used, all settling included in sediment deposition - KDW
         !if (rchdep > 0.001) rs5_k = Theta(ch_nut(jnut)%rs5,thrs5,wtmp) / rchdep 
         !factk=-rs5_k
         factk = 0.
         alg_m_orgp = alg_md * ch_nut(jnut)%ai2 ! KDW
         factm =bc4_m + alg_m_orgp ! added alg_m_orgp - KDW
         ch(jrch)%organicp = wq_semianalyt (tday, rt_delt, factm, factk, orgpcon, orgpin)
         if (ch(jrch)%organicp < 1.e-6) ch(jrch)%organicp = 0. 

         !Water column adsorption-desorption - KDW 03/03/22
         wtclm_mixed_minpa = (minpacon * vol_0 + minpain * wtrin) ! KDW 02/24/22
         wtclm_mixed_minpa_conc = wtclm_mixed_minpa / qdin! KDW 02/24/22
         wtclm_mixed_minps = (minpscon * vol_0 + minpsin * wtrin) ! KDW 02/24/22
         wtclm_mixed_minps_conc = wtclm_mixed_minps / qdin! KDW 02/24/22
         wtclm_mixed_solp = (solpcon * vol_0 + dispin * wtrin)! KDW 02/24/22
         wtclm_mixed_solp_conc = wtclm_mixed_solp / qdin! KDW 02/24/22
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
         P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc
         xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
         yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
         if(xx > 0.)then ! KDW 04/18/22
            P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
               4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
         else
            P_min_new = 0.
         endif
         if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
         if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22 
         if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
         P_adsorbed = P_min_new * xx - wtclm_mixed_minpa_conc! KDW 04/15/22         
         ! wtclm_por = 1 - (sedin / qdin) / 1.2 
         ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 02/26/22
         ! P_min_new = (adsts_wtclm + P_tot + yyy - ((adsts_wtclm + P_tot + yyy)**2 - 4*P_tot*adsts_wtclm) **0.5) / 2
         ! P_adsorbed = P_min_new - wtclm_mixed_minpa_conc    
         P_adsorbed = MIN(P_adsorbed, wtclm_mixed_solp_conc)
         P_adsorbed = MAX(P_adsorbed, -wtclm_mixed_minpa_conc)     

         !! calculate dissolved phosphorus concentration at end of day QUAL2E section 3.4.2 equation III-25
         factk = 0.
         alg_m_disp = alg_mg * ch_nut(jnut)%ai2 ! KDW
         factm = -bc4_m - alg_m_disp - P_adsorbed ! changed from ads_m to P_adsorbed, KDW 03/03/22
         ch(jrch)%disolvp = wq_semianalyt (tday, rt_delt, factm, 0., solpcon, dispin)
         if (ch(jrch)%disolvp < 1.e-6) ch(jrch)%disolvp = 0.         
         
         !! calculate active mineral phosphorus concentration at end of day
         if (vol_0 + wtrin > 0) then  
            ssd_m1_max = wtclm_mixed_minpa_conc ! KDW 02/24/22
            ssd_m2_max = wtclm_mixed_minps_conc ! KDW 02/24/22
         else
            ssd_m1_max = 0
            ssd_m2_max = 0            
         endif
         factk = 7. * bk ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
         ssd_m1 = -wq_k2m (tday, rt_delt, -factk, minpacon, minpain) !Made RHS and factk negative on 05/24, KDW
         ssd_m1 = Min(ssd_m1, ssd_m1_max)
         factk = bk
         ssd_m2 = -wq_k2m (tday, rt_delt, -factk, minpscon, minpsin) !Made RHS and factk negative on 05/24, KDW
         ssd_m2 = Min(ssd_m2, ssd_m2_max)
         if (ssd_m1 > ssd_m2) then
            ssd_m = ssd_m1 - ssd_m2
            !ssd_m = Min(ssd_m, minpacon) ! is this needed??
         else 
            factk = 7. * bk * 0.1 ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
            ssd_m1 = -wq_k2m (tday, rt_delt, -factk, minpacon, minpain) !Made RHS and factk negative on 05/24, KDW
            ssd_m1 = Min(ssd_m1, ssd_m1_max)
            factk = bk  * 0.1
            ssd_m2 = -wq_k2m (tday, rt_delt, -factk, minpscon, minpsin) !Made RHS and factk negative on 05/24, KDW
            ssd_m2 = Min(ssd_m1, ssd_m2_max)
            ssd_m = ssd_m1 - ssd_m2
            !ssd_m = Max(ssd_m, -minpscon) ! is this needed??
         endif

         factk = 0.
         factm = P_adsorbed - ssd_m ! changed from ads_m to P_adsorbed, KDW 03/03/22
         ch(jrch)%minpa = wq_semianalyt (tday, rt_delt, factm, factk, minpacon, minpain)
         if (ch(jrch)%minpa < 1.e-6) ch(jrch)%minpa = 0. 
         !!! Check to see if minpa exceeds adsts_wtclm !!!
         xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/18/22
         if (ch(jrch)%minpa > (ch_nut(jnut)%srpm2 * xx)) then ! KDW 11/30/21, had typo, wtclmmm, fixed 04/12/22 ! KDW 04/18/22 replaced adsts with ! KDW 04/18/22 here and next 3 lines
            xyz = ch(jrch)%minpa - (ch_nut(jnut)%srpm2 * xx)
            ch(jrch)%minpa = (ch_nut(jnut)%srpm2 * xx)
            ch(jrch)%disolvp = ch(jrch)%disolvp + xyz
         endif

         !! calculate stable mineral phosphorus concentration at end of day
         factk = 0.
         factm = ssd_m
         ch(jrch)%minps = wq_semianalyt (tday, rt_delt, factm, factk, minpscon, minpsin)
         if (ch(jrch)%minps < 1.e-6) ch(jrch)%minps = 0. 

         !! end phosphorus calculations

        else
        !! all water quality variables set to zero when no flow
        !! Note: transformations within the aerobic and anaerobic layers are being foregone as well - KDW
         algin = 0.0
         chlin = 0.0
         orgnin = 0.0
         ammoin = 0.0
         nitritin = 0.0
         nitratin = 0.0
         orgpin = 0.0
         minpain = 0.0 ! - KDW
         minpsin = 0.0 ! - KDW
         dispin = 0.0
         cbodin = 0.0
         disoxin = 0.0
         ch(jrch)%algae = 0.0
         ch(jrch)%chlora = 0.0
         ch(jrch)%organicn = 0.0
         ch(jrch)%ammonian = 0.0
         ch(jrch)%nitriten = 0.0
         ch(jrch)%nitraten = 0.0
         ch(jrch)%organicp = 0.0
         ch(jrch)%minpa = 0.0 ! - KDW
         ch(jrch)%minps = 0.0 ! - KDW
         ch(jrch)%disolvp = 0.0
         ch(jrch)%rch_cbod = 0.0
         ch(jrch)%rch_dox = 0.0
         soxy = 0.0
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

        algcon2 = 1000. * ch(jrch)%chlora / ch_nut(jnut)%ai0 ! changed to match watqual4
        orgncon2 = ch(jrch)%organicn
        nh3con2 = ch(jrch)%ammonian 
        no2con2 =  ch(jrch)%nitriten
        no3con2 = ch(jrch)%nitraten 
        orgpcon2 =ch(jrch)%organicp 
        minpacon2 =ch(jrch)%minpa ! - KDW
        minpscon2 =ch(jrch)%minps ! - KDW
        solpcon2 = ch(jrch)%disolvp 
        cbodcon2 = ch(jrch)%rch_cbod 
        o2con2 = ch(jrch)%rch_dox 
        aer_disp2 = ch(jrch)%aer_disp
        anaer_disp2 = ch(jrch)%anaer_disp
        aer_orgp2 = ch(jrch)%aer_orgp
        anaer_orgp2 = ch(jrch)%anaer_orgp
        aer_minpa2 = ch(jrch)%aer_minpa
        anaer_minpa2 = ch(jrch)%anaer_minpa
        aer_minps2 = ch(jrch)%aer_minps
        anaer_minps2 = ch(jrch)%anaer_minps

        ch_d(jrch)%aer_disp = ch(jrch)%aer_disp ! KDW, 06/04
        ch_d(jrch)%anaer_disp = ch(jrch)%anaer_disp ! KDW, 06/04
        ch_d(jrch)%aer_orgp = ch(jrch)%aer_orgp ! KDW, 06/04
        ch_d(jrch)%anaer_orgp = ch(jrch)%anaer_orgp ! KDW, 06/04
        ch_d(jrch)%aer_minpa = ch(jrch)%aer_minpa ! KDW, 06/04
        ch_d(jrch)%anaer_minpa = ch(jrch)%anaer_minpa ! KDW, 06/04
        ch_d(jrch)%aer_minps = ch(jrch)%aer_minps ! KDW, 06/04
        ch_d(jrch)%anaer_minps = ch(jrch)%anaer_minps ! KDW, 06/04
        ch_d(jrch)%phos_diff = (diff_aer2wtclm + diff_aer2wtclm2 + diff_anaer2wtclm)/ 1000 ! KDW, 06/04, 04/18/22 added diff_aer2wtclm2
      return
      end subroutine ch_watqual5
      

