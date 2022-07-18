      subroutine res_nutrient2 (jres, inut, iob, mxsa, iswet)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine performs reservoir nutrient transformations and water
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
    use hydrograph_module      
    use climate_module
    use reservoir_data_module
    use time_module
    use reservoir_module
    use channel_data_module
    use channel_module
    use basin_module
    use water_body_module

    integer, intent (in) :: iob
    integer, intent (in) :: jres
    integer, intent (in) :: inut
    real :: tday, wtmp, fll, gra
    real :: lambda, fnn, fpp, algi, fl_1, xx, yy, zz, ww, qq, xyz, zyx, yyy
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
    real :: alg_m, alg_md, alg_mg, alg_set, con_out
    real :: alg_m_o2, alg_m_disp, alg_m_orgp, alg_nh4_m, alg_no3_m
    real :: fines, aer_mass, anaer_mass, aer_mass_remain, anaer_mass_remain
    real :: sedin, sanin, silin, clain, sagin
    real :: topcla, topsil, topsag, topsan
    real :: cla_aer, san_aer, sil_aer, sag_aer
    real :: cla_anaer, san_anaer, sil_anaer, sag_anaer
    real :: andep_cla, andep_sil, andep_sag, andep_san, andep_tot
    real :: w12, diff_aer2anaer, diff_aer2wtclm
    real :: dox_factor, cbod_factor, btrb_orgp, btrb_minpa, btrb_minps
    real :: mnrlz_aer, mnrlz_anaer
    real :: conv_mass2vol, conv_conc, wtclm_por
    !real :: adsts_wtclm, adsts_aer, adsts_anaer
    !real :: frsts_wtclm, frsts_aer, frsts_anaer
    real :: pai_wtclm, pai_aer, pai_anaer
    real :: rto_wtclm, rto_aer, rto_anaer
    real :: sfarea_anaer, sfarea_aer, sfarea_wtclm
    real :: aer_vol, anaer_vol
    real :: rmp_aer, rmp_anaer, roc_aer, roc_anaer
    real :: ads_m, ads_m1, ads_m2
    real :: ssd_m, ssd_m1, ssd_m2
    real :: qdin
    real :: vol_ave, resdep, vol_init 
    real :: tstep
    real, intent (in) :: mxsa ! - KDW
    integer, intent (in) :: iswet ! KDW 03/30
    real :: equil_conc
    real :: disox_ave ! KDW 02/11/22
    real :: diff, max_diff
    real :: wtclm_mixed_minpa, wtclm_mixed_minpa_conc, wtclm_mixed_solp, wtclm_mixed_solp_conc ! KDW 02/23/22
    real :: diffremain! KDW 02/23/2
    real :: P_adsorbed, P_min_new, P_tot, P_adsorbed_mass, P_desorbed_mass, ads_inflow! KDW 03/03/22
    !real :: max_diff2, B_aer, B_anaer, B_wtclm ! KDW 03/05/22
    real :: max_diff_minpa, max_diff_disp, max_disp_conc, max_Ptot_conc ! KDW 03/07/22, 03/21/22
    real :: min_disp_conc, min_Ptot_conc, max_diff1, max_diff2, max_diff3 ! KDw 03/22/22
    real :: rs2_ratio, diff_aer2wtclm2, Ptot_conc_balanced,aer_balanced ! KDw 03/29/22, 04/01/22
    real :: aer_disp_temp, aer_minpa_temp, anaer_disp_temp, anaer_minpa_temp ! KDw 04/01/22
    real :: wtclm_mixed_solp_conc_temp, wtclm_mixed_minpa_conc_temp ! KDw 04/01/22

    !! initialize new variables - KDW
    rs3_s = 0. ; rk4_s = 0. ; rk1_k = 0. ; rk1_m = 0. ; rk3_k = 0. ; rk2_m = 0. ; rk2_k = 0. 
    bc1_k = 0. ; bc2_m  = 0. ; bc3_k = 0. ; bc3_m = 0. ; rs2_s = 0. ; bc4_k = 0. ; bc4_m = 0.
    factk = 0. ; factm = 0.
    alg_m = 0. ; alg_md = 0. ; alg_mg = 0. ; alg_set = 0. ; con_out = 0.
    alg_m_o2 = 0. ; alg_m_disp = 0. ; alg_m_orgp = 0. ; alg_nh4_m = 0. ; alg_no3_m = 0.
    fines = 0. ; aer_mass = 0. ; anaer_mass = 0. ; aer_mass_remain = 0. ; anaer_mass_remain = 0.
    sedin = 0. ; sanin = 0. ; silin = 0. ; clain = 0. ; sagin = 0.
    topcla = 0. ; topsil = 0. ; topsag = 0. ; topsan = 0.
    cla_aer = 0. ; san_aer = 0. ; sil_aer = 0. ; sag_aer = 0.
    cla_anaer = 0. ; san_anaer = 0. ; sil_anaer = 0. ; sag_anaer = 0.
    andep_cla = 0. ; andep_sil = 0. ; andep_sag = 0. ; andep_san = 0. ; andep_tot = 0.
    w12 = 0. ; rs2_s = 0. ; diff_aer2anaer = 0. ; diff_aer2wtclm = 0.
    dox_factor = 0. ; cbod_factor = 0. ; btrb_orgp = 0. ; btrb_minpa = 0. ; btrb_minps = 0.
    mnrlz_aer = 0. ; mnrlz_anaer = 0.
    conv_mass2vol = 0. ; conv_conc = 0. ; wtclm_por = 0.
    !adsts_wtclm = 0. ; adsts_aer = 0. ; adsts_anaer = 0.
    !frsts_wtclm = 0. ; frsts_aer = 0. ; frsts_anaer = 0.
    pai_wtclm = 0. ; pai_aer = 0. ; pai_anaer = 0.
    rto_wtclm = 0. ; rto_aer = 0. ; rto_anaer = 0.
    sfarea_anaer = 0. ; sfarea_aer = 0. ; sfarea_wtclm = 0.
    rmp_aer = 0. ; rmp_anaer = 0. ; roc_aer = 0. ; roc_anaer = 0.
    ads_m = 0. ; ads_m1 = 0. ; ads_m2 = 0.
    ssd_m = 0. ; ssd_m1 = 0. ; ssd_m2 = 0.
    aer_vol = 0. ; anaer_vol = 0. ; tstep = 0.
    vol_ave = 0.; vol_init = 0.; resdep = 0.
    equil_conc = 0. 
    xyz = 0. ; zyx = 0. ; diff = 0. ; max_diff = 0. ! KDW 02/23/22
    disp_aer_init = 0. ! KDW 02/23/22
    wtclm_mixed_minpa = 0. ; wtclm_mixed_minpa_conc = 0. ! KDW 02/23/22
    wtclm_mixed_solp = 0. ; wtclm_mixed_solp_conc = 0. ! KDW 02/23/22
    diffremain = 0. ! KDW 02/23/22
    yyy = 0.! KDW 02/28/22
    P_adsorbed = 0.; P_min_new = 0.; P_tot = 0.; P_adsorbed_mass = 0.; P_desorbed_mass = 0.; ads_inflow = 0.  ! KDW 03/03/22
    max_diff2 = 0.; B_aer = 0.; B_anaer = 0.; B_wtclm = 0. ! KDW 03/05/22
    max_diff_minpa = 0.; max_diff_disp = 0.; max_disp_conc = 0.; max_Ptot_conc = 0. ! KDW 03/07/22, 03/21/22
    min_disp_conc = 0.; min_Ptot_conc = 0.; max_diff1 = 0.; max_diff2 = 0.; max_diff3 = 0. ! KDw 03/22/22
    rs2_ratio = 0.1 ! KDW ratio of rs2_s b/w wtclm and aerobic layer vs. b/w aerobic and anaerobic layer, 03/29/22
    diff_aer2wtclm2 = 0.; Ptot_conc_balanced = 0.; aer_balanced = 0.; ! KDw 03/29/22, 04/01/22
    aer_disp_temp = 0.; aer_minpa_temp = 0.; anaer_disp_temp = 0.; anaer_minpa_temp = 0. ! KDw 04/01/22
    wtclm_mixed_solp_conc_temp = 0.; wtclm_mixed_minpa_conc_temp = 0. ! KDw 04/01/22
    SA_VOL_wtclm = 0. ; SA_VOL_anaer = 0. ; SA_VOL_aer = 0. ! KDW 04/15/22

    iwst = ob(iob)%wst
    !! calculate flow duration
    if (iswet == 0) then ! KDW 03/30
      vol_ave = (2 * wbody%flo - ht1%flo - ht2%flo) / 2 !flo includes inflow but outflow has not been subtracted, see res_control - KDW
      vol_init = wbody%flo - ht1%flo
    else
      vol_ave = ((wbody%flo-ht2%flo) + &
         (wbody%flo - ht1%flo - wbody_wb%precip + wbody_wb%evap + &
         wbody_wb%seep)) / 2 ! KDW 03/30
      !vol_init = wbody%flo - ht1%flo - wbody_wb%precip + &
      !   wbody_wb%evap + wbody_wb%seep ! KDW 03/30
      vol_init = wbody%flo - ht1%flo !Note - this isn't, for wetlands, the volume at beginning of day...
      !since precip, evap and seepage accounted for already. Need to include those here so no missing water in routine below - KDW 02/25/22
    endif
    vol_ave = max(0., vol_ave)
    vol_init = max(0., vol_init) ! KDW 03/30
    ! if (wbody_wb%area_ha > 0.01) then
      resdep = vol_ave / (wbody_wb%area_ha * 10000) !for algae calculations - KDW
    ! else
    !   resdep = vol_ave / (0.01 * 10000)
    ! endif
    if (ht1%flo > 0) then ! KDW 08/21
      tday = vol_init / ht1%flo
    else
      tday = 10000.
    endif
    tday = min(10000., tday) ! sets arbitratily large maximum residence time of water
    tday = max(.0001, tday) ! sets arbitrarily small minimum residence time of water
    !inflow_adj = vol_init/tday ! for extremely low inflows, an adjusted inflow for benthic source calculations 
    !vol_adj = tday * ht1%flo ! for extremely low initial volumes, an adjusted vol_init 

    if (time%step == 0) then
      tstep = 1.
    else
      tstep = 1. / time%step  
    endif

    !! calculate temperature in stream Stefan and Preudhomme. 1993.  Stream temperature estimation 
    !! from air temperature.  Water Res. Bull. p. 27-45 SWAT manual equation 2.3.13
    wtmp = 5.0 + 0.75 * wst(iwst)%weat%tave
    if (wtmp <= 0.) wtmp = 0.1

    !! benthic sources/losses in mg   
    
    rs3_s =  Theta(ch_nut(inut)%rs3,thrs3,wtmp) * 10 * wbody_wb%area_ha * tstep
    rk4_s =  Theta(ch_nut(inut)%rk4,thrk4,wtmp) * 10 * wbody_wb%area_ha * tstep

    !! initialize concentration of nutrient in reservoir
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
    wbody%cbod = amax1(1.e-6, wbody%cbod)

    !wtrin = ht1%flo
    qdin = wbody%flo
    if (ht1%flo > 0.0001) then ! change from > 0.0001 , very small 
     chlin = 1000. * ht1%chla * tstep / ht1%flo 
     algin = 1000. * chlin / ch_nut(inut)%ai0        !! QUAL2E equation III-1
     orgnin = 1000. * ht1%orgn  * tstep/ ht1%flo  
     ammoin = 1000. * ht1%nh3 * tstep / ht1%flo  
     nitritin = 1000. * ht1%no2 * tstep / ht1%flo
     nitratin = 1000. * ht1%no3 * tstep  / ht1%flo
     orgpin = 1000. * ht1%orgp * tstep / ht1%flo
     minpain = 1000. * ht1%minpa *tstep / ht1%flo
     minpsin = 1000. * ht1%minps *tstep / ht1%flo
     dispin = 1000. * ht1%solp* tstep  / ht1%flo
     cbodin = 1000. * ht1%cbod * tstep / ht1%flo
     disoxin = 1000. * ht1%dox  * tstep / ht1%flo
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
     orgpin = 0
     minpain = 0
     minpsin = 0
     dispin = 0
     cbodin = 0
     disoxin = 0
     cinn = 0
    end if    
  
   if (wbody%flo > 0.0001 .and. resdep > 0.0001 .and. wbody_wb%area_ha > 0.01) then ! KDW 03/30, added resdep condition, 10/15/21 added area_ha condition
    if (vol_init > 0.0001) then
    !algcon =ch(jrch)%algae 
        algcon = 1000. * (1000. * wbody%chla / vol_init) / ch_nut(inut)%ai0 ! changed to match watqual4
        orgncon = 1000. * wbody%orgn / vol_init
        nh3con = 1000. * wbody%nh3 / vol_init 
        no2con =  1000. * wbody%no2 / vol_init 
        no3con = 1000. * wbody%no3 / vol_init  
        orgpcon = 1000. * wbody%orgp / vol_init 
        minpacon = 1000. * wbody%minpa / vol_init 
        minpscon = 1000. * wbody%minps / vol_init 
        solpcon = 1000. * wbody%solp / vol_init 
        cbodcon= 1000. * wbody%cbod / vol_init 
        o2con = 1000. * wbody%dox / vol_init 
    else
        algcon = 0
        orgncon = 0
        nh3con = 0 
        no2con =  0
        no3con = 0
        orgpcon = 0
        minpacon = 0
        minpscon = 0
        solpcon = 0
        cbodcon= 0
        o2con = 0
    endif
    if (orgncon < 1.e-6) orgncon = 0.0
    if (nh3con < 1.e-6) nh3con = 0.0
    if (no2con < 1.e-6) no2con = 0.0
    if (no3con < 1.e-6) no3con = 0.0
    if (orgpcon < 1.e-6) orgpcon = 0.0
    if (minpacon < 1.e-6) minpacon = 0.0 
    if (minpscon < 1.e-6) minpscon = 0.0 
    if (solpcon < 1.e-6) solpcon = 0.0
    if (cbodcon < 1.e-6) cbodcon = 0.0
    if (o2con < 1.e-6) o2con = 0.0

    !!!! add to closing else statement, so when no inflow, benthic source injected straight to initial conc. - KDWnew!!!!!!!! 
    if (tday < 10000) then ! changed from ht1%flo > 0.0001 - KDW
     !! Update dispin here for diffusion between water column and aerobic layer
     disoxin = disoxin - rk4_s /ht1%flo 
     if (disoxin < 0) then   ! subtract from o2con - KDW ! KDW 02/11/22, copied from watqual5
        if (vol_init  > 0) then
            o2con = ((o2con * vol_init ) + (disoxin * ht1%flo)) / vol_init 
        else
            o2con= 0
        endif
        disoxin = 0
     endif

     disoxin = max(0., disoxin)
     !dispin = dispin + rs2_s / wtrin ! new formulation in P section below - KDW
     ammoin = ammoin + rs3_s /ht1%flo
     !! calculate effective concentration of available nitrogen QUAL2E equation III-15
     !cinn = nh3con + no3con !commented out by KDW and replaced below, linne 439, 12/13/21
    else
      if (vol_init  > 0) then
        o2con = o2con - rk4_s /vol_init
        o2con = max(0., o2con)
          !dispin = dispin + rs2_s / wtrin ! new formulation in P section below - KDW
        nh3con = nh3con + rs3_s /vol_init
          !! calculate effective concentration of available nitrogen QUAL2E equation III-15
        !cinn = nh3con + no3con !commented out by KDW and replaced below, linne 439, 12/13/21
      else !if/else added 03/04/22, KDW
        o2con = 0.
        nh3con = 0.
      endif 
    endif

    cinn = ((nh3con + no3con) * vol_init + (ammoin +nitratin) * ht1%flo)/(ht1%flo + vol_init) ! KDW 12/13/21 in case no water in reach to start day, this formulation prevents instability in k2m
    cinp = (solpcon * vol_init + dispin * ht1%flo)/(ht1%flo + vol_init)

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
    bc1mod = ch_nut(inut)%bc1 * cordo !!! not in watqual4
    bc2mod = ch_nut(inut)%bc2 * cordo !!! not in watqual4
    !! end O2 impact calculations

    !! algal growth
    !! calculate light extinction coefficient 
    !! (algal self shading) QUAL2E equation III-12
    if (ch_nut(inut)%ai0 * algcon > 1.e-6) then
      lambda = ch_nut(inut)%lambda0 + (ch_nut(inut)%lambda1 *      &
          ch_nut(inut)%ai0 * algcon)                                &
          + ch_nut(inut)%lambda2 * (ch_nut(inut)%ai0 *              & 
          algcon) ** (.66667)
    else
      lambda = ch_nut(inut)%lambda0
    endif

    if (lambda > ch_nut(inut)%lambda0) lambda = ch_nut(inut)%lambda0
    !! calculate algal growth limitation factors for nitrogen
    !! and phosphorus QUAL2E equations III-13 & III-14
    fnn = cinn / (cinn + ch_nut(inut)%k_n)
    !fpp = solpcon / (solpcon + ch_nut(inut)%k_p)
    fpp = cinp / (cinp + ch_nut(inut)%k_p)  ! replaced with below by KDW in case no water in reach at beginning of time step, 12/13/21

    !! calculate daylight average, photosynthetically active,
    !! light intensity QUAL2E equation III-8
    !! Light Averaging Option # 2
    iwgn = wst(iwst)%wco%wgn
    if (wst(iwst)%weat%daylength > 0.) then ! Corrected from wgn_pms(iwgn)%daylth (the dormancy threshold) to wst(iwst)%weat%daylength - KDW, 02/12/22
      algi = wst(iwst)%weat%solrad * ch_nut(inut)%tfact /  wst(iwst)%weat%daylength  ! Corrected from wgn_pms(iwgn)%daylth (the dormancy threshold) to wst(iwst)%weat%daylength - KDW, 02/12/22
    else
      algi = 0.00001
    end if

    !! calculate growth attenuation factor for light, based on
    !! daylight average light intensity QUAL2E equation III-7b
    fl_1 = (1. / (lambda * resdep)) *                               &                             
         Log((ch_nut(inut)%k_l + algi) / (ch_nut(inut)%k_l + algi *  &
         (Exp(-lambda * resdep))))
         fll = 0.92 * (wst(iwst)%weat%daylength / 24.) * fl_1 ! Corrected from wgn_pms(iwgn)%daylth (the dormancy threshold) to wst(iwst)%weat%daylength - KDW, 05/24

    !! calculcate local algal growth rate
    if (algcon < 5000.) then
      select case (ch_nut(inut)%igropt)
       case (1)
         !! multiplicative QUAL2E equation III-3a
         gra = ch_nut(inut)%mumax * fll * fnn * fpp
       case (2)
         !! limiting nutrient QUAL2E equation III-3b
         gra = ch_nut(inut)%mumax * fll * Min(fnn, fpp)
       case (3)
         !! harmonic mean QUAL2E equation III-3c
         if (fnn > 1.e-6 .and. fpp > 1.e-6) then
           gra = ch_nut(inut)%mumax * fll * 2. / ((1. / fnn) + (1. / fpp))
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
    alg_mg = wq_k2m (tday, tstep, factk, algcon, algin) ! the growth rate ; changed from semianalyt to k2m - KDW
    if (alg_mg < 1.e-6) alg_mg = 0 ! KDW 12/13/21
    if ((alg_mg * ch_nut(jnut)%ai2) > (cinp*(ht1%flo + vol_init))) then !KDW 02/12/22, so growth doesn't require more disp than available
      alg_mg  = (cinp*(ht1%flo + vol_init)) / ch_nut(jnut)%ai2
    endif
    if ((alg_mg * ch_nut(jnut)%ai1) > (cinn*(ht1%flo + vol_init))) then!KDW 02/12/22, so growth doesn't require more N than available
      alg_mg  = (cinn*(ht1%flo + vol_init)) / ch_nut(jnut)%ai1
    endif
    if (alg_mg  > 500) then  !KDW 02/12/22, so algal growth isn't unstable, KDW debug 04/20/22, 2500 was 500
      alg_mg  = 500
    endif
    factk = Theta(gra,thgra,wtmp) - Theta(ch_nut(inut)%rhoq, thrho, wtmp)      
    alg_m = wq_k2m(tday, tstep, factk, algcon, algin) ! net growth rate; changed from semianalyt to k2m - KDW
    alg_md = alg_mg-alg_m ! the death rate ; edited by negating RHS, was the negative of the death rate - KDW
    if (alg_md < 1.e-6) alg_md = 0 ! KDW 12/13/21
    if (ABS(alg_m) < 1.e-6) alg_m = 0 ! KDW 12/13/21
    alg_set = 0.
    if (resdep > 0.001) alg_set = Theta (ch_nut(inut)%rs1, thrs1, wtmp) / resdep 
     
    con_out = wq_semianalyt (tday, tstep, alg_m, -alg_set, algcon, algin)           

    wbody_bed%algae = con_out
    if (wbody_bed%algae < 1.e-6) wbody_bed%algae = 0.
    if (wbody_bed%algae > 5000.) wbody_bed%algae = 5000.

    !! calculate chlorophyll-a concentration at end of day QUAL2E equation III-1
    wbody%chla = wbody_bed%algae * (ch_nut(inut)%ai0 / 1000.) * (wbody%flo / 1000.)
    con_out = 0.
    !! end algal growth 

    !! oxygen calculations
    !! calculate carbonaceous biological oxygen demand at end of day QUAL2E section 3.5 equation III-26
    !! adjust rk1 to m-term and BOD & O2 mass availability
     
    cbodocon = min (cbodcon,o2con)
    cbodoin = min (cbodin,disoxin)
    rk1_k = -Theta (ch_nut(inut)%rk1, thrk1,wtmp)
    rk1_m = wq_k2m (tday, tstep, rk1_k, cbodocon, cbodoin)
    ! calculate corresponding m-term
    rk3_k=0.
    if (resdep > 0.001)  rk3_k = -Theta (ch_nut(inut)%rk3, thrk3, wtmp) / resdep
    factm = rk1_m
    factk = rk3_k
    wbody%cbod = 0.
    con_out = wq_semianalyt (tday, tstep, factm, factk, cbodcon, cbodin)
    wbody%cbod = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%cbod < 1.e-6) wbody%cbod = 0.

    !! calculate dissolved oxygen concentration if reach at 
    !! end of day QUAL2E section 3.6 equation III-28

    rk2_m = Theta (ch_nut(inut)%rk2, thrk2, wtmp) * soxy ! this is the piece multiplied by saturated oxy level - KDW
    rk2_k = -Theta (ch_nut(inut)%rk2, thrk2, wtmp) 

    bc1_k = -Theta(ch_nut(inut)%bc1,thbc1,wtmp)
    bc2_k = -Theta(ch_nut(inut)%bc2,thbc2,wtmp)
    bc1_m = wq_k2m (tday, tstep, bc1_k, nh3con, ammoin) ! changed from rk2_k to bc1_k - KDW
    bc2_m = wq_k2m (tday, tstep, bc2_k, no2con, nitritin)

    alg_m_o2 = -ch_nut(inut)%ai4 * alg_md + ch_nut(inut)%ai3 * alg_mg
    disox_ave = (o2con * vol_init + disoxin * ht1%flo) / (ht1%flo + vol_init) ! KDW 02/11/22 so o2 consumed in respiration doesn't exceed available o2
    if(alg_m_o2 + disox_ave < 0) then
        alg_m_o2 = -disox_ave  !KDW 02/11/22 so o2 consumed in respiration doesn't exceed available o2 (including from photosynthesis in current step)
    endif
     
    factm = rk1_m + rk2_m + alg_m_o2 + bc1_m * ch_nut(inut)%ai5 + bc2_m * ch_nut(inut)%ai6 ! added alg_m_02, removed rk4_m b/c not defined and already accounted for in doxin - KDW
    con_out = wq_semianalyt (tday, tstep, factm, rk2_k, o2con, disoxin)

    wbody%dox = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (vol_init > 0.) o2con = 1000. * wbody%dox / vol_init ! if added on 03/04/22, KDW
    if (o2con.le.0.001) o2con=0.001 ! KDW 02/11/22, repeated from lines 465-66 so o2con doesn't grow unreasonably
    if (o2con.gt.30.) o2con=30. ! KDW 02/11/22, repeated from lines 465-66 so o2con doesn't grow unreasonably
    wbody%dox = vol_init * o2con / 1000. ! KDW 02/11/22, repeated from lines 465-66 so o2con doesn't grow unreasonably
    if (wbody%dox <0.) wbody%dox = 0.
    !! end oxygen calculations 

    !! nitrogen calculations
    !! calculate fraction of algal nitrogen uptake from ammonia pool QUAL2E equation III-18
    f1 = ch_nut(inut)%p_n * nh3con / (ch_nut(inut)%p_n * nh3con +     &
    (1. - ch_nut(inut)%p_n) * no3con + 1.e-6)

    !! calculate organic N concentration at end of day
    !! QUAL2E section 3.3.1 equation III-16
    bc3_k = Theta(ch_nut(inut)%bc3,thbc3,wtmp) 
    rs4_k=0.
    if (resdep > 0.001)  rs4_k = Theta (ch_nut(inut)%rs4, thrs4, wtmp) / resdep   
    bc3_m = wq_k2m (tday, tstep, -bc3_k, orgncon, orgnin)
    factk =-rs4_k
    factm = bc3_m ! this formulation omits organic N source from dying algea (only to nitrate/nitrite) - KDW
    con_out = wq_semianalyt (tday, tstep, factm, factk, orgncon, orgnin)
    wbody%orgn = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%orgn < 1.e-6) wbody%orgn = 0.

    !! calculate ammonia nitrogen concentration at end of day QUAL2E section 3.3.2 equation III-17
    factk= 0. ! changed from -bc1_k to 0 for clarity since not used as k factor below - KDW
    alg_nh4_m = alg_mg * f1 * ch_nut(inut)%ai1
    factm= bc1_m - bc3_m - alg_nh4_m
    con_out = wq_semianalyt(tday,tstep,factm,0.,nh3con,ammoin) 
    wbody%nh3 = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%nh3 < 1.e-6) wbody%nh3 = 0.

    !! calculate concentration of nitrite at end of day QUAL2E section 3.3.3 equation III-19
    factm=-bc1_m + bc2_m
    con_out = wq_semianalyt(tday,tstep,factm,0.,no2con,nitritin)
    wbody%no2 = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%no2 < 1.e-6) wbody%no2 = 0.

    !! calculate nitrate concentration at end of day QUAL2E section 3.3.4 equation III-20
    factk = 0.
    alg_no3_m = alg_md * (1. - f1) * ch_nut(inut)%ai1
    factm = -bc2_m - alg_no3_m
    
    con_out = wq_semianalyt (tday, tstep, factm,0., no3con, nitratin)
    wbody%no3 = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%no3 < 1.e-6) wbody%no3 = 0.

    !! end nitrogen calculations

    !! phosphorus calculations

    !! Determine surface area of pools
    sedin = ht1%sed  + wbody%sed
    sanin = ht1%san  + wbody%san
    silin = ht1%sil  + wbody%sil
    clain = ht1%cla  + wbody%cla
    sagin = ht1%sag  + wbody%sag

   !  sfarea_anaer = wbody_bed%anaer_srfarea
   !  sfarea_aer = wbody_bed%aer_srfarea
   !  sfarea_wtclm = clain + silin * .2 + sagin * (1./15.) + sanin * .01

   !  aer_mass = .001 * 1 * (mxsa * wbody_bed%bed_bd) ! depth of aer layer is 1mm
   !  anaer_mass = .1 * (mxsa * wbody_bed%bed_bd)
   !  aer_vol = .001 * mxsa * 1000 ! liters
   !  anaer_vol = .1 * mxsa * 1000 ! liters
    aer_mass = .002 * 1 * (mxsa * wbody_bed%bed_bd) ! depth of aer layer is 1mm! KDW debug 04/28/22, test thicker layers, replacing above
    anaer_mass = .2 * (mxsa * wbody_bed%bed_bd)! KDW debug 04/28/22, test thicker layers, replacing above
    aer_vol = .002 * mxsa * 1000 ! liters! KDW debug 04/28/22, test thicker layers, replacing above
    anaer_vol = .2 * mxsa * 1000 ! liters! KDW debug 04/28/22, test thicker layers, replacing above
    vol_h2o_aer = wbody_bed%bed_por * aer_vol ! KDW 03/21/22
    vol_h2o_anaer = wbody_bed%bed_por * anaer_vol ! KDW 03/21/22
    
    SA_VOL_anaer = wbody_bed%anaer_srfarea / anaer_mass ! KDW 04/14/22
    SA_VOL_aer = wbody_bed%aer_srfarea / aer_mass ! KDW 04/14/22
    SA_VOL_wtclm = 3. * ((clain/0.000002) + (silin/0.00001) + (sagin/0.00003) + (sanin/0.0002)) / sedin ! KDW 04/14/22


    !! Update aerobic and anaerobic layers for mineralization
    bc4_k = Theta(ch_nut(inut)%bc4,thbc4,wtmp)
    mnrlz_aer = bc4_k * wbody_bed%aer_orgp * tstep
    mnrlz_anaer = bc4_k * wbody_bed%anaer_orgp * tstep
    wbody_bed%aer_orgp = wbody_bed%aer_orgp - mnrlz_aer
    wbody_bed%anaer_orgp = wbody_bed%anaer_orgp - mnrlz_anaer
    wbody_bed%aer_disp = wbody_bed%aer_disp + mnrlz_aer * wbody_bed%bed_bd / wbody_bed%bed_por !  "* wbody_bed%bed_bd" added 02/28/22, KDW
    wbody_bed%anaer_disp = wbody_bed%anaer_disp + mnrlz_anaer * wbody_bed%bed_bd / wbody_bed%bed_por !  "* wbody_bed%bed_bd" added 02/28/22, KDW

    !! Update aerobic and anaerobic layers for diffusion between
       rs2_s =  Theta(ch_nut(inut)%rs2,thrs2,wtmp)* 10 * wbody_wb%area_ha *1000 ! repurposed rs2, units [(g_P/day)/(m^2 * mg_P/L)], rs2_s units [(g_p/day)/(mg_P/L)]
    ! KDW 2/14/22 added "* 1000" and changed unit of rs2 to [(g_P/day)/(m^2 * mg_P/L)] because issue with calibration max value
    !!!!KDW 02/21/22 !!!!
   ! dox_factor = 1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half)))
   ! adsts_wtclm = 1000000 * ch_nut(jnut)%spkgc * (sfarea_wtclm / (wbody%flo * 1000)) * &
   !    (0.75 + 0.25 * dox_factor)! mg P/L KDW 07/23
   ! adsts_aer = 1000000 * ch_nut(jnut)%spkgc * (sfarea_aer / aer_vol) * &
   !    (0.75 + 0.25 * dox_factor) ! mg P/L kg P/Mg sed  *  Mg clay equivilent / Liters = mg / L
   ! adsts_anaer = 1000000 * ch_nut(jnut)%spkgc * (sfarea_anaer / anaer_vol) * 0.75 
   !! Before diffusion, re-equilibrate dissolved v adsorbed because sediment erosion and deposition - KDW 03/03/22
   !Aerobic layer
   dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
   P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
   xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
   yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
   P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
      4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
   if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
   if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
   if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
   wbody_bed%aer_minpa = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
   wbody_bed%aer_disp = (P_tot - (wbody_bed%aer_minpa * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
   if (wbody_bed%aer_disp < .0) wbody_bed%aer_disp = 0.! kDW 04/28/22 
   ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
   ! P_tot = wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por
   ! P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
   ! P_adsorbed = P_min_new - wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd
   ! P_adsorbed = MIN(P_adsorbed, wbody_bed%aer_disp * wbody_bed%bed_por)
   ! P_adsorbed = MAX(P_adsorbed, (-wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd))
   ! wbody_bed%aer_minpa =  wbody_bed%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
   ! wbody_bed%aer_disp = wbody_bed%aer_disp - P_adsorbed / wbody_bed%bed_por
   !Anaerobic layer
   dox_factor = 0.5
   P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
   xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
   yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
   P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
      4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
   if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
   if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
   if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
   wbody_bed%anaer_minpa = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
   wbody_bed%anaer_disp = (P_tot - (wbody_bed%anaer_minpa * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
   if (wbody_bed%anaer_disp < .0) wbody_bed%anaer_disp = 0.! kDW 04/28/22 
   ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
   ! P_tot = wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por
   ! P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
   ! P_adsorbed = P_min_new - wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd
   ! P_adsorbed = MIN(P_adsorbed, wbody_bed%anaer_disp * wbody_bed%bed_por)
   ! P_adsorbed = MAX(P_adsorbed, (-wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd))
   ! wbody_bed%anaer_minpa =  wbody_bed%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
   ! wbody_bed%anaer_disp = wbody_bed%anaer_disp - P_adsorbed / wbody_bed%bed_por

   !! Before diffusion, re-equilibrate dissolved v adsorbed because inflow mixing
   wtclm_mixed_minpa = (minpacon * vol_init + minpain * ht1%flo) ! KDW 02/24/22
   wtclm_mixed_minpa_conc = wtclm_mixed_minpa / wbody%flo! KDW 02/24/22
   wtclm_mixed_minps = (minpscon * vol_init + minpsin * ht1%flo) ! KDW 02/24/22
   wtclm_mixed_minps_conc = wtclm_mixed_minps / wbody%flo! KDW 02/24/22
   wtclm_mixed_solp = (solpcon * vol_init + dispin * ht1%flo)! KDW 02/24/22
   wtclm_mixed_solp_conc = wtclm_mixed_solp / wbody%flo! KDW 02/24/22
   !wtclm_por = 1 - (sedin / wbody%flo) / 1.2 , KDW 04/18/22

   dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
   P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc ! wtclm_por approximately 1, so leave out here and below (in contrast to bed, for per TOTAL volume calc)
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
   P_adsorbed_mass = (P_min_new * xx - wtclm_mixed_minpa_conc)* wbody%flo ! KDW 04/15/22 
   wtclm_mixed_minpa_conc = P_min_new * xx ! KDW 04/15/22, back to g/ kg
   wtclm_mixed_solp_conc = P_tot - wtclm_mixed_minpa_conc ! KDW 04/15/22, back to g/m^3 pore water volume
   if (wtclm_mixed_solp_conc < .0) wtclm_mixed_solp_conc = 0.! kDW 04/28/22 
   wtclm_mixed_minpa = wtclm_mixed_minpa_conc * wbody%flo ! KDW 04/15/22
   wtclm_mixed_solp = wtclm_mixed_solp_conc * wbody%flo ! KDW 04/15/22 
   ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 02/26/22
   ! P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc
   ! P_min_new = (adsts_wtclm + P_tot + yyy - ((adsts_wtclm + P_tot + yyy)**2 - 4*P_tot*adsts_wtclm) **0.5) / 2
   ! P_adsorbed = P_min_new - wtclm_mixed_minpa_conc    
   ! P_adsorbed = MIN(P_adsorbed, wtclm_mixed_solp_conc)
   ! P_adsorbed = MAX(P_adsorbed, -wtclm_mixed_minpa_conc)
   ! P_adsorbed_mass = P_adsorbed * wbody%flo  
   ! wtclm_mixed_minpa_conc = P_min_new ! KDW 03/03/22
   ! wtclm_mixed_minpa = wtclm_mixed_minpa_conc * wbody%flo! KDW 03/03/22
   ! wtclm_mixed_solp_conc = wtclm_mixed_solp_conc - P_adsorbed ! KDW 03/03/22
   ! wtclm_mixed_solp = wtclm_mixed_solp_conc * wbody%flo ! KDW 03/03/22
   if (P_adsorbed_mass >= 0.) then
      if (P_adsorbed_mass < solpcon * vol_init) then !all adsorption/desorption from stored dissolved and mineral P
         solpcon = (solpcon * vol_init - P_adsorbed_mass) / (vol_init + 0.0001)
         minpacon = (minpacon * vol_init + P_adsorbed_mass) / (vol_init + 0.0001)
      else
         ads_inflow = P_adsorbed_mass - solpcon * vol_init
         solpcon = 0.
         minpacon = (minpacon * vol_init + (P_adsorbed_mass - ads_inflow)) / (vol_init + 0.0001)
         dispin = (dispin * ht1%flo - ads_inflow) / ht1%flo
         minpain = (minpain * ht1%flo + ads_inflow) / ht1%flo
      endif
   else
      P_desorbed_mass = -P_adsorbed_mass
      if (P_desorbed_mass < minpacon * vol_init) then
         minpacon = (minpacon * vol_init - P_desorbed_mass) / (vol_init + 0.0001)
         solpcon = (solpcon * vol_init + P_desorbed_mass) / (vol_init + 0.0001)
      else
         ads_inflow = P_desorbed_mass - minpacon * vol_init
         minpacon = 0.
         solpcon = (solpcon * vol_init + (P_desorbed_mass - ads_inflow)) / (vol_init + 0.0001)
         minpain = (minpain * ht1%flo - ads_inflow) / ht1%flo
         dispin = (dispin * ht1%flo + ads_inflow) / ht1%flo
      endif
   endif
  !!Exchangable P via diffusion!!
  ! Initial PAI used to estimate how much diffusion sourced from dissolved versus labile/adsorbed P
   ! frsts_wtclm = adsts_wtclm - wtclm_mixed_minpa_conc
   ! frsts_aer = adsts_aer - wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd ! mg P/L - aer_minpa is mgP/kg and is multiplied by Mg/m^3 (i.e. mg/L)
   ! frsts_anaer = adsts_anaer - wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd
   ! if (frsts_wtclm > 0) then
   !    pai_wtclm = 1 - (1 / (1 + (wtclm_por/((ch_nut(jnut)%adskeq / 100) * frsts_wtclm)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
   ! else
   !    frsts_wtclm = 0
   !    pai_wtclm = 1
   ! endif
   ! if (frsts_aer > 0) then
   !    pai_aer = 1 - (1 / (1 + (wbody_bed%bed_por/((ch_nut(jnut)%adskeq / 100) * frsts_aer)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
   ! else
   !    frsts_aer = 0
   !    pai_aer = 1
   ! endif
   ! if (frsts_anaer > 0) then
   !    pai_anaer = 1 - (1 / (1 + (wbody_bed%bed_por/((ch_nut(jnut)%adskeq / 100) * frsts_anaer)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
   ! else
   !    frsts_anaer = 0
   !    pai_anaer = 1
   ! endif
   if ((wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) > 0.) then ! KDW 04/18/22
      pai_wtclm = wtclm_mixed_solp_conc / (wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) !KDW 04/15/22
   else
      pai_wtclm = 1.
   endif
   if ((wbody_bed%aer_disp * wbody_bed%bed_por + wbody_bed%aer_minpa * wbody_bed%bed_bd) > 0.) then ! KDW 04/18/22
      pai_aer = wbody_bed%aer_disp * wbody_bed%bed_por / (wbody_bed%aer_disp * wbody_bed%bed_por + wbody_bed%aer_minpa * wbody_bed%bed_bd) !KDW 04/15/22
   else
      pai_aer = 0.01
   endif
   if((wbody_bed%anaer_disp * wbody_bed%bed_por + wbody_bed%anaer_minpa * wbody_bed%bed_bd) > 0.) then ! KDW 04/18/22
      pai_anaer = wbody_bed%anaer_disp * wbody_bed%bed_por / (wbody_bed%anaer_disp * wbody_bed%bed_por + wbody_bed%anaer_minpa * wbody_bed%bed_bd) !KDW 04/15/22
   else
      pai_anaer = 0.01
   endif

   !!Diffusion between the aerobic layer and the water column
   disp_aer_init = wbody_bed%aer_disp
   if (wbody%flo > 0) then !change from vol_init to wbody%flo, KDW 03/24/22
      diff_aer2wtclm = rs2_s * (wbody_bed%aer_disp - wtclm_mixed_solp_conc) * tstep ! unit grams P, by convention positive is diffusion from aerobic layer to water column
      !Check for overshooting equilibrium
      if (diff_aer2wtclm <= 0) then  
         diff = abs(diff_aer2wtclm)
         max_diff = (wtclm_mixed_solp_conc - wbody_bed%aer_disp) * & !labile/adsorbed P can compensate for strictly dissolved P
         (1/ ( (pai_aer/(vol_h2o_aer / 1000)) + (pai_wtclm/wbody%flo) ) )
         max_diff1 = (wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) * wbody%flo

         max_diff_disp = pai_wtclm * max_diff ! KDW 03/22/22
         P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc ! KDW 03/22/22
         min_disp_conc = wtclm_mixed_solp_conc - (max_diff_disp/wbody%flo) ! KDW 03/22/22
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
         xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, 04/18/22
         yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
         min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
         max_diff_minpa = (P_tot - min_Ptot_conc)*wbody%flo - max_diff_disp  ! KDW 04/15/22, grams
         max_diff_minpa = MIN(max_diff_minpa, ((1-pai_wtclm) * max_diff)) ! KDW 04/15/22
         max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
         ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
         ! min_Ptot_conc = min_disp_conc *(1+(adsts_wtclm/(yyy + min_disp_conc))) ! KDW 03/22/22
         ! max_diff_minpa = (P_tot - min_Ptot_conc)*wbody%flo - max_diff_disp  ! KDW 03/22/22
         ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_wtclm) * max_diff)) ! KDW 03/22/22
         ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

         max_diff_disp = pai_aer * max_diff ! KDW 03/07/22
         P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! KDW 03/21/22
         max_disp_conc = wbody_bed%aer_disp * wbody_bed%bed_por + (max_diff_disp/(aer_vol/1000)) ! KDW 03/21/22
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
         xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
         yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
         max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
         ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
         ! max_Ptot_conc = max_disp_conc *(1+(adsts_aer/(yyy + max_disp_conc))) ! KDW 03/21/22
         max_diff_minpa = (max_Ptot_conc - P_tot)*(aer_vol/1000) - max_diff_disp  ! KDW 03/21/22
         max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 03/21/22
         max_diff3 =  max_diff_disp + max_diff_minpa ! KDW 03/07/22

         max_diff = MIN(max_diff, max_diff1, max_diff2, max_diff3) ! KDW 03/22/22
         max_diff = MAX(max_diff, 0.)
         if(diff > max_diff) diff = max_diff
         if(diff < 0.) diff = 0. ! KDW 04/28/22
         diff_aer2wtclm = -diff
      endif
      if (diff_aer2wtclm >= 0) then  
         diff = abs(diff_aer2wtclm)
         max_diff = (wbody_bed%aer_disp - wtclm_mixed_solp_conc) * & !labile/adsorbed P can compensate for strictly dissolved P
         (1/ ( (pai_wtclm/wbody%flo) + (pai_aer/(vol_h2o_aer / 1000)) ) )
         max_diff1 = (wbody_bed%aer_disp * wbody_bed%bed_por + wbody_bed%aer_minpa * wbody_bed%bed_bd) * aer_vol/1000

         max_diff_disp = pai_aer * max_diff ! KDW 03/22/22
         P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! KDW 03/22/22
         min_disp_conc = wbody_bed%aer_disp * wbody_bed%bed_por - (max_diff_disp/(aer_vol/1000)) ! KDW 03/22/22
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
         xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
         yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
         min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
         max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
         max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 04/15/22
         max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
         ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
         ! min_Ptot_conc = min_disp_conc *(1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/22/22
         ! max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 03/22/22
         ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 03/22/22
         ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

         max_diff_disp = pai_wtclm * max_diff ! KDW 03/07/22
         P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc! KDW 03/21/22
         max_disp_conc = wtclm_mixed_solp_conc + (max_diff_disp/wbody%flo) ! KDW 03/21/22
         dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
         xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
         yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, KDW 04/18/22
         max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
         ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
         ! max_Ptot_conc = max_disp_conc *(1+(adsts_wtclm/(yyy + max_disp_conc))) ! KDW 03/21/22
         max_diff_minpa = (max_Ptot_conc - P_tot)*wbody%flo - max_diff_disp  ! KDW 03/21/22
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
      wbody_bed%aer_disp = wbody_bed%aer_disp + (diff*pai_aer/(vol_h2o_aer/1000))
      wbody_bed%aer_minpa = wbody_bed%aer_minpa + (diff*(1-pai_aer)/aer_mass)
      if((diff*pai_wtclm)<(dispin*ht1%flo))then
         dispin = dispin - (diff*pai_wtclm/ht1%flo)
         diffremain = diff - diff*pai_wtclm
      else
         diffremain = diff - (dispin*ht1%flo)
         dispin = 0
      endif
      if ((diff*(1-pai_wtclm))<(minpain*ht1%flo))then
         minpain = minpain - (diff*(1-pai_wtclm)/ ht1%flo)
         diffremain = diffremain - diff*(1-pai_wtclm)
      else
         diffremain = diffremain - (minpain*ht1%flo)
         minpain = 0
      endif
      diff = diffremain
      if(diff > 0.)then
         if((diff*pai_wtclm)<(solpcon*vol_init))then
            solpcon = solpcon - (diff*pai_wtclm/(vol_init + 0.0001))
            diffremain = diff - diff*pai_wtclm
         else
            diffremain = diff - (solpcon*vol_init)
            solpcon = 0
         endif
         if ((diff*(1-pai_wtclm))<(minpacon*vol_init))then
            minpacon = minpacon - (diff*(1-pai_wtclm)/ (vol_init + 0.0001))
            diffremain = diffremain - diff*(1-pai_wtclm)
         else
            diffremain = diffremain - (minpacon*vol_init)
            minpacon = 0
         endif
      endif
      if(diffremain > 0.) then !necessary because don't know fraction of P from inflow vs storage
         if(diffremain < (dispin*ht1%flo))then
            dispin = dispin - (diffremain/ht1%flo)
            diffremain = 0
         else
            diffremain = diffremain - (dispin*ht1%flo)
            dispin = 0
            if(diffremain < (solpcon*vol_init))then
               solpcon = solpcon - (diffremain/(vol_init + 0.0001))
               diffremain = 0
            else
               diffremain = diffremain - (solpcon*vol_init)
               solpcon = 0
               if(diffremain < (minpain*ht1%flo))then
                  minpain = minpain - (diffremain/ht1%flo)
                  diffremain = 0
               else
                  diffremain = diffremain - (minpain*ht1%flo)
                  minpain = 0
                  if(diffremain < (minpacon*vol_init))then
                     minpacon = minpacon - (diffremain/(vol_init + 0.0001))
                     diffremain = 0
                  else
                     diffremain = diffremain - (minpacon*vol_init)
                     minpacon = 0
                  endif
               endif
            endif
         endif
      endif 
   endif
   if (diff_aer2wtclm > 0) then
      diff = abs(diff_aer2wtclm)
      if ((diff*pai_wtclm/ht1%flo) < 10.) then ! in case ht1%flo is low, don't want excessive dispin causing instability
         dispin = dispin + (diff*pai_wtclm/ht1%flo)
         minpain = minpain + (diff*(1-pai_wtclm)/ht1%flo)
      else
         dispin = dispin + 10
         minpain = minpain + 10*(1-pai_wtclm)/pai_wtclm
         diffremain = diff - 10*(1/pai_wtclm) * ht1%flo
         solpcon = solpcon + (diffremain*pai_wtclm/(vol_init + 0.0001))
         minpacon = minpacon + (diffremain*(1-pai_wtclm)/(vol_init + 0.0001))
      endif
      wbody_bed%aer_disp = wbody_bed%aer_disp - (diff*pai_aer/(vol_h2o_aer/1000))
      wbody_bed%aer_minpa = wbody_bed%aer_minpa - (diff*(1-pai_aer)/aer_mass)
   endif  

   !! Bring water column and aerobic layer to adsorption equilibrium before anaerobic diffusion - KDW 03/22/22
   !Aerobic layer
   dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
   P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
   xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
   yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
   P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
      4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
   if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
   if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
   if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
   wbody_bed%aer_minpa = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
   wbody_bed%aer_disp = (P_tot - (wbody_bed%aer_minpa * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
   if (wbody_bed%aer_disp < .0) wbody_bed%aer_disp = 0.! kDW 04/28/22 
   ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
   ! P_tot = wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por
   ! P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
   ! P_adsorbed = P_min_new - wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd
   ! P_adsorbed = MIN(P_adsorbed, wbody_bed%aer_disp * wbody_bed%bed_por)
   ! P_adsorbed = MAX(P_adsorbed, (-wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd))
   ! wbody_bed%aer_minpa =  wbody_bed%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
   ! wbody_bed%aer_disp = wbody_bed%aer_disp - P_adsorbed / wbody_bed%bed_por
   !Water column
   wtclm_mixed_minpa = (minpacon * vol_init + minpain * ht1%flo) 
   wtclm_mixed_minpa_conc = wtclm_mixed_minpa / wbody%flo
   wtclm_mixed_solp = (solpcon * vol_init + dispin * ht1%flo)
   wtclm_mixed_solp_conc = wtclm_mixed_solp / wbody%flo
   ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) 
   P_tot = wtclm_mixed_minpa_conc + wtclm_mixed_solp_conc
   xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) ! KDW 04/15/22, KDW 04/18/22
   yy = (1 - (sedin / qdin) / 1.2) / ch_nut(jnut)%adskeq ! KDW 04/15/22, 04/18/22
   if(xx > 0.)then ! KDW 04/18/22
      P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
         4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
   else
      P_min_new = 0.
   endif
   if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
   if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22 
   if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
   P_adsorbed_mass = (P_min_new * xx - wtclm_mixed_minpa_conc)* wbody%flo ! KDW 04/15/22 
   wtclm_mixed_minpa_conc = P_min_new * xx ! KDW 04/15/22, back to g/ kg
   wtclm_mixed_solp_conc = P_tot - wtclm_mixed_minpa_conc ! KDW 04/15/22, back to g/m^3 pore water volume
   if (wtclm_mixed_solp_conc < .0) wtclm_mixed_solp_conc = 0.! kDW 04/28/22 
   wtclm_mixed_minpa = wtclm_mixed_minpa_conc * wbody%flo ! KDW 04/15/22
   wtclm_mixed_solp = wtclm_mixed_solp_conc * wbody%flo ! KDW 04/15/22 
   ! P_min_new = (adsts_wtclm + P_tot + yyy - ((adsts_wtclm + P_tot + yyy)**2 - 4*P_tot*adsts_wtclm) **0.5) / 2
   ! P_adsorbed = P_min_new - wtclm_mixed_minpa_conc    
   ! P_adsorbed = MIN(P_adsorbed, wtclm_mixed_solp_conc)
   ! P_adsorbed = MAX(P_adsorbed, -wtclm_mixed_minpa_conc)
   ! P_adsorbed_mass = P_adsorbed * wbody%flo  
   ! wtclm_mixed_minpa_conc = P_min_new
   ! wtclm_mixed_minpa = wtclm_mixed_minpa_conc * wbody%flo
   ! wtclm_mixed_solp_conc = wtclm_mixed_solp_conc - P_adsorbed
   ! wtclm_mixed_solp = wtclm_mixed_solp_conc * wbody%flo 
   if (P_adsorbed_mass >= 0.) then
      if (P_adsorbed_mass < solpcon * vol_init) then !all adsorption/desorption from stored dissolved and mineral P
         solpcon = (solpcon * vol_init - P_adsorbed_mass) / (vol_init + 0.0001)
         minpacon = (minpacon * vol_init + P_adsorbed_mass) / (vol_init + 0.0001)
      else
         ads_inflow = P_adsorbed_mass - solpcon * vol_init
         solpcon = 0.
         minpacon = (minpacon * vol_init + (P_adsorbed_mass - ads_inflow)) / (vol_init + 0.0001)
         dispin = (dispin * ht1%flo - ads_inflow) / ht1%flo
         minpain = (minpain * ht1%flo + ads_inflow) / ht1%flo
      endif
   else
      P_desorbed_mass = -P_adsorbed_mass
      if (P_desorbed_mass < minpacon * vol_init) then
         minpacon = (minpacon * vol_init - P_desorbed_mass) / (vol_init + 0.0001)
         solpcon = (solpcon * vol_init + P_desorbed_mass) / (vol_init + 0.0001)
      else
         ads_inflow = P_desorbed_mass - minpacon * vol_init
         minpacon = 0.
         solpcon = (solpcon * vol_init + (P_desorbed_mass - ads_inflow)) / (vol_init + 0.0001)
         minpain = (minpain * ht1%flo - ads_inflow) / ht1%flo
         dispin = (dispin * ht1%flo + ads_inflow) / ht1%flo
      endif
   endif

   !!Diffusion between the aerobic and anaerobic layers (or anaerobic through aerobic to water column)
   diff_aer2anaer = rs2_s * (wbody_bed%aer_disp - wbody_bed%anaer_disp) * tstep / 10. ! unit grams P, / 10 added 03/26/22, KDW
   !Update mixed wtclm concentrations before calculating diffusion with anaerobic
  !  wtclm_mixed_minpa = (minpacon * vol_init + minpain * ht1%flo)
  !  wtclm_mixed_minpa_conc = wtclm_mixed_minpa / wbody%flo
  !  wtclm_mixed_solp = (solpcon * vol_init + dispin * ht1%flo)
  !  wtclm_mixed_solp_conc = wtclm_mixed_solp / wbody%flo
   ! frsts_wtclm = adsts_wtclm - wtclm_mixed_minpa_conc
   ! frsts_aer = adsts_aer - wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd
   ! if (frsts_wtclm > 0) then
   !    pai_wtclm = 1 - (1 / (1 + (wtclm_por/((ch_nut(jnut)%adskeq / 100) * frsts_wtclm)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
   ! else
   !    frsts_wtclm = 0
   !    pai_wtclm = 1
   ! endif
   ! if (frsts_aer > 0) then
   !    pai_aer = 1 - (1 / (1 + (wbody_bed%bed_por/((ch_nut(jnut)%adskeq / 100) * frsts_aer)))) !denominator is [ha*mm/kg P] * [mg P/L] * [1/100] = [10 m^3/kg P] * [kg P/10^3 m^3] * [1/100] = unitless
   ! else
   !    frsts_aer = 0
   !    pai_aer = 1
   ! endif
   if ((wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) >  0.) then ! KDW 04/18/22
      pai_wtclm = wtclm_mixed_solp_conc / (wtclm_mixed_solp_conc + wtclm_mixed_minpa_conc) !KDW 04/15/22
   else
      pai_wtclm = 1.
   endif
   if ((wbody_bed%aer_disp * wbody_bed%bed_por + wbody_bed%aer_minpa * wbody_bed%bed_bd) > 0.) then !KDW 04/18/22
      pai_aer = wbody_bed%aer_disp * wbody_bed%bed_por / (wbody_bed%aer_disp * wbody_bed%bed_por + wbody_bed%aer_minpa * wbody_bed%bed_bd) !KDW 04/15/22
   else
      pai_aer = 0.01
   endif

         !! Calculate aer conc where flows to/from wtclm and anaer are balance 
   aer_balanced = (wtclm_mixed_solp_conc + wbody_bed%anaer_disp * rs2_ratio) / (1 +rs2_ratio) ! KDW 03/30/22
   !! if anaer conc. is between aer and wtclm, allow aer to rise/fall to level of anaer (if sufficient diffusion potential)
   if ((wbody_bed%aer_disp < wbody_bed%anaer_disp) .and. (wbody_bed%anaer_disp <= wtclm_mixed_solp_conc)) then
      diff = abs(diff_aer2anaer)
      max_diff = (wbody_bed%anaer_disp - wbody_bed%aer_disp) * &
      (1/ ( (pai_aer/(vol_h2o_aer / 1000)) + (pai_aer/(vol_h2o_anaer / 1000)) ) )
      max_diff1 = (wbody_bed%anaer_disp * wbody_bed%bed_por + wbody_bed%anaer_minpa * wbody_bed%bed_bd) * anaer_vol / 1000 

      max_diff_disp = pai_anaer * max_diff ! KDW 03/22/22
      P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por ! KDW 03/22/22
      min_disp_conc = wbody_bed%anaer_disp * wbody_bed%bed_por - (max_diff_disp/(anaer_vol/1000)) ! KDW 03/22/22
      dox_factor = 0.5 
      xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
      max_diff_minpa = (P_tot - min_Ptot_conc)*(anaer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
      max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 04/15/22
      max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
      ! min_Ptot_conc = min_disp_conc *(1+(adsts_anaer/(yyy + min_disp_conc))) ! KDW 03/22/22
      ! max_diff_minpa = (P_tot - min_Ptot_conc)*(anaer_vol/1000) - max_diff_disp  ! KDW 03/22/22
      ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 03/22/22
      ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

      max_diff_disp = pai_aer * max_diff ! KDW 03/07/22
      P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! KDW 03/21/22
      max_disp_conc = wbody_bed%aer_disp * wbody_bed%bed_por + (max_diff_disp/(aer_vol/1000)) ! KDW 03/21/22
      dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
      xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
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
   elseif ((wbody_bed%aer_disp > wbody_bed%anaer_disp) .and. (wbody_bed%anaer_disp >= wtclm_mixed_solp_conc)) then
      diff = abs(diff_aer2anaer)
      max_diff = (wbody_bed%aer_disp - wbody_bed%anaer_disp) * &
      (1/ ( (pai_aer/(vol_h2o_aer / 1000)) + (pai_anaer/(vol_h2o_anaer / 1000)) ) )
      max_diff1 = (wbody_bed%aer_disp * wbody_bed%bed_por + wbody_bed%aer_minpa * wbody_bed%bed_bd) * aer_vol / 1000

      max_diff_disp = pai_aer * max_diff ! KDW 03/22/22
      P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! KDW 03/22/22
      min_disp_conc = wbody_bed%aer_disp * wbody_bed%bed_por - (max_diff_disp/(aer_vol/1000)) ! KDW 03/22/22
      dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
      xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
      max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
      max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 04/15/22
      max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
      ! min_Ptot_conc = min_disp_conc *(1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/22/22
      ! max_diff_minpa = (P_tot - min_Ptot_conc)*(aer_vol/1000) - max_diff_disp  ! KDW 03/22/22
      ! max_diff_minpa = MIN(max_diff_minpa, ((1-pai_aer) * max_diff)) ! KDW 03/22/22
      ! max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 03/22/22

      max_diff_disp = pai_anaer * max_diff ! KDW 03/07/22
      P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por ! KDW 03/21/22
      max_disp_conc = wbody_bed%anaer_disp * wbody_bed%bed_por + (max_diff_disp/(anaer_vol/1000)) ! KDW 03/21/22
      dox_factor = 0.5
      xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
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
      if (wbody_bed%anaer_disp < wtclm_mixed_solp_conc) then
         if (wbody_bed%aer_disp < aer_balanced) then
         !! disP will move through aer between anaer and wtclm
            diff_anaer2wtclm = -diff_aer2anaer 
            diff_aer2anaer  = 0.
         else
         !! allow aer_disp to be brought down as low as a aer_balanced
            ! how much diffusion to bring aer to aer_balanced
            P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! KDW 03/21/22
            min_disp_conc = aer_balanced * wbody_bed%bed_por ! KDW 03/21/22, 04/02/22
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
            yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
            ! min_Ptot_conc = min_disp_conc *(1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/21/22
            max_diff = (P_tot - min_Ptot_conc)*(aer_vol/1000)  ! KDW 03/21/22
            ! make sure doesn't cause P in anaer to exceed aer_balance
            P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por ! KDW 04/04/22
            max_disp_conc = aer_balanced * wbody_bed%bed_por ! KDW 04/04/22
            dox_factor = 0.5
            xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
            yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
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
      else !i.e. wbody_bed%anaer_disp >= wtclm_mixed_solp_conc
         if (wbody_bed%aer_disp > aer_balanced) then
         !! disP will move through aer between anaer and wtclm
            diff_anaer2wtclm = -diff_aer2anaer 
            diff_aer2anaer  = 0.
         else
         !! allow aer_disp to be brought up as high as a aer_balanced
            ! how much diffusion to bring aer to aer_balanced
            P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! KDW 03/21/22
            max_disp_conc = aer_balanced * wbody_bed%bed_por ! KDW 03/21/22, 04/02/22
            dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
            xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
            yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
            ! max_Ptot_conc = max_disp_conc *(1+(adsts_aer/(yyy + max_disp_conc))) ! KDW 03/21/22
            max_diff = (max_Ptot_conc - P_tot)*(aer_vol/1000)  ! KDW 03/21/22
            ! make sure enough P in anaer to supply diffusion
            P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por ! KDW 04/04/22
            min_disp_conc = aer_balanced * wbody_bed%bed_por ! KDW 04/04/22
            dox_factor = 0.5 
            xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
            yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
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
      aer_disp_temp = wbody_bed%aer_disp - diff_aer2anaer / (vol_h2o_aer/1000)
      dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
      P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + aer_disp_temp * wbody_bed%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22, 04/27/22
      xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
         4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
      if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
      if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
      if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
      aer_minpa_temp = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
      aer_disp_temp = (P_tot - (aer_minpa_temp * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume 
      if (aer_disp_temp < .0) aer_disp_temp = 0.! kDW 04/28/22                 
      ! aer_disp_temp = wbody_bed%aer_disp - diff_aer2anaer / (vol_h2o_aer/1000)
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
      ! P_tot = wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd + aer_disp_temp * wbody_bed%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
      ! P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
      ! P_adsorbed = P_min_new - wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
      ! P_adsorbed = MIN(P_adsorbed, aer_disp_temp * wbody_bed%bed_por) !  "* wbody_bed%bed_por" added 03/21/22
      ! P_adsorbed = MAX(P_adsorbed, (-wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd))
      ! aer_minpa_temp =  wbody_bed%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
      ! aer_disp_temp = aer_disp_temp - P_adsorbed / wbody_bed%bed_por !  "/ wbody_bed%bed_por" added 03/21/22

      anaer_disp_temp = wbody_bed%anaer_disp + diff_aer2anaer / (vol_h2o_anaer/1000)
      dox_factor = 0.5 
      P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + anaer_disp_temp * wbody_bed%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22, 04/27/22
      xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
         4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
      if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
      if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
      if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
      anaer_minpa_temp = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
      anaer_disp_temp = (P_tot - (anaer_minpa_temp * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume 
      if (anaer_disp_temp < .0) anaer_disp_temp = 0.! kDW 04/28/22           
      ! anaer_disp_temp = wbody_bed%anaer_disp + diff_aer2anaer / (vol_h2o_anaer/1000)
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
      ! P_tot = wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd + anaer_disp_temp * wbody_bed%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
      ! P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
      ! P_adsorbed = P_min_new - wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
      ! P_adsorbed = MIN(P_adsorbed, anaer_disp_temp * wbody_bed%bed_por) !  "* wbody_bed%bed_por" added 03/21/22
      ! P_adsorbed = MAX(P_adsorbed, (-wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd))
      ! anaer_minpa_temp =  wbody_bed%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
      ! anaer_disp_temp = anaer_disp_temp - P_adsorbed / wbody_bed%bed_por !  "/ wbody_bed%bed_por" added 03/21/22
   endif
   
   !! allow exchange between anaer and wtclm, if diffusion potential available
   if(diff_anaer2wtclm > 0) then
      !! check diffusion for overshooting equilibrium
      diff = abs(diff_anaer2wtclm)
      max_diff = (anaer_disp_temp - wtclm_mixed_solp_conc) * &
      (1/ ( (pai_wtclm/qdin) + (pai_anaer/(vol_h2o_anaer / 1000)) ) )
      max_diff1 = (anaer_disp_temp * wbody_bed%bed_por + anaer_minpa_temp * wbody_bed%bed_bd) * anaer_vol / 1000

      max_diff_disp = pai_anaer * max_diff ! KDW 03/22/22
      P_tot = anaer_minpa_temp * wbody_bed%bed_bd + anaer_disp_temp * wbody_bed%bed_por ! KDW 03/22/22
      min_disp_conc = anaer_disp_temp * wbody_bed%bed_por - (max_diff_disp/(anaer_vol/1000)) ! KDW 03/22/22
      dox_factor = 0.5 
      xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
      max_diff_minpa = (P_tot - min_Ptot_conc)*(anaer_vol/1000) - max_diff_disp  ! KDW 04/15/22, grams
      max_diff_minpa = MIN(max_diff_minpa, ((1-pai_anaer) * max_diff)) ! KDW 04/15/22
      max_diff2 =  max_diff_disp + max_diff_minpa ! KDW 04/15/22
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/22/22
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
      max_diff = (wtclm_mixed_solp_conc - wbody_bed%anaer_disp) * &
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
      P_tot = anaer_minpa_temp * wbody_bed%bed_bd + anaer_disp_temp * wbody_bed%bed_por ! KDW 03/21/22
      max_disp_conc = anaer_disp_temp * wbody_bed%bed_por + (max_diff_disp/(anaer_vol/1000)) ! KDW 03/21/22
      dox_factor = 0.5
      xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
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

   !! update wtclm_mixed_solp_conc and wbody_bed%anaer_disp (use temp variables)
   wtclm_mixed_solp_conc_temp = wtclm_mixed_solp_conc + diff_anaer2wtclm / wbody%flo
   ! yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) 
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
   P_adsorbed_mass = (P_min_new * xx - wtclm_mixed_minpa_conc)* wbody%flo ! KDW 04/15/22 
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
   ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
   P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + anaer_disp_temp * wbody_bed%bed_por ! unit here is mg/L TOTAL volume (not just pore volume), KDW , " / ch(jrch)%bed_por" added 03/21/22
   dox_factor = 0.5
   xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
   yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
   P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
      4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
   if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
   if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
   if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
   anaer_minpa_temp = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
   anaer_disp_temp = (P_tot - (anaer_minpa_temp * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume   
   if ( anaer_disp_temp < .0)  anaer_disp_temp = 0.! kDW 04/28/22     
   ! P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
   ! P_adsorbed = P_min_new - wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd ! unit mg/L TOTAL volume
   ! P_adsorbed = MIN(P_adsorbed, anaer_disp_temp * wbody_bed%bed_por) !  "* wbody_bed%bed_por" added 03/21/22
   ! P_adsorbed = MAX(P_adsorbed, (-wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd))
   ! anaer_minpa_temp =  wbody_bed%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
   ! anaer_disp_temp = anaer_disp_temp - P_adsorbed / wbody_bed%bed_por !  "/ wbody_bed%bed_por" added 03/21/22


   if(abs(diff_anaer2wtclm) > 0) then
      !! find new rs2_balanced aer conc
      aer_balanced = (wtclm_mixed_solp_conc_temp + anaer_disp_temp * rs2_ratio) / (1 +rs2_ratio) ! KDW 03/30/22

      !! allow aer_disp to reach aer_balanced again
      P_tot = aer_minpa_temp * wbody_bed%bed_bd + aer_disp_temp * wbody_bed%bed_por ! KDW 03/21/22
      min_disp_conc = aer_balanced * wbody_bed%bed_por ! KDW 03/21/22, 04/02/22
      dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
      xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
      yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
      Ptot_conc_balanced = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
      ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/21/22
      ! Ptot_conc_balanced = min_disp_conc * (1+(adsts_aer/(yyy + min_disp_conc))) ! KDW 03/21/22
      diff = (P_tot - Ptot_conc_balanced)*(aer_vol/1000)  ! KDW 03/21/22
      !! assume that change in aer_disp (to aer_balanced) is due to flow with whichever pool is in correct direction (anaerobic or watercolumn)
      !!!!!!!!! need to check if sufficient P available in source pool!!!!!!!!!!!!
      if (diff > 0) then ! i.e. P leaving aer
         if (wtclm_mixed_solp_conc_temp > anaer_disp_temp) then
            ! make sure enough anaer doesn't overshoot aer_balanced
            P_tot = anaer_minpa_temp * wbody_bed%bed_bd + anaer_disp_temp * wbody_bed%bed_por ! KDW 04/04/22
            max_disp_conc = aer_balanced * wbody_bed%bed_por ! KDW 04/04/22
            dox_factor = 0.5
            xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
            yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            max_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+max_disp_conc))+1)*max_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
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
            P_tot = anaer_minpa_temp * wbody_bed%bed_bd + anaer_disp_temp * wbody_bed%bed_por ! KDW 04/04/22
            min_disp_conc = aer_balanced * wbody_bed%bed_por ! KDW 04/04/22
            dox_factor = 0.5 
            xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
            yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
            min_Ptot_conc = ((xx*ch_nut(jnut)%srpm2*dox_factor/(yy+min_disp_conc))+1)*min_disp_conc! KDW 04/15/22, g/m^3 bulk
            ! yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 04/04/22
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
   if (diff_anaer2wtclm + diff_aer2wtclm2 < 0) then  ! KDW 04/01/22, added diff_aer2wtclm2
      diff = abs(diff_anaer2wtclm + diff_aer2wtclm2)  ! KDW 04/01/22, added diff_aer2wtclm2
      wbody_bed%anaer_disp = wbody_bed%anaer_disp + (-diff_anaer2wtclm*pai_anaer/(vol_h2o_anaer/1000)) ! KDW 04/01/22, changed from "diff" to "-diff_anaer2wtclm"
      wbody_bed%anaer_minpa = wbody_bed%anaer_minpa + (-diff_anaer2wtclm*(1-pai_anaer)/anaer_mass) ! KDW 04/01/22, changed from "diff" to "-diff_anaer2wtclm"
      wbody_bed%aer_disp = wbody_bed%aer_disp + (-diff_aer2wtclm2*pai_aer/(vol_h2o_aer/1000))! KDW 04/01/22
      wbody_bed%aer_minpa = wbody_bed%aer_minpa + (-diff_aer2wtclm2*(1-pai_aer)/aer_mass)! KDW 04/01/22

      if((diff*pai_wtclm)<(dispin*ht1%flo))then
         dispin = dispin - (diff*pai_wtclm/ht1%flo)
         diffremain = diff - diff*pai_wtclm
      else
         diffremain = diff - (dispin*ht1%flo)
         dispin = 0
      endif
      if ((diff*(1-pai_wtclm))<(minpain*ht1%flo))then
         minpain = minpain - (diff*(1-pai_wtclm)/ ht1%flo)
         diffremain = diffremain - diff*(1-pai_wtclm)
      else
         diffremain = diffremain - (minpain*ht1%flo)
         minpain = 0
      endif
      diff = diffremain
      if(diff > 0.)then
         if((diff*pai_wtclm)<(solpcon*vol_init))then
            solpcon = solpcon - (diff*pai_wtclm/(vol_init + 0.0001))
            diffremain = diff - diff*pai_wtclm
         else
            diffremain = diff - (solpcon*vol_init)
            solpcon = 0
         endif
         if ((diff*(1-pai_wtclm))<(minpacon*vol_init))then
            minpacon = minpacon - (diff*(1-pai_wtclm)/ (vol_init + 0.0001))
            diffremain = diffremain - diff*(1-pai_wtclm)
         else
            diffremain = diffremain - (minpacon*vol_init)
            minpacon = 0
         endif
      endif
      if(diffremain > 0.) then !necessary because don't know fraction of P from inflow vs storage
         if(diffremain < (dispin*ht1%flo))then
            dispin = dispin - (diffremain/ht1%flo)
            diffremain = 0
         else
            diffremain = diffremain - (dispin*ht1%flo)
            dispin = 0
            if(diffremain < (solpcon*vol_init))then
               solpcon = solpcon - (diffremain/(vol_init + 0.0001))
               diffremain = 0
            else
               diffremain = diffremain - (solpcon*vol_init)
               solpcon = 0
               if(diffremain < (minpain*ht1%flo))then
                  minpain = minpain - (diffremain/ht1%flo)
                  diffremain = 0
               else
                  diffremain = diffremain - (minpain*ht1%flo)
                  minpain = 0
                  if(diffremain < (minpacon*vol_init))then
                     minpacon = minpacon - (diffremain/(vol_init + 0.0001))
                     diffremain = 0
                  else
                     diffremain = diffremain - (minpacon*vol_init)
                     minpacon = 0
                  endif
               endif
            endif
         endif
      endif 
   endif
   if (diff_anaer2wtclm + diff_aer2wtclm2 > 0) then ! KDW 04/01/22, added diff_aer2wtclm2
      diff = abs(diff_anaer2wtclm + diff_aer2wtclm2) ! KDW 04/01/22, added diff_aer2wtclm2
      if ((diff*pai_wtclm/ht1%flo) < 10.) then ! in case ht1%flo is low, don't want excessive dispin causing instability
         dispin = dispin + (diff*pai_wtclm/ht1%flo)
         minpain = minpain + (diff*(1-pai_wtclm)/ht1%flo)
      else
         dispin = dispin + 10
         minpain = minpain + 10*(1-pai_wtclm)/pai_wtclm
         diffremain = diff - 10*(1/pai_wtclm) * ht1%flo
         solpcon = solpcon + (diffremain*pai_wtclm/(vol_init + 0.0001))
         minpacon = minpacon + (diffremain*(1-pai_wtclm)/(vol_init + 0.0001))
      endif
      wbody_bed%anaer_disp = wbody_bed%anaer_disp - (diff_anaer2wtclm*pai_anaer/(vol_h2o_anaer/1000)) ! KDW 04/01/22, changed from "diff" to "diff_anaer2wtclm"
      wbody_bed%anaer_minpa = wbody_bed%anaer_minpa - (diff_anaer2wtclm*(1-pai_anaer)/anaer_mass) ! KDW 04/01/22, changed from "diff" to "diff_anaer2wtclm"\
      wbody_bed%aer_disp = wbody_bed%aer_disp - (diff_aer2wtclm2*pai_aer/(vol_h2o_aer/1000))! KDW 04/01/22
      wbody_bed%aer_minpa = wbody_bed%aer_minpa - (diff_aer2wtclm2*(1-pai_aer)/aer_mass)! KDW 04/01/22
   endif
   if (diff_aer2anaer > 0) then
      diff = abs(diff_aer2anaer)
      wbody_bed%anaer_disp = wbody_bed%anaer_disp + (diff*pai_anaer/(vol_h2o_anaer/1000))
      wbody_bed%anaer_minpa = wbody_bed%anaer_minpa + (diff*(1-pai_anaer)/anaer_mass)
      wbody_bed%aer_disp = wbody_bed%aer_disp - (diff*pai_aer/(vol_h2o_aer/1000))
      wbody_bed%aer_minpa = wbody_bed%aer_minpa - (diff*(1-pai_aer)/aer_mass)
   endif
   if (diff_aer2anaer < 0) then
      diff = abs(diff_aer2anaer)
      wbody_bed%anaer_disp = wbody_bed%anaer_disp - (diff*pai_anaer/(vol_h2o_anaer/1000))
      wbody_bed%anaer_minpa = wbody_bed%anaer_minpa - (diff*(1-pai_anaer)/anaer_mass)
      wbody_bed%aer_disp = wbody_bed%aer_disp + (diff*pai_aer/(vol_h2o_aer/1000))
      wbody_bed%aer_minpa = wbody_bed%aer_minpa + (diff*(1-pai_aer)/aer_mass)
   endif
   dispin = MAX(0.,dispin)
   solpcon = MAX(0.,solpcon)
   minpain = MAX(0.,minpain)
   minpacon = MAX(0.,minpacon)
   wbody_bed%aer_disp = MAX(0.,wbody_bed%aer_disp)
   wbody_bed%aer_minpa = MAX(0.,wbody_bed%aer_minpa)
   wbody_bed%anaer_disp = MAX(0.,wbody_bed%anaer_disp)
   wbody_bed%anaer_minpa = MAX(0.,wbody_bed%anaer_minpa)

   wbody_wb%phos_diff = (diff_aer2wtclm + diff_aer2wtclm2 + diff_anaer2wtclm) / 1000 ! KDW, 06/07, 04/18/22 added diff_aer2wtclm2
    
    w12 = Theta(ch_nut(inut)%w12,thw12,wtmp) * 10 * wbody_wb%area_ha
    !dox_factor = wbody%dox / (wbody%dox + ch_nut(inut)%dox_half) ! KDW 10/06
    dox_factor = 1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))) ! KDW 10/06/21, replace above to be more sensistive near DO threshold and use external DO meas
    cbod_factor = wbody%cbod / (wbody%cbod + ch_nut(inut)%cbod_half)
    w12 = w12 * dox_factor * cbod_factor
    !if (w12 > 1) w12 = 1 ! commented out 06/01, KDW
    btrb_orgp = w12 * (wbody_bed%aer_orgp - wbody_bed%anaer_orgp)
    !for section below, had "aer_mass / 1000" and "anaer_mass / 1000" when /1000 doesn't belong, fixed 02/24/22, KDW
    equil_conc = ((wbody_bed%aer_orgp * aer_mass) + (wbody_bed%anaer_orgp * anaer_mass)) / & ! KDW 06/02, and all equil_conc terms below
      ((aer_mass) + (anaer_mass))
    if (btrb_orgp > (wbody_bed%aer_orgp - equil_conc) * aer_mass .and. btrb_orgp > 0) then ! added ".and. btrb_orgp >= 0" KDW 11/04/21, and similar below
      btrb_orgp = (wbody_bed%aer_orgp - equil_conc) * aer_mass ! added "* anaer_vol / 1000" on 06/01, KDW
    endif
    if (btrb_orgp < -(wbody_bed%anaer_orgp - equil_conc) * anaer_mass.and. btrb_orgp < 0) then
      btrb_orgp = -(wbody_bed%anaer_orgp - equil_conc) * anaer_mass ! added "* aer_vol / 1000" on 06/01, KDW
    endif
    btrb_minpa = w12 * (wbody_bed%aer_minpa - wbody_bed%anaer_minpa)
    equil_conc = ((wbody_bed%aer_minpa * aer_mass) + (wbody_bed%anaer_minpa * anaer_mass)) / & ! KDW 06/02, and all equil_conc terms below
      ((aer_mass) + (anaer_mass))
    if (btrb_minpa > (wbody_bed%aer_minpa - equil_conc) * aer_mass .and. btrb_minpa > 0) then
      btrb_minpa = (wbody_bed%aer_minpa - equil_conc) * aer_mass ! added "* anaer_vol / 1000" on 06/01, KDW
    endif
    if (btrb_minpa < -(wbody_bed%anaer_minpa - equil_conc) * anaer_mass .and. btrb_minpa  < 0) then
      btrb_minpa = -(wbody_bed%anaer_minpa - equil_conc) * anaer_mass ! added "* aer_vol / 1000" on 06/01, KDW
    endif
    btrb_minps = w12 * (wbody_bed%aer_minps - wbody_bed%anaer_minps)
    equil_conc = ((wbody_bed%aer_minps * aer_mass) + (wbody_bed%anaer_minps * anaer_mass)) / & ! KDW 06/02, and all equil_conc terms below
      ((aer_mass) + (anaer_mass))
    if (btrb_minps > (wbody_bed%aer_minps - equil_conc) * aer_mass .and. btrb_minps > 0) then
      btrb_minps = (wbody_bed%aer_minps - equil_conc) * aer_mass ! added "* anaer_vol / 1000" on 06/01, KDW
    endif
    if (btrb_minps < -(wbody_bed%anaer_minps - equil_conc) * anaer_mass .and. btrb_minps < 0) then
      btrb_minps = -(wbody_bed%anaer_minps - equil_conc) * anaer_mass ! added "* aer_vol / 1000" on 06/01, KDW
    endif
    
    ! wbody_bed%aer_disp = wbody_bed%aer_disp - (diff_aer2anaer/(aer_vol/1000)) - (diff_aer2wtclm/(aer_vol/1000))
    ! wbody_bed%anaer_disp = wbody_bed%anaer_disp + (diff_aer2anaer/(anaer_vol/1000)) - (diff_anaer2wtclm/(anaer_vol/1000)) ! KDW 06/15
    !Commented out above with below because diffusion mass transfer accounted for above now, KDW 02/25/22
    wbody_bed%aer_orgp = wbody_bed%aer_orgp - btrb_orgp/aer_mass
    wbody_bed%anaer_orgp = wbody_bed%anaer_orgp + btrb_orgp/anaer_mass
    wbody_bed%aer_minpa = wbody_bed%aer_minpa - btrb_minpa/aer_mass
    wbody_bed%anaer_minpa = wbody_bed%anaer_minpa + btrb_minpa/anaer_mass
    wbody_bed%aer_minps = wbody_bed%aer_minps - btrb_minps/aer_mass
    wbody_bed%anaer_minps = wbody_bed%anaer_minps + btrb_minps/anaer_mass

    !! Update aerobic and anaerobic layers for adsorption
    !!Assume adsorbed and dissolved reach equilibrium - KDW 03/03/22
    !Aerobic layer
    dox_factor = 0.5 + 0.5 * (1 / ( 1 + exp(-(dox_meas%ts(time%day,time%yrs) - ch_nut(jnut)%dox_half))))
    P_tot = wbody_bed%aer_minpa * wbody_bed%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
    xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/15/22
    yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
    P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
       4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2
    if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
    if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
    if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
    wbody_bed%aer_minpa = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
    wbody_bed%aer_disp = (P_tot - (wbody_bed%aer_minpa * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
    if ( wbody_bed%aer_disp < .0)  wbody_bed%aer_disp = 0.! kDW 04/28/22 
   !  yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
   !  P_tot = wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd + wbody_bed%aer_disp * wbody_bed%bed_por
   !  P_min_new = (adsts_aer + P_tot + yyy - ((adsts_aer + P_tot + yyy)**2 - 4*P_tot*adsts_aer) **0.5) / 2
   !  P_adsorbed = P_min_new - wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd
   !  P_adsorbed = MIN(P_adsorbed, wbody_bed%aer_disp * wbody_bed%bed_por)
   !  P_adsorbed = MAX(P_adsorbed, (-wbody_bed%aer_minpa * ch_sed(jsed)%bed_bd))
   !  wbody_bed%aer_minpa =  wbody_bed%aer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
   !  wbody_bed%aer_disp = wbody_bed%aer_disp - P_adsorbed / wbody_bed%bed_por

    !Anaerobic layer
    dox_factor = 0.5
    P_tot = wbody_bed%anaer_minpa * wbody_bed%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por ! bulk concentrations, g/m^3! KDW 04/15/22
    xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/15/22
    yy = wbody_bed%bed_por / ch_nut(jnut)%adskeq ! KDW 04/15/22
    P_min_new = (xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot - ((xx*ch_nut(jnut)%srpm2*dox_factor + yy + P_tot)**2 - &
       4*xx*ch_nut(jnut)%srpm2*dox_factor*P_tot) **0.5) / (2*xx) ! KDW 04/15/22, g/m^2, dox_factor only for aerobic and wtclm
    if (P_min_new < 0.) P_min_new = 0. ! KDW 04/15/22
    if (P_min_new > ch_nut(jnut)%srpm2) P_min_new = ch_nut(jnut)%srpm2 ! KDW 04/15/22
    if (P_min_new * xx > P_tot) P_min_new = P_tot / xx ! KDW 04/28/22
    wbody_bed%anaer_minpa = P_min_new * xx / wbody_bed%bed_bd ! KDW 04/15/22, back to g/ kg
    wbody_bed%anaer_disp = (P_tot - (wbody_bed%anaer_minpa * wbody_bed%bed_bd)) / wbody_bed%bed_por ! KDW 04/15/22, back to g/m^3 pore water volume
    if ( wbody_bed%anaer_disp < .0)  wbody_bed%anaer_disp = 0.! kDW 04/28/22 
   !  yyy = wbody_bed%bed_por/(ch_nut(jnut)%adskeq / 100) ! KDW 03/03/22
   !  P_tot = wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd + wbody_bed%anaer_disp * wbody_bed%bed_por
   !  P_min_new = (adsts_anaer + P_tot + yyy - ((adsts_anaer + P_tot + yyy)**2 - 4*P_tot*adsts_anaer) **0.5) / 2
   !  P_adsorbed = P_min_new - wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd
   !  P_adsorbed = MIN(P_adsorbed, wbody_bed%anaer_disp * wbody_bed%bed_por)
   !  P_adsorbed = MAX(P_adsorbed, (-wbody_bed%anaer_minpa * ch_sed(jsed)%bed_bd))
   !  wbody_bed%anaer_minpa =  wbody_bed%anaer_minpa + P_adsorbed / ch_sed(jsed)%bed_bd
   !  wbody_bed%anaer_disp = wbody_bed%anaer_disp - P_adsorbed / wbody_bed%bed_por

    roc_aer = bk * (7. * wbody_bed%aer_minpa - wbody_bed%aer_minps) ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
    if (roc_aer > 0.) then
        roc_aer = Min(roc_aer, wbody_bed%aer_minpa)
    endif
    if (roc_aer < 0.) then
        roc_aer = roc_aer * .1
        roc_aer = Max(roc_aer, -wbody_bed%aer_minps)
    endif
    roc_anaer = bk * (7. * wbody_bed%anaer_minpa - wbody_bed%anaer_minps) ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
    if (roc_anaer > 0.) then
        roc_anaer = Min(roc_anaer, wbody_bed%anaer_minpa)
    endif
    if (roc_anaer < 0.) then
        roc_anaer = roc_anaer * .1
        roc_anaer = Max(roc_anaer, -wbody_bed%anaer_minps)
    endif

    if (wbody_bed%aer_minpa - roc_aer < 0)  roc_aer = wbody_bed%aer_minpa ! KDW 03/03/22
    xx = SA_VOL_aer * (1-wbody_bed%bed_por) ! KDW 04/18/22
    if ((wbody_bed%aer_minpa - roc_aer) > (ch_nut(jnut)%srpm2* xx / wbody_bed%bed_bd)) then! KDW 04/18/22
        roc_aer = wbody_bed%aer_minpa - (ch_nut(jnut)%srpm2* xx / wbody_bed%bed_bd) ! KDW 03/03/22! KDW 04/18/22
    endif
    if (wbody_bed%anaer_minpa - roc_anaer < 0)  roc_anaer = wbody_bed%anaer_minpa ! KDW 03/03/22
    xx = SA_VOL_anaer * (1-wbody_bed%bed_por) ! KDW 04/18/22
    if ((wbody_bed%anaer_minpa - roc_anaer) > (ch_nut(jnut)%srpm2* xx / wbody_bed%bed_bd)) then! KDW 04/18/22
        roc_anaer = wbody_bed%anaer_minpa - (ch_nut(jnut)%srpm2* xx / wbody_bed%bed_bd)! KDW 03/03/22! KDW 04/18/22
    endif

    wbody_bed%aer_minps = wbody_bed%aer_minps + roc_aer
    if (wbody_bed%aer_minps < 0.) wbody_bed%aer_minps = 0.
    wbody_bed%anaer_minps = wbody_bed%anaer_minps + roc_anaer
    if (wbody_bed%anaer_minps < 0.) wbody_bed%anaer_minps = 0.

    wbody_bed%aer_minpa = wbody_bed%aer_minpa - roc_aer  ! KDW 03/03/22
    if (wbody_bed%aer_minpa < 0.) wbody_bed%aer_minpa = 0.
    wbody_bed%anaer_minpa = wbody_bed%anaer_minpa - roc_anaer  ! KDW 03/03/22
    if (wbody_bed%anaer_minpa < 0.) wbody_bed%anaer_minpa = 0.

    !! Update water column pools, adding new m factors from above - KDW

    !! calculate organic phosphorus concentration at end of day QUAL2E section 3.3.6 equation III-24
    bc4_m = wq_k2m(tday,tstep,-bc4_k,orgpcon,orgpin) 
    !rs5_k=0. ! rs5 no longer used, all settling included in sediment deposition - KDW
    !if (rchdep > 0.001) rs5_k = Theta(ch_nut(jnut)%rs5,thrs5,wtmp) / rchdep 
    !factk=-rs5_k
    factk = 0.
    alg_m_orgp = alg_md * ch_nut(inut)%ai2 ! KDW
    factm =bc4_m + alg_m_orgp ! added alg_m_orgp - KDW
    con_out = wq_semianalyt (tday, tstep, factm, factk, orgpcon, orgpin)
    wbody%orgp = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%orgp < 1.e-6) wbody%orgp = 0.

    !Water column adsorption-desorption - KDW 03/03/22
    wtclm_mixed_minpa = (minpacon * vol_init + minpain * ht1%flo) ! KDW 02/24/22
    wtclm_mixed_minpa_conc = wtclm_mixed_minpa / wbody%flo! KDW 02/24/22
    wtclm_mixed_minps = (minpscon * vol_init + minpsin * ht1%flo) ! KDW 02/24/22
    wtclm_mixed_minps_conc = wtclm_mixed_minps / wbody%flo! KDW 02/24/22
    wtclm_mixed_solp = (solpcon * vol_init + dispin * ht1%flo)! KDW 02/24/22
    wtclm_mixed_solp_conc = wtclm_mixed_solp / wbody%flo! KDW 02/24/22
    !wtclm_por = 1 - (sedin / wbody%flo) / 1.2 , KDW 04/18/22
   !  yyy = wtclm_por/(ch_nut(jnut)%adskeq / 100) ! KDW 02/26/22
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
   !  P_min_new = (adsts_wtclm + P_tot + yyy - ((adsts_wtclm + P_tot + yyy)**2 - 4*P_tot*adsts_wtclm) **0.5) / 2
   !  P_adsorbed = P_min_new - wtclm_mixed_minpa_conc 
   !  P_adsorbed = MIN(P_adsorbed, wtclm_mixed_solp_conc)
   !  P_adsorbed = MAX(P_adsorbed, -wtclm_mixed_minpa_conc)

    !! calculate dissolved phosphorus concentration at end of day QUAL2E section 3.4.2 equation III-25
    factk = 0.
    alg_m_disp = alg_mg * ch_nut(inut)%ai2 ! KDW
    factm = -bc4_m - alg_m_disp - P_adsorbed ! ads_m replaced with P_adsorbed, KDW 03/03/22
    con_out = wq_semianalyt (tday, tstep, factm, 0., solpcon, dispin)
    wbody%solp = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%solp < 1.e-6) wbody%solp = 0.

    !! calculate active mineral phosphorus concentration at end of day
    !minpain = minpain + P_adsorbed * tstep! ads_m replaced with P_adsorbed, KDW 03/03/22
    if (wbody%flo > 0) then  
      ssd_m1_max = wtclm_mixed_minpa_conc ! KDW 02/24/22
      ssd_m2_max = wtclm_mixed_minps_conc ! KDW 02/24/22
   else
      ssd_m1_max = 0
      ssd_m2_max = 0            
   endif
    factk = 7. * bk ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
    ssd_m1 = -wq_k2m (tday, tstep, -factk, minpacon, minpain)!Made RHS and factk negative on 05/24, KDW 
    ssd_m1 = Min(ssd_m1, ssd_m1_max)  
    factk = bk
    ssd_m2 = -wq_k2m (tday, tstep, -factk, minpscon, minpsin)!Made RHS and factk negative on 05/24, KDW   
    ssd_m2 = Min(ssd_m2, ssd_m2_max)
    if (ssd_m1 > ssd_m2) then
       ssd_m = ssd_m1 - ssd_m2
       ssd_m = Min(ssd_m, minpacon)
    else 
       factk = 7. * bk * 0.1 ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, 02/17/22
       ssd_m1 = -wq_k2m (tday, tstep, -factk, minpacon, minpain)!Made RHS and factk negative on 05/24, KDW 
       ssd_m1 = Min(ssd_m1, ssd_m1_max)  
       factk = bk * 0.1
       ssd_m2 = -wq_k2m (tday, tstep, -factk, minpscon, minpsin)!Made RHS and factk negative on 05/24, KDW  
       ssd_m2 = Min(ssd_m2, ssd_m2_max) 
       ssd_m = ssd_m1 - ssd_m2
       ssd_m = Max(ssd_m, -minpscon)
    endif

    factk = 0.
    factm = -ssd_m + P_adsorbed ! ads_m/P_adsorbed used here rather than line 1288 KDW 03/03/22
    con_out = wq_semianalyt (tday, tstep, factm, factk, minpacon, minpain)
    wbody%minpa = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%minpa < 1.e-6) wbody%minpa = 0.
    xx = SA_VOL_wtclm * ((sedin / qdin) / 1.2) !KDW 04/18/22
    if (wbody%minpa > (ch_nut(jnut)%srpm2 * xx)) then ! KDW 11/30/21  !KDW 04/18/22 replaced adsts with (ch_nut(jnut)%srpm2 * xx) here and next 3 lines
      xyz = wbody%minpa - (ch_nut(jnut)%srpm2 * xx)
      wbody%minpa = (ch_nut(jnut)%srpm2 * xx)
      wbody%solp = wbody%solp + xyz
    endif

    !! calculate stable mineral phosphorus concentration at end of day
    factk = 0.
    factm = ssd_m
    con_out = wq_semianalyt (tday, tstep, factm, factk, minpscon, minpsin)
    wbody%minps = con_out * (wbody%flo / 1000.)
    con_out = 0.
    if (wbody%minps < 1.e-6) wbody%minps = 0.

    !! end phosphorus calculations

  else
    !! all water quality variables set to zero when no water in reservoir
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
    !ch(jrch)%algae = 0.0
    wbody%chla = 0.0
    wbody%orgn = 0.0
    wbody%nh3 = 0.0
    wbody%no2 = 0.0
    wbody%no3 = 0.0
    wbody%orgp = 0.0
    wbody%minpa = 0.0 ! - KDW
    wbody%minps = 0.0 ! - KDW
    wbody%solp = 0.0
    wbody%cbod = 0.0
    wbody%dox = 0.0
    soxy = 0.0
    wbody_wb%phos_diff = 0.0 ! KDW 04/18/22
  endif
  
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


  !! calculate amount of nutrients leaving reservoir - commented out because outflow calculated after sediment settling
  ! ht2%no3 = wbody%no3 * ht2%flo / wbody%flo
  ! ht2%orgn = wbody%orgn * ht2%flo / wbody%flo
  ! ht2%orgp = wbody%orgp * ht2%flo / wbody%flo
  ! ht2%minpa = wbody%minpa * ht2%flo / wbody%flo
  ! ht2%minps = wbody%minps * ht2%flo / wbody%flo
  ! ht2%solp = wbody%solp * ht2%flo / wbody%flo
  ! ht2%chla = wbody%chla * ht2%flo / wbody%flo
  ! ht2%nh3 = wbody%nh3 * ht2%flo / wbody%flo
  ! ht2%no2 = wbody%no2 * ht2%flo / wbody%flo

  !! Some old res_nutrient code, in case need during debugging - KDW
  !!      integer, intent (in) :: iob
  !!      integer :: jres            !none          |reservoir number
  !!      integer :: inut            !none          |counter
  !!      integer :: iwst            !none          |weather station number
  !!      iwst = ob(iob)%wst
  !!      tpco = 1.e+6 * (wbody%orgp + wbody%minpa + wbody%minps + wbody%solp) / &
  !!      (wbody%flo + ht2%flo) ! - KDW
  !!      if (tpco > 1.e-4) then
  !!      !! equation 29.1.6 in SWAT manual
  !!      chlaco = res_nut(inut)%chlar * 0.551 * (tpco**0.76)
  !!      wbody%chla = chlaco * (wbody%flo + ht2%flo) * 1.e-6
  !!      endif
  return
  end subroutine res_nutrient2