      subroutine ch_watqual3

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

      real :: tday, wtmp, fll, gra
      real :: lambda, fnn, fpp, algi, fl_1, xx, yy, zz, ww, qq
      real :: uu, vv, cordo, f1, algcon, orgncon, nh3con, no2con, no3con
      real :: orgpcon, minpacon, minpscon, solpcon, cbodcon, o2con, wtrtot, bc1mod, bc2mod
      real :: sedpcon, sedpcon2, sedpin, totp ! KDW - all sediment-associated P, which was called "orgp" in old formulation
      real :: thgra = 1.047, thrho = 1.047, thrs1 = 1.024
      real :: thrs2 = 1.074, thrs3 = 1.074, thrs4 = 1.024, thrs5 = 1.024
      real :: thbc1 = 1.083, thbc2 = 1.047, thbc3 = 1.047, thbc4 = 1.047
      real :: thrk1 = 1.047, thrk2 = 1.024, thrk3 = 1.024, thrk4 = 1.060
      real, parameter :: bk = .0006     ! - KDW
      real :: rs3_s, rk4_s, rk1_k, rk1_m, rk3_k, rk2_m, rk2_k ! - KDW
      real :: bc1_k, bc2_m, bc3_k, bc3_m, rs2_s, bc4_k, bc4_m ! - KDW
      real :: chlin            !mg chl-a/L    |chlorophyll-a concentration in inflow - KDW
      real :: algin            !mg alg/L      |algal biomass concentration in inflow- KDW
      real :: orgnin           !mg N/L        |organic N concentration in inflow - KDW
      real :: ammoin           !mg N/L        |ammonium N concentration in inflow - KDW
      real :: nitratin         !mg N/L        |nitrate concentration in inflow - KDW
      real :: nitritin         !mg N/L        |nitrite concentration in inflow - KDW
      real :: orgpin           !mg P/L        |organic P concentration in inflow - KDW
      real :: minpain          !mg P/L        |active mineral P concentration in inflow - KDW
      real :: minpsin          !mg P/L        |stable mineral P concentration in inflow - KDW
      real :: dispin           !mg P/L        |soluble P concentration in inflow - KDW
      real :: cbodin           !mg/L          |carbonaceous biological oxygen demand - KDW
      real :: disoxin          !mg O2/L       |dissolved oxygen concentration in inflow - KDW
      real :: soxy             !mg O2/L       |saturation concetration of dissolved oxygen - KDW
      real :: factk, factm ! - KDW
      real :: disp_init, solp_init ! KDW 05/12/22

      rs3_s = 0. ; rk4_s = 0. ; rk1_k = 0. ; rk1_m = 0. ; rk3_k = 0. ; rk2_m = 0. ; rk2_k = 0. ! - KDW
      bc1_k = 0. ; bc2_m  = 0. ; bc3_k = 0. ; bc3_m = 0. ; rs2_s = 0. ; bc4_k = 0. ; bc4_m = 0. ! - KDW
      factk = 0. ; factm = 0. ! - KDW
      alg_m = 0. ; alg_md = 0. ; alg_mg = 0. ; alg_set = 0. ; algcon_out = 0. ! - KDW
      alg_m_o2 = 0. ; alg_m_disp = 0. ; alg_m_orgp = 0. ; alg_nh4_m = 0. ; alg_no3_m = 0. ! - KDW
      sedpcon2 = 0. ; totp = 0. ! KDW

!!    rnum1        |none          |fraction of overland flow

         !! calculate flow duration
         ! tday = rttime / 24.0 ! replaced with below - KDW
        !! below added for new "residence time" by KDW
         if (wtrin > 0) then
            tday = ch(jrch)%rchstor / wtrin ! KDW - this code really uses this ratio, not truly residence time
         else
            tday = 1000
         endif
         tday = min(1000., tday) ! sets arbitratily large maximum residence time of water
         tday = max(.001, tday) ! sets arbitrarily small minimum residence time of water
        !! above added for new "residence time" by KDW

         !! calculate temperature in stream Stefan and Preudhomme. 1993.  Stream temperature estimation 
         !! from air temperature.  Water Res. Bull. p. 27-45 SWAT manual equation 2.3.13
         wtmp = 5.0 + 0.75 * wst(iwst)%weat%tave ! moved up to precede benthic source/loss calculations - KDW
         if (wtmp <= 0.) wtmp = 0.1 ! moved up to precede benthic source/loss calculations - KDW

!        benthic sources/losses in mg   ! actually, in g - KDW
         rs2_s = ch_nut(jnut)%rs2 * 1000 ! rs2_s parameter is now in g/m^2/day, but formulation below assumes in mg/m^2/day, KDW 05/12/22
         rs2_s =  Theta(rs2_s,thrs2,wtmp)*ch_hyd(jhyd)%l *ch_hyd(jhyd)%w * rt_delt ! rs2_s replaces ch_nut(jnut)%rs2, see above, KDW 05/12/22
         rs3_s =  Theta(ch_nut(jnut)%rs3,thrs3,wtmp)*ch_hyd(jhyd)%l *ch_hyd(jhyd)%w * rt_delt
         rk4_s =  Theta(ch_nut(jnut)%rk4,thrk4,wtmp)*ch_hyd(jhyd)%l *ch_hyd(jhyd)%w * rt_delt

            !! initialize concentration of nutrient in reach
         wtrtot = 0.
         algcon = 0.
         orgncon = 0.
         nh3con = 0.
         no2con = 0.
         no3con = 0.
         orgpcon = 0.
         sedpcon = 0.
         minpacon = 0. ! - KDW
         minpscon = 0. ! - KDW
         solpcon = 0.
         cbodcon = 0.
         o2con = 0.
         ch(jrch)%rch_cbod = amax1(1.e-6, ch(jrch)%rch_cbod)
         
         !! below added by KDW, to calculate incoming concentration
         qdin = wtrin + ch(jrch)%rchstor
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
          sedpin = orgpin + minpain + minpsin ! KDW
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
          if (sedpin < 1.e-6) sedpin = 0.0 ! - KDW, 03/26
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
          sedpin = 0
          minpain = 0
          minpsin = 0
          dispin = 0
          cbodin = 0
          disoxin = 0
         end if    
         !! above added by KDW, to calculate incoming concentration


         !algcon =ch(jrch)%algae !replaced with below by KDW
         algcon = 1000. * ch(jrch)%chlora / ch_nut(jnut)%ai0 ! changed for consistency with channel_control line 145 and watqual4 line 64 - KDW
         orgncon = ch(jrch)%organicn
         nh3con = ch(jrch)%ammonian 
         no2con =  ch(jrch)%nitriten
         no3con = ch(jrch)%nitraten 
         orgpcon =ch(jrch)%organicp 
         minpacon =ch(jrch)%minpa ! - KDW
         minpscon =ch(jrch)%minps ! - KDW
         sedpcon = orgpcon + minpacon + minpscon ! KDW
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
     if (sedpcon < 1.e-6) sedpcon = 0.0 ! - KDW
	   if (solpcon < 1.e-6) solpcon = 0.0
	   if (cbodcon < 1.e-6) cbodcon = 0.0
	   if (o2con < 1.e-6) o2con = 0.0
         
         if (wtrin> 0.0001) then !changed from > 0 - KDW
            disoxin = disoxin - rk4_s /wtrin
            if (disoxin < 0) then   ! loop added to subtract "excess sink" from o2con - KDW
              if (ch(jrch)%rchstor > 0) then
                  o2con = ((o2con * ch(jrch)%rchstor) + (disoxin * wtrin)) / ch(jrch)%rchstor
              else
                  o2con= 0
              endif
              disoxin = 0
            endif ! end added loop, KDW
            if (o2con < 1.e-6) o2con = 0.0
            !disoxin = max(0., disoxin) ! not need b/c above loop covers - KDW
            ch_d(jrch)%phos_diff = rs2_s/ 1000 ! KDW 05/12/22, kg
            disp_init = dispin! KDW 05/12/22, for outputting diffusion
            solp_init = solpcon! KDW 05/12/22, for outputting diffusion
            dispin = dispin + rs2_s / wtrin 
            if (dispin< 0) then   ! loop added to subtract "excess sink" from solpcon - KDW
              if (ch(jrch)%rchstor > 0) then
                  solpcon = ((solpcon * ch(jrch)%rchstor) + (dispin * wtrin)) / ch(jrch)%rchstor
              else
                  solpcon= 0
              endif
              dispin = 0
              ch_d(jrch)%phos_diff = (((solp_init - solpcon) * ch(jrch)%rchstor) + ((disp_init - dispin) * wtrin))/ 1000 ! KDW 05/12/22, kg
            endif ! end added loop, KDW
            if (solpcon < 1.e-6) solpcon = 0.0
            ammoin = ammoin + rs3_s / wtrin
            if (ammoin < 0) then   ! loop added to subtract "excess sink" from nh3con - KDW
              if (ch(jrch)%rchstor > 0) then   
                nh3con = ((nh3con * ch(jrch)%rchstor) + (ammoin * wtrin)) / ch(jrch)%rchstor
              else
                  nh3con = 0
              endif
            ammoin = 0
            endif ! end added loop, KDW
            if (nh3con < 1.e-6) nh3con = 0.0


         !! calculate effective concentration of available nitrogen QUAL2E equation III-15
         cinn = nh3con + no3con

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
        !bc1mod = ch_nut(jnut)%bc1 * cordo !commented out for clarity b/c is not used anywhere - KDW
        !bc2mod = ch_nut(jnut)%bc2 * cordo !commented out for clarity b/c is not used anywhere - KDW
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
         fpp = solpcon / (solpcon + ch_nut(jnut)%k_p)

         !! calculate daylight average, photosynthetically active,
         !! light intensity QUAL2E equation III-8
         !! Light Averaging Option # 2
         iwgn = wst(iwst)%wco%wgn
         if (wgn_pms(iwgn)%daylth > 0.) then
           algi = wst(iwst)%weat%solrad * ch_nut(jnut)%tfact /  wgn_pms(iwgn)%daylth
         else
           algi = 0.00001
         end if

         !! calculate growth attenuation factor for light, based on
         !! daylight average light intensity QUAL2E equation III-7b
         if (rchdep > 0) then !pul inside loop by KDW to avoid divide by zero
            fl_1 = (1. / (lambda * rchdep)) *                               &                             
                Log((ch_nut(jnut)%k_l + algi) / (ch_nut(jnut)%k_l + algi *  &
                (Exp(-lambda * rchdep))))
         else
            fl_1 = 1
         endif
         fll = 0.92 * (wst(iwst)%weat%daylength / 24.) * fl_1 ! Corrected from wgn_pms(iwgn)%daylth (the dormancy threshold) to wst(iwst)%weat%daylength - KDW, 05/24

         !! calculcate local algal growth rate
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
         !ch(jrch)%algae = 0. ! commented out by KDW
         !factm = 0. ! commented out by KDW
         factk = Theta(gra,thgra,wtmp)   ! temperature correction for specific growth rate - added by KDW         
         !alg_m1 = wq_semianalyt (tday, rt_delt, 0., factk, algcon, algin) !commented out, should be k2m, see below - KDW
         alg_m1 = wq_k2m(tday, rt_delt, factk, algcon, algin) !algae growth rate - added by KDW
         factk = Theta(gra,thgra,wtmp) - Theta(ch_nut(jnut)%rhoq, thrho, wtmp) ! this provides the temperature corrected specific net growth rate - note by KDW
         !alg_m = wq_semianalyt(tday, rt_delt, 0., factk, algcon, algin) !commented out, should be k2m, see below - KDW
         alg_m = wq_k2m(tday, rt_delt, factk, algcon, algin) !algea net growth rate - added by KDW
         alg_m2 = alg_m-alg_m1 ! this provides the (negative of) algae death rate - note by KDW
         !alg_no3_m = -alg_m * (1. - f1) * ch_nut(jnut)%ai1 !replaced by line 555 by KDW
         !alg_nh4_m = -alg_m * f1 * ch_nut(jnut)%ai1 !replaced by line 537 by KDW
         !alg_P_m = -alg_m * ch_nut(jnut)%ai2 !replaced by lines 571 and 597 by KDW
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
         
         cbodocon = min (cbodcon,o2con) ! corrected to cbodocon, was cbodo, KDW
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
         if (ch(jrch)%rch_cbod < 1.e-6) ch(jrch)%rch_cbod = 0. ! KDW

         !! calculate dissolved oxygen concentration if reach at 
         !! end of day QUAL2E section 3.6 equation III-28

         rk2_m = Theta (ch_nut(jnut)%rk2, thrk2, wtmp) * soxy ! this is the piece multiplied by saturated oxy level - KDW
         rk2_k = -Theta (ch_nut(jnut)%rk2, thrk2, wtmp) ! made negative, eliminated line 499, and replaced factk with rk2_k in 510 - KDW
         bc1_k = -Theta(ch_nut(jnut)%bc1,thbc1,wtmp) ! moved from line 382 by KDW, and negated to match eqn 7:3.5.1
         bc1_m = wq_k2m (tday, rt_delt, bc1_k, nh3con, ammoin) ! added by KDW ! changed from rk2_k to bc1_k - KDW
         bc2_k = -Theta(ch_nut(jnut)%bc2,thbc2,wtmp) ! moved from line 426 by KDW
         bc2_m= wq_k2m(tday,rt_delt,bc2_k,no2con,nitritin) ! moved from line 427 by KDW

         alg_m_o2 = ch_nut(jnut)%ai4 * alg_m2 + ch_nut(jnut)%ai3 * alg_m1
     
         !factk = - rk2_k
         !factm = rk1_m + rk2_m - rk4_m + bc1_m * ch_nut(jnut)%ai5 + bc2_m * ch_nut(jnut)%ai6  !removing rk4_m b/c already accounted for in line 337 and adding alg_m_o2, see below - KDW
         factm = rk1_m + rk2_m + alg_m_o2 + bc1_m * ch_nut(jnut)%ai5 + bc2_m * ch_nut(jnut)%ai6 ! added by KDW
         ch(jrch)%rch_dox = wq_semianalyt (tday, rt_delt, factm, rk2_k, o2con, disoxin) !factk replaced with rk2_k, KDW
         if (ch(jrch)%rch_dox <0.) ch(jrch)%rch_dox = 0.
!! end oxygen calculations   


!! nitrogen calculations
        !! calculate fraction of algal nitrogen uptake from ammonia pool QUAL2E equation III-18
         f1 = ch_nut(jnut)%p_n * nh3con / (ch_nut(jnut)%p_n * nh3con +     &
         (1. - ch_nut(jnut)%p_n) * no3con + 1.e-6) !! moved up, had been after organic N calc, KDW

         !! calculate organic N concentration at end of day
         !! QUAL2E section 3.3.1 equation III-16
         !bc1_k = Theta(ch_nut(jnut)%bc1,thbc1,wtmp) ! moved to line 500 by KDW
         bc3_k = Theta(ch_nut(jnut)%bc3,thbc3,wtmp) 
         rs4_k=0.
         if (rchdep > 0.001)  rs4_k = Theta (ch_nut(jnut)%rs4, thrs4, wtmp) / rchdep   

         bc3_m = wq_k2m (tday, rt_delt, -bc3_k, orgncon, orgnin)
         factk =-rs4_k
         alg_orgN_m = -alg_m2 * ch_nut(jnut)%ai1 ! algae death contribution - added by KDW
         !factm = bc3_m  ! replaced with below by KDW
         factm = bc3_m + alg_orgN_m ! added by KDW
         ch(jrch)%organicn = wq_semianalyt (tday, rt_delt, factm, factk, orgncon, orgnin)

          if (ch(jrch)%organicn < 1.e-6) ch(jrch)%organicn = 0.

        !! calculate ammonia nitrogen concentration at end of day QUAL2E section 3.3.2 equation III-17
         alg_nh4_m = alg_m1 * f1 * ch_nut(jnut)%ai1 !moved down from line 350, -alg_m changed to +alg_m1 by KDW
         !factk=-bc1_k ! folded in to lines 500-501, KDW
         !bc1_m=wq_k2m(tday,rt_delt,factk,nh3con,ammoin) ! folded in to lines 500-501, KDW
         !factm=bc1_m-bc3_m  !replaced by below by KDW
         factm= bc1_m - bc3_m - alg_nh4_m ! added by KDW
         ch(jrch)%ammonian= wq_semianalyt(tday,rt_delt,factm,0.,nh3con,ammoin)
         if (ch(jrch)%ammonian < 1.e-6) ch(jrch)%ammonian = 0.
  
        !! calculate concentration of nitrite at end of day QUAL2E section 3.3.3 equation III-19
         !bc2_k = -Theta(ch_nut(jnut)%bc2,thbc2,wtmp) ! moved to line 403 by KDW
         !bc2_m= wq_k2m(tday,rt_delt,bc2_k,no2con,nitritin) ! moved to line 404 by KDW
 
        factm=-bc1_m+bc2_m
        ch(jrch)%nitriten = wq_semianalyt(tday,rt_delt,factm,0.,no2con,nitritin)
        if (ch(jrch)%nitriten < 1.e-6) ch(jrch)%nitriten = 0.

        !! calculate nitrate concentration at end of day QUAL2E section 3.3.4 equation III-20
        factk = 0.
        alg_no3_m = alg_m1 * (1. - f1) * ch_nut(jnut)%ai1 !moved down from line 349, -alg_m changed to +alg_m1 by KDW
        !factm = -bc2_m+alg_m_no3 !replaced by below by KDW
        factm = -bc2_m - alg_no3_m !change sign for alg_m_no3 and change to alg_no3_m - added by KDW
        
        ch(jrch)%nitraten = wq_semianalyt (tday, rt_delt, factm,0., no3con, nitratin)
        if (ch(jrch)%nitraten < 1.e-6) ch(jrch)%nitraten = 0.
!! end nitrogen calculations

!! phosphorus calculations
        !! calculate organic phosphorus concentration at end of day QUAL2E section 3.3.6 equation III-24
         bc4_k = Theta(ch_nut(jnut)%bc4,thbc4,wtmp)
         !ch(jrch)%nitriten = 0.  ! commented out - KDW
         bc4_m = wq_k2m(tday,rt_delt,-bc4_k,sedpcon,sedpin) 
         ch_d(jrch)%phos_dep = - bc4_m * qdin / 1000 ! KDW 05/12/22
         rs5_k=0.
         if (rchdep > 0.001) rs5_k = Theta(ch_nut(jnut)%rs5,thrs5,wtmp) / rchdep 
         factk=-rs5_k
        alg_m_sedp = -alg_m2 * ch_nut(jnut)%ai2 ! KDW, replaces line 351

        !factm =bc4_m !replace by below by KDW
        factm =bc4_m + alg_m_sedp ! added alg_m_sedp term - KDW

        sedpcon2 = wq_semianalyt (tday, rt_delt, factm, factk, sedpcon, sedpin)! - KDW
        if (sedpcon2 < 1.e-6) sedpcon2 = 0.! - KDW
        totp = sedpcon * ch(jrch)%rchstor + sedpin * wtrin ! KDW
        if (totp >  1.e-6) then
            ch(jrch)%organicp =  ((orgpcon * ch(jrch)%rchstor + orgpin * wtrin) &
            / totp) * sedpcon2 ! - KDW
            ch(jrch)%minpa = ((minpacon * ch(jrch)%rchstor + minpain * wtrin) &
              / totp) * sedpcon2   ! - KDW
            ch(jrch)%minps = ((minpscon * ch(jrch)%rchstor + minpsin * wtrin) &
            / totp) * sedpcon2    ! - KDW
        else
            ch(jrch)%organicp =  sedpcon2 / 3. ! - KDW
            ch(jrch)%minpa = sedpcon2 / 3.   ! - KDW
            ch(jrch)%minps = sedpcon2 / 3.    ! - KDW
        endif
        if (ch(jrch)%organicp < 1.e-6) ch(jrch)%organicp = 0. ! - KDW
        if (ch(jrch)%minpa < 1.e-6) ch(jrch)%minpa = 0. ! - KDW
        if (ch(jrch)%minps < 1.e-6) ch(jrch)%minps = 0. ! - KDW
    
        !! calculate dissolved phosphorus concentration at end of day QUAL2E section 3.4.2 equation III-25
        factk = 0.
        alg_m_disp = alg_m1 * ch_nut(jnut)%ai2 !added by KDW
        !factm = -bc4_m+ch_nut(jnut)%ai2 *alg  !replaced by below by KDW
        factm = -bc4_m - alg_m_disp !following form of other algae death/growth lines - added by KDW
        ch(jrch)%disolvp = wq_semianalyt (tday, rt_delt, factm, 0., solpcon, dispin)
        if (ch(jrch)%disolvp < 1.e-6) ch(jrch)%disolvp = 0.
!! end phosphorus calculations

      else
        !! all water quality variables set to zero when no flow
        algin = 0.0
        chlin = 0.0
        orgnin = 0.0
        ammoin = 0.0
        nitritin = 0.0
        nitratin = 0.0
        orgpin = 0.0
        minpain = 0.0 ! - KDW
        minpsin = 0.0 ! - KDW
        sedpin = 0.0 ! - KDW
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
      return
      end subroutine ch_watqual3
      

