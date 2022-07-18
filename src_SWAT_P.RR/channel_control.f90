      subroutine channel_control
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine simulates channel routing     

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition  
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bankst(:)   |m^3 H2O       |bank storage
!!    inum1       |none          |reach number
!!    inum2       |none          |inflow hydrograph storage location number
!!    pet_ch      |mm H2O        |potential evapotranspiration on day
!!    rchdep      |m             |depth of flow on day
!!    rttlc       |m^3 H2O       |transmission losses from reach on day
!!    rtwtr       |m^3 H2O       |water leaving reach on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    revapday    |m^3 H2O       |amount of water moving from bank storage
!!                               |into the soil profile or being taken
!!                               |up by plant roots in the bank storage zone
!!    rtwtr       |m^3 H2O       |water leaving reach on day
!!    qdbank      |m^3 H2O       |streamflow contribution from bank storage
!!    sedrch      |metric tons   |sediment transported out of reach on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter
!!    jrch        |none          |reach number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min
!!    SWAT: rchinit, rtover, rtday, rtmusk, rthourly, rtsed, rthsed, watqual
!!    SWAT: noqual, hhwatqual, hhnoqual, rtpest, rthpest
!!    SWAT: rchuse, reachout

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use climate_module
      use hydrograph_module
      use basin_module
      use channel_data_module
      use time_module
      use channel_module
      use constituent_mass_module
      
      implicit none

      integer :: ii       !units         |description
      !real :: qdbank      !m^3 H2O       |streamflow contribution from bank storage
      real :: revapday    !m^3 H2O       |amount of water moving from bank storage
                          !              |into the soil profile or being taken
                          !              |up by plant roots in the bank storage zone
      real :: sedcon      !g/m^3         |sediment concentration
      real :: reactw      !mg pst        |amount of pesticide in reach that is lost
                          !              |through reactions  
      real :: volatpst    !mg pst        |amount of pesticide lost from reach by
                          !              |volatilization
      real :: setlpst     !mg pst        |amount of pesticide moving from water to
                          !              |sediment due to settling
      real :: resuspst    !mg pst        |amount of pesticide moving from sediment to
                          !              |reach due to resuspension
      real :: difus       !mg pst        |diffusion of pesticide from sediment to reach
      real :: reactb      !mg pst        |amount of pesticide in sediment that is lost
                          !              |through reactions
                          !              |up by plant roots in the bank storage zone
      real :: bury        !mg pst        |loss of pesticide from active sediment layer
                          !              |by burial
      real :: sedpest     !mg pst        |pesticide in river bed sediment 
      real :: soxy             !mg O2/L       |saturation concetration of dissolved oxygen
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
      real :: cinn, cinp       !mg N/L        |effective available nitrogen concentration - cinp added by KDW 12/13/21
      real :: rttlc_d
      !real :: rtdelt ! commented out by KDW
      real :: rtevp_d
      integer :: inum1

      ch_d(jrch)%flo_out = 0.
      ch_d(jrch)%evap = 0.
      jhyd = ch_dat(jrch)%hyd
      jsed = ch_dat(jrch)%sed
      jnut = ch_dat(jrch)%nut

      iwst = ob(icmd)%wst
      pet_ch = wst(iwst)%weat%pet
      qdbank = 0.
      revapday = 0.
      
      !! initialize variables for route command loop
      call ch_rchinit

      !! zero flow out variables
      ob(icmd)%hd(1) = hz
      if (time%step > 0) then
        do ii = 1, time%step
          ob(icmd)%ts(1,ii) = hz
        end do
      end if  

      ch(jrch)%vel_chan = 0.
      ch(jrch)%dep_chan = 0.
	  sedrch = 0.
	  rch_san = 0.
	  rch_sil = 0.
	  rch_cla = 0.
	  rch_sag = 0.
	  rch_lag = 0.
	  rch_gra = 0.
      wtrin = 0.
      chlin = 0.
         algin = 0.
         orgnin = 0.
         ammoin = 0.
         nitritin = 0.
         nitratin = 0.
         orgpin = 0.
         minpain = 0. ! - KDW
         minpsin = 0. ! - KDW
         dispin = 0.
         cbodin = 0.
         disoxin = 0.
         cinn = 0.
!! route water through reach
        rtwtr_d=0.
        rttlc_d=0.
        rtevp_d =0.

!! route water through reach
      if (time%step == 0) then
          rt_delt = 1.
        wtrin = ob(icmd)%hin%flo
        if ( ob(icmd)%hin%flo > 0.001) then
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
        end if            
!! initialize inflow concentrations
        
        if (bsn_cc%rte == 0) call ch_rtday
        if (bsn_cc%rte == 1) call ch_rtmusk

        rtwtr_d = rtwtr ! moved up by KDW from line 253
        rttlc_d = rttlc
        rtevp_d = rtevp
! lines 340-360 moved up by KDW        
!! add transmission losses to bank storage/deep aquifer in subbasin
        if (rttlc > 0.) then
          ch(jrch)%bankst = ch(jrch)%bankst + rttlc * (1.-bsn_prm%trnsrch)
  !!!! Jeff add to hydrograph -----------------------------
          rchsep(jrch) = rttlc * bsn_prm%trnsrch
        end if
   
  !! compute revap from bank storage
        revapday = 0.6 * pet_ch *ch_hyd(jhyd)%l *ch_hyd(jhyd)%w
        revapday = Min(revapday, ch(jrch)%bankst)
        ch(jrch)%bankst = ch(jrch)%bankst - revapday
  
  !! compute contribution of water in bank storage to streamflow
        qdbank = ch(jrch)%bankst * (1. - ch_hyd(jhyd)%alpha_bnk)
        ch(jrch)%bankst = ch(jrch)%bankst - qdbank
        rtwtr_d = rtwtr_d + qdbank
        if (time%step > 0) then
          do ii = 1, time%step
            hrtwtr(ii) = hrtwtr(ii) + qdbank / real(time%step)
          end do
        end if

        ben_area = ch_hyd(jhyd)%l *ch_hyd(jhyd)%w 
        if (bsn_cc%wq == 1 .or. bsn_cc%wq == 0) then
          if (bsn_cc%wq == 1) call ch_watqual3
          if (ch_sed(jsed)%eqn == 0) call ch_rtsed
          if (ch_sed(jsed)%eqn == 1) call ch_rtsed_bagnold
          if (ch_sed(jsed)%eqn == 2) call ch_rtsed_kodatie
          if (ch_sed(jsed)%eqn == 3) call ch_rtsed_Molinas_Wu
          if (ch_sed(jsed)%eqn == 4) call ch_rtsed_yangsand
        endif
        if (bsn_cc%wq == 2) then ! SWAT+P water quality processes and coupling with sediment transport
          call ch_watqual5
          if (ch_sed(jsed)%eqn == 0) call ch_rtsed2
          if (ch_sed(jsed)%eqn == 1) call ch_rtsed_bagnold2 
          if (ch_sed(jsed)%eqn == 2) call ch_rtsed_kodatie2
          if (ch_sed(jsed)%eqn == 3) call ch_rtsed_Molinas_Wu2
          if (ch_sed(jsed)%eqn == 4) call ch_rtsed_yangsand2
        endif

        !! Below moved from prior location (immediately before sediment routing) - KDW
        ob(icmd)%hd(1)%orgn = ch(jrch)%organicn * rtwtr_d /1000. !changed from rtwtr to rtwtr_d, here and next 13 lines, KDW 04/04/22
        !ob(icmd)%hd(1)%sedp = ch(jrch)%organicp * rtwtr_d  /1000.
        ob(icmd)%hd(1)%orgp = ch(jrch)%organicp * rtwtr_d /1000. ! - KDW
        ob(icmd)%hd(1)%minpa = ch(jrch)%minpa * rtwtr_d /1000. ! - KDW
        ob(icmd)%hd(1)%minps = ch(jrch)%minps * rtwtr_d /1000. ! - KDW
        ob(icmd)%hd(1)%no3 = ch(jrch)%nitraten * rtwtr_d / 1000.
        ob(icmd)%hd(1)%solp = ch(jrch)%disolvp * rtwtr_d / 1000.
        ob(icmd)%hd(1)%chla = ch(jrch)%chlora * rtwtr_d / 1000.
        ob(icmd)%hd(1)%nh3 = ch(jrch)%ammonian * rtwtr_d / 1000.
        ob(icmd)%hd(1)%no2 = ch(jrch)%nitriten * rtwtr_d / 1000.
        ob(icmd)%hd(1)%cbod = ch(jrch)%rch_cbod *  rtwtr_d / 1000.
        ob(icmd)%hd(1)%dox = ch(jrch)%rch_dox *  rtwtr_d / 1000.
        !! Above moved from prior location (immediately before sediment routing) - KDW

!        if (cs_db%num_pests > 0) then
!          call ch_rtpest 
!        end if
        
        ! rtwtr_d = rtwtr moved up to line 206, KDW 04/04/22
        ! rttlc_d = rttlc
        ! rtevp_d = rtevp
      else
        rt_delt = 1. / time%step  
        wtrin = ob(icmd)%hin%flo * rt_delt
        if (ob(icmd)%hin%flo > 0.001) then
         chlin = 1000. * ob(icmd)%hin%chla * rt_delt / wtrin 
         algin = 1000. * chlin / ch_nut(jnut)%ai0        !! QUAL2E equation III-1
         orgnin = 1000. * ob(icmd)%hin%orgn  * rt_delt/ wtrin  
         ammoin = 1000. * ob(icmd)%hin%nh3 * rt_delt / wtrin  
         nitritin = 1000. * ob(icmd)%hin%no2 * rt_delt / wtrin
         nitratin = 1000. * ob(icmd)%hin%no3 * rt_delt  / wtrin
         !orgpin = 1000. * ob(icmd)%hin%sedp *rt_delt / wtrin 
         orgpin = 1000. * ob(icmd)%hin%orgp *rt_delt / wtrin ! - KDW
         minpain = 1000. * ob(icmd)%hin%minpa *rt_delt / wtrin ! - KDW
         minpsin = 1000. * ob(icmd)%hin%minps *rt_delt / wtrin ! - KDW
         dispin = 1000. * ob(icmd)%hin%solp * rt_delt  / wtrin
         cbodin = 1000. * ob(icmd)%hin%cbod * rt_delt / wtrin
         disoxin = 1000. * ob(icmd)%hin%dox  * rt_delt / wtrin
        end if
        
        rtwtr = 0.
        do ii = 1, time%step
        ! wtrin = ob(icmd)%ts(1,ii)%flo  
        call ch_rtday
        if (bsn_cc%wq == 1)  call ch_watqual5
!        call ch_rthpest
               rtwtr_d=rtwtr_d+rtwtr
               rttlc_d=rttlc_d+ rttlc
               rtevp_d =rtevp_d+ rtevp

               
        ob(icmd)%hd(1)%orgn = ob(icmd)%hd(1)%orgn+ ch(jrch)%organicn  * rtwtr_d /1000.
        !ob(icmd)%hd(1)%sedp = ob(icmd)%hd(1)%sedp + ch(jrch)%organicp * rtwtr_d /1000.
        ob(icmd)%hd(1)%orgp = ob(icmd)%hd(1)%orgp + ch(jrch)%organicp * rtwtr_d /1000. ! - KDW
        ob(icmd)%hd(1)%minpa = ob(icmd)%hd(1)%minpa + ch(jrch)%minpa * rtwtr_d /1000. ! - KDW
        ob(icmd)%hd(1)%minps = ob(icmd)%hd(1)%minps + ch(jrch)%minps * rtwtr_d /1000. ! - KDW
       ob(icmd)%hd(1)%no3 = ob(icmd)%hd(1)%no3+ ch(jrch)%nitraten * rtwtr_d / 1000.
       ob(icmd)%hd(1)%solp = ob(icmd)%hd(1)%solp+ ch(jrch)%disolvp * rtwtr_d / 1000.
       ob(icmd)%hd(1)%chla = ob(icmd)%hd(1)%chla + ch(jrch)%chlora * rtwtr_d / 1000.
       ob(icmd)%hd(1)%nh3 =  ob(icmd)%hd(1)%nh3+ ch(jrch)%ammonian * rtwtr_d / 1000.
       ob(icmd)%hd(1)%no2 = ob(icmd)%hd(1)%no2 + ch(jrch)%nitriten * rtwtr_d / 1000.
       ob(icmd)%hd(1)%cbod = ob(icmd)%hd(1)%cbod + ch(jrch)%rch_cbod *  rtwtr_d / 1000.
       ob(icmd)%hd(1)%dox = ob(icmd)%hd(1)%dox + ch(jrch)%rch_dox *  rtwtr_d / 1000.    

!       hsedyld(ii) = ob(icmd)%ts(1,ii)%sed 
!       sedrch = sedrch + hsedyld(ii)
        rch_san = 0.
!       rch_sil = rch_sil + hsedyld(ii)  !!All are assumed to be silt type particles
        rch_cla = 0.
        rch_sag = 0.
        rch_lag = 0.
        rch_gra = 0.
        call ch_rthsed

        end do
      
! lines 340-360 moved up by KDW        
!! add transmission losses to bank storage/deep aquifer in subbasin
        if (rttlc > 0.) then
          ch(jrch)%bankst = ch(jrch)%bankst + rttlc * (1.-bsn_prm%trnsrch)
  !!!! Jeff add to hydrograph -----------------------------
          rchsep(jrch) = rttlc * bsn_prm%trnsrch
        end if
   
  !! compute revap from bank storage
        revapday = 0.6 * pet_ch *ch_hyd(jhyd)%l *ch_hyd(jhyd)%w
        revapday = Min(revapday, ch(jrch)%bankst)
        ch(jrch)%bankst = ch(jrch)%bankst - revapday
  
  !! compute contribution of water in bank storage to streamflow
        qdbank = ch(jrch)%bankst * (1. - ch_hyd(jhyd)%alpha_bnk)
        ch(jrch)%bankst = ch(jrch)%bankst - qdbank
        rtwtr_d = rtwtr_d + qdbank
        if (time%step > 0) then
          do ii = 1, time%step
            hrtwtr(ii) = hrtwtr(ii) + qdbank / real(time%step)
          end do
        end if

      endif

!! average daily water depth for sandi doty 09/26/07
      ch(jrch)%dep_chan = rchdep

! moved up to follow ch_rtday by KDW
! !! add transmission losses to bank storage/deep aquifer in subbasin
!       if (rttlc > 0.) then
!         ch(jrch)%bankst = ch(jrch)%bankst + rttlc * (1.-bsn_prm%trnsrch)
! !!!! Jeff add to hydrograph -----------------------------
!         rchsep(jrch) = rttlc * bsn_prm%trnsrch
!       end if
 
! !! compute revap from bank storage
!       revapday = 0.6 * pet_ch *ch_hyd(jhyd)%l *ch_hyd(jhyd)%w
!       revapday = Min(revapday, ch(jrch)%bankst)
!       ch(jrch)%bankst = ch(jrch)%bankst - revapday

! !! compute contribution of water in bank storage to streamflow
!       qdbank = ch(jrch)%bankst * (1. - ch_hyd(jhyd)%alpha_bnk)
!       ch(jrch)%bankst = ch(jrch)%bankst - qdbank
!       rtwtr_d = rtwtr_d + qdbank
!       if (time%step > 0) then
!         do ii = 1, time%step
!           hrtwtr(ii) = hrtwtr(ii) + qdbank / real(time%step)
!         end do
!       end if

!! perform in-stream sediment calculations

!! perform in-stream pesticide calculations
!!      call ch_biofilm
      
!! perform in-stream pathogen calculations
      ! call ch_rtbact  dont merge

!! remove water from reach for consumptive water use
      call ch_rchuse

!! summarize output/determine loadings to next routing unit
      call ch_rtout
      
!   output_channel
      ch_d(jrch)%flo_in = ob(icmd)%hin%flo / 86400.     !m^3 -> m^3/s - KDW
      ch_d(jrch)%flo_out = rtwtr_d / 86400.     !m^3 -> m^3/s - KDW
      ch_d(jrch)%evap =rtevp_d / 10000.                  !m^3 -> ha-m
      ch_d(jrch)%tloss = rttlc_d / 10000.  
      ch_d(jrch)%sed_in = ob(icmd)%hin%sed   
      ch_d(jrch)%sed_out = sedrch              
      ch_d(jrch)%sed_conc = sedrch / (ch_d(jrch)%flo_in + .01)
      ch_d(jrch)%orgn_in = ob(icmd)%hin%orgn   
      ch_d(jrch)%orgn_out = ob(icmd)%hd(1)%orgn              
      !ch_d(jrch)%orgp_in = ob(icmd)%hin%sedp    
      !ch_d(jrch)%orgp_out = ob(icmd)%hd(1)%sedp         
      ch_d(jrch)%orgp_in = ob(icmd)%hin%orgp ! - KDW
      ch_d(jrch)%orgp_out = ob(icmd)%hd(1)%orgp ! - KDW
      ch_d(jrch)%minpa_in = ob(icmd)%hin%minpa ! - KDW
      ch_d(jrch)%minpa_out = ob(icmd)%hd(1)%minpa ! - KDW 
      ch_d(jrch)%minps_in = ob(icmd)%hin%minps ! - KDW   
      ch_d(jrch)%minps_out = ob(icmd)%hd(1)%minps ! - KDW
      ch_d(jrch)%no3_in = ob(icmd)%hin%no3  
      ch_d(jrch)%no3_out = ob(icmd)%hd(1)%no3                     
      ch_d(jrch)%nh4_in = ob(icmd)%hin%nh3    
      ch_d(jrch)%nh4_out = ob(icmd)%hd(1)%nh3               
      ch_d(jrch)%no2_in = ob(icmd)%hin%no2 
      ch_d(jrch)%no2_out = ob(icmd)%hd(1)%no2                       
      ch_d(jrch)%solp_in = ob(icmd)%hin%solp         
      ch_d(jrch)%solp_out = ob(icmd)%hd(1)%solp     
      ch_d(jrch)%chla_in = ob(icmd)%hin%chla     
      ch_d(jrch)%chla_out = ob(icmd)%hd(1)%chla                   
      ch_d(jrch)%cbod_in = ob(icmd)%hin%cbod    
      ch_d(jrch)%cbod_out = ob(icmd)%hd(1)%cbod                    
      ch_d(jrch)%dis_in = ob(icmd)%hin%dox         
      ch_d(jrch)%dis_out = ob(icmd)%hd(1)%dox                                    
      ch_d(jrch)%sand_in = ob(icmd)%hin%san 
      ch_d(jrch)%sand_out = ob(icmd)%hd(1)%san                         
      ch_d(jrch)%silt_in = ob(icmd)%hin%sil         
      ch_d(jrch)%silt_out = ob(icmd)%hd(1)%sil                       
      ch_d(jrch)%clay_in = ob(icmd)%hin%cla            
      ch_d(jrch)%clay_out = ob(icmd)%hd(1)%cla                         
      ch_d(jrch)%smag_in = ob(icmd)%hin%sag             
      ch_d(jrch)%smag_out = ob(icmd)%hd(1)%sag                       
      ch_d(jrch)%lag_in = ob(icmd)%hin%lag           
      ch_d(jrch)%lag_out = ob(icmd)%hd(1)%lag                        
      ch_d(jrch)%grvl_in = ob(icmd)%hin%grv          
      ch_d(jrch)%grvl_out = ob(icmd)%hd(1)%grv                      
      ch_d(jrch)%bnk_ero = bnkrte
      ch_d(jrch)%ch_deg = degrte
!!    Channel Deposition (Only new deposits during the current time step)
      if (ch(jrch)%depch >= ch(jrch)%depprch) then
	    ch_d(jrch)%ch_dep = ch(jrch)%depch - ch(jrch)%depprch
	  else
	    ch_d(jrch)%ch_dep = 0.
	  end if
!!    Floodplain Deposition (Only new deposits during the current time step)
      if (ch(jrch)%depfp >= ch(jrch)%depprfp) then
	    ch_d(jrch)%fp_dep = ch(jrch)%depfp - ch(jrch)%depprfp
	  else
	    ch_d(jrch)%fp_dep = 0.
	  end if
!!    Total suspended sediments (only silt and clay)
      if (ch_sed(jsed)%eqn == 0) then
        ch_d(jrch)%tot_ssed = sedrch
      else
        ch_d(jrch)%tot_ssed = rch_sil + rch_cla
      endif
      
      return

      end subroutine channel_control