      subroutine soiltest_init (isol, isolt)
    
      use soil_module  
      use soil_data_module
      use organic_mineral_mass_module
      use basin_module, only : bsn_prm, bsn_cc
      use landuse_data_module, only: lum
      
      implicit none 
      
      integer :: ly          !none         |counter
      integer, intent (in) :: isol        !none         |counter
      integer, intent (in) :: isolt       !             |
      real :: dep_frac       !             |  
      real :: wt1, zdst, solp, actp, ssp
      real :: temp
      real :: dep_prev ! KDW 11/22/21
      real :: adsts ! KDW 11/30/21
      real :: aa, bb, cc, bray_P ! KDW 12/14/21
      real :: xx, yy, SA_VOL, wt2, P_min_new, P_tot ! KDW 04/14/22

      do ly = 1, soil(isol)%nly !isol is actually ihru - KDW

!!Below added by KDW (from soil_nutcarb_init)
          wt1 = soil(isol)%phys(ly)%bd * soil(isol)%phys(ly)%thick / 100.    ! mg/kg => kg/ha
          wt2 = 0.1 / (soil(isol)%phys(ly)%thick / 1000.)  ! kg/ha => g/m^3 ! KDW 04/14/22

          !set initial no3 pools
          if (solt_db(isolt)%inorgn <= 0.) then
            !zdst = Exp(-soil(isol)%phys(ly)%d / 1000.)
            !soil1(isol)%mn(ly)%no3 = 10. * zdst * .7
            soil1(isol)%mn(ly)%no3 = 7 ! KDW , 05/27
          else
            soil1(isol)%mn(ly)%no3 = solt_db(isolt)%inorgn ! KDW , 05/27
          end if
          soil1(isol)%mn(ly)%no3 =  soil1(isol)%mn(ly)%no3 * wt1      !! mg/kg => kg/ha
          
          if (bsn_cc%sol_P_model == 0 .or. bsn_cc%sol_P_model == 1) then  
              
            ! commented out by KDW 05/11/22
              ! !set initial labile P pool  
              ! if (solt_db(isolt)%watersol_p > 0.0001) then
              !   soil1(isol)%mp(ly)%lab = solt_db(isolt)%watersol_p * wt1   !! mg/kg => kg/ha
              ! else
              ! !! assume initial concentration of 5 mg/kg
              !   soil1(isol)%mp(ly)%lab = 5. * wt1
              ! end if
              ! !! set active P pool based on dynamic PSP MJW
              ! !! if (bsn_cc%sol_P_model == 0) then !! commented out so as to run no matter the model option (also, should have run if == 1 originally) - KDW
              !   !! Allow Dynamic PSP Ratio
              ! !! convert to concentration
              ! solp = soil1(isol)%mp(ly)%lab / wt1
              ! !! PSP = -0.045*log (% clay) + 0.001*(Solution P, mg kg-1) - 0.035*(% Organic C) + 0.43
              ! if (soil(isol)%phys(ly)%clay > 0.) then
              !   bsn_prm%psp = -0.045 * log(soil(isol)%phys(ly)%clay) + (0.001 * solp)
              !   bsn_prm%psp = bsn_prm%psp - (0.035 * soil1(isol)%tot(ly)%c / soil1(isol)%tot(ly)%m) + 0.43 
              ! else
              !   bsn_prm%psp = 0.4
              ! endif   		

                  !! Limit PSP range
              if (bsn_prm%psp < .05) then
                bsn_prm%psp = 0.05
              else if (bsn_prm%psp > 0.9) then
                bsn_prm%psp = 0.9
              end if
              !! end if
              ! Bray P based initialization - KDW 05/11/22
              if (solt_db(isolt)%bray_strong_p > 0.0001) then
                soil1(isol)%mp(ly)%lab = 1.5 * solt_db(isolt)%bray_strong_p * wt1 * bsn_prm%psp 
              else
                soil1(isol)%mp(ly)%lab = 5.
              endif
              soil1(isol)%mp(ly)%act = soil1(isol)%mp(ly)%lab * (1. - bsn_prm%psp) / bsn_prm%psp
          else
          !!! Calculate initial labile P according to Bray P, KDW 02/14/22 (replaces formulation from KDW on 11/30/21), updated for errors on 04/06/22 KDW
              ! adsts = 10 * soil(isol)%phys(ly)%thick * bsn_prm%spkgc * soil(isol)%phys(ly)%bd *    &
              !   (soil(isol)%phys(ly)%clay + .2 * soil(isol)%phys(ly)%silt +                 &
              !   .01 * soil(isol)%phys(1)%sand) / 100 ! kg P / ha
              ! if (1.5 * solt_db(isolt)%bray_strong_p * wt1 < adsts) then
              !   bray_P = 1.5 * solt_db(isolt)%bray_strong_p * wt1
              ! else
              !   bray_P = adsts
              ! endif
              !   aa = bsn_prm%adskeq
              ! if (solt_db(isolt)%bray_strong_p > 0.0001) then
              !   bb = adsts * bsn_prm%adskeq - (bray_P*bsn_prm%adskeq) + soil(isol)%phys(ly)%fc !! wts for conversion mg/kg => kg/ha, assume at field capacity, KDW 04/06/22
              !   cc = -bray_P*soil(isol)%phys(ly)%fc
              ! else
              !   bb = adsts * bsn_prm%adskeq - 30 + soil(isol)%phys(ly)%fc
              !   cc = -30.
              ! endif
              ! bb = Max(1.0, bb)
              ! soil1(isol)%mp(ly)%lab = (-bb + (((bb**2)-4*aa*cc)**0.5))/(2*aa) !quadratic formula, taking positive solution (bb negative)
              ! soil1(isol)%mp(ly)%act = bray_P - soil1(isol)%mp(ly)%lab 
              ! xx = soil1(isol)%mp(ly)%act ! KDW debug 04/13/22
              ! yy = soil(isol)%phys(ly)%fc / bsn_prm%adskeq ! KDW debug 04/13/22
              ! P_min_new = (adsts + bray_P + yy - ((adsts + bray_P + yy)**2 - 4*bray_P*adsts) **0.5) / 2 ! KDW debug 04/13/22
            if (solt_db(isolt)%bray_strong_p > 0.0001) then
              P_tot = 1.5 * solt_db(isolt)%bray_strong_p * wt1 * wt2 ! g/m^3 ! KDW 04/14/22
            else
              P_tot = 30.
            endif
            SA_VOL = 3. * ((soil(isol)%phys(ly)%clay/0.000002) + (soil(isol)%phys(ly)%silt/0.00001) + (soil(isol)%phys(1)%sand/0.0002)) / 100. ! KDW 04/14/22
            xx = SA_VOL * (1-soil(isol)%phys(ly)%por) ! KDW 04/14/22
            yy = soil(isol)%phys(ly)%up / bsn_prm%adskeq ! KDW 04/14/22
            P_min_new = (xx*bsn_prm%srpm2 + yy + P_tot - ((xx*bsn_prm%srpm2 + yy + P_tot)**2 - 4*xx*bsn_prm%srpm2*P_tot) **0.5) / (2*xx) ! KDW 04/14/22, g/m^2
            if (P_min_new < 0.) P_min_new = 0.
            if (P_min_new > bsn_prm%srpm2) P_min_new = bsn_prm%srpm2 ! KDW 04/15/22
            soil1(isol)%mp(ly)%act = P_min_new * xx / wt2 ! kg/ha
            soil1(isol)%mp(ly)%lab = (P_tot / wt2) - soil1(isol)%mp(ly)%act ! kg/ha 
          endif
  
          !! Set Stable pool based on dynamic coefficient
            !! convert to concentration for ssp calculation
            !! Commented out lines 64-70 and replaced with 71 b/c dividing by wt1 not properly handled (units don't match, relationship changes with layer depth...), KDW 10/25/21
            ! actp = soil1(isol)%mp(ly)%act / wt1 
            ! solp = soil1(isol)%mp(ly)%lab / wt1 
            ! !! estimate Total Mineral P in this soil based on data from sharpley 2004
            ! ssp = 25.044 * (actp + solp)** -0.3833
            ! !!limit SSP Range
            ! if (ssp > 7.) ssp = 7.
            ! if (ssp < 1.) ssp = 1.  
            ssp = 7. ! Replaces above, KDW 10/25/21! changed from "4. " to "7. " by KDW based on ratios in Callister and Logan, 02/17/22
            soil1(isol)%mp(ly)%sta = ssp * (soil1(isol)%mp(ly)%act + soil1(isol)%mp(ly)%lab)
!! End addition by KDW (from soil_nutcarb_init)
          ! new depth fraction (average in layer, not according to bottom) by KDW to have more consistent inital nutrient profiles, 11/22/21
          if (ly == 1) then
            dep_prev = 0
          else
            dep_prev = soil(isol)%phys(ly-1)%d
          endif
          dep_frac=-(1/(solt_db(isolt)%exp_co*(soil(isol)%phys(ly)%d-dep_prev))) * &
              (Exp(-solt_db(isolt)%exp_co * soil(isol)%phys(ly)%d) - Exp(-solt_db(isolt)%exp_co * dep_prev))
          ! end new depth fraction by KDW
          !dep_frac=Exp(-solt_db(isolt)%exp_co * soil(isol)%phys(ly)%d) ! commented out and replace with above by KDW, 11/22/21
          !soil1(isol)%mn(ly)%no3 = solt_db(isolt)%inorgn * dep_frac ! replaced with below by KDW
          soil1(isol)%mn(ly)%no3 = soil1(isol)%mn(ly)%no3 * dep_frac ! KDW
          soil1(isol)%mp(ly)%lab = soil1(isol)%mp(ly)%lab * dep_frac ! KDW
          soil1(isol)%mp(ly)%act = soil1(isol)%mp(ly)%act * dep_frac ! KDW
          soil1(isol)%mp(ly)%sta = soil1(isol)%mp(ly)%sta * dep_frac ! KDW
          !soil1(isol)%mp(ly)%lab = solt_db(isolt)%inorgp * dep_frac
          !soil1(isol)%hp(ly)%n = solt_db(isolt)%orgn * dep_frac
          !soil1(isol)%hp(ly)%p = solt_db(isolt)%orgp * dep_frac
          !soil(j)%ly(ly)%watersol_p = solt_db(isolt)%watersol_p* dep_frac
          !soil(j)%ly(ly)%h3a_p = solt_db(isolt)%h3a_p * dep_frac
          !soil(j)%ly(ly)%mehlich_p = solt_db(isolt)%mehlich_p * dep_frac
          !soil(j)%ly(ly)%bray_strong_p = solt_db(isolt)%bray_strong_p    
          !   &                                                      * dep_frac
      end do
      
      return
      end subroutine soiltest_init