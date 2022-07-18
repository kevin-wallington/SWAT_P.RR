      subroutine nut_pminrl3
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine computes p flux between the labile, active mineral
!!    and stable mineral p pools.     
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use basin_module
      use organic_mineral_mass_module
      use hru_module, only : ihru
      use soil_module
      use output_landscape_module, only : hnb_d
      
      implicit none

      real, parameter :: bk = .0006     !              |
      integer :: j                      !none          |HRU number
      integer :: l                      !none          |counter 
      real :: rto                       !              |
      real :: frsts                     !              |free adsorption sites - KDW
      real :: adsts                     !              |adsorption sites - KDW
      real :: pai                       !              |phosphorus availability index - KDW
      real :: rmp1                      !kg P/ha       |amount of phosphorus moving from the solution
                                        !              |mineral to the active mineral pool in the
                                        !              |soil layer
      real :: roc                       !kg P/ha       |amount of phosphorus moving from the active
                                        !              |mineral to the stable mineral pool in the 
                                        !              |soil layer
      real :: xx,yy,P_tot, P_min_new, SA_VOL, wt2 ! KDW 04/15/22

      j = ihru

      hnb_d(j)%lab_min_p = 0.
      hnb_d(j)%act_sta_p = 0.

      !!rto = bsn_prm%psp / (1. - bsn_prm%psp) - moved/editted by KDW
      
      do l = 1, soil(j)%nly
        !! Below added by KDW
        if (soil1(j)%mp(l)%sta < 1.e-6) soil1(j)%mp(l)%sta = 0. ! KDW, 03/24
        if (soil1(j)%mp(l)%act < 1.e-6) soil1(j)%mp(l)%act = 0. ! KDW, 03/24
        if (soil1(j)%mp(l)%lab <1.e-6) soil1(j)%mp(l)%lab = 0. ! KDW, 03/24

        wt2 = 0.1 / (soil(j)%phys(l)%thick / 1000.)  ! kg/ha => g/m^3 ! KDW 04/14/22

        P_tot = (soil1(j)%mp(l)%lab + soil1(j)%mp(l)%act) * wt2
        SA_VOL = 3. * ((soil(j)%phys(l)%clay/0.000002) + (soil(j)%phys(l)%silt/0.00001) + (soil(j)%phys(1)%sand/0.0002)) / 100. ! KDW 04/14/22
        xx = SA_VOL * (1-soil(j)%phys(l)%por) ! KDW 04/14/22
        yy = ((soil(j)%phys(l)%st + soil(j)%phys(l)%wpmm) / soil(j)%phys(l)%thick) / bsn_prm%adskeq ! KDW 04/14/22
        P_min_new = (xx*bsn_prm%srpm2 + yy + P_tot - ((xx*bsn_prm%srpm2 + yy + P_tot)**2 - 4*xx*bsn_prm%srpm2*P_tot) **0.5) / (2*xx) ! KDW 04/14/22, g/m^2
        if (P_min_new < 0.) P_min_new = 0.
        if (P_min_new > bsn_prm%srpm2) P_min_new = bsn_prm%srpm2 ! KDW 04/15/22
        rmp1 =  P_min_new * xx / wt2 - soil1(j)%mp(l)%act
        soil1(j)%mp(l)%act = P_min_new * xx / wt2 ! kg/ha
        soil1(j)%mp(l)%lab = (P_tot / wt2) - soil1(j)%mp(l)%act ! kg/ha


        ! adsts = 10 * soil(j)%phys(l)%thick * bsn_prm%spkgc * soil(j)%phys(l)%bd *    &
        !         (soil(j)%phys(l)%clay + .2 * soil(j)%phys(l)%silt +                 &
        !         .01 * soil(j)%phys(1)%sand) / 100 ! kg P / ha
        ! frsts = adsts - soil1(j)%mp(l)%act
        ! if (frsts > 0) then
        !     !pai = 1 - (1 / (1 + ((soil(j)%phys(1)%st + soil(j)%phys(1)%wp * soil(j)%phys(1)%thick)/(bsn_prm%adskeq * frsts))))
        !     pai = 1 - (1 / (1 + (soil(j)%phys(l)%st + soil(j)%phys(l)%wpmm)/(bsn_prm%adskeq * frsts))) !updated from wp to wpmm on 11/12/21, KDW
        !  else
        !     frsts = 0
        !     pai = 1
        ! endif
        
        ! rto = pai / (1.0001 - pai) ! kept for correspondence with old code; simplifies to rto = (soil(j)%phys(l)%st + soil(j)%phys(l)%wpmm)/(bsn_prm%adskeq * frsts), KDW
        ! if (rto > 100) rto = 100
        !! End additions by KDW
        ! rmp1 = (soil1(j)%mp(l)%lab - soil1(j)%mp(l)%act * rto)
        !!! mike changed/added per isabelle beaudin"s email from 01/21/09
        !if (rmp1 > 0.) rmp1 = rmp1 * 0.1
        !if (rmp1 < 0.) rmp1 = rmp1 * 0.6
        !! mike changed/added per isabelle beaudin"s email from 01/21/09
        !rmp1 = Min(rmp1, soil1(j)%mp(l)%lab)
        !roc = bk * (4. * soil1(j)%mp(l)%act - soil1(j)%mp(l)%sta)
        !if (roc < 0.) roc = roc * .1
        !roc = Min(roc, soil1(j)%mp(l)%act)
      ! Above commented out and replaced with below by KDW
        ! if (rmp1 > 0.) then 
        !     rmp1 = rmp1 * 0.1 
        !     rmp1 = Min(rmp1, soil1(j)%mp(l)%lab)
        ! endif
        ! if (rmp1 < 0.) then
        !     rmp1 = rmp1 * 0.6
        !     rmp1 = Max(rmp1, -soil1(j)%mp(l)%act)
        ! endif
        roc = bk * (7. * soil1(j)%mp(l)%act - soil1(j)%mp(l)%sta) ! changed from "4. *" to "7. *" by KDW based on ratios in Callister and Logan, and Pizzeghello, 02/17/22
        if (roc > 0.) then
            roc = Min(roc, soil1(j)%mp(l)%act)
        endif
        if (roc < 0.) then
            roc = roc * .1
            roc = Max(roc, -soil1(j)%mp(l)%sta)
        endif
      ! End new rmp/roc code by KDW

        ! if (soil1(j)%mp(l)%act - roc + rmp1 < 0) then ! KDW 11/30/21
        !   if (soil1(j)%mp(l)%act + rmp1 < 0) rmp1 = -soil1(j)%mp(l)%act
        !   if (soil1(j)%mp(l)%act - roc + rmp1 < 0) roc = soil1(j)%mp(l)%act + rmp1
        ! endif
        ! if (soil1(j)%mp(l)%act - roc + rmp1 > adsts) then ! KDW 11/30/21
        !   if (soil1(j)%mp(l)%act + rmp1 > adsts) rmp1 = adsts - soil1(j)%mp(l)%act
        !   if (soil1(j)%mp(l)%act - roc + rmp1 > adsts) roc = soil1(j)%mp(l)%act + rmp1 - adsts
        ! endif
        if (soil1(j)%mp(l)%act - roc < 0) roc = soil1(j)%mp(l)%act ! KDW 04/15/22, replaces above b/c no longer using rmp1
        if (soil1(j)%mp(l)%act - roc > (bsn_prm%srpm2 * xx / wt2)) roc = soil1(j)%mp(l)%act - (bsn_prm%srpm2 * xx / wt2) ! KDW 04/15/22, replaces above b/c no longer using rmp1, 04/18/22, replace adsts with srpms calc

        soil1(j)%mp(l)%sta = soil1(j)%mp(l)%sta + roc
        if (soil1(j)%mp(l)%sta < 1.e-6) soil1(j)%mp(l)%sta = 0. ! changed from < 0. - KDW, 03/24

        soil1(j)%mp(l)%act = soil1(j)%mp(l)%act - roc ! + rmp1, KDW commented out b/c rmp1 just used to hold value, pools already adjusted
        if (soil1(j)%mp(l)%act < 1.e-6) soil1(j)%mp(l)%act = 0. ! changed from < 0. - KDW, 03/24

        ! soil1(j)%mp(l)%lab = soil1(j)%mp(l)%lab - rmp1
        ! if (soil1(j)%mp(l)%lab <1.e-6) soil1(j)%mp(l)%lab = 0. ! changed from < 0. - KDW, 03/24
        
        hnb_d(j)%lab_min_p = hnb_d(j)%lab_min_p + rmp1
        hnb_d(j)%act_sta_p = hnb_d(j)%act_sta_p + roc

      end do

      return
      end subroutine nut_pminrl3