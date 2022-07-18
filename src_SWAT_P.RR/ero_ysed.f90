      subroutine ero_ysed
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine predicts daily soil loss caused by water erosion
!!    using the modified universal soil loss equation

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cvm(:)      |none          |natural log of USLE_C (the minimum value
!!                               |of the USLE C factor for the land cover)
!!    hru_km(:)   |km**2         |area of HRU in square kilometers
!!    peakr       |m^3/s         |peak runoff rate
!!    surfq(:)    |mm H2O        |surface runoff for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cklsp(:)    |
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bio_frcov   |              |fraction of cover by biomass - adjusted for
!!                                  canopy height
!!    grcov_fr    |              |fraction of cover by biomass as function of lai
!!    rsd_frcov   |              |fraction of cover by residue
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use hru_module, only : hru, usle_cfac, cklsp, surfq, sedyld, sanyld, silyld, clayld, lagyld, sagyld,  &
         ihru, peakr, usle_ei, surfq_pervious
      use soil_module
      use urban_data_module ! KDW 11/22/21
      
      implicit none

      integer :: j           !none                   |HRU number
      real :: c              !                       |
      real :: usle           !!metric tons/ha        | daily soil loss predicted with USLE equation
      real :: ulu! KDW 

      j = ihru
      ulu = hru(j)%luse%urb_lu
      !! initialize variables
      cklsp(j) = usle_cfac(j) * hru(j)%lumv%usle_mult
      
      !! compute sediment yield with musle
      if (hru(j)%luse%urb_lu > 0) then ! KDW 11/22/21, corrected accounting for urban HRU
        sedyld(j) = (surfq_pervious(j) * peakr * 1000. * hru(j)%km * (1. - urbdb(ulu)%fimp)) ** .56 * cklsp(j)
      else
        sedyld(j) = (surfq(j) * peakr * 1000. * hru(j)%km) ** .56 * cklsp(j)
      endif
      
      if (sedyld(j) < 0.) sedyld(j) = 0.

      !!adjust sediment yield for protection of snow cover
      if (hru(j)%sno_mm > 0.) then
        if (sedyld(j) < 1.e-6) sedyld(j) = 0.0
      else if (hru(j)%sno_mm > 100.) then
        sedyld(j) = 0.
      else
        sedyld(j) = sedyld(j) / Exp(hru(j)%sno_mm * 3. / 25.4)
      end if
      

	!!Particle size distribution of sediment yield
	  sanyld(j) = sedyld(j) * soil(j)%det_san    !! Sand yield
	  silyld(j) = sedyld(j) * soil(j)%det_sil    !! Silt yield
	  clayld(j) = sedyld(j) * soil(j)%det_cla    !! Clay yield
	  sagyld(j) = sedyld(j) * soil(j)%det_sag    !! Small Aggregate yield
	  lagyld(j) = sedyld(j) * soil(j)%det_lag    !! Large Aggregate yield

      !! compute erosion with usle (written to output for comparison)
      usle = 1.292 * usle_ei * cklsp(j) / 11.8

      return
      end subroutine ero_ysed