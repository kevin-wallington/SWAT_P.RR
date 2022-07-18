      subroutine pest_enrsb

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates the enrichment ratio for nutrient and
!!    pesticide transport with runoff

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    surfq(:)    |mm H2O        |surface runoff generated on day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    enratio     |none          |enrichment ratio
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use hru_module, only : hru, surfq, sedyld, sanyld, silyld, clayld, sagyld, lagyld, ihru, enratio, surfq_pervious
      use urban_data_module
      
      implicit none       
      
      integer :: j           !none          |HRU number
      real :: cy             !              |
      real :: ulu ! KDW 01/10/22

      j = ihru

      if (sedyld(j) < 1.e-4) then
        sedyld(j) = 0.0
        sanyld(j) = 0.0
        silyld(j) = 0.0
        clayld(j) = 0.0
        sagyld(j) = 0.0
        lagyld(j) = 0.0
      endif

      !! CREAMS method for calculating enrichment ratio
      if (hru(j)%luse%urb_lu > 0) then ! if statement and code therein (if = true) added by KDW 01/10/22
        ulu = hru(j)%luse%urb_lu
        cy = .1 * sedyld(j) / ((1. - urbdb(ulu)%fimp) * hru(j)%area_ha * surfq_pervious(j) + 1.e-6)
      else
        cy = .1 * sedyld(j) / (hru(j)%area_ha * surfq(j) + 1.e-6)
      endif
      if (cy > 1.e-6) then
        enratio = .78 * cy ** (-.2468)
      else
        enratio = 0.
      endif

      if (enratio > 3.0) enratio = 3.0

      return
      end subroutine pest_enrsb