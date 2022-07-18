      subroutine ero_pkq
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine computes the peak runoff rate for each HRU
!!    and the entire subbasin using a modification of the rational formula

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hru_km(:)    |km^2          |area of HRU in square kilometers
!!    ihru         |none          |HRU number
!!    tconc(:)     |hr            |time of concentration for HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    peakr        |m^3/s         |peak runoff rate
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log, Expo

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use hru_module, only: hru, tconc, ihru, peakr, qday, surfq, surfq_pervious
      use hydrograph_module
      use climate_module
      use urban_data_module ! KDW 11/22/21
      
      implicit none

      integer :: j      !none          |HRU number
      real :: altc      !              |
      real :: expo      !              | 
      integer :: iob    !              | 
      real :: ulu ! KDW 
      
      j = ihru
      ulu = hru(j)%luse%urb_lu
      iob = hru(j)%obj_no
      iwst = ob(iob)%wst
      altc = 1. - expo(2. * tconc(j) * Log(1. - wst(iwst)%weat%precip_half_hr))
      if (hru(j)%luse%urb_lu > 0) then ! first part of if statement added by KDW to handle urban sediment yield peak rate, KDW 11/22/21
            peakr = altc * qday * (surfq_pervious(j)/surfq(j)) / tconc(j)           !! mm/h ! should qday be surfq? - KDW
            peakr = peakr * hru(j)%km / 3.6 ! "* (1. - urbdb(ulu)%fimp) / 3.6" ? - KDW          !! m^3/s
      else
            peakr = altc * qday / tconc(j)           !! mm/h 
            peakr = peakr * hru(j)%km / 3.6          !! m^3/s
      endif
    
      return
      end subroutine ero_pkq