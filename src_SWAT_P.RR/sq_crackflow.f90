      subroutine sq_crackflow
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this surboutine modifies surface runoff to account for crack flow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhqday(:)   |mm H2O        |surface runoff for the hour in HRUS
 
!!    surfq(:)    |mm H2O        |surface runoff in the HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhqday(:)   |mm H2O        |surface runoff for the hour in HRU
!!    surfq(:)    |mm H2O        |surface runoff in the HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use basin_module
      use hru_module, only : hru, surfq, hhqday, ihru, voltot, surfq_pervious
      use time_module
      use urban_data_module ! KDW 01/10/22
      
      implicit none

      integer :: j      !none          |HRU number
      real :: voli      !none          |volume available for crack flow
      integer :: ii     !none          |counter
      real :: ulu ! KDW 01/10/22

      j = ihru

      !! subtract crack flow from surface runoff
      if (hru(j)%luse%urb_lu > 0) then !if statement and first code therein added by KDW 01/10/22
        ulu = hru(j)%luse%urb_lu
        if (surfq_pervious(j) > voltot) then
          surfq(j) = surfq(j) - voltot * (1. - urbdb(ulu)%fimp)  ! subtracts pervious surface runoff lost to cracks, scaled by fraction impervious
          surfq_pervious(j) = surfq_pervious(j) - voltot
        else
          surfq(j) = surfq(j) - surfq_pervious(j) * (1. - urbdb(ulu)%fimp) 
          surfq_pervious(j) = 0.
        endif  
      else
        if (surfq(j) > voltot) then
          surfq(j) = surfq(j) - voltot
        else
          surfq(j) = 0.
        endif
      endif 

      if (time%step > 0) then
        voli = 0.
        voli = voltot
        do ii = 1, time%step  !j.jeong 4/24/2009
          if (hhqday(j,ii) > voli) then
            hhqday(j,ii) = hhqday(j,ii) - voli
            voli = 0.
          else
            voli = voli - hhqday(j,ii)
            hhqday(j,ii) = 0.
          endif
        end do
      end if

      return
      end subroutine sq_crackflow