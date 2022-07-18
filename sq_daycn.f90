      subroutine sq_daycn

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Predicts daily runoff given daily precipitation and snow melt
!!    using a modified SCS curve number approach

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cnday(:)    |none          |curve number for current day, HRU and at 
!!                               |current soil moisture
!!    fcimp(:)    |fraction      |fraction of HRU area that is classified
!!                               |as directly connected impervious
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use urban_data_module
      use hru_module, only : hru, cnday, surfq, ihru, precip_eff, surfq_pervious, precip_eff_imp ! KDW added precip_eff_imp 01/10/22
      
      implicit none
       
      real :: bb        !none          |variable used to store intermediate 
                        !              |calculation result
      real :: cnimp     !none          |curve number for impervious areas
      integer :: j      !none          |HRU number
      real :: pb        !none          |variable used to store intermediate
                        !              |calculation result
      real :: r2        !none          |retention parameter in CN equation
      real :: surfqimp  !mm H2O        |surface runoff from impervious area
      integer :: ulu    !              |

      j = ihru

      r2 = 25400. / cnday(j) - 254.
      bb = .2 * r2
      pb = precip_eff - bb

      surfq(j) = 0.
      if (pb > 0.) then
        surfq(j) = pb * pb / (precip_eff + .8 * r2)
      end if

      if (hru(j)%luse%urb_lu > 0) then
        ulu = hru(j)%luse%urb_lu
        surfq_pervious(j) = surfq(j) !* (1. - urbdb(ulu)%fimp) ! KDW 11/22/21
        surfqimp = 0.
        cnimp = 98.
        r2 = 25400. / cnimp - 254.
        bb = .2 * r2
        pb = precip_eff_imp - bb ! changed from precip_eff to precip_eff_imp by KDW 01/10/22
        if (pb > 0.) then
          surfqimp = pb * pb / (precip_eff_imp + .8 * r2) ! changed from precip_eff to precip_eff_imp by KDW 01/10/22
        end if
        surfq(j) = surfq(j) * (1. - urbdb(ulu)%fimp) +                 &    ! KDW 12/20/21 editted   "surfq(j) * (1. - urbdb(ulu)%fcimp)" to 'surfq(j) * (1. - urbdb(ulu)%fimp)'           
                                          surfqimp * urbdb(ulu)%fcimp
      end if            
          
      return
      end subroutine sq_daycn