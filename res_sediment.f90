subroutine res_sediment (jres, ihyd, ised)

  use reservoir_data_module
  use reservoir_module
  use conditional_module
  use climate_module
  use time_module
  use hydrograph_module
  use water_body_module
  
  implicit none

  integer, intent (in) :: jres          !none          |reservoir number
  integer, intent (in) :: ihyd          !none          |res hydrologic data pointer
  integer, intent (in) :: ised          !none          |res sediment data pointer
  real :: trapres                       !              |
  real :: susp                          !              |
  real :: velofl                        !              |  
  real :: sed_ppm, sil_ppm, cla_ppm 
  real :: xx, sag_ppm, san_ppm, lag_ppm, grv_ppm
  real :: sed_settle, sil_settle, cla_settle, fines_settle, settle_remain            ! tons - KDW
  real :: san_settle, sag_settle, lag_settle, grv_settle    ! tons - KDW

  sil_settle = 0. ; cla_settle = 0. ; sed_settle = 0. ; fines_settle = 0.
  san_settle = 0. ; sag_settle = 0. ; lag_settle = 0. ; grv_settle = 0.

  if (wbody%flo < 1.e-6) then
    ! reservoir is empty
    wbody%sed = 0.
  else

    !! compute new sediment concentration in reservoir
  if (ht1%sed < 1.e-6) ht1%sed = 0.0      
    !! velsetl = 1.35 for clay particle m/d
  if (wbody_wb%area_ha > 1.e-6 .and. ht2%flo > 1.e-6) then
    velofl = (ht2%flo / wbody_wb%area_ha) / 10000.  ! m3/d / ha * 10000. = m/d - changed from wbody%flo to ht2%flo to match old formulation - KDW
  trapres = res_sed(ised)%velsetlr / velofl 
  else
      trapres = 1 ! - KDW 05/12
  end if
  if (trapres > 1.) trapres = 1.
  susp =  1. - res_sed(ised)%trapeff * trapres !! editted so not overly trapping,(was just "1. - trapres") - KDW

    !! compute concentrations
  if (wbody%flo > 0.) then
      sed_ppm = 1000000. * (ht1%sed * susp + wbody%sed) / wbody%flo
      sed_ppm = Max(1.e-6, sed_ppm)
      sil_ppm = 1000000. * (ht1%sil * susp + wbody%sil) / wbody%flo
      sil_ppm = Max(1.e-6, sil_ppm)
      cla_ppm = 1000000. * (ht1%cla * susp + wbody%cla) / wbody%flo
      cla_ppm = Max(1.e-6, cla_ppm)
  else
      sed_ppm = 1.e-6
      sil_ppm = 1.e-6
      cla_ppm = 1.e-6
  endif
    
    !! compute change in sediment concentration due to settling 
    if (sed_ppm > res_sed(ised)%nsed) then
      sed_ppm = (sed_ppm - res_sed(ised)%nsed) * res_sed(ised)%sed_stlr + res_sed(ised)%nsed
      wbody%sed = sed_ppm * wbody%flo / 1000000.      ! ppm -> t
      
      sil_ppm = (sil_ppm - res_sed(ised)%nsed) * res_sed(ised)%sed_stlr + res_sed(ised)%nsed
      wbody%sil = sil_ppm * wbody%flo / 1000000.      ! ppm -> t
      
      cla_ppm = (cla_ppm - res_sed(ised)%nsed) * res_sed(ised)%sed_stlr + res_sed(ised)%nsed
      wbody%cla = cla_ppm * wbody%flo / 1000000.      ! ppm -> t

      !! assume all sand aggregates and gravel settles
      wbody%san = 0.
      wbody%sag = 0.
      wbody%lag = 0.
      wbody%grv = 0.
    end if

    !! compute sediment leaving reservoir - ppm -> t
    ht2%sed = sed_ppm * ht2%flo / 1000000.
    ht2%sil = sil_ppm * ht2%flo / 1000000.
    ht2%cla = cla_ppm * ht2%flo / 1000000.

  end if

  return
  end subroutine res_sediment