      subroutine rls_routesurf (iob)
      
!!    ~ ~ ~ PURPOSE ~ ~ ~

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!
!!    ht1== deposition: write to deposition.out
!!    ht2== outflow from inflow: added to hru generated flows

      use hru_module, only : hru, ihru, usle_cfac, ls_overq, precip_eff, rls_solp, rls_no3
      use hydrograph_module
      use organic_mineral_mass_module, only : soil1 ! - KDW, updated 04/07
      !use topography_data_module
      
      implicit none
      
      integer :: j              !            |
      integer :: iob            !            |
      integer :: ifield         !            |
      real :: sed               !            |
      real :: trancap           !            |
      real :: tranperc          ! - KDW
      real, parameter :: rtof=0.5         !none          |weighting factor used to partition the runon organic nutrients, copied from pl_fert - KDW

      j = ihru
      ifield = hru(j)%dbs%field

!!    compute infiltration from surface runon to next landscape unit
      ls_overq = ob(iob)%hin_sur%flo        !/ (10. * hru(j)%area_ha)   ! m3/10*ha = mm
      precip_eff = precip_eff + ls_overq
      !rls_solp = ob(iob)%hin_sur%no3 ! - KDW, commented out 04/07
      !rls_no3 = ob(iob)%hin_sur%solp ! - KDW, commented out 04/07
      tranperc = 0 ! - KDW
      
!!    compute infiltration from surface runon to next landscape unit
      !ls_overq = ob(iob)%hin%flo    !/ (10. * hru(j)%area_ha) 
      !if (ls_overq > 1.e-6) then
      !  qs = ls_overq / 24.   !mm/hr
      !  vs = (qs ** .4) * (hru(j)%topo%slope ** .3) / (hru(j)%luse%ovn ** .6)
      !  if (vs > 1.e-6) then
      !    trt = hru(j)%topo%slope_len / (3600. * vs)
      !    inflrout = soil(j)%phys(1)%k * trt
      !    inflrout = Min (inflrout, ls_overq)
      !  else
      !    inflrout = ls_overq
      !  end if
      !  ht1%flo = inflrout
      !  ht2%flo = ls_overq - inflrout
      !end if
      
!!    sediment deposition across the landscape
      sed = ob(iob)%hin_sur%sed / hru(j)%area_ha
      !! use surface runoff (mm) for eiq - m3/(10 * ha) = mm
      trancap = hru(j)%topo%dep_co * usle_cfac(j) * ls_overq *        &
                hru(j)%topo%slope**1.4 * hru(j)%field%wid**1.4
      if (sed > trancap) then
        !ht1%sed = (sed - trancap) * hru(j)%area_ha ! commented out by KDW
        !ht2%sed = trancap * hru(j)%area_ha ! commented out by KDW
        tranperc = trancap/sed
        ht1%sed = (1 - tranperc) * sed * hru(j)%area_ha ! - KDW
        ht2%sed = tranperc * sed * hru(j)%area_ha ! - KDW
        ht1%san = (1 - tranperc) * ob(iob)%hin_sur%san  ! - KDW
        ht2%san = tranperc * ob(iob)%hin_sur%san ! - KDW
        ht1%sil = (1 - tranperc) * ob(iob)%hin_sur%sil  ! - KDW
        ht2%sil = tranperc * ob(iob)%hin_sur%sil ! - KDW
        ht1%cla = (1 - tranperc) * ob(iob)%hin_sur%cla  ! - KDW
        ht2%cla = tranperc * ob(iob)%hin_sur%cla ! - KDW
        ht1%sag = (1 - tranperc) * ob(iob)%hin_sur%sag  ! - KDW
        ht2%sag = tranperc * ob(iob)%hin_sur%sag ! - KDW
        ht1%lag = (1 - tranperc) * ob(iob)%hin_sur%lag  ! - KDW
        ht2%lag = tranperc * ob(iob)%hin_sur%lag ! - KDW
        ! for lines below, add to soil top layers pools instead of ht1
        soil1(j)%microb(1)%n = soil1(j)%microb(1)%n + rtof * (1 - tranperc) * ob(iob)%hin_sur%orgn / hru(j)%area_ha  ! - KDW, also changed from tot to microb
        soil1(j)%hs(1)%n = soil1(j)%hs(1)%n + (1. - rtof) * (1 - tranperc) * ob(iob)%hin_sur%orgn / hru(j)%area_ha  ! - KDW
        ht2%orgn = tranperc * ob(iob)%hin_sur%orgn ! - KDW
        soil1(j)%microb(1)%p = soil1(j)%microb(1)%p + rtof * (1 - tranperc) * ob(iob)%hin_sur%orgp / hru(j)%area_ha  ! - KDW, also changed from tot to microb
        soil1(j)%hp(1)%p = soil1(j)%hp(1)%p + (1. - rtof) * (1 - tranperc) * ob(iob)%hin_sur%orgp / hru(j)%area_ha  ! - KDW
        ht2%orgp = tranperc * ob(iob)%hin_sur%orgp ! - KDW
        soil1(j)%mp(1)%act = soil1(j)%mp(1)%act + (1 - tranperc) * ob(iob)%hin_sur%minpa / hru(j)%area_ha  ! - KDW
        ht2%minpa = tranperc * ob(iob)%hin_sur%minpa ! - KDW
        soil1(j)%mp(1)%sta = soil1(j)%mp(1)%sta + (1 - tranperc) * ob(iob)%hin_sur%minps / hru(j)%area_ha  ! - KDW
        ht2%minps = tranperc * ob(iob)%hin_sur%minps ! - KDW
      else
        ht1%sed = 0.
        ht1%san = 0.
        ht1%sil = 0.
        ht1%cla = 0.
        ht1%sag = 0.
        ht1%lag = 0.
        ht2%sed = sed * hru(j)%area_ha
        ht2%san = ob(iob)%hin_sur%san
        ht2%sil = ob(iob)%hin_sur%sil
        ht2%cla = ob(iob)%hin_sur%cla
        ht2%sag = ob(iob)%hin_sur%sag
        ht2%lag = ob(iob)%hin_sur%lag
        ht2%orgn = ob(iob)%hin_sur%orgn
        ht2%orgp = ob(iob)%hin_sur%orgp
        ht2%minpa = ob(iob)%hin_sur%minpa
        ht2%minps = ob(iob)%hin_sur%minps
      end if
      ht1%no3 = ob(iob)%hin_sur%no3
      ht1%solp = ob(iob)%hin_sur%solp
      
!!    organic nitrogen
!      orgn = hd(ihout)%orgn 
!      cy = hd(ihout)%sed  / (hd(ihout)%flo  + 1.e-6)
!      if (cy > .01) then
!        enratio = .78 * cy ** (-.2468)
!      else
!        enratio = 3.
!      end if
!      enratio = Min(enratio,3.)
!      dr_er = drls * enratio
!      dr_er = Min(dr_er,1.)
!      hd(ihout)%orgn = orgn * dr_er
      
!!    organic phosphorus
!      orgp = hd(ihout)%sedp 
!      hd(ihout)%sedp = orgp * dr_er
      
!!    nitrate (& nitrite)
!      no3 = (hd(ihout)%no3 + hd(ihout)%no2)
!!    soluble phosphorus
!      slbp = hd(ihout)%solp 
!!    soluble pesticide not routed
!!    sorbed pesticide not routed
!!    chlorophyll-a not routed
!!    ammonium
!      nh3 = hd(ihout)%nh3 
!!    CBOD not routed
!!    dissolved oxygen not routed
!!    persistent bacteria not routed
!!    less persistent bacteria not routed

      return
      end subroutine rls_routesurf