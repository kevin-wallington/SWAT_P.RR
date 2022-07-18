      subroutine hru_hyds
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine summarizes data for subbasins with multiple HRUs and
!!    prints the daily output.hru file

      use hru_module, only : cbodu, chl_a, clayld, doxq, hhsurfq, hru, ihru, itb, lagyld, latq, peakr, percn, qday,  &
         sagyld, sanyld, silyld, sedminpa, sedminps, sedorgn, sedorgp, sepbtm, surqno3, surqsolp, tileno3,     &
         tmpav, uh, sedyld, latno3, qtile, latsed, latsan, latsil, latcla, latsag, latlag,  &
         tilesed, tilesan, tilesil, tilecla, latorgp, tileorgp, latorgn, tileorgn, &
         latminpa, tileminpa, latminps, tileminps, tilesag, tilelag, tilesolp, latsolp, percp
      use hydrograph_module
      use basin_module
      use time_module
      use constituent_mass_module
      use output_ls_pesticide_module
      use urban_data_module ! KDW 01/10/22
      
      implicit none

      integer :: j                   !none          |same as ihru (hru number)
      integer :: kk                  !none          |counter
      integer :: ii                  !none          |counter
      real :: cnv                    !none          |conversion factor (mm/ha => m^3)
      real :: wtmp                   !deg C         |temperature of water in reach
      ! hqd, hsd locally defined. J.Jeong 4/26/2009
      !real :: hqd(4*time%step)
      !real :: hsd(4*time%step)
      !real :: hqdtst(time%step)  
      real :: cnv_m3                 !              |
      real :: cnv_kg                 !              |
      integer :: ihdmx               !              | 
      integer :: iob                 !              |
      integer :: ihyd                !none          |counter  
      integer :: iday                !              | 
      integer :: iday_prev           !              |
      real :: ssq                    !              |
      real :: sumflo                 !              |
      real :: sumflo_day             !              |
      integer :: itot                !              |
      integer :: ib                  !none          |counter
      real :: rto                    !none          |cloud cover factor
      real :: ratio                  !              |   
      integer :: iadj                !none          |counter
      integer :: istep               !none          |counter 
      integer :: ipest               !none          |counter
      integer :: ipath               !none          |counter 
      real :: ulu ! KDW 01/10/22     

      j = ihru
      cnv_m3 = hru(j)%area_ha * 10.
      cnv_kg = hru(j)%area_ha
      ihdmx = 2

      !! assign reach loadings for subbasin
      !! zero out hydrograph storage locations
      iob = icmd 
      ob(icmd)%hd(3) = hz

      !! surface runoff hydrograph (3)
      ob(icmd)%peakrate = peakr
      ob(icmd)%hd(3)%temp = 5. + .75 * tmpav(j)       !!wtmp
      ob(icmd)%hd(3)%flo = qday * cnv_m3              !!qdr m3/d
      ob(icmd)%hd(3)%sed = sedyld(j)                  !!sedyld
      ob(icmd)%hd(3)%orgn = sedorgn(j) * cnv_kg       !!sedorgn
      !ob(icmd)%hd(3)%sedp = (sedorgp(j) + sedminpa(j) +                 &
      !                  sedminps(j)) * cnv_kg         !!sedorgp & sedminps
      ob(icmd)%hd(3)%orgp = sedorgp(j) * cnv_kg       !!sedorgp  - KDW 
      ob(icmd)%hd(3)%minpa = sedminpa(j) * cnv_kg       !!sedminpa  - KDW 
      ob(icmd)%hd(3)%minps = sedminps(j) * cnv_kg       !!sedminps  - KDW 

      ob(icmd)%hd(3)%no3 = surqno3(j) * cnv_kg        !!surqno3 & latno3 & no3gw
      ob(icmd)%hd(3)%solp = surqsolp(j) * cnv_kg      !!surqsolp & sedminpa
      ob(icmd)%hd(3)%chla = chl_a(j) * ob(icmd)%hd(3)%flo / 1000        !!chl_a, changed from "* cnv_kg", KDW 07/15
      ob(icmd)%hd(3)%nh3 = 0.                         !! NH3
      ob(icmd)%hd(3)%no2 = 0.                         !! NO2
      ob(icmd)%hd(3)%cbod = cbodu(j) * ob(icmd)%hd(3)%flo / 1000        !!cbodu, changed from "* cnv_kg", KDW 07/15
      ob(icmd)%hd(3)%dox = doxq(j) * ob(icmd)%hd(3)%flo / 1000        !!dissolved oxygen, changed from "* cnv_kg", KDW 07/15
      ob(icmd)%hd(3)%san = sanyld(j)                  !! detached sand
      ob(icmd)%hd(3)%sil = silyld(j)                  !! detached silt
      ob(icmd)%hd(3)%cla = clayld(j)                  !! detached clay
      ob(icmd)%hd(3)%sag = sagyld(j)                  !! detached small aggregates
      ob(icmd)%hd(3)%lag = lagyld(j)                  !! detached large aggregates


      !set constituents
      do ipest = 1, cs_db%num_pests
        obcs(icmd)%hd(3)%pest(ipest) = (hpestb_d(j)%pest(ipest)%surq + hpestb_d(j)%pest(ipest)%sed) * cnv_kg
      end do
      do ipath = 1, cs_db%num_paths
        obcs(icmd)%hd(3)%path(ipath) = 0
      end do
      
      !recharge hydrograph (2)
      if (hru(j)%luse%urb_lu > 0) then ! if statement and code therein (if = true) added by KDW 01/10/22
        ulu = hru(j)%luse%urb_lu
        ob(icmd)%hd(2)%flo = (1. - urbdb(ulu)%fimp) * sepbtm(j) * cnv_m3           !! recharge flow
        ob(icmd)%hd(2)%no3 = (1. - urbdb(ulu)%fimp) * percn(j) * cnv_kg            !! recharge nitrate
        ob(icmd)%hd(2)%no3 = (1. - urbdb(ulu)%fimp) * percp(j) * cnv_kg            !! recharge dissolved P - KDW
      else
        ob(icmd)%hd(2)%flo = sepbtm(j) * cnv_m3           !! recharge flow
        ob(icmd)%hd(2)%no3 = percn(j) * cnv_kg            !! recharge nitrate
        ob(icmd)%hd(2)%no3 = percp(j) * cnv_kg            !! recharge dissolved P - KDW
      endif
      !set constituents
      do ipest = 1, cs_db%num_pests
        obcs(icmd)%hd(2)%pest(ipest) = hpestb_d(j)%pest(ipest)%perc * cnv_kg
      end do
      do ipath = 1, cs_db%num_paths
        obcs(icmd)%hd(2)%path(ipath) = 0
      end do
      
      !lateral soil flow hydrograph (4)
      if (hru(j)%luse%urb_lu > 0) then ! if statement and code therein (if = true) added by KDW 01/10/22
          ulu = hru(j)%luse%urb_lu
          ob(icmd)%hd(4)%flo = (1. - urbdb(ulu)%fimp) * latq(j) * cnv_m3             !! lateral flow
          ob(icmd)%hd(4)%no3 = (1. - urbdb(ulu)%fimp) * latno3(j) * cnv_kg
          ob(icmd)%hd(4)%solp = (1. - urbdb(ulu)%fimp) * latsolp(j) * cnv_kg           !! lat flow soluble phos -KDW
          ob(icmd)%hd(4)%orgn = (1. - urbdb(ulu)%fimp) * latorgn(j) * cnv_kg           !! lat flow organic n -KDW
          ob(icmd)%hd(4)%orgp = (1. - urbdb(ulu)%fimp) * latorgp(j) * cnv_kg           !! lat flow organic p -KDW
          ob(icmd)%hd(4)%minpa = (1. - urbdb(ulu)%fimp) * latminpa(j) * cnv_kg           !! lat flow active mineral p -KDW
          ob(icmd)%hd(4)%minps = (1. - urbdb(ulu)%fimp) * latminps(j) * cnv_kg           !! lat flow stable mineral p -KDW
          ob(icmd)%hd(4)%chla = (1. - urbdb(ulu)%fimp) * chl_a(j) * ob(icmd)%hd(4)%flo / 1000        !!chl_a, added subsurface chla, cbod, and dox, KDW 12/13/21
          ob(icmd)%hd(4)%cbod = (1. - urbdb(ulu)%fimp) * cbodu(j) * ob(icmd)%hd(4)%flo / 1000        !!cbodu, added subsurface chla, cbod, and dox, KDW 12/13/21
          ob(icmd)%hd(4)%dox = (1. - urbdb(ulu)%fimp) * doxq(j) * ob(icmd)%hd(4)%flo / 1000        !!dissolved oxygen, added subsurface chla, cbod, and dox, KDW 12/13/21
          ob(icmd)%hd(4)%sed = (1. - urbdb(ulu)%fimp) * latsed(j)                  !! detached sand -KDW
          ob(icmd)%hd(4)%san = (1. - urbdb(ulu)%fimp) * latsan(j)                  !! detached sand -KDW
          ob(icmd)%hd(4)%sil = (1. - urbdb(ulu)%fimp) * latsil(j)                  !! detached silt -KDW
          ob(icmd)%hd(4)%cla = (1. - urbdb(ulu)%fimp) * latcla(j)                  !! detached clay
          ob(icmd)%hd(4)%sag = (1. - urbdb(ulu)%fimp) * latsag(j)                  !! detached small aggregates -KDW
          ob(icmd)%hd(4)%lag = (1. - urbdb(ulu)%fimp) * latlag(j)                  !! detached large aggregates -KDW
      else
          ob(icmd)%hd(4)%flo = latq(j) * cnv_m3             !! lateral flow
          ob(icmd)%hd(4)%no3 = latno3(j) * cnv_kg
          ob(icmd)%hd(4)%solp = latsolp(j) * cnv_kg           !! lat flow soluble phos -KDW
          ob(icmd)%hd(4)%orgn = latorgn(j) * cnv_kg           !! lat flow organic n -KDW
          ob(icmd)%hd(4)%orgp = latorgp(j) * cnv_kg           !! lat flow organic p -KDW
          ob(icmd)%hd(4)%minpa = latminpa(j) * cnv_kg           !! lat flow active mineral p -KDW
          ob(icmd)%hd(4)%minps = latminps(j) * cnv_kg           !! lat flow stable mineral p -KDW
          ob(icmd)%hd(4)%chla = chl_a(j) * ob(icmd)%hd(4)%flo / 1000        !!chl_a, added subsurface chla, cbod, and dox, KDW 12/13/21
          ob(icmd)%hd(4)%cbod = cbodu(j) * ob(icmd)%hd(4)%flo / 1000        !!cbodu, added subsurface chla, cbod, and dox, KDW 12/13/21
          ob(icmd)%hd(4)%dox = doxq(j) * ob(icmd)%hd(4)%flo / 1000        !!dissolved oxygen, added subsurface chla, cbod, and dox, KDW 12/13/21
          ob(icmd)%hd(4)%sed = latsed(j)                  !! detached sand -KDW
          ob(icmd)%hd(4)%san = latsan(j)                  !! detached sand -KDW
          ob(icmd)%hd(4)%sil = latsil(j)                  !! detached silt -KDW
          ob(icmd)%hd(4)%cla = latcla(j)                  !! detached clay
          ob(icmd)%hd(4)%sag = latsag(j)                  !! detached small aggregates -KDW
          ob(icmd)%hd(4)%lag = latlag(j)                  !! detached large aggregates -KDW
      endif

      !set constituents
      do ipest = 1, cs_db%num_pests
        obcs(icmd)%hd(4)%pest(ipest) = hpestb_d(j)%pest(ipest)%latq * cnv_kg
      end do
      do ipath = 1, cs_db%num_paths
        obcs(icmd)%hd(4)%path(ipath) = 0
      end do
      
      !tile flow hydrograph (5)
      if (hru(j)%luse%urb_lu > 0) then ! if statement and code therein (if = true) added by KDW 01/10/22
        ulu = hru(j)%luse%urb_lu
        ob(icmd)%hd(5)%flo = (1. - urbdb(ulu)%fimp) * qtile * cnv_m3               !! tile flow
        ob(icmd)%hd(5)%no3 = (1. - urbdb(ulu)%fimp) * tileno3(j) * cnv_kg          !! tile flow nitrate 
        ob(icmd)%hd(5)%solp = (1. - urbdb(ulu)%fimp) * tilesolp(j) * cnv_kg          !! tile flow soluble phos - KDW
        ob(icmd)%hd(5)%orgn = (1. - urbdb(ulu)%fimp) * tileorgn(j) * cnv_kg           !! lat flow organic n -KDW
        ob(icmd)%hd(5)%orgp = (1. - urbdb(ulu)%fimp) * tileorgp(j) * cnv_kg           !! lat flow organic p -KDW
        ob(icmd)%hd(5)%minpa = (1. - urbdb(ulu)%fimp) * tileminpa(j) * cnv_kg           !! lat flow active mineral p -KDW
        ob(icmd)%hd(5)%minps = (1. - urbdb(ulu)%fimp) * tileminps(j) * cnv_kg           !! lat flow stable mineral p -KDW
        ob(icmd)%hd(5)%chla = (1. - urbdb(ulu)%fimp) * chl_a(j) * ob(icmd)%hd(5)%flo / 1000        !!chl_a, added subsurface chla, cbod, and dox, KDW 12/13/21
        ob(icmd)%hd(5)%cbod = (1. - urbdb(ulu)%fimp) * cbodu(j) * ob(icmd)%hd(5)%flo / 1000        !!cbodu, added subsurface chla, cbod, and dox, KDW 12/13/21
        ob(icmd)%hd(5)%dox = (1. - urbdb(ulu)%fimp) * doxq(j) * ob(icmd)%hd(5)%flo / 1000        !!dissolved oxygen, added subsurface chla, cbod, and dox, KDW 12/13/21
        ob(icmd)%hd(5)%sed = (1. - urbdb(ulu)%fimp) * tilesed(j)                  !! detached sand -KDW
        ob(icmd)%hd(5)%san = (1. - urbdb(ulu)%fimp) * tilesan(j)                  !! detached sand -KDW
        ob(icmd)%hd(5)%sil = (1. - urbdb(ulu)%fimp) * tilesil(j)                  !! detached silt -KDW
        ob(icmd)%hd(5)%cla = (1. - urbdb(ulu)%fimp) * tilecla(j)                  !! detached clay -KDW
        ob(icmd)%hd(5)%sag = (1. - urbdb(ulu)%fimp) * tilesag(j)                  !! detached small aggregates -KDW
        ob(icmd)%hd(5)%lag = (1. - urbdb(ulu)%fimp) * tilelag(j)                  !! detached small aggregates -KDW
      else
        ob(icmd)%hd(5)%flo = qtile * cnv_m3               !! tile flow
        ob(icmd)%hd(5)%no3 = tileno3(j) * cnv_kg          !! tile flow nitrate 
        ob(icmd)%hd(5)%solp = tilesolp(j) * cnv_kg          !! tile flow soluble phos - KDW
        ob(icmd)%hd(5)%orgn = tileorgn(j) * cnv_kg           !! lat flow organic n -KDW
        ob(icmd)%hd(5)%orgp = tileorgp(j) * cnv_kg           !! lat flow organic p -KDW
        ob(icmd)%hd(5)%minpa = tileminpa(j) * cnv_kg           !! lat flow active mineral p -KDW
        ob(icmd)%hd(5)%minps = tileminps(j) * cnv_kg           !! lat flow stable mineral p -KDW
        ob(icmd)%hd(5)%chla = chl_a(j) * ob(icmd)%hd(5)%flo / 1000        !!chl_a, added subsurface chla, cbod, and dox, KDW 12/13/21
        ob(icmd)%hd(5)%cbod = cbodu(j) * ob(icmd)%hd(5)%flo / 1000        !!cbodu, added subsurface chla, cbod, and dox, KDW 12/13/21
        ob(icmd)%hd(5)%dox = doxq(j) * ob(icmd)%hd(5)%flo / 1000        !!dissolved oxygen, added subsurface chla, cbod, and dox, KDW 12/13/21
        ob(icmd)%hd(5)%sed = tilesed(j)                  !! detached sand -KDW
        ob(icmd)%hd(5)%san = tilesan(j)                  !! detached sand -KDW
        ob(icmd)%hd(5)%sil = tilesil(j)                  !! detached silt -KDW
        ob(icmd)%hd(5)%cla = tilecla(j)                  !! detached clay -KDW
        ob(icmd)%hd(5)%sag = tilesag(j)                  !! detached small aggregates -KDW
        ob(icmd)%hd(5)%lag = tilelag(j)                  !! detached small aggregates -KDW
      endif


      !set constituents
      do ipest = 1, cs_db%num_pests
        obcs(icmd)%hd(5)%pest(ipest) = hpestb_d(j)%pest(ipest)%tileq * cnv_kg
      end do
      do ipath = 1, cs_db%num_paths
        obcs(icmd)%hd(5)%path(ipath) = 0.
      end do
      
      !sum to obtain the total outflow hydrograph (1)
      ob(icmd)%hd(1) = hz
      do ihyd = 3, 5
        ob(icmd)%hd(1) = ob(icmd)%hd(1) + ob(icmd)%hd(ihyd)
      end do
      
      !set constituents
      do ipest = 1, cs_db%num_pests
        obcs(icmd)%hd(1)%pest(ipest) = obcs(icmd)%hd(3)%pest(ipest) + obcs(icmd)%hd(4)%pest(ipest) +    &
                                                                      obcs(icmd)%hd(5)%pest(ipest)
      end do
      do ipath = 1, cs_db%num_paths
        obcs(icmd)%hd(1)%path(ipath) = 0
      end do
      
      !! set subdaily hydrographs
      if (time%step > 0) then
      iday = ob(icmd)%day_cur
      iday_prev = iday - 1
      if (iday_prev < 1) iday_prev = 2
        
      !! subsurface flow = lateral + tile
      ssq = (ob(icmd)%hd(4)%flo + ob(icmd)%hd(5)%flo) * cnv_m3  / time%step
        
      !! zero previous days hyds - current day is the hyd from yesterday so its set
      do kk = 1, time%step
        ob(icmd)%ts(iday_prev,kk) = hz
      end do

      if (qday > 1.e-9) then
          
        !! use unit hydrograph to compute subdaily flow hydrographs
        sumflo = 0.  !sum flow in case hydrograph exceeds max days 
        
        do ii = 1, time%step !loop for total time steps in a day
          itot = ii
          do ib = 1, itb(j)  !loop for number of steps in the unit hydrograph base time
            itot = itot + ib - 1
            if (itot > time%step) then
              iday = iday + 1
              if (iday > ihdmx) iday = 1
              itot = 1
            end if

            !! check to see if day has gone past the max allocated days- uh > 1 day
            if (iday <= ihdmx) then
              ob(icmd)%ts(iday,itot)%flo = ob(icmd)%ts(iday,itot)%flo + hhsurfq(j,ii) * uh(j,ib) * cnv_m3
              sumflo = sumflo + ob(icmd)%ts(iday,itot)%flo
            else
              !! adjust if flow exceeded max days
              rto = Max (1., ob(icmd)%hd(3)%flo / sumflo)
              do iadj = 1, itot - 1
                iday = iadj / time%step + 1
                istep = iadj - (iday - 1) * time%step
                ob(icmd)%ts(iday,itot)%flo = ob(icmd)%ts(iday,itot)%flo * rto
              end do
            end if
          end do
        end do
        
        sumflo_day = 0.
        iday = ob(icmd)%day_cur
        do istep = 1, time%step
          ob(icmd)%ts(iday,istep)%flo = ob(icmd)%ts(iday,istep)%flo + ssq
          sumflo_day = sumflo_day + ob(icmd)%ts(iday,istep)%flo
        end do

        !! set values for other routing variables - assume constant concentration
        !! storage locations set to zero are not currently used
        do ii = 1, time%step
          ratio = ob(icmd)%ts(iday,ii)%flo / sumflo_day
            if (ob(icmd)%hd(1)%flo > 0.) then
              ob(icmd)%ts(iday,ii)%temp = wtmp                                !!wtmp
              ob(icmd)%ts(iday,ii)%sed = ob(icmd)%hd(1)%sed * ratio           !!sedyld
              ob(icmd)%ts(iday,ii)%orgn = ob(icmd)%hd(1)%orgn * ratio         !!sedorgn
              !ob(icmd)%ts(iday,ii)%sedp = ob(icmd)%hd(1)%sedp * ratio         !!sedorgp
              ob(icmd)%ts(iday,ii)%orgp = ob(icmd)%hd(1)%orgp * ratio         !!sedorgp - KDW
              ob(icmd)%ts(iday,ii)%minpa = ob(icmd)%hd(1)%minpa * ratio         !!sedminpa - KDW
              ob(icmd)%ts(iday,ii)%minps = ob(icmd)%hd(1)%minps * ratio         !!sedminps - KDW
              ob(icmd)%ts(iday,ii)%no3 = ob(icmd)%hd(1)%no3 * ratio           !!no3
              ob(icmd)%ts(iday,ii)%solp = ob(icmd)%hd(1)%solp * ratio         !!minp
              !ob(icmd)%ts(iday,ii)%psol = ob(icmd)%hd(1)%psol * ratio         !!sol pst
              !ob(icmd)%ts(iday,ii)%psor = ob(icmd)%hd(1)%psor * ratio         !!sorb pst
              ob(icmd)%ts(iday,ii)%chla = ob(icmd)%hd(1)%chla * ratio         !!chl_a
              ob(icmd)%ts(iday,ii)%nh3 = 0.                                   !! NH3
              ob(icmd)%ts(iday,ii)%no2 = 0.                                   !! NO2
              ob(icmd)%ts(iday,ii)%cbod = ob(icmd)%hd(1)%cbod * ratio         !!cbodu
              ob(icmd)%ts(iday,ii)%dox = ob(icmd)%hd(1)%dox * ratio           !!doxq & soxy 
            end if
          end do
        else
          !! no surface runoff on current day so zero hyds
          do istep = 1, time%step
            ob(icmd)%ts(iday,istep)%flo = ssq
          end do
        end if  ! qday > 0
      end if  ! time%step  > 0

      return   
      end subroutine hru_hyds