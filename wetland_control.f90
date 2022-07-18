      subroutine wetland_control
    
      use reservoir_data_module
      use reservoir_module
      use hru_module, only : hru, sedyld, sanyld, silyld, clayld, sagyld, lagyld, grayld, sedminps, sedminpa,   &
        surqno3, sedorgn, sedorgp, qdr, ihru, pet_day, qday, precipday, surqsolp
      use conditional_module
      use climate_module
      use hydrograph_module
      use time_module
      use basin_module
      use channel_module
      use water_body_module
      
      implicit none
     
      real :: bypass                  !              | 
      real :: fracwet                 !              | 
      integer :: j                    !none          |counter
      integer :: iprop                !              |  
      integer :: iac                  !none          |counter
      character(len=1) :: action           !         |
      integer :: ial                  !none          |counter
      real :: b_lo                    !              |
      real :: res_h                   !              |
      real :: x1                      !              |
      real :: wet_h                   !              |
      real :: wet_h1                  !              |
      real :: flwi                    !m^3 H2O       |water entering pothole on day  
      real :: flwo                    !              |
      real :: sedi                    !metric tons   |sediment entering pothole on day
      real :: sedo                    !metric tons   |sed leaving res 
      integer :: k                    !              | 
      integer :: ii                   !none          |counter 
      integer :: jres                 !none          |reservoir number
      integer :: idat                 !              |
      integer :: ihyd                 !none          |counter
      integer :: ised                 !none          |counter
      integer :: irel                 !              |
      integer :: inut                 !none          |counter
      integer :: ipst                 !none          |counter
      integer :: irch                 !none          |counter
      integer :: ires = 0
      integer :: iob = 0
      real :: wet_fr = 0.
      real :: pvol_m3
      real :: evol_m3
      integer :: ch_par = 0
      real :: mxsa                    !m^2              ! maximum surface area, i.e. esa - KDW
      integer :: iswet                !binary         |is this a wetland (1) or a reservoir (0), KDW 03/30

      iswet = 1 ! KDW 03/30

      j = ihru
      ires= hru(j)%dbs%surf_stor
      ihyd = wet_dat(ires)%hyd
      ised = wet_dat(ires)%sed
      irel = wet_dat(ires)%release
      iob = hru(j)%obj_no ! - KDW
      hru(j)%water_fr = 0.
      mxsa = wet_ob(ihru)%esa * 10000  ! - KDW

      !! initialize variables for reservoir daily simulation
      hru(ihru)%water_seep = 0.

      bypass = 1. - wet_hyd(ihyd)%frac
      fracwet = 1. - bypass 
      fracwet = max (fracwet,0.)

      !! set incoming flow, sediment and nutrients
      ht1%flo = fracwet * qday * 10 * hru(ihru)%area_ha ! all mutliplied by fracwet * 10  * area_ha by KDW, 06/02
      ht1%sed = fracwet * sedyld(ihru)! all sediment yields mutliplied by fracwetby KDW; likewise for constituents below, 03/30
      ht1%san = fracwet * sanyld(ihru) 
      ht1%sil = fracwet * silyld(ihru)
	    ht1%cla = fracwet * clayld(ihru) 
	    ht1%sag = fracwet * sagyld(ihru)
	    ht1%lag = fracwet * lagyld(ihru)
	    ht1%grv = fracwet * grayld(ihru)
      ht1%orgn = fracwet * sedorgn(ihru) * hru(ihru)%area_ha ! all mutliplied by fracwet * area_ha by KDW; likewise for constituents below, 06/02
      !ht1%sedp = sedorgp(ihru)
      ht1%orgp = fracwet * sedorgp(ihru) * hru(ihru)%area_ha ! - KDW
      ht1%no3 = fracwet * surqno3(ihru) * hru(ihru)%area_ha
      ht1%nh3 = 0. 
      ht1%no2 = 0.
      ht1%cbod = 0. ! - KDW
      ht1%dox = 0. ! - KDW
      !ht1%solp = sedminps(ihru) + sedminpa(ihru)
      ht1%solp = fracwet * surqsolp(ihru) * hru(ihru)%area_ha ! - KDW
      ht1%minpa = fracwet * sedminpa(ihru) * hru(ihru)%area_ha ! - KDW
      ht1%minps = fracwet * sedminps(ihru) * hru(ihru)%area_ha ! - KDW

      !! add inflow - KDW
      if (bsn_cc%wq == 1 .or. bsn_cc%wq == 0) then
          wet(ihru) = wet(ihru) + ht1 ! all inflows
      else
          wet(ihru)%flo =  wet(ihru)%flo + ht1%flo ! flow only - KDW
      endif

      !! add precipitation - mm*ha*10.=m3 (used same area for infiltration and soil evap)
      wet_wat_d(ihru)%precip = precipday * wet_wat_d(ihru)%area_ha * 10.
      wet(ihru)%flo =  wet(ihru)%flo + wet_wat_d(ihru)%precip
      
      !! subtract evaporation and seepage - mm*ha*10.=m3
      wet_wat_d(ihru)%evap = pet_day * wet_hyd(ihyd)%evrsv * wet_wat_d(ihru)%area_ha * 10.
      wet_wat_d(ihru)%evap = min(wet_wat_d(ihru)%evap, wet(ihru)%flo)
      wet(ihru)%flo =  wet(ihru)%flo - wet_wat_d(ihru)%evap
      hru(ihru)%water_evap = wet_wat_d(ihru)%evap / (10. * hru(ihru)%area_ha)
        
      !! save hru(ihru)%water_seep to add to infiltration on next day
      wet_wat_d(ihru)%seep = wet_wat_d(ihru)%area_ha * wet_hyd(ihyd)%k * 10.* 24.
      wet_wat_d(ihru)%seep = min(wet(ihru)%flo, wet_wat_d(ihru)%seep)
      wet(ihru)%flo = wet(ihru)%flo - hru(ihru)%water_seep
      hru(ihru)%water_seep = wet_wat_d(ihru)%seep / (10. * hru(ihru)%area_ha)
        
      !! calc release from decision table
      d_tbl => dtbl_res(irel)
      wbody => wet(ihru)
      wbody_wb => wet_wat_d(ihru)
      wbody_bed => wet_bed(ihru)
      pvol_m3 = wet_ob(ihru)%pvol
      evol_m3 = wet_ob(ihru)%evol
      call conditions (ihru, irel)
      call res_hydro (ihru, irel, ihyd, pvol_m3, evol_m3)

      !! update surface area - solve quadratic to find new depth - copied to preceeed nutrient routine which uses area (still re-calc below though) - KDW 10/21/15
      wet_wat_d(ihru)%area_ha = 0.
      if (wet(ihru)%flo > 0.) then
        x1 = wet_hyd(ihyd)%bcoef ** 2 + 4. * wet_hyd(ihyd)%ccoef * (1. - wet(ihru)%flo / wet_ob(ihru)%pvol)
        if (x1 < 1.e-6) then
          wet_h = 0.
        else
          wet_h1 = (-wet_hyd(ihyd)%bcoef - sqrt(x1)) / (2. * wet_hyd(ihyd)%ccoef)
          wet_h = wet_h1 + wet_hyd(ihyd)%bcoef
        end if
        wet_fr = (1. + wet_hyd(ihyd)%acoef * wet_h)
        wet_fr = min(wet_fr,1.)
        wet_wat_d(ihru)%area_ha = hru(ihru)%area_ha * wet_fr
                
        hru(ihru)%water_fr =  wet_wat_d(ihru)%area_ha / hru(ihru)%area_ha

      end if 

      if (bsn_cc%wq == 1 .or. bsn_cc%wq == 0) then
        call res_sediment (ihru, ihyd, ised)
        !! below lines moved up from lines ~125 - KDW
        inut = wet_dat(ires)%nut 
        ch_par = res_nut(inut)%ch_params !this indicates which *channel* nutcha to use - KDW
        call res_nutrient (ires, inut, ihru) !KDW 08/18        
        ipst = wet_dat(ires)%pst
        !call res_pest (ires)
        !! above lines moved up from lines ~125 - KDW
      endif
      if (bsn_cc%wq == 2) then
        call res_trap (ihru, ihyd, ised, mxsa) ! - KDW
        inut = wet_dat(ires)%nut 
        ch_par = res_nut(inut)%ch_params !this indicates which *channel* nutcha to use - KDW
        call res_nutrient2 (ihru, ch_par, iob, mxsa, iswet) ! first input changed from ires, third changed from ihru - KDW 03/30
        ipst = wet_dat(ires)%pst
        !call res_pest (ires)
        call res_sediment2 (ihru, ihyd, ised, mxsa)
      endif

      !! subtract outflow from storage
      !wet(ihru)%flo =  wet(ihru)%flo - ht2%flo
      !subtract all nutrients, etc. not just flo, replaced above, KDW 10/15/21
      wet(ihru)%flo = wet(ihru)%flo - ht2%flo
      wet(ihru)%sed = wet(ihru)%sed - ht2%sed        
      wet(ihru)%orgn = wet(ihru)%orgn - ht2%orgn        
      wet(ihru)%orgp = wet(ihru)%orgp - ht2%orgp
      wet(ihru)%minpa = wet(ihru)%minpa - ht2%minpa
      wet(ihru)%minps = wet(ihru)%minps - ht2%minps
      wet(ihru)%no3 = wet(ihru)%no3 - ht2%no3
      wet(ihru)%solp = wet(ihru)%solp - ht2%solp
      wet(ihru)%chla = wet(ihru)%chla - ht2%chla
      wet(ihru)%nh3 = wet(ihru)%nh3 - ht2%nh3
      wet(ihru)%no2 = wet(ihru)%no2 - ht2%no2
      wet(ihru)%cbod = wet(ihru)%cbod - ht2%cbod
      wet(ihru)%dox = wet(ihru)%dox - ht2%dox
      wet(ihru)%san = wet(ihru)%san - ht2%san
      wet(ihru)%sil = wet(ihru)%sil - ht2%sil
      wet(ihru)%cla = wet(ihru)%cla - ht2%cla
      wet(ihru)%sag = wet(ihru)%sag - ht2%sag
      wet(ihru)%lag = wet(ihru)%lag - ht2%lag
      wet(ihru)%grv = wet(ihru)%grv - ht2%grv



      !! update surface area - solve quadratic to find new depth
      wet_wat_d(ihru)%area_ha = 0.
      if (wet(ihru)%flo > 0.) then
        x1 = wet_hyd(ihyd)%bcoef ** 2 + 4. * wet_hyd(ihyd)%ccoef * (1. - wet(ihru)%flo / wet_ob(ihru)%pvol)
        if (x1 < 1.e-6) then
          wet_h = 0.
        else
          wet_h1 = (-wet_hyd(ihyd)%bcoef - sqrt(x1)) / (2. * wet_hyd(ihyd)%ccoef)
          wet_h = wet_h1 + wet_hyd(ihyd)%bcoef
        end if
        wet_fr = (1. + wet_hyd(ihyd)%acoef * wet_h)
        wet_fr = min(wet_fr,1.)
        wet_wat_d(ihru)%area_ha = hru(ihru)%area_ha * wet_fr
                
        hru(ihru)%water_fr =  wet_wat_d(ihru)%area_ha / hru(ihru)%area_ha

      end if 
 
      !! subtract sediment leaving from reservoir
      !wet(ihru)%sed = wet(ihru)%sed - ht2%sed
      !wet(ihru)%sil = wet(ihru)%sil - ht2%sil
      !wet(ihru)%cla = wet(ihru)%cla - ht2%cla
          

      !! set values for routing variables
      ob(icmd)%hd(1)%temp = 0.                  !!undefined

      qdr(ihru) = ht2%flo / (10. * hru(ihru)%area_ha) + ht1%flo * bypass
      sedyld(ihru) = ht2%sed / hru(ihru)%area_ha + sedyld(ihru) * bypass
      sanyld(ihru) = ht2%san / hru(ihru)%area_ha + sanyld(ihru) * bypass
      silyld(ihru) = ht2%sil / hru(ihru)%area_ha + silyld(ihru) * bypass
	    clayld(ihru) = ht2%cla / hru(ihru)%area_ha + clayld(ihru) * bypass 
	    sagyld(ihru) = ht2%sag / hru(ihru)%area_ha + sagyld(ihru) * bypass
	    lagyld(ihru) = ht2%lag / hru(ihru)%area_ha + lagyld(ihru) * bypass
	    grayld(ihru) = ht2%grv / hru(ihru)%area_ha + grayld(ihru) * bypass
      sedorgn(ihru) = ht2%orgn / hru(ihru)%area_ha + sedorgn(ihru) * bypass
      !sedorgp(ihru) = ht2%sedp / hru(ihru)%area_ha + sedorgp(ihru) * bypass
      sedorgp(ihru) = ht2%orgp / hru(ihru)%area_ha + sedorgp(ihru) * bypass ! - KDW
      surqno3(ihru) = ht2%no3/ hru(ihru)%area_ha  + surqno3(ihru) * bypass
      !nh3 = resnh3o + 0.  !add ammonium 
      !no2  = resno2o + 0.  !add no2
      !sedminps(ihru) = ht2%solp / hru(ihru)%area_ha / 2. + sedminps(ihru) * bypass
      !sedminpa(ihru) = ht2%solp / hru(ihru)%area_ha / 2. + sedminpa(ihru) * bypass
      sedminps(ihru) = ht2%minps / hru(ihru)%area_ha + sedminps(ihru) * bypass ! - KDW
      sedminpa(ihru) = ht2%minpa / hru(ihru)%area_ha + sedminpa(ihru) * bypass ! - KDW
      surqsolp(ihru) = ht2%solp / hru(ihru)%area_ha  + surqsolp(ihru) * bypass
      ! surqsolp not added because routine only deals with sediment attached nutrients
      
      !! set inflow and outflow variables for reservoir_output
      if (time%yrs > pco%nyskip) then
        wet_in_d(ihru) = ht1 
        wet_out_d(ihru) = ht2
        !wet_in_d(ihru)%flo = wet(ihru)%flo / 10000.   !m^3 -> ha-m
        !wet_out_d(ihru)%flo = wet(ihru)%flo / 10000.  !m^3 -> ha-m
      end if  
      
      return
      end subroutine wetland_control