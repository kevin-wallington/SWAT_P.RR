      subroutine nut_psed

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates the amount of organic and mineral phosphorus
!!    attached to sediment in surface runoff

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    enratio       |none         |enrichment ratio calculated for day in HRU
!!    erorgp(:)     |none         |organic P enrichment ratio, if left blank
!!                                |the model will calculate for every event
!!    ihru          |none         |HRU number
!!    sedyld(:)     |metric tons  |daily soil loss caused by water erosion in
!!                                |HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sedminpa(:)  |kg P/ha       |amount of active mineral phosphorus sorbed to
!!                                |sediment in surface runoff in HRU for day
!!    sedminps(:)  |kg P/ha       |amount of stable mineral phosphorus sorbed to
!!                                |sediment in surface runoff in HRU for day
!!    sedorgp(:)   |kg P/ha       |amount of organic phosphorus in surface
!!                                |runoff in HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

        use hru_module, only : hru, sedyld, sedorgp, sedminpa, sedminps, ihru, enratio,  &
          ihru, ipl 
        use soil_module
        use plant_module
        use organic_mineral_mass_module
      
        implicit none       

        integer :: j                !none           |HRU number
        integer :: sb               !none           |subbasin number
        real :: sedp_attach         !kg P/ha        |amount of phosphorus attached to sediment 
                                    !               |in soil
        real :: wt1                 !none           |conversion factor (mg/kg => kg/ha)
        real :: er                  !none           |enrichment ratio
        real :: conc                !               |concentration of organic N in soil
        real :: sedp                !kg P/ha        |total amount of P removed in sediment erosion 
        real :: sed_orgp            !kg P/ha        |total amount of P in organic pools
        real :: sed_hump            !kg P/ha        |amount of P in humus pool
        real :: sed_microbp            !kg P/ha        |amount of P in manure soil pool
        real :: sed_rsdp       !kg P/ha        |maount of P in residue manure pool
        real :: fr_orgp             !kg P/ha        |fraction of organic phosphorus in soil (humus + manure in soil + manure in residue)
        real :: fr_actmin           !kg P/ha        |fraction of active mineral phosphorus in soil
        real :: fr_stamin           !kg P/ha        |fraction of stable mineral phosphorus in soil
        real :: rsd                 !kg P/ha in plant residue - KDW
        real :: rsd_perc

        j = ihru
        rsd = 0
        !! HRU sediment calculations
        do ipl = 1, pcom(j)%npl ! - KDW
          rsd = rsd + rsd1(j)%tot(ipl)%p ! - KDW
        end do ! - KDW
        !sedp_attach = soil1(j)%hp(1)%p + soil1(j)%man(1)%p + rsd1(j)%man%p + soil1(j)%mp(1)%sta + soil1(j)%mp(1)%act !  man%p only exists when cswat=1, otherwise need tot%p - KDW 
        ! above, rsd%man doesn't exist but plant residue does - KDW
        sedp_attach = soil1(j)%hp(1)%p + soil1(j)%microb(1)%p + rsd + soil1(j)%mp(1)%sta + soil1(j)%mp(1)%act ! - KDW , add soil1(j)%man(1)%p if want to include cswat=1 application
        if (sedp_attach > 1.e-6) then !KDW 03/29, was 1.e-9
          fr_orgp = (soil1(j)%hp(1)%p + soil1(j)%microb(1)%p  + rsd) / sedp_attach !kDW 10/21/21 changed soil1(j)tot(1)%p to soil1(j)microb(1)p (also in line 64)
          fr_actmin = soil1(j)%mp(1)%act / sedp_attach ! KDW corrected 02/18/22, had flipped active and stable on right hand side
          fr_stamin = soil1(j)%mp(1)%sta / sedp_attach
        else !KDW 03/29
          sedp_attach = 0
          fr_orgp = 0
          fr_actmin = 0
          fr_stamin = 0
        end if

        wt1 = soil(j)%phys(1)%bd * soil(j)%phys(1)%d / 100.

        if (hru(j)%hyd%erorgp > .001) then
          er = hru(j)%hyd%erorgp
        else
          er = enratio
        end if
      
        conc = sedp_attach * er / wt1
        sedp = .001 * conc * sedyld(j) / hru(j)%area_ha
        
        if (sedp > 1.e-9) then
          sedorgp(j) = sedp * fr_orgp
          sedminpa(j) = sedp * fr_actmin
          sedminpa(j) = amin1 (sedminpa(j), soil1(j)%mp(1)%act)
          soil1(j)%mp(1)%act = soil1(j)%mp(1)%act - sedminpa(j)
          if (soil1(j)%mp(1)%act < 1.e-6) soil1(j)%mp(1)%act = 0. ! KDW 03/31
          sedminps(j) = sedp * fr_stamin
          sedminps(j) = amin1 (sedminps(j), soil1(j)%mp(1)%sta)
          soil1(j)%mp(1)%sta = soil1(j)%mp(1)%sta - sedminps(j)
          if (soil1(j)%mp(1)%sta < 1.e-6) soil1(j)%mp(1)%sta = 0. ! KDW 03/31
        
          ! sed_orgp = soil1(j)%hp(1)%p + soil1(j)%man(1)%p  + rsd1(j)%man%p 
          sed_orgp = soil1(j)%hp(1)%p + soil1(j)%microb(1)%p  + rsd ! KDW 10/21/21, replaces above
          if (sed_orgp > 1.e-6) then
            sed_hump = sedorgp(j) * (soil1(j)%hp(1)%p / sed_orgp)
            sed_hump = amin1 (sed_hump, soil1(j)%hp(1)%p)
            soil1(j)%hp(1)%p = soil1(j)%hp(1)%p - sed_hump
            if (soil1(j)%hp(1)%p < 1.e-6) soil1(j)%hp(1)%p = 0. ! KDW 03/31
        
            sed_microbp = sedorgp(j) * (soil1(j)%microb(1)%p / sed_orgp) ! KDW 10/21/21, replaced all instances of "tot" with "microb" in this stanza
            sed_microbp = amin1 (sed_microbp, soil1(j)%microb(1)%p)
            soil1(j)%microb(1)%p = soil1(j)%microb(1)%p - sed_microbp
            if (soil1(j)%microb(1)%p < 1.e-6) soil1(j)%microb(1)%p = 0. ! KDW 03/31, 10/21/21 chaanged soil1(j)%man(1)%p to soil1(j)%microb(1)%p
        
            sed_rsdp = sedorgp(j) * (rsd / sed_orgp)
            sed_rsdp = amin1 (sed_rsdp, rsd)
            do ipl = 1, pcom(j)%npl ! - KDW
              if (rsd > 1.e-6) then ! KDW 04/06
                rsd_perc = rsd1(j)%tot(ipl)%p / rsd ! - KDW
              else! KDW 04/06
                rsd_perc = 0! KDW 04/06
              endif
              rsd1(j)%tot(ipl)%p = rsd1(j)%tot(ipl)%p - rsd_perc * sed_rsdp
              if (rsd1(j)%tot(ipl)%p < 1.e-6) rsd1(j)%tot(ipl)%p = 0. ! KDW 03/31
            end do ! - KDW
          end if
        else
          sedorgp(j) = 0.
          sedminpa(j) = 0.
          sedminps(j) = 0.
        end if 

      return
      end subroutine nut_psed