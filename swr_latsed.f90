      subroutine swr_latsed

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates the sediment load contributed in lateral flow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hru_km(:)   |km^2          |area of HRU in square kilometers
!!    lat_sed(:)  |g/L           |sediment concentration in lateral flow
!!    latq(:)     |mm H2O        |total lateral flow in soil profile for the
!!                               |day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use hru_module, only : hru, ihru, latsed, latsan, latsil, latcla, latsag, latlag,  &
         tilesed, tilesan, tilesil, tilecla, latorgp, tileorgp, latorgn, tileorgn, latq, &
         latminpa, tileminpa, latminps, tileminps, tilesag, tilelag, qtile, enratio
      use soil_module
      use organic_mineral_mass_module
      use basin_module
      
      implicit none

      integer :: j              !none          |HRU number
      integer :: jj             !none          |counter for soil layer loop
      real :: sed               !metric tons   |sediment load in lateral flow
      real :: lyrorgp              !kg P          |organic P attached to sediment
      real :: lyrorgn              !kg P          |organic N attached to sediment
      real :: lyrminpa             !kg P          |active mineral P attached to sediment
      real :: lyrminps             !kg P          |stable mineral P attached to sediment 
      real :: lyrsedp              !kg P          |for old version (sol_p model = 0 or 1), to maintain mass balance in layers, KDW 08/13
      real :: lyrsedn              !kg P          |for old version (sol_p model = 0 or 1), to maintain mass balance in layers, KDW 08/13
      real :: fracminpa            !kg P          |for old version (sol_p model = 0 or 1), to maintain mass balance in layers, KDW 08/13
      real :: fracminps            !kg P          |for old version (sol_p model = 0 or 1), to maintain mass balance in layers, KDW 08/13
      real :: totsedp              !kg P          |for old version (sol_p model = 0 or 1), to maintain mass balance in layers, KDW 08/13
      real :: totorgp           !kg P/ha       |organic P in soil layer, all forms
      real :: totorgn          !kg P/ha       |organic N in soil layer, all forms
      real :: frachp            !kg P/ha       |humus N/P in soil layer
      real :: frachs            !kg P/ha       |humis N in soil layer
      real :: fracmicrobp          !kg P/ha       |"total" fresh P in soil layer - changed from fractotp to fracmicrobp by KDW 10/21/21
      real :: er                ! --            | enrichment ratio for sediment-associated P - KDW 10/16/21
      real :: clay, silt, sag !! KDW temporary 01/30/22

      j = 0
      j = ihru

      !! KDW temporary 01/30/22 no large sediment in tile/lateral flow !!
      clay = soil(j)%det_cla /(soil(j)%det_cla + soil(j)%det_sil + soil(j)%det_sag)
      silt = soil(j)%det_sil /(soil(j)%det_cla + soil(j)%det_sil + soil(j)%det_sag)
      sag = soil(j)%det_sag /(soil(j)%det_cla + soil(j)%det_sil + soil(j)%det_sag)
      latsil(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * silt
      latcla(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * clay
      latsag(j) = latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * sag
      tilesil(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed * silt
      tilecla(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed * clay
      tilesag(j) = qtile * hru(j)%km * hru(j)%hyd%tile_sed * sag
      latsan(j) = 0.
      latlag(j) = 0.
      tilesan(j) = 0.
      tilelag(j) = 0.
      !! End KDW no large sediment !!
      !! update sediment yield for sediment in lateral flow
!     sedyld(j)=sedyld(j) + latq(j) * hru(j)%km * hru(j)%hyd%lat_sed
!     sanyld(j)=sanyld(j)+latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_san
!     silyld(j)=silyld(j)+latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_sil
!     clayld(j)=clayld(j)+latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_cla
!     sagyld(j)=sagyld(j)+latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_sag
!     lagyld(j)=lagyld(j)+latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_lag
      !! update to seperate surface and subsurface sediment
      latsed(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed ! mm * 10^6 m^2 * g/L - comment by KDW
      ! latsan(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_san
      ! latsil(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_sil
      ! latcla(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_cla
      ! latsag(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_sag
      ! latlag(j)=latq(j) * hru(j)%km * hru(j)%hyd%lat_sed * soil(j)%det_lag
      tilesed(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed ! KDW, 06/23, changed from lat_sed to tile_sed
      ! tilesan(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed * soil(j)%det_san
      ! tilesil(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed * soil(j)%det_sil
      ! tilecla(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed * soil(j)%det_cla
      ! tilesag(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed * soil(j)%det_sag
      ! tilelag(j)=qtile * hru(j)%km * hru(j)%hyd%tile_sed * soil(j)%det_lag
      if (latsed(j) < 1.e-9) latsed(j) = 0.0 ! KDW, 03/26
      if (latsan(j) < 1.e-9) latsan(j) = 0.0 ! KDW, 03/26
      if (latsil(j) < 1.e-9) latsil(j) = 0.0 ! KDW, 03/26
      if (latcla(j) < 1.e-9) latcla(j) = 0.0 ! KDW, 03/26
      if (latsag(j) < 1.e-9) latsag(j) = 0.0 ! KDW, 03/26
      if (latlag(j) < 1.e-9) latlag(j) = 0.0 ! KDW, 03/26
      if (tilesed(j) < 1.e-9) tilesed(j) = 0.0 ! KDW, 03/26
      if (tilesan(j) < 1.e-9) tilesan(j) = 0.0 ! KDW, 03/26
      if (tilesil(j) < 1.e-9) tilesil(j) = 0.0 ! KDW, 03/26
      if (tilecla(j) < 1.e-9)  tilecla(j) = 0.0 ! KDW, 03/26
      if (tilesag(j) < 1.e-9)  tilesag(j) = 0.0 ! KDW, 03/26
      if (tilelag(j) < 1.e-9)  tilelag(j) = 0.0 ! KDW, 03/26


      !! organic n and p in the lateral flow     - by J.Jeong BREC 2011  
      !mm*mg/L*1000L/m3*kg/1000000mg*m3/10ha-mm=kg/ha
      !sedorgn(j) = sedorgn(j) + latq(j) * hru(j)%hyd%lat_orgn /10000.
      !!sedorgp(j) = sedorgp(j) + latq(j) * hru(j)%hyd%lat_orgp /10000. - eliminated by KDW

      !! Below added by KDW

      lyrorgp = 0
      lyrorgn = 0
      lyrminpa = 0
      lyrminps = 0
      lyrsedp = 0
      lyrsedn = 0
      totsedp = 0
      tileorgp(j) = 0
      tileorgn(j) = 0
      tileminpa(j) = 0
      tileminps(j) = 0
      latorgp(j) = 0
      latminpa(j) = 0
      latminps(j) = 0
      latorgn(j) = 0

      if (hru(j)%hyd%erorgp > .001) then
            er = hru(j)%hyd%erorgp
          else
            er = enratio
      end if

   if (bsn_cc%sol_P_model == 0 .or. bsn_cc%sol_P_model == 1) then  
      do jj = 1, soil(j)%nly
            totsedp = soil1(j)%hp(jj)%p + soil1(j)%man(jj)%p + soil1(j)%mp(jj)%act + soil1(j)%mp(jj)%sta + 0.001
            lyrsedp = soil(j)%ly(jj)%flat * hru(j)%hyd%lat_orgp /10000.
            if (lyrsedp < 1.e-6)  lyrsedp = 0.0 ! KDW, 03/26  
            frachp = soil1(j)%hp(jj)%p / totsedp
            fracmicrobp = soil1(j)%microb(jj)%p / totsedp ! KDW 10/21/21, changed tot(l)% to microb(l)%
            fracminpa = soil1(j)%mp(jj)%act / totsedp
            fracminpa = soil1(j)%mp(jj)%sta / totsedp
            latorgp(j) = latorgp(j) + lyrsedp * (frachp + fracmicrobp)
            soil1(j)%hp(jj)%p = soil1(j)%hp(jj)%p - lyrsedp * frachp
            soil1(j)%microb(jj)%p = soil1(j)%microb(jj)%p - lyrsedp * fracmicrobp ! KDW 10/21/21, changed tot(l)% to microb(l)%
            latminpa(j) = latminpa(j) + lyrsedp * fracminpa
            soil1(j)%mp(jj)%act = soil1(j)%mp(jj)%act - lyrsedp * fracminpa
            latminps(j) = latminps(j) + lyrsedp * fracminps
            soil1(j)%mp(jj)%sta = soil1(j)%mp(jj)%sta - lyrsedp * fracminps
     
            totorgn = soil1(j)%hp(jj)%n + soil1(jj)%hs(1)%n + 0.001
            lyrorgn = soil(j)%ly(jj)%flat * hru(j)%hyd%lat_orgn /10000.
            if (lyrsedn < 1.e-6)  lyrorgn = 0.0 ! KDW, 03/26
            latorgn(j) = latorgn(j) + lyrorgn  
            frachp = soil1(j)%hp(jj)%n / totorgn
            frachs = soil1(j)%hs(jj)%n / totorgn
            soil1(j)%hp(jj)%n = soil1(j)%hp(jj)%n - lyrorgn * frachp
            soil1(j)%hs(jj)%n = soil1(j)%hs(jj)%n - lyrorgn * frachs
     
            if (hru(j)%lumv%ldrain == jj) then
                  lyrsedp = qtile * hru(j)%hyd%lat_orgp /10000.
                  if (lyrsedp < 1.e-6)  lyrsedp = 0.0 ! KDW, 03/26  
                  latorgp(j) = latorgp(j) + lyrsedp * (frachp + fracmicrobp) ! KDW 10/21/21, changed tot(l)% to microb(l)%
                  soil1(j)%hp(jj)%p = soil1(j)%hp(jj)%p - lyrsedp * frachp
                  soil1(j)%microb(jj)%p = soil1(j)%microb(jj)%p - lyrsedp * fracmicrobp ! KDW 10/21/21, changed tot(l)% to microb(l)%
                  latminpa(j) = latminpa(j) + lyrsedp * fracminpa
                  soil1(j)%mp(jj)%act = soil1(j)%mp(jj)%act - lyrsedp * fracminpa
                  latminps(j) = latminps(j) + lyrsedp * fracminps
                  soil1(j)%mp(jj)%sta = soil1(j)%mp(jj)%sta - lyrsedp * fracminps
                  lyrorgn = qtile * hru(j)%hyd%lat_orgn /10000.
                  if (lyrorgn < 1.e-6)  lyrorgn = 0.0 ! KDW, 03/26
                  latorgn(j) = latorgn(j) + lyrorgn  
                  soil1(j)%hp(jj)%n = soil1(j)%hp(jj)%n - lyrorgn * frachp
                  soil1(j)%hs(jj)%n = soil1(j)%hs(jj)%n - lyrorgn * frachs
            end if   
            if (soil1(j)%hp(jj)%p < 1.e-6) soil1(j)%hp(jj)%p = 0.0 ! Changed all below from < 0. - KDW 03/24, moved inside do loop, 04/06
            if (soil1(j)%microb(jj)%p < 1.e-6) soil1(j)%microb(jj)%p = 0.0 ! KDW 10/21/21, changed tot(l)% to microb(l)%
            if (soil1(j)%hp(jj)%n < 1.e-6) soil1(j)%hp(jj)%n = 0.0
            if (soil1(j)%hs(jj)%n < 1.e-6) soil1(j)%hs(jj)%n = 0.0
            if (soil1(j)%mp(jj)%act < 1.e-6) soil1(j)%mp(jj)%act = 0.0
            if (soil1(j)%mp(jj)%sta < 1.e-6) soil1(j)%mp(jj)%sta = 0.0   
           end do




  else !SWAT+P calculates nutrients lost from each layer based on concentration in layer, rather than parameters lat_orgp and lat_orgn, KDW
      do jj = 1, soil(j)%nly
       sed = soil(j)%ly(jj)%flat * hru(j)%km * hru(j)%hyd%lat_sed
       if (sed < 1.e-6)  sed = 0.0 ! KDW, 03/26
       totorgp = soil1(j)%hp(jj)%p + soil1(j)%microb(jj)%p + 0.001 ! KDW 10/21/21 changed man(jj)%p to microb(jj)%p
       lyrorgp = er *totorgp * sed / (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * hru(j)%area_ha)    ! added er * , the enrichment ratio, KDW 10/16/21
       if (lyrorgp < 1.e-6)  lyrorgp = 0.0 ! KDW, 03/26
       latorgp(j) = latorgp(j) + lyrorgp  
       frachp = soil1(j)%hp(jj)%p / totorgp
       fracmicrobp = soil1(j)%microb(jj)%p / totorgp
       soil1(j)%hp(jj)%p = soil1(j)%hp(jj)%p - lyrorgp * frachp
       soil1(j)%microb(jj)%p = soil1(j)%microb(jj)%p - lyrorgp * fracmicrobp

       totorgn = soil1(j)%hp(jj)%n + soil1(jj)%hs(1)%n + 0.001
       lyrorgn = totorgn * sed / (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * hru(j)%area_ha)
       if (lyrorgn < 1.e-6)  lyrorgn = 0.0 ! KDW, 03/26
       latorgn(j) = latorgn(j) + lyrorgn  
       frachp = soil1(j)%hp(jj)%n / totorgn
       frachs = soil1(j)%hs(jj)%n / totorgn
       soil1(j)%hp(jj)%n = soil1(j)%hp(jj)%n - lyrorgn * frachp
       soil1(j)%hs(jj)%n = soil1(j)%hs(jj)%n - lyrorgn * frachs

       lyrminpa = er * soil1(j)%mp(jj)%act * sed /                                      &
                   (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * hru(j)%area_ha) ! added er * , the enrichment ratio, KDW 10/16/21
       if (lyrminpa < 1.e-6)  lyrminpa = 0.0 ! KDW, 03/26
       latminpa(j) = latminpa(j) + lyrminpa  
       soil1(j)%mp(jj)%act = soil1(j)%mp(jj)%act - lyrminpa
       lyrminps = er * soil1(j)%mp(jj)%sta * sed /                                      &
                   (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * hru(j)%area_ha) ! added er * , the enrichment ratio, KDW 10/16/21
       if (lyrminps < 1.e-6)  lyrminps = 0.0 ! KDW, 03/26
       latminps(j) = latminps(j) + lyrminps 
       soil1(j)%mp(jj)%sta = soil1(j)%mp(jj)%sta - lyrminps
       if (hru(j)%lumv%ldrain == jj) then
            sed = qtile * hru(j)%km * hru(j)%hyd%tile_sed ! KDW 06/23, was lat_sed
            if (sed < 1.e-6)  sed = 0.0 ! KDW, 03/26
            tileorgp(j) = er * totorgp * sed / (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * & ! added er * , the enrichment ratio, KDW 10/16/21
                  hru(j)%area_ha) ! KDW, 03/26 - from lyrorgp to totorgp 
            soil1(j)%hp(jj)%p = soil1(j)%hp(jj)%p - tileorgp(j) * frachp
            soil1(j)%microb(jj)%p = soil1(j)%microb(jj)%p - tileorgp(j) * fracmicrobp
            tileorgn(j) = totorgn * sed / (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * & 
                  hru(j)%area_ha)! KDW, 03/26 - from lyrorgn to totorgn
            soil1(j)%hp(jj)%n = soil1(j)%hp(jj)%n - tileorgn(j) * frachp
            soil1(j)%hs(jj)%n = soil1(j)%hs(jj)%n - tileorgn(j) * frachs
            tileminpa(j) = er * soil1(j)%mp(jj)%act * sed /                                      &
                        (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * hru(j)%area_ha) ! added er * , the enrichment ratio, KDW 10/16/21
            soil1(j)%mp(jj)%act = soil1(j)%mp(jj)%act - tileminpa(j)
            tileminps(j) = er * soil1(j)%mp(jj)%sta * sed /                                      &
                        (10*soil(j)%phys(jj)%bd * soil(j)%phys(jj)%thick * hru(j)%area_ha) ! added er * , the enrichment ratio, KDW 10/16/21
            soil1(j)%mp(jj)%sta = soil1(j)%mp(jj)%sta - tileminps(j)
       end if   
       if (soil1(j)%hp(jj)%p < 1.e-6) soil1(j)%hp(jj)%p = 0.0 ! Changed all below from < 0. - KDW 03/24, moved inside do loop, 04/06
       if (soil1(j)%microb(jj)%p < 1.e-6) soil1(j)%microb(jj)%p = 0.0
       if (soil1(j)%hp(jj)%n < 1.e-6) soil1(j)%hp(jj)%n = 0.0
       if (soil1(j)%hs(jj)%n < 1.e-6) soil1(j)%hs(jj)%n = 0.0
       if (soil1(j)%mp(jj)%act < 1.e-6) soil1(j)%mp(jj)%act = 0.0
       if (soil1(j)%mp(jj)%sta < 1.e-6) soil1(j)%mp(jj)%sta = 0.0   
      end do

  endif
      !! Above added by KDW

      !if (sedyld(j) < 0.) sedyld(j) = 0.
      !if (sanyld(j) < 0.) sanyld(j) = 0.0
      !if (silyld(j) < 0.) silyld(j) = 0.0
      !if (clayld(j) < 0.) clayld(j) = 0.0
      !if (sagyld(j) < 0.) sagyld(j) = 0.0
      !if (lagyld(j) < 0.) lagyld(j) = 0.0
      !if (sedorgn(j) < 0.) sedorgn(j) = 0.0
      
      !! Below added by KDW 
      if (latorgp(j) < 1.e-6) latorgp(j) = 0.0
      if (latorgn(j) < 1.e-6) latorgn(j) = 0.0
      if (latminpa(j) < 1.e-6) latminpa(j) = 0.0
      if (latminps(j) < 1.e-6) latminps(j) = 0.0
      if (tileorgp(j) < 1.e-6) tileorgp(j) = 0.0
      if (tileorgn(j) < 1.e-6) tileorgn(j) = 0.0
      if (tileminpa(j) < 1.e-6) tileminpa(j) = 0.0
      if (tileminps(j) < 1.e-6) tileminps(j) = 0.0
      !! Above added by KDW

      return
      end subroutine swr_latsed