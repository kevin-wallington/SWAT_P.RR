      subroutine swr_substor
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine stores and lags lateral soil flow and nitrate

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bss(1,:)      |mm H2O       |amount of lateral flow lagged
!!    bss(2,:)      |kg N/ha      |amount of nitrate in lateral flow lagged

!!    bss(3,:)      |mm           |amount of tile flow lagged
!!    bss(4,:)      |kg N/ha      |amount of nitrate in tile flow lagged
!!    ihru          |none         |HRU number
!!    lat_pst(:)    |kg pst/ha    |amount of pesticide in lateral flow in HRU
!!                                |for the day
!!    lat_ttime(:)  |none         |Exponential of the lateral flow travel time
!!    latq(:)       |mm H2O       |amount of water in lateral flow in HRU for
!!                                |the day
!!    qtile(:)      |mm H2O       |amount of water in tile flow in HRU for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bss(1,:)      |mm H2O        |amount of lateral flow lagged
!!    bss(2,:)      |kg N/ha       |amount of nitrate in lateral flow lagged
!!    bss(3,:)      |mm            |amount of tile flow lagged
!!    bss(4,:)      |kg N/ha       |amount of nitrate in tile flow lagged
!!    lat_pst(:)    |kg pst/ha     |amount of pesticide in lateral flow in HRU
!!                                 |for the day
!!    latq(:)       |mm H2O        |amount of water in lateral flow in HRU for
!!                                 |the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use pesticide_data_module
      use hru_module, only : hru, ihru, latq, latno3, qtile, bss, tileno3, &
      latsed, latsan, latsil, latcla, latsag, latlag,  &
      tilesed, tilesan, tilesil, tilecla, latorgp, tileorgp, latorgn, tileorgn, &
      latminpa, tileminpa, latminps, tileminps, tilesag, tilelag, tilesolp, latsolp
      use constituent_mass_module
      use output_ls_pesticide_module
      
      implicit none      
      
      integer :: j           !none          |HRU number
      integer :: k           !none          |counter
      integer :: icmd        !              |  
      integer :: ipst_db     !              |  


      j = ihru

      bss(1,j) = bss(1,j) + latq(j)
      bss(2,j) = bss(2,j) + latno3(j)
      bss(3,j) = bss(3,j) + qtile
      bss(4,j) = bss(4,j) + tileno3(j)
      !! Below added by KDW
      bss(5,j) = bss(5,j) + latorgn(j)
      bss(6,j) = bss(6,j) + tileorgn(j)
      bss(7,j) = bss(7,j) + latsolp(j)
      bss(8,j) = bss(8,j) + tilesolp(j)
      bss(9,j) = bss(9,j) + latorgp(j)
      bss(10,j) = bss(10,j) + tileorgp(j)
      bss(11,j) = bss(11,j) + latminpa(j)
      bss(12,j) = bss(12,j) + tileminpa(j)
      bss(13,j) = bss(13,j) + latminps(j)
      bss(14,j) = bss(14,j) + tileminps(j)
      bss(15,j) = bss(15,j) + latsed(j)
      bss(16,j) = bss(16,j) + tilesed(j)
      bss(17,j) = bss(17,j) + latsan(j)
      bss(18,j) = bss(18,j) + tilesan(j)
      bss(19,j) = bss(19,j) + latsil(j)
      bss(20,j) = bss(20,j) + tilesil(j)
      bss(21,j) = bss(21,j) + latcla(j)
      bss(22,j) = bss(22,j) + tilecla(j)
      bss(23,j) = bss(23,j) + latsag(j)
      bss(24,j) = bss(24,j) + tilesag(j)
      bss(25,j) = bss(25,j) + latlag(j)
      bss(26,j) = bss(26,j) + tilelag(j)
      !! Above added by KDW

      if (bss(1,j) < 1.e-6) bss(1,j) = 0.0
      if (bss(2,j) < 1.e-6) bss(2,j) = 0.0
      if (bss(3,j) < 1.e-6) bss(3,j) = 0.0
      if (bss(4,j) < 1.e-6) bss(4,j) = 0.0
      !! Below added by KDW
      if (bss(5,j) < 1.e-6) bss(5,j) = 0.0
      if (bss(6,j) < 1.e-6) bss(6,j) = 0.0
      if (bss(7,j) < 1.e-6) bss(7,j) = 0.0
      if (bss(8,j) < 1.e-6) bss(8,j) = 0.0
      if (bss(9,j) < 1.e-6) bss(9,j) = 0.0
      if (bss(10,j) < 1.e-6) bss(10,j) = 0.0
      if (bss(11,j) < 1.e-6) bss(11,j) = 0.0
      if (bss(12,j) < 1.e-6) bss(12,j) = 0.0
      if (bss(13,j) < 1.e-6) bss(13,j) = 0.0
      if (bss(14,j) < 1.e-6) bss(14,j) = 0.0
      if (bss(15,j) < 1.e-6) bss(15,j) = 0.0
      if (bss(16,j) < 1.e-6) bss(16,j) = 0.0
      if (bss(17,j) < 1.e-6) bss(17,j) = 0.0
      if (bss(18,j) < 1.e-6) bss(18,j) = 0.0
      if (bss(19,j) < 1.e-6) bss(19,j) = 0.0
      if (bss(20,j) < 1.e-6) bss(20,j) = 0.0
      if (bss(21,j) < 1.e-6) bss(21,j) = 0.0
      if (bss(22,j) < 1.e-6) bss(22,j) = 0.0
      if (bss(23,j) < 1.e-6) bss(23,j) = 0.0
      if (bss(24,j) < 1.e-6) bss(24,j) = 0.0
      if (bss(25,j) < 1.e-6) bss(25,j) = 0.0
      if (bss(26,j) < 1.e-6) bss(26,j) = 0.0
      !! Above added by KDW

      latq(j) = bss(1,j) * hru(j)%hyd%lat_ttime
      latno3(j) = bss(2,j) * hru(j)%hyd%lat_ttime
      qtile = bss(3,j) * hru(j)%lumv%tile_ttime
      tileno3(j) = bss(4,j) * hru(j)%lumv%tile_ttime
      !! Below added by KDW
      latorgn(j) = bss(5,j) * hru(j)%hyd%lat_ttime
      tileorgn(j) = bss(6,j) * hru(j)%lumv%tile_ttime
      latsolp(j) = bss(7,j) * hru(j)%hyd%lat_ttime
      tilesolp(j) = bss(8,j) * hru(j)%lumv%tile_ttime
      latorgp(j) = bss(9,j) * hru(j)%hyd%lat_ttime
      tileorgp(j) = bss(10,j) * hru(j)%lumv%tile_ttime
      latminpa(j) = bss(11,j) * hru(j)%hyd%lat_ttime 
      tileminpa(j) = bss(12,j) * hru(j)%lumv%tile_ttime
      latminps(j) = bss(13,j) * hru(j)%hyd%lat_ttime
      tileminps(j) = bss(14,j) * hru(j)%lumv%tile_ttime
      latsed(j) = bss(15,j) * hru(j)%hyd%lat_ttime
      tilesed(j) = bss(16,j) * hru(j)%lumv%tile_ttime
      latsan(j) = bss(17,j) * hru(j)%hyd%lat_ttime
      tilesan(j) = bss(18,j) * hru(j)%lumv%tile_ttime
      latsil(j) = bss(19,j) * hru(j)%hyd%lat_ttime 
      tilesil(j) = bss(20,j) * hru(j)%lumv%tile_ttime
      latcla(j) = bss(21,j) * hru(j)%hyd%lat_ttime
      tilecla(j) = bss(22,j) * hru(j)%lumv%tile_ttime 
      latsag(j) = bss(23,j) * hru(j)%hyd%lat_ttime
      tilesag(j) = bss(24,j) * hru(j)%lumv%tile_ttime
      latlag(j) = bss(25,j) * hru(j)%hyd%lat_ttime
      tilelag(j) = bss(26,j) * hru(j)%lumv%tile_ttime
      !! Above added by KDW

      if (latq(j) < 1.e-6) latq(j) = 0.
      if (latno3(j) < 1.e-6) latno3(j) = 0.
      if (qtile < 1.e-6) qtile = 0.
      if (tileno3(j) < 1.e-6) tileno3(j) = 0.
      !! Below added by KDW
      if (latorgn(j) < 1.e-6) latorgn(j) = 0.
      if (tileorgn(j) < 1.e-6) tileorgn(j) = 0.
      if (latsolp(j) < 1.e-6) latsolp(j) = 0.
      if (tilesolp(j) < 1.e-6) tilesolp(j) = 0.
      if (latorgp(j) < 1.e-6) latorgp(j) = 0.
      if (tileorgp(j) < 1.e-6) tileorgp(j) = 0.
      if (latminpa(j) < 1.e-6) latminpa(j) = 0.
      if (tileminpa(j) < 1.e-6) tileminpa(j) = 0.
      if (latminps(j) < 1.e-6) latminps(j) = 0.
      if (tileminps(j) < 1.e-6) tileminps(j) = 0.
      if (latsed(j) < 1.e-6) latsed(j) = 0.
      if (tilesed(j) < 1.e-6) tilesed(j) = 0.
      if (latsan(j) < 1.e-6) latsan(j) = 0.
      if (tilesan(j) < 1.e-6) tilesan(j) = 0.
      if (latsil(j) < 1.e-6) latsil(j) = 0.
      if (tilesil(j) < 1.e-6) tilesil(j) = 0.
      if (latcla(j) < 1.e-6) latcla(j) = 0.
      if (tilecla(j) < 1.e-6) tilecla(j) = 0.
      if (latsag(j) < 1.e-6) latsag(j) = 0.
      if (tilesag(j) < 1.e-6) tilesag(j) = 0.
      if (latlag(j) < 1.e-6) latlag(j) = 0.
      if (tilelag(j) < 1.e-6) tilelag(j) = 0.
      !! Above added by KDW

      bss(1,j) = bss(1,j) - latq(j)
      bss(2,j) = bss(2,j) - latno3(j)
      bss(3,j) = bss(3,j) - qtile
      bss(4,j) = bss(4,j) - tileno3(j)
      !! Below added by KDW
      bss(5,j) = bss(5,j) - latorgn(j)
      bss(6,j) = bss(6,j) - tileorgn(j)
      bss(7,j) = bss(7,j) - latsolp(j)
      bss(8,j) = bss(8,j) - tilesolp(j)
      bss(9,j) = bss(9,j) - latorgp(j)
      bss(10,j) = bss(10,j) - tileorgp(j)
      bss(11,j) = bss(11,j) - latminpa(j)
      bss(12,j) = bss(12,j) - tileminpa(j)
      bss(13,j) = bss(13,j) - latminps(j)
      bss(14,j) = bss(14,j) - tileminps(j)
      bss(15,j) = bss(15,j) - latsed(j)
      bss(16,j) = bss(16,j) - tilesed(j)
      bss(17,j) = bss(17,j) - latsan(j)
      bss(18,j) = bss(18,j) - tilesan(j)
      bss(19,j) = bss(19,j) - latsil(j)
      bss(20,j) = bss(20,j) - tilesil(j)
      bss(21,j) = bss(21,j) - latcla(j)
      bss(22,j) = bss(22,j) - tilecla(j)
      bss(23,j) = bss(23,j) - latsag(j)
      bss(24,j) = bss(24,j) - tilesag(j)
      bss(25,j) = bss(25,j) - latlag(j)
      bss(26,j) = bss(26,j) - tilelag(j)
      !! Above added by KDW

      return
      end subroutine swr_substor