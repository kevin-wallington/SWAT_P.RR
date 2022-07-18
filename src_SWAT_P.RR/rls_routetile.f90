      subroutine rls_routetile (iob)
      
!!    ~ ~ ~ PURPOSE ~ ~ ~

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!
      use hru_module, only : ihru, hru
      use soil_module
      use hydrograph_module
      use organic_mineral_mass_module
      
      implicit none
 
      integer, intent (in) :: iob   !           |object number
      integer :: j                  !           |hru number
      integer :: lyr                !           |tile soil layer 

      j = ihru

      !! add tile inflow and nitrate to soil layer the tile is in
      !! if exceeds saturation, it will be redistributed in swr_satexcess
      lyr = hru(j)%lumv%ldrain
      soil(j)%phys(lyr)%st = soil(j)%phys(lyr)%st + ob(iob)%hin_til%flo
      soil1(j)%mn(lyr)%no3 = soil1(j)%mn(lyr)%no3 + ob(iob)%hin_til%no3
      soil1(j)%mp(lyr)%lab = soil1(j)%mp(lyr)%lab + ob(iob)%hin_til%solp ! - KDW
      soil1(j)%mp(lyr)%act = soil1(j)%mp(lyr)%act + ob(iob)%hin_til%minpa ! - KDW
      soil1(j)%mp(lyr)%sta = soil1(j)%mp(lyr)%sta + ob(iob)%hin_til%minps ! - KDW
      soil1(j)%hp(lyr)%p = soil1(j)%hp(lyr)%p + ob(iob)%hin_til%orgp ! - KDW
      soil1(j)%hs(lyr)%n = soil1(j)%hs(lyr)%n + ob(iob)%hin_til%orgn ! - KDW

      return
      end subroutine rls_routetile