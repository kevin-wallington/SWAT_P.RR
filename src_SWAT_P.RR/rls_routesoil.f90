      subroutine rls_routesoil (iob)
      
!!    ~ ~ ~ PURPOSE ~ ~ ~

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!
      use hru_module, only : ihru, latqrunon
      use soil_module
      use hydrograph_module
      use organic_mineral_mass_module

      
      implicit none

      real :: latqlyr               !           |
      real :: xslat                 !           | 
      integer :: j                  !           |
      integer, intent (in) :: iob   !           | 
      integer :: dep                !           |  
      integer :: lyr                !none       |counter
      real :: percent_lyr

      j = ihru
      
      latqrunon = ob(iob)%hin_lat%flo
      if (latqrunon > 1.e-9) then
      !!put in soil layers - weighted by depth of soil layer
        dep = 0.
        xslat = 0.
        do lyr = 1, soil(j)%nly
          !latqlyr = ((soil(j)%phys(lyr)%d - dep) / soil(j)%phys(soil(j)%nly)%d) * latqrunon ! replace below by KDW
          percent_lyr = ((soil(j)%phys(lyr)%d - dep) / soil(j)%phys(soil(j)%nly)%d) ! - KDW
          latqlyr = percent_lyr * latqrunon ! - KDW
          dep = soil(j)%phys(lyr)%d
          soil(j)%phys(lyr)%st = soil(j)%phys(lyr)%st + latqlyr
          soil1(j)%mn(lyr)%no3 = soil1(j)%mn(lyr)%no3 + percent_lyr * ob(iob)%hin_lat%no3
          soil1(j)%mp(lyr)%lab = soil1(j)%mp(lyr)%lab + percent_lyr * ob(iob)%hin_lat%solp ! - KDW
          soil1(j)%mp(lyr)%act = soil1(j)%mp(lyr)%act + percent_lyr * ob(iob)%hin_lat%minpa ! - KDW
          soil1(j)%mp(lyr)%sta = soil1(j)%mp(lyr)%sta + percent_lyr * ob(iob)%hin_lat%minps ! - KDW
          soil1(j)%hp(lyr)%p = soil1(j)%hp(lyr)%p + percent_lyr * ob(iob)%hin_lat%orgp ! - KDW
          soil1(j)%hs(lyr)%n = soil1(j)%hs(lyr)%n + percent_lyr * ob(iob)%hin_lat%orgn ! - KDW

          !if (soil(j)%phys(lyr)%st > soil(j)%phys(lyr)%ul) then
          !  xslat = xslat + (soil(j)%phys(lyr)%st - soil(j)%phys(lyr)%ul)
          !  soil(j)%phys(lyr)%st = soil(j)%phys(lyr)%ul
          !end if
        end do
        !add excess to surface storage
      end if

      return
      end subroutine rls_routesoil