      module water_body_module
    
      !! water body (reservoir, wetland, and channel) output not uncluded in hyd_output object

      type water_body
        real :: area_ha = 0.            !ha         |water body surface area
        real :: precip = 0.             !ha-m       |precip on the water body
        real :: evap = 0.               !ha-m       |evaporation from the water surface
        real :: seep = 0.               !ha-m       |seepage from bottom of water body
        real :: phos_dep = 0.           !kg         |deposition of phosphorus in water body bed - KDW, 06/07
        real :: phos_diff = 0.          !kg         |diffusion of phosphorus from bed of water body to water column - KDW, 06/07
      end type water_body
      type (water_body) :: wbodz

      !! Added by KDW
      type water_body_bed
      real :: bed_por = 0.3 ! fraction - KDW
      real :: bed_bd = 1.50     ! (g/cc)        |bulk density of channel bed sediment (1.1-1.9) - KDW
      real :: algae = 0.     ! mg alg/L      |algal biomass concentration in reservoir
      real :: aer_srfarea = 0. !
      real :: aer_orgp = 5. ! mg/kg or g/Mg
      real :: aer_minpa = 29.  ! mg/kg or g/Mg; - updated IC 08/10, KDW ,
      real :: aer_minps = 200.  ! mg/kg or g/Mg; - updated IC 08/10, KDW
      real :: aer_disp = 0.1  ! mg/L - was 3? KDW 08/10
      real :: anaer_srfarea = 0. !
      real :: anaer_orgp = 5. ! mg/kg or g/Mg
      real :: anaer_minpa = 14. ! mg/kg or g/Mg; - updated IC 08/10, KDW
      real :: anaer_minps = 100. ! mg/kg or g/Mg; - updated IC 08/10, KDW
      real :: anaer_disp = 0.2 ! mg/L  - was 3? KDW 08/10
      real :: ext_orgp = 5. ! mg/kg or g/Mg 
      real :: ext_minpa = 10. ! mg/kg or g/Mg; - updated IC 08/10, KDW 
      real :: ext_minps = 70. ! mg/kg or g/Mg; - updated IC 08/10, KDW
      real :: ext_disp = 0.3 ! mg/L  - was 3? KDW 08/10
      end type water_body_bed
      type (water_body_bed), dimension(:), allocatable, target :: res_bed
      type (water_body_bed), dimension(:), allocatable, target :: wet_bed
      type (water_body_bed), pointer :: wbody_bed       !! used for reservoir and wetlands
      type (water_body_bed) :: wbody_bed_init ! KDW 03/26
      !! Added by KDW
      
      type (water_body), dimension(:), allocatable, target :: res_wat_d
      type (water_body), dimension(:), allocatable :: res_wat_m
      type (water_body), dimension(:), allocatable :: res_wat_y
      type (water_body), dimension(:), allocatable :: res_wat_a
      type (water_body), dimension(:), allocatable, target :: wet_wat_d
      type (water_body), dimension(:), allocatable :: wet_wat_m
      type (water_body), dimension(:), allocatable :: wet_wat_y
      type (water_body), dimension(:), allocatable :: wet_wat_a
      type (water_body), dimension(:), allocatable :: ch_wat_d
      type (water_body), dimension(:), allocatable :: ch_wat_m
      type (water_body), dimension(:), allocatable :: ch_wat_y
      type (water_body), dimension(:), allocatable :: ch_wat_a
      type (water_body) :: bch_wat_d
      type (water_body) :: bch_wat_m
      type (water_body) :: bch_wat_y
      type (water_body) :: bch_wat_a
      type (water_body) :: bres_wat_d
      type (water_body) :: bres_wat_m
      type (water_body) :: bres_wat_y
      type (water_body) :: bres_wat_a
      type (water_body), pointer :: wbody_wb       !! used for reservoir and wetlands


      ! KDW, 06/07 below
      type water_body_header  
      character (len=15) :: bed_por           = "        bed_por"       
      character (len=15) :: bed_bd            = "         bed_bd"
      character (len=15) :: algae        = "          algae"
      character (len=15) :: aer_srfarea           = "    aer_srfarea"       
      character (len=15) :: aer_orgp           = "       aer_orgp"        
      character (len=15) :: aer_minpa            = "      aer_minpa"       
      character (len=15) :: aer_minps         = "      aer_minps"   
      character (len=15) :: aer_disp      = '       aer_disp'
      character (len=15) :: anaer_srfarea       = '  anaer_srfarea'
      character (len=15) :: anaer_orgp         = '     anaer_orgp'
      character (len=15) :: anaer_minpa         = '    anaer_minpa'
      character (len=15) :: anaer_minps         = '    anaer_minps'
      character (len=15) :: anaer_disp         = '     anaer_disp'
      character (len=15) :: ext_orgp         = '     anaer_orgp'
      character (len=15) :: ext_minpa         = '    anaer_minpa'
      character (len=15) :: ext_minps         = '    anaer_minps'
      character (len=15) :: ext_disp         = '     anaer_disp'
    end type water_body_header 
    type (water_body_header) :: water_body_hdr
    
    type water_body_header_units
    character (len=15) :: bed_por           = "           %   "       
    character (len=15) :: bed_bd            = "         Mg/m^3"
    character (len=15) :: algae        = "          kg   "
    character (len=15) :: aer_srfarea           = "    kg clay    "       
    character (len=15) :: aer_orgp           = "       kg      "        
    character (len=15) :: aer_minpa            = "      kg       "       
    character (len=15) :: aer_minps         = "      kg       "   
    character (len=15) :: aer_disp      = '       kg      '
    character (len=15) :: anaer_srfarea       = '    kg clay    '
    character (len=15) :: anaer_orgp         = '     kg        '
    character (len=15) :: anaer_minpa         = '       kg      '
    character (len=15) :: anaer_minps         = '        kg     '
    character (len=15) :: anaer_disp         = '         kg    '
    character (len=15) :: ext_orgp         = '         kg    '
    character (len=15) :: ext_minpa         = '         kg    '
    character (len=15) :: ext_minps         = '         kg    '
    character (len=15) :: ext_disp         = '         kg    '
    end type water_body_header_units
    type (water_body_header_units) :: water_body_hdr_units
   ! KDW, 06/07 above


       interface operator (+)
        module procedure watbod_add
      end interface

      interface operator (/)
        module procedure watbod_div
      end interface   
      
      interface operator (//)
        module procedure watbod_ave
      end interface   
      
      contains
      
     !! routines for reservoir module
      function watbod_add (wbod1, wbod2) result (wbod3)
        type (water_body), intent (in) :: wbod1
        type (water_body), intent (in) :: wbod2
        type (water_body) :: wbod3
        wbod3%area_ha = wbod1%area_ha + wbod2%area_ha
        wbod3%precip = wbod1%precip + wbod2%precip
        wbod3%evap = wbod1%evap + wbod2%evap        
        wbod3%seep = wbod1%seep + wbod2%seep   
        wbod3%phos_dep = wbod1%phos_dep + wbod2%phos_dep
        wbod3%phos_diff = wbod1%phos_diff + wbod2%phos_diff
      end function watbod_add
      
      function watbod_div (wbod1,const) result (wbod2)
        type (water_body), intent (in) :: wbod1
        real, intent (in) :: const
        type (water_body) :: wbod2
        wbod2%area_ha = wbod1%area_ha
        wbod2%precip = wbod1%precip / const
        wbod2%evap = wbod1%evap / const
        wbod2%seep = wbod1%seep / const
        wbod2%phos_dep = wbod1%phos_dep / const
        wbod2%phos_diff = wbod1%phos_diff / const
      end function watbod_div
            
      function watbod_ave (wbod1,const) result (wbod2)
        type (water_body), intent (in) :: wbod1
        real, intent (in) :: const
        type (water_body) :: wbod2
        wbod2%area_ha = wbod1%area_ha / const
        wbod2%precip = wbod1%precip
        wbod2%evap = wbod1%evap
        wbod2%seep = wbod1%seep
        wbod2%phos_dep = wbod1%phos_dep
        wbod2%phos_diff = wbod1%phos_diff
      end function watbod_ave
      
      end module water_body_module