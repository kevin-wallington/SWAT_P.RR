      subroutine ch_initial (idat, irch)

      use channel_data_module
      use channel_module
	  use hydrograph_module, only : om_init_water, om_init_name
	  use maximum_data_module, only : db_mx
      
      implicit none
      
      integer :: irch     !none      |counter
      integer :: ised     !none      |counter
      integer :: idat     !units     |description
      real :: bnksize     !units     |description
      real :: bedsize     !units     |description
      real :: sc          !units     |description   
	  integer :: i ! KDW
	  integer :: init ! KDW 
	  integer :: ics ! KDW
      
      ised = ch_dat(idat)%sed
      bnksize = ch_sed(ised)%bnk_d50 / 1000.  !! Units conversion Micrometer to Millimeters
!!    Channel sediment particle size distribution
!!    Clayey bank
	if (bnksize <= 0.005) then
	  ch(irch)%bnk_cla = 0.65
        ch(irch)%bnk_sil = 0.15
	  ch(irch)%bnk_san = 0.15
	  ch(irch)%bnk_gra = 0.05
	end if

!!    Silty bank
	if (bnksize > 0.005 .and. bnksize <= 0.05) then
        ch(irch)%bnk_sil = 0.65
	  ch(irch)%bnk_cla = 0.15
	  ch(irch)%bnk_san = 0.15
	  ch(irch)%bnk_gra = 0.05
	end if

!!    Sandy bank
	if (bnksize > 0.05 .and. bnksize <= 2.) then
	  ch(irch)%bnk_san = 0.65
        ch(irch)%bnk_sil = 0.15
	  ch(irch)%bnk_cla = 0.15
	  ch(irch)%bnk_gra = 0.05
	end if
      
!!    Gravel bank
	if (bnksize > 2.) then
	  ch(irch)%bnk_gra = 0.65
	  ch(irch)%bnk_san = 0.15
        ch(irch)%bnk_sil = 0.15
	  ch(irch)%bnk_cla = 0.05
	end if

!!    Channel sediment particle size distribution
!!    Clayey bed
      bedsize = ch_sed(ised)%bed_d50 / 1000.  !! Units conversion Micrometer to Millimeters
	if (bedsize <= 0.005) then
	  ch(irch)%bed_cla = 0.65
        ch(irch)%bed_sil = 0.15
	  ch(irch)%bed_san = 0.15
	  ch(irch)%bed_gra = 0.05
	end if

!!    Silty bed
	if (bedsize > 0.005 .and. bedsize <= 0.05) then
        ch(irch)%bed_sil = 0.65
	  ch(irch)%bed_cla = 0.15
	  ch(irch)%bed_san = 0.15
	  ch(irch)%bed_gra = 0.05
	end if

!!    Sandy bed
	if (bedsize > 0.05 .and. bedsize <= 2.) then
	  ch(irch)%bed_san = 0.65
        ch(irch)%bed_sil = 0.15
	  ch(irch)%bed_cla = 0.15
	  ch(irch)%bed_gra = 0.05
	end if
      
!!    Gravel bed
	if (bedsize > 2.) then
	  ch(irch)%bed_gra = 0.65
	  ch(irch)%bed_san = 0.15
        ch(irch)%bed_sil = 0.15
	  ch(irch)%bed_cla = 0.05
      end if
      
!!    An estimate of Critical shear stress if it is not given (N/m^2)
!!	Critical shear stress based on silt and clay %
!!	Critical Shear Stress based on Julian and Torres (2005)
!!    Units of critical shear stress (N/m^2)
	SC = 0.
	if  (ch_sed(ised)%tc_bnk <= 1.e-6) then 
	  SC = (ch(irch)%bnk_sil + ch(irch)%bnk_cla) * 100.
        ch_sed(ised)%tc_bnk = (0.1 + (0.1779*SC) + (0.0028*(SC)**2)       &
                         - ((2.34E-05)*(SC)**3)) * ch_sed(ised)%cov1
      end if

	if  (ch_sed(ised)%tc_bed <= 1.e-6) then
	  SC = (ch(irch)%bed_sil + ch(irch)%bed_cla) * 100.
        ch_sed(ised)%tc_bed = (0.1 + (0.1779*SC) + (0.0028*(SC)**2)       &
                         - ((2.34E-05)*(SC)**3)) * ch_sed(ised)%cov2
	  end if
	  
!! Copied from ch_read_sed to recalculate bnk_kd and bed_kd - KDW
!!  An estimate of channel bank erodibility coefficient from jet test if it is not available
!!  Units of kd is (cm^3/N/s)
!!  Base on Hanson and Simon, 2001
	if (ch_sed(ised)%bnk_kd <= 1.e-6) then ! if loop added 08/09, KDW
		if (ch_sed(ised)%tc_bnk <= 1.e-6) then
		  ch_sed(ised)%bnk_kd = 0.2
		else 
			ch_sed(ised)%bnk_kd = 0.2 / sqrt(ch_sed(ised)%tc_bnk)
		end if
	end if
  !!  An estimate of channel bed erodibility coefficient from jet test if it is not available
  !!  Units of kd is (cm^3/N/s)
  !!  Base on Hanson and Simon, 2001
	if (ch_sed(ised)%bed_kd <= 1.e-6) then
		if (ch_sed(ised)%tc_bed <= 1.e-6) then
		  ch_sed(ised)%bed_kd = 0.2
		else 
		ch_sed(ised)%bed_kd = 0.2 / sqrt(ch_sed(ised)%tc_bed)
		end if
	end if
!! Copied from ch_read_sed to recalculate bnk_kd and bed_kd - KDW

!! Below from KDW, based on res_initial and sd_hydsed_init, to initialize channel storage and water quality
	if (ch_hyd(irch)%l > 1.e-3) then
		i = ch_dat(idat)%init
		do ics = 1, db_mx%om_water_init
			if (ch_init(i)%org_min == om_init_name(ics)) then
				init = ics
			endif
		enddo
		!init = ch_init(i)%org_min
		ch(irch)%rchstor = om_init_water(init)%flo * ch_hyd(irch)%w * ch_hyd(irch)%d * ch_hyd(irch)%l * 1000 ! to convert from % full to m^3
		ch(irch)%chlora = om_init_water(init)%chla
		ch(irch)%disolvp = om_init_water(init)%solp
		ch(irch)%nitraten = om_init_water(init)%no3
		ch(irch)%nitriten = om_init_water(init)%no2
		ch(irch)%organicn = om_init_water(init)%orgn
		ch(irch)%organicp = om_init_water(init)%orgp
		ch(irch)%minpa = om_init_water(init)%minpa
		ch(irch)%minps = om_init_water(init)%minps
		ch(irch)%rch_cbod = om_init_water(init)%cbod
		ch(irch)%rch_dox = om_init_water(init)%dox
		ch(irch)%sedst = om_init_water(init)%sed * ch(irch)%rchstor / 1000000 ! to convert from mg/L to Mg
		ch(irch)%sanst = om_init_water(init)%san * ch(irch)%rchstor / 1000000
		ch(irch)%silst = om_init_water(init)%sil * ch(irch)%rchstor / 1000000
		ch(irch)%clast = om_init_water(init)%cla * ch(irch)%rchstor / 1000000
		ch(irch)%sagst = om_init_water(init)%sag * ch(irch)%rchstor / 1000000
		ch(irch)%lagst = om_init_water(init)%lag * ch(irch)%rchstor / 1000000
		ch(irch)%grast = om_init_water(init)%grv * ch(irch)%rchstor / 1000000
		ch(irch)%ammonian = om_init_water(init)%nh3
		!! Below add by KDW on 11/04/21
		do ics = 1, db_mx%streambed_init
			if (ch_init(i)%bed == streambed_init_name(ics)) then
				init = ics
			endif
		enddo
		ch(irch)%aer_orgp = streambed_init(init)%aer_orgp
		ch(irch)%aer_minpa = streambed_init(init)%aer_minpa
		ch(irch)%aer_minps = streambed_init(init)%aer_minps
		ch(irch)%aer_disp = streambed_init(init)%aer_disp
		ch(irch)%anaer_orgp = streambed_init(init)%anaer_orgp
		ch(irch)%anaer_minpa = streambed_init(init)%anaer_minpa
		ch(irch)%anaer_minps = streambed_init(init)%anaer_minps
		ch(irch)%anaer_disp = streambed_init(init)%anaer_disp
		! ch(irch)%ext_orgp = streambed_init(init)%anaer_orgp ! KDW 02/15/22, added external assignment here and below.
		! ch(irch)%ext_minpa = streambed_init(init)%anaer_minpa
		! ch(irch)%ext_minps = streambed_init(init)%anaer_minps
		! ch(irch)%ext_disp = streambed_init(init)%anaer_disp

		if(irch == 25) then ! KDW debug 04/07/22, test/temp
			ch(irch)%aer_orgp = 0
			ch(irch)%aer_minpa = 0
			ch(irch)%aer_minps = 0
			ch(irch)%aer_disp = 0
			ch(irch)%anaer_orgp = 0
			ch(irch)%anaer_minpa = 0
			ch(irch)%anaer_minps = 0
			ch(irch)%anaer_disp = 0
			ch(irch)%ext_orgp = 0
			ch(irch)%ext_minpa = 0
			ch(irch)%ext_minps = 0
			ch(irch)%ext_disp = 0
		endif
		
	end if

      return
      end subroutine ch_initial