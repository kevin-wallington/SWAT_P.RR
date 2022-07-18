      subroutine hrudb_init
    
      use hydrograph_module, only : sp_ob, sp_ob1, ob
      use hru_module, only : hru, hru_db, ilu
      use urban_data_module
      use landuse_data_module ! KDW 12/21/21
      
      implicit none

      integer :: eof                  !           |end of file
      integer :: imp                  !           |
      integer :: ihru                 !none       |counter 
      integer :: iob                  !           |
      integer :: ihru_db              !           | 

      !!assign database pointers for the hru
      imp = 0
      do ihru = 1, sp_ob%hru
        iob = sp_ob1%hru + ihru - 1
        ihru_db = ob(iob)%props    !points to hru.dat
        hru(ihru)%dbs = hru_db(ihru_db)%dbs
        hru(ihru)%dbsc = hru_db(ihru_db)%dbsc
        hru(ihru)%parms = hru_db(ihru_db)%parms
        hru(ihru)%obj_no = sp_ob1%hru + ihru - 1
        hru(ihru)%area_ha = ob(iob)%area_ha
        hru(ihru)%km = ob(iob)%area_ha / 100.
        !Added by KDW 12/20/21, to fix yields from pervious and impervious areas in Urban HRUs
        ! if (hru(ihru)%luse%urb_lu > 0) then
        !   ulu = hru(ihru)%luse%urb_lu
        !   hru(ihru)%area_ha = hru(ihru)%area_ha * (1 - urbdb(ulu)%fimp)
        !   hru(ihru)%km_tot = hru(ihru)%km
        !   hru(ihru)%km = hru(ihru)%km * (1 - urbdb(ulu)%fimp)
        !End addition by KDW
        hru(ihru)%land_use_mgt_c = hru_db(ihru_db)%dbsc%land_use_mgt

      !!assign land use pointers for the hru - moved from plant_init by KDW, 12/21/21, so that occurs before topohyd_init
        ilu = hru(ihru)%dbs%land_use_mgt
        hru(ihru)%land_use_mgt = ilu
        hru(ihru)%plant_cov = lum_str(ilu)%plant_cov
        hru(ihru)%lum_group_c = lum(ilu)%cal_group
        hru(ihru)%mgt_ops = lum_str(ilu)%mgt_ops
        hru(ihru)%tiledrain = lum_str(ilu)%tiledrain
        hru(ihru)%septic = lum_str(ilu)%septic
        hru(ihru)%fstrip = lum_str(ilu)%fstrip
        hru(ihru)%grassww = lum_str(ilu)%grassww
        hru(ihru)%bmpuser = lum_str(ilu)%bmpuser
        hru(ihru)%luse%cn_lu = lum_str(ilu)%cn_lu
        hru(ihru)%luse%cons_prac = lum_str(ilu)%cons_prac
        !!assign land use pointers for the hru - moved from plant_init by KDW

      end do

      return
      end subroutine hrudb_init