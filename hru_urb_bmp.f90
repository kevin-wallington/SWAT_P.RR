      subroutine hru_urb_bmp
!!    ~ ~ ~ PURPOSE ~ ~ ~

    use hru_module, only : hru, ihru, sedyld, sedorgn, sedorgp, sedminpa, sedminps, sed_con, soln_con, solp_con,  &
       surqno3, latno3, surqsolp, qdr, orgn_con, orgp_con, &
       clayld, silyld, sanyld, sagyld, lagyld ! KDW 08/02
    
    implicit none
    
    integer :: j          !none          |HRU number 
    real :: xx         !              | changed from integer to real - KDW 08/02
    real :: sedppm        !              |
    real :: solnppm       !              |
    real :: solpppm       !              | 
    real :: sednppm       !              |
    real :: sedpppm       !              |
    real :: percent_remain
    real :: yy

	j = 0
	j = ihru

!! convert to ppm -> (kg/ha)*100./mm = ppm
      if (qdr(j) > 0.1) then
	  xx = 100. / qdr(j)
        sedppm = 1000. * xx * sedyld(j) / hru(j)%area_ha
        solnppm = xx * (surqno3(j) + latno3(j))
        solpppm = xx * surqsolp(j)
        sednppm = xx * sedorgn(j)
        sedpppm = xx * (sedorgp(j) + sedminpa(j) + sedminps(j))
      
        if (sedppm > sed_con (j)) then
          yy = sedyld(j) ! KDW 08/02
          sedyld(j) = sed_con(j) * hru(j)%area_ha / xx / 1000.
          percent_remain = sedyld(j) / yy ! KDW 08/02
          clayld(j) = clayld(j) * percent_remain ! KDW 08/02
          silyld(j) = silyld(j) * percent_remain ! KDW 08/02
          sanyld(j) = sanyld(j) * percent_remain ! KDW 08/02
          sagyld(j) = sagyld(j) * percent_remain ! KDW 08/02
          lagyld(j) = lagyld(j) * percent_remain ! KDW 08/02
        endif

        if (solnppm > soln_con(j)) then
	    surqno3(j) = soln_con(j) / xx
          latno3(j) = soln_con(j) / xx
        endif

        if (solpppm > solp_con(j)) then
          surqsolp(j) = solp_con(j) / xx
        endif

        if (sednppm > orgn_con(j)) then
          sedorgn(j) = orgn_con(j) / xx
        endif

        if (sedpppm > orgp_con(j)) then 
	    sedorgp(j)= orgp_con(j) / xx !corrected from sedorgn(j) = to sedorgp(j) = , KDW 10/21/21
          sedminpa(j)= orgp_con(j) / xx
	    sedminps(j)= orgp_con(j) / xx
        endif

	endif

      return
      end subroutine hru_urb_bmp