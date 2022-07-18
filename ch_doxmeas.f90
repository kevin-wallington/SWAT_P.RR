      subroutine ch_doxmeas
      
      use climate_module
      use maximum_data_module
      use basin_module
      use input_file_module
      use time_module
      
      implicit none
            
      character (len=80) :: titldum   !           |title of file
      character (len=80) :: header    !           |header of file
      integer :: eof                  !           |end of file
      integer :: iyr                  !none       |number of years 
      logical :: i_exist              !none       |check to determine if file exists 
      integer :: istep                !           |
      integer :: iyr_prev             !none       |previous year
      integer :: iyrs                 !           |
            
        eof = 0
        
      ! weather path code
        open (108,file = "wyck.dox")   
        read (108,*,iostat=eof) titldum
        read (108,*,iostat=eof) header
        read (108,*,iostat=eof) dox_meas%nbyr, dox_meas%tstep
        allocate (dox_meas%ts(366,dox_meas%nbyr))
        read (108,*,iostat=eof) iyr, istep
        backspace (108)
        iyr_prev = iyr
        iyrs = 1
       
        do  !!!!!!!this is where data is read!!!!!!!!!!
            read (108,*,iostat=eof)iyr, istep, dox_meas%ts(istep,iyrs)
            if (eof < 0) exit
            !check to see when next year
            if (istep == 365 .or. istep == 366) then
                read (108,*,iostat=eof) iyr, istep
                if (eof < 0) exit
                backspace (108)
                if (iyr /= iyr_prev) then
                    iyr_prev = iyr
                    iyrs = iyrs + 1
                end if
            end if
        end do
        close (108)

        return
      end subroutine ch_doxmeas