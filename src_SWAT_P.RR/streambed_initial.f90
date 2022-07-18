      subroutine streambed_initial
      
      use basin_module
      use input_file_module
      use maximum_data_module
      use channel_data_module
      use hydrograph_module
      use sd_channel_module
      use constituent_mass_module
      use channel_module, only : streambed_init, streambed_init_name

      implicit none
      
      character (len=80) :: titldum     !          |title of file
      character (len=80) :: header      !          |header of file
      integer :: eof                    !          |end of file
      integer :: imax                   !units     |description
      logical :: i_exist                !          |check to determine if file exists
      integer :: ichi                   !none      |counter
      
      eof = 0
      imax = 0
      
      inquire (file=in_init%streambed, exist=i_exist)
      if (.not. i_exist .or. in_init%streambed == "null") then
        allocate (streambed_init(0:0))
        allocate (streambed_init_name(0:0))
      else   
      do
       open (105,file=in_init%streambed)
       read (105,*,iostat=eof) titldum
       if (eof < 0) exit
       read (105,*,iostat=eof) header
       if (eof < 0) exit
        do while (eof == 0)
          read (105,*,iostat=eof) titldum
          if (eof < 0) exit
          imax = imax + 1
        end do
        
      db_mx%streambed_init = imax
      
      allocate (streambed_init(0:imax))
      allocate (streambed_init_name(0:imax))
      rewind (105)
      read (105,*,iostat=eof) titldum
      if (eof < 0) exit
      read (105,*,iostat=eof) header
      if (eof < 0) exit
      
       do ichi = 1, db_mx%streambed_init
         read (105,*,iostat=eof) titldum
         if (eof < 0) exit
         backspace (105)
         read (105,*,iostat=eof) streambed_init_name(ichi), streambed_init(ichi)
         if (eof < 0) exit
       end do
       close (105)
      exit
      enddo
      endif

      return    
      end subroutine streambed_initial