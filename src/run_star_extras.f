! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def

      implicit none

      integer :: time0, time1, clock_rate

      ! these routines are called by the standard run_star check_model
      contains

        integer function how_many_extra_history_columns(id, id_extra)
          integer, intent(in) :: id, id_extra
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_history_columns = 3
        end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)

   ! MESA provides routines for taking logarithms which will handle
   ! cases where you pass them zero, etc.  Take advantage of this!
      use crlibm_lib, only: safe_log10_cr

      integer, intent(in) :: id, id_extra, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      real(dp), parameter :: frac = 0.90, Rsun = 6.69d10
      real(dp), parameter :: pi = 3.1415927
      integer :: i
      real(dp) :: B_max, B_ave, B_region
      integer :: k, j

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

   ! calculate the total nuclear energy release rate by integrating
   ! the specific rate (eps_nuc) over the star.  using the dot_product
   ! intrinsic is a common idiom for calculating integrated quantities.
   ! note that one needs to explicitly limit the range of the arrays.
   ! NEVER assume that the array size is the same as the number of zones.
   k = s% nz
   do while ( s% conv_vel(k) /= 0 )
    k = k - 1
  end do
!all these values are to adjust these constants for the general case
   i = k
   j = k

   do while ((s% r(i) .ge. (s% r(k) - s% mlt_mixing_length(k))) .and. s% conv_vel(1) /= 0)
    !write(*,*) 'R(i)= ', s% r(i)
    !write(*,*) 'MLT LengtH=  ', s% mlt_mixing_length(k)
    !write(*,*) 'r(k) - MLT Lenght=  ', (s% r(k) - s% mlt_mixing_length(k))
    i = i + 1
  end do
  do while (s% r(k) >= s% scale_height(j))
    j = j + 1
  end do

  !real(dp), dimension(k-s% nz) :: B_val

  !do while (s% r(k) .ge. s% scale_height(j)) then
   !B_val(j) = (sqrt(4 * pi * (s% conv_vel(j) ** 2) * s% rho(j)))
   !j = j + 1
 !end do
  !do while (s% r(k) .le. s% scale_height(j))
   !B_val(j) = 0
   !j = j + 1
 !end do
   B_max = maxval(sqrt(4 * pi * (s% conv_vel(k:j) ** 2) * s% rho(k:j)))
   B_ave = sum(sqrt(4 * pi * (s% conv_vel(k:j) ** 2) * s% rho(k:j)))/size(s% rho(k:j))
  B_region = sum(sqrt(4 * pi * (s% conv_vel(k:i) ** 2) * s% rho(k:i)))/size(s% rho(k:i))

    !edot = dot_product(s% dm(1:s% nz), s% eps_nuc(1:s% nz))

   ! the center of the star is at i = s% nz and the surface at i = 1 .
   ! so go from the center outward until 90% of the integrated eps_nuc
   ! is enclosed.  exit and then i will contain the desired cell index

   ! note: do NOT add these names to history_columns.list
   ! the history_columns.list is only for the built-in log column options.
   ! it must not include the new column names you are adding here.

   ! column 1
   names(1) = "B_max"
   vals(1) = B_max  ! in solar masses

   ! column 2
   names(2) = "B_ave"
   vals(2) = B_ave ! in solar radii

   names(3) = "B_region"
   vals(3) = B_region

   ierr = 0
end subroutine data_for_extra_history_columns

integer function how_many_extra_profile_columns(id, id_extra)
   use star_def, only: star_info
   integer, intent(in) :: id, id_extra
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_profile_columns = 1
end function how_many_extra_profile_columns


subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
   use star_def, only: star_info, maxlen_profile_column_name
   use const_def, only: dp
   integer, intent(in) :: id, id_extra, n, nz
   character (len=maxlen_profile_column_name) :: names(n)
   real(dp) :: vals(nz,n)
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   integer :: k, j, l
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   ! note: do NOT add the extra names to profile_columns.list
   ! the profile_columns.list is only for the built-in profile column options.
   ! it must not include the new column names you are adding here.

   ! here is an example for adding a profile column
   k = s% nz
   do while ( s% conv_vel(k) /= 0 )
    k = k - 1
  end do
  j = k

   if (n /= 1) stop 'data_for_extra_profile_columns'
   names(1) = 'B_equipartition'

   iloop: do l = 1,s% nz
      if (l >= k) then
      if (s% r(k) .ge. s% scale_height(j)) then
      vals(l,1) = sqrt(4 * pi * (s% conv_vel(l) ** 2) * s% rho(l))
      j = j + 1
    else if (s% r(k) .le. s% scale_height(j)) then
      vals(l,1) = vals(l-1,1)
      j = j+1
      end if
  end if
   if (l < k) then
      vals(l,1) = sqrt(4 * pi * (s% conv_vel(l) ** 2) * s% rho(l))
    end if
   if (vals(l,1) < 0) then
      print *, 'An error has occured in '
      print *, l
      stop
    end if
  end do iloop

end subroutine data_for_extra_profile_columns

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
      end subroutine extras_controls


      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         call system_clock(time0,clock_rate)
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup


      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call system_clock(time1,clock_rate)
         dt = dble(time1 - time0) / clock_rate / 60
         write(*,'(/,a50,f12.2,99i10/)') 'runtime (minutes), retries, backups, steps', &
            dt, s% num_retries, s% num_backups, s% model_number
         ierr = 0
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
      end function extras_check_model


    !  integer function how_many_extra_history_columns(id, id_extra)
    !     integer, intent(in) :: id, id_extra
    !     integer :: ierr
    !     type (star_info), pointer :: s
    !     ierr = 0
    !     call star_ptr(id, s, ierr)
    !     if (ierr /= 0) return
    !     how_many_extra_history_columns = 0
    !  end function how_many_extra_history_columns


    !  subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
    !     integer, intent(in) :: id, id_extra, n
    !     integer, intent(out) :: ierr
  !       type (star_info), pointer :: s
  !       ierr = 0
  !       call star_ptr(id, s, ierr)
  !       if (ierr /= 0) return
  !    end subroutine data_for_extra_history_columns


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_finish_step(id, id_extra)
        !implicit none
         use chem_def
         integer, intent(in) :: id, id_extra
         integer :: ierr
         integer :: a = 0, b = 0, c = 0, d = 0, e = 0 , f = 0, g = 0, h=0, i=0, j=0, k=0, l=0
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         !s% need_to_save_profiles_now = .false.
         !Lets try this with a for loop later

         if (a == 0 .and. s% center_h1 <= 0.6) then
           s% need_to_save_profiles_now = .true.
           print *, 'center h1 <= 0.6'
           a = 1
         end if
         if (b == 0 .and. s% center_h1 <= 0.5) then
           s% need_to_save_profiles_now = .true.
           print *, 'center h1 <= 0.5'
           b = 1
         end if
         if (c == 0 .and. s% center_h1 <= 0.4) then
           s% need_to_save_profiles_now = .true.
           print *, 'center h1 <= 0.4'
           c = 1
         end if
         if (d == 0 .and. s% center_h1 <= 0.3) then
           s% need_to_save_profiles_now = .true.
           print *, 'center h1 <= 0.3'
           d = 1
         end if
         if (e == 0 .and. s% center_h1 <= 0.2) then
           s% need_to_save_profiles_now = .true.
           print *, 'center h1 <= 0.2'
           e= 1
         end if
         if (f == 0 .and. s% center_h1 <= 0.1) then
           s% need_to_save_profiles_now = .true.
           print *, 'center h1 <= 0.1'
           f = 1
         end if
         if (g == 0 .and. s% center_h1 <= 0.01) then
           s% need_to_save_profiles_now = .true.
           print *, 'center h1 <= 0.01'
           g = 1
         end if
         if (h == 0 .and. s% center_he4 >= 0.98 ) then
           s% need_to_save_profiles_now = .true.
           print *, 'center he4 <= 0.98'
           h = 1
           i = 1
         end if
         if (i == 1 .and. s% center_he4 <= 0.5 ) then
           s% need_to_save_profiles_now = .true.
           print *, 'center he4 <= 0.5'
           i = 2
           j = 1
         end if
         if (j == 1 .and. s% center_he4 <= 0.01) then
           s% need_to_save_profiles_now = .true.
           print *, 'center he4 <= 0.01'
           j = 2
         end if
         if (k == 0 .and. s% center_c12 >= 0.5) then
           s% need_to_save_profiles_now = .true.
           print *, 'center c12 <= 0.5'
           k = 1
           l=1
         end if
         if (l == 1 .and. s% center_c12 <= 0.3) then
           s% need_to_save_profiles_now = .true.
           print *, 'center c12 <= 0.3'
           l = 2
         end if

      end function extras_finish_step
      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
