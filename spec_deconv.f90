!-----------------------------------------------------------------------
! Ad-hoc wiener filter method to calculate RFs
!-----------------------------------------------------------------------
module const
  implicit none 
  integer, parameter :: ncmp = 3, iz = 1, ir = 2, it = 3
  
  integer, parameter :: clen_max = 200, ntrc_max = 500, nsmp_max = 20000
  real(8), parameter :: pi = 3.1415926535897932
  real(8), parameter :: pi2 = pi * 2.d0
  
  integer, parameter :: max_eig = 2, nshift = 5
  real(8), parameter :: over = 0.75d0, pro = 4.0d0
  

end module const

!-----------------------------------------------------------------------

program main
  use const
  implicit none 
  integer :: ntrc, nsmp
  character(clen_max) :: list, out_dir, fnames(ntrc_max,ncmp)
  character(clen_max) :: stack_files(ncmp)
  real(8) :: a_gus, t0, tlen, delta, pcnt, tpre
  real(8) :: xs(nsmp_max, ntrc_max, ncmp),xn(nsmp_max, ntrc_max, ncmp)
  real(8) :: rf(nsmp_max, ntrc_max, ncmp), srf(nsmp_max, ncmp)

  complex(kind(0d0)) :: rff(nsmp_max, ntrc_max, ncmp)
  logical :: stack_flag, freq_flag, sp_flag
  
  call get_param(list, a_gus, t0, tlen, pcnt, out_dir, tpre, stack_flag, &
       & freq_flag, sp_flag)
  
  call read_list(10, list, ntrc, fnames)

  
  call read_sac(20, ntrc, t0, tlen, fnames, xs, xn, nsmp, delta, sp_flag)
  
  call calc_rf(ntrc, nsmp, a_gus, delta, pcnt, tpre, &
       & xs(1:nsmp, 1:ntrc, 1:ncmp), xn(1:nsmp, 1:ntrc, 1:ncmp), &
       & rf(1:nsmp, 1:ntrc, 1:ncmp), rff(1:nsmp, 1:ntrc, 1:ncmp), &
       & sp_flag)
  
  if (stack_flag) then
     call stack_rf(ntrc, nsmp, rf(1:nsmp, 1:ntrc, 1:ncmp), &
          & srf(1:nsmp, 1:ncmp), stack_files(1:ncmp), list)
  end if
  
  call write_sac(30, ntrc, nsmp, a_gus, delta, tpre, pcnt, &
       & rf(1:nsmp, 1:ntrc, 1:ncmp), rff(1:nsmp, 1:ntrc, 1:ncmp), &
       & fnames(1:ntrc, 1:ncmp), out_dir, &
       & stack_flag, srf(1:nsmp, 1:ncmp), stack_files(1:ncmp), &
       & freq_flag)
  
  write(*,*)"END SUCCESSFULLY"
  
  stop
end program main

!-----------------------------------------------------------------------
subroutine stack_rf(ntrc, nsmp, rf, srf, stack_files, list)
  use const
  implicit none 
  integer, intent(in) :: ntrc, nsmp
  real(8), intent(in) :: rf(nsmp, ntrc, ncmp)
  real(8), intent(out) :: srf(nsmp, ncmp)
  character(clen_max), intent(in)  :: list
  character(clen_max), intent(out) :: stack_files(ncmp)
  integer :: j, itrc, icmp
  
  ! Stack 
  srf(1:nsmp, 1:ncmp) = 0.d0
  do icmp = 1, ncmp
     do itrc = 1, ntrc
        do j = 1, nsmp
           srf(j, icmp) = srf(j, icmp) + rf(j, itrc, icmp)
        end do
     end do
     srf(1:nsmp, icmp) = srf(1:nsmp, icmp) / ntrc
  end do
  
  ! Get stack file names
  do icmp = 1, ncmp
     if (icmp == iz) then
        stack_files(icmp) = trim(list) // ".z"
     else if (icmp == ir) then
        stack_files(icmp) = trim(list) // ".r"
     else if (icmp == it) then
        stack_files(icmp) = trim(list) // ".t"
     end if
  end do

  return 
end subroutine stack_rf

!-----------------------------------------------------------------------

subroutine get_file_name(fname, suffix, out_dir, outfile)
  use const
  implicit none 
  character(clen_max), intent(in) :: fname, out_dir
  character(*), intent(in) :: suffix
  character(clen_max), intent(out) :: outfile
  integer :: ndirname

  ndirname = index(fname, "/", back = .true.)
  if (ndirname /= 0) then
     outfile = trim(out_dir) // "/" // &
          & trim(fname(ndirname+1:clen_max)) // trim(suffix)
  else
     outfile = trim(out_dir) // trim(fname) // trim(suffix)
  end if
  
  return 
end subroutine get_file_name

!-----------------------------------------------------------------------

subroutine write_sac(io_unit, ntrc, nsmp, a_gus, delta, tpre, pcnt, &
     & rf, rff, fnames, out_dir, stack_flag, srf, &
     & stack_files, freq_flag)
  use const
  implicit none 
  integer, intent(in) :: io_unit, ntrc, nsmp
  character(clen_max), intent(in) :: fnames(ntrc, ncmp), out_dir
  character(clen_max), intent(in) :: stack_files(ncmp)
  real(8), intent(in) :: a_gus, delta, pcnt, rf(nsmp, ntrc, ncmp)
  real(8), intent(in) :: srf(nsmp, ncmp), tpre
  complex(kind(0d0)), intent(in) :: rff(nsmp, ntrc, ncmp)
  logical, intent(in) :: stack_flag, freq_flag
  character(clen_max) :: outfile
  integer :: itrc, icmp, ierr, ihdr(158, ntrc, ncmp), i, nf
  
  ! read original header
  do icmp = 1, ncmp
     do itrc = 1, ntrc
        open(unit = io_unit, file = fnames(itrc, icmp), &
             & status = 'old', access = 'direct', &
             & iostat = ierr, recl = 4)
        if (ierr /= 0) then
           write(*,*)"ERROR: cannot open ", trim(fnames(itrc, icmp))
           stop
        end if
        do i = 1, 158 
           read(unit = io_unit, rec = i) ihdr(i, itrc, icmp)
        end do
        close(io_unit)
     end do
  end do
  
  ! output
  write(*,*)
  write(*,*)"--- output time-series SAC files ---"
  do icmp = 1, ncmp
     do itrc = 1, ntrc
        call get_file_name(fnames(itrc, icmp), ".rft", &
             & out_dir, outfile)
        open(unit = io_unit, file = outfile, status = "unknown", &
             & access = "direct", recl = 4, iostat = ierr)
        if (ierr /= 0) then
           write(*,*)"ERROR: cannot create ", trim(outfile)
           stop
        end if
        write(*,*)trim(outfile)    
        
        do i = 1, 158
           write(io_unit, rec = i)ihdr(i, itrc, icmp)
        end do
        write(io_unit, rec = 1)real(delta)
        write(io_unit, rec = 6)real(-tpre)
        write(io_unit, rec = 45)real(a_gus)
        write(io_unit, rec = 46)real(pcnt)
        write(io_unit, rec = 80)nsmp
        do i = 1, nsmp
           write(io_unit, rec = 158 + i)real(rf(i, itrc, icmp))
        end do
        close(io_unit)
     end do
  end do
  write(*,*)
  
  ! output spectrum in SAC real-imaginary format
  if (freq_flag) then
     nf = nsmp / 2 + 1
     do icmp = 1, ncmp
        do itrc = 1, ntrc
           call get_file_name(fnames(itrc, icmp), ".rff", &
                & out_dir, outfile)
           open(unit = io_unit, file = outfile, status = "unknown", &
                & access = "direct", recl = 4, iostat = ierr)
           if (ierr /= 0) then
              write(*,*)"ERROR: cannot create ", trim(outfile)
              stop
           end if
           write(*,*)trim(outfile)  
           do i = 1, 158
              write(io_unit, rec = i)ihdr(i, itrc, icmp)
           end do
           write(io_unit, rec = 1) real(1.0 / (nsmp * delta))
           write(io_unit, rec = 6) 0.0
           write(io_unit, rec = 86) 2 ! IFTYPE=real-imaginary
           write(io_unit, rec = 45)real(a_gus)
           write(io_unit, rec = 46)real(pcnt)
           write(io_unit, rec = 80)nf
           do i = 1, nf
              write(io_unit, rec = 158 + i) &
                   & real(real(rff(i, itrc, icmp)))
           end do
           do i = 1, nf
              write(io_unit, rec = 158 + nf + i) &
                   & real(aimag(rff(i, itrc, icmp)))
           end do
           close(io_unit)
        end do
     end do
  end if

  if (stack_flag) then
     do icmp = 1, ncmp
        call get_file_name(stack_files(icmp), ".rft", &
             & out_dir, outfile)
        open(unit = io_unit, file = outfile, status = "unknown", &
             & access = "direct", recl = 4, iostat = ierr)
        if (ierr /= 0) then
           write(*,*)"ERROR: cannot create ", trim(outfile)
           stop
        end if
        write(*,*)trim(outfile)    
        
        do i = 1, 158
           write(io_unit, rec = i)ihdr(i, 1, icmp)
        end do
        write(io_unit, rec = 1)real(delta)
        write(io_unit, rec = 6)real(-tpre)
        write(io_unit, rec = 45)real(a_gus)
        write(io_unit, rec = 46)real(pcnt)
        write(io_unit, rec = 80)nsmp
        do i = 1, nsmp
           write(io_unit, rec = 158 + i)real(srf(i, icmp))
        end do
        close(io_unit)
     end do
  end if

  return 
end subroutine write_sac

!-----------------------------------------------------------------------

subroutine calc_rf(ntrc, nsmp, a_gus, delta, pcnt, tpre, xs, xn, rf, &
     & rff, sp_flag)
  use const
  implicit none 
  include 'fftw3.f'
  integer, intent(in) :: ntrc, nsmp
  real(8), intent(inout) :: xs(nsmp, ntrc, ncmp), xn(nsmp, ntrc, ncmp)
  real(8), intent(in) :: a_gus, delta, pcnt, tpre
  logical, intent(in) :: sp_flag
  real(8), intent(out) :: rf(nsmp, ntrc, ncmp)
  complex(kind(0d0)) :: rff(nsmp, ntrc, ncmp), cfft(nsmp)
  complex(kind(0d0)) :: fs(nsmp, ntrc, ncmp)
  complex(kind(0d0)) :: fn(nsmp, ntrc, ncmp)
  real(8) :: xfft(nsmp), gauss(nsmp/2+1)
  integer(8) :: ifft
  integer :: itrc, icmp, i, j, npre
  real(8) :: df, fac_normal
  


  ! apply coine taper
  do icmp = 1, ncmp
     do itrc = 1, ntrc
        call cos_taper(xs(1:nsmp, itrc, icmp), nsmp, pcnt)
        call cos_taper(xn(1:nsmp, itrc, icmp), nsmp, pcnt)
     end do
  end do

  
  ! FFT to frequency domain
  xfft(1:nsmp) = 0.d0
  call dfftw_plan_dft_r2c_1d(ifft, nsmp, xfft, cfft, FFTW_ESTIMATE)
  do icmp = 1, ncmp
     do itrc = 1, ntrc
        ! signal
        xfft(1:nsmp) = xs(1:nsmp, itrc, icmp)
        call dfftw_execute(ifft)
        fs(1:nsmp, itrc, icmp) = cfft(1:nsmp)
        ! noise
        xfft(1:nsmp) = xn(1:nsmp, itrc, icmp)
        call dfftw_execute(ifft)
        fn(1:nsmp, itrc, icmp) =  cfft(1:nsmp)
     end do
  end do
  call dfftw_destroy_plan(ifft)

  ! deconvolution
  do icmp = 1, ncmp
     do itrc = 1, ntrc
        call decon(nsmp, fs(1:nsmp, itrc, icmp), fs(1:nsmp, itrc, iz), &
             & fn(1:nsmp, itrc, icmp), fn(1:nsmp, itrc, iz), &
             & rff(1:nsmp, itrc, icmp))
     end do
  end do
  
  ! apply Gaussian filter
  df = 1.d0 / (nsmp * delta)
  if (a_gus > 0.d0) then
     do j = 1, nsmp / 2 + 1
        gauss(j) = exp(-(pi2 * (j - 1) * df / (2.d0 * a_gus))**2)
     end do
     do icmp = 1, ncmp
        do itrc = 1, ntrc
           do j = 1, nsmp / 2 + 1
              rff(j, itrc, icmp) = rff(j, itrc, icmp) * gauss(j)
           end do
        end do
     end do
  end if

  ! FFT to time domain
  npre = nint(tpre / delta)
  if (a_gus > 0.d0) then
     fac_normal = nsmp * a_gus * delta / sqrt(pi)
  else
     fac_normal = 1.d0
  end if
  
  call dfftw_plan_dft_c2r_1d(ifft, nsmp, cfft(1:nsmp), &
       & xfft(1:nsmp), FFTW_ESTIMATE)
  
  do icmp = 1, ncmp
     do itrc = 1, ntrc
        cfft(1:nsmp) = rff(1:nsmp, itrc, icmp)
        call dfftw_execute(ifft)
        
        if (.not. sp_flag) then
           do i = 1, nsmp
              j = mod(nsmp - npre + i, nsmp)
              if (j == 0) j = nsmp
              rf(i, itrc, icmp) = xfft(j) / fac_normal
           end do
        else
           do i = 1, nsmp
              j = mod(nsmp + npre - i, nsmp)
              if (j == 0) j = nsmp
              rf(i, itrc, icmp) = -xfft(j) / fac_normal
           end do
        end if
     end do
  end do
  
  call dfftw_destroy_plan(ifft)

  return 
end subroutine calc_rf

!-----------------------------------------------------------------------

subroutine decon(n, fs_child, fs_parent, fn_child, fn_parent, rff)
  use const
  implicit none 
  integer, intent(in) ::n
  complex(kind(0d0)), intent(in) :: fs_child(n), fs_parent(n)
  complex(kind(0d0)), intent(in) :: fn_child(n)
  complex(kind(0d0)), intent(in) :: fn_parent(n)
  complex(kind(0d0)) :: rff(n), num, din, fsp, fsc, fnp, fnc
  integer :: i
  
  do i = 1, n/2 + 1
     fsc = fs_child(i)
     fsp = fs_parent(i)
     fnc = fn_child(i)
     fnp = fn_parent(i)
     num = fsc * conjg(fsp)
     din = fsp * conjg(fsp)
     din = din + fnp * conjg(fnp)
     rff(i) = num / din
  end do
  
  return 
end subroutine decon
!-----------------------------------------------------------------------

subroutine cos_taper(x,n,pcnt)
  use const, only:  pi
  implicit none 
  integer, intent(in) :: n
  real(8), intent(in) :: pcnt
  real(8), intent(inout) :: x(n)
  integer :: i, nleng
  real(8) :: fac
  
  nleng = int(n*pcnt)
  do i = 1, nleng
     fac = 0.5d0 * (1.d0 - cos((i-1)*pi/nleng))
     x(i) = x(i) * fac
     x(n-i+1) = x(n-i+1)*fac
  end do
  return 
end subroutine cos_taper

!-----------------------------------------------------------------------

subroutine read_sac(io_unit, ntrc, t0, tlen, fnames, xs, xn, nsmp, &
     & delta, sp_flag)
  use const
  implicit none 
  integer, intent(in) :: ntrc, io_unit
  logical, intent(in) :: sp_flag
  integer, intent(out) :: nsmp
  real(8), intent(in)  :: t0, tlen
  real(8), intent(out) :: delta
  real(8), intent(out) :: xs(nsmp_max,ntrc_max,ncmp)
  real(8), intent(out) :: xn(nsmp_max,ntrc_max,ncmp)
  character(clen_max), intent(in) :: fnames(ntrc_max,ncmp)
  integer :: itrc, ierr, icmp, npts
  integer :: j, istart, iend, istart2, iend2
  real :: prev_del, tmp_del, beg, tmpx(nsmp_max)

  prev_del = 0.00
  do icmp = 1, ncmp
     do itrc = 1, ntrc

        ! open file
        open(unit = io_unit, file = fnames(itrc,icmp), &
             & access = 'direct', &
             & status = 'old', recl = 4, iostat = ierr)
        if (ierr /= 0) then
           write(*,*)"ERROR: Cannot open ", trim(fnames(itrc,icmp))
           stop
        end if

        ! get delta
        read(io_unit, rec = 1)tmp_del
        if (icmp /= 1 .or. itrc /= 1) then
           if (abs(prev_del - tmp_del) > 1.0e-8) then
              write(*,*)"ERROR: invalid sampling interval"
              stop
           end if
        end if
        prev_del = tmp_del
        
        ! get file begin time
        read(io_unit, rec = 6)beg
        
        ! get npts
        read(io_unit, rec = 80)npts
        
        ! get data
        istart = nint((t0 - beg) / tmp_del) + 1
        nsmp   = nint(tlen / tmp_del)
        iend   = istart + nsmp - 1
        istart2 = istart - nsmp
        iend2   = istart2 + nsmp - 1
        if (nsmp <= 0) then
           write(*,*)"ERROR: invalid time length"
           write(*,*)"       (nsmp <= 0)"
           stop
        end if
        if (iend > npts) then
           write(*,*)"ERROR: insufficient signal length"
           stop
        end if
        if (istart2 < 1) then
           write(*,*)"ERROR: insufficient noise length"
           stop
        end if
        
        ! noise 
        do j = 1, nsmp
           read(io_unit, rec = 158 + istart2 + j - 1) tmpx(j)
        end do
        if (icmp == it .or. .not. sp_flag) then
           xn(1:nsmp,itrc,icmp) = dble(tmpx(1:nsmp))
        else if (icmp == ir .and. sp_flag) then
           xn(1:nsmp,itrc,iz) = dble(tmpx(1:nsmp))
        else if (icmp == iz .and. sp_flag) then
           xn(1:nsmp,itrc,ir) = dble(tmpx(1:nsmp))
        end if
        
        ! signal
        do j = 1, nsmp
           read(io_unit, rec = 158 + istart + j - 1) tmpx(j)
        end do
        if (icmp == it .or. .not. sp_flag) then
           xs(1:nsmp,itrc,icmp) = dble(tmpx(1:nsmp))
        else if (icmp == ir .and. sp_flag) then
           xs(1:nsmp,itrc,iz) = dble(tmpx(1:nsmp))
        else if (icmp == iz .and. sp_flag) then
           xs(1:nsmp,itrc,ir) = dble(tmpx(1:nsmp))
        end if
        
        close(io_unit)
     end do
  end do

  delta = dble(tmp_del)
  
  write(*,*)
  write(*,*)"--- read from SAC files ---"
  write(*,*)"Sampling interval (s): ", delta
  write(*,*)"Number of data in time window: ", nsmp
  write(*,*)
  
  return 
end subroutine read_sac

!-----------------------------------------------------------------------

subroutine read_list(io_unit, list, ntrc, fnames)
  use const
  implicit none 
  integer, intent(in)  :: io_unit
  integer, intent(out) :: ntrc
  character(clen_max), intent(in) :: list
  character(clen_max), intent(out) :: fnames(ntrc_max,ncmp)
  character(clen_max) :: line, zfile, rfile, tfile
  integer :: ierr, nchar, itrc
  logical :: lfile
  
  open(unit = io_unit, file = list, status = 'old', iostat = ierr)
  if (ierr /= 0) then
     write(*,*)"ERROR: cannot open ", trim(list)
     stop
  end if
  
  ! get file names
  ntrc = 0
  do 
     call get_line(io_unit, line, ierr)
     if (ierr /= 0) then 
        exit 
     end if
     ! check file extention
     nchar = len_trim(line)
     if (nchar < 2) then
        cycle 
     else if (line(nchar-1:nchar) == ".z") then
        zfile = line
        inquire(file=zfile,exist=lfile)
        if (.not. lfile) then
           write(*,*)"WARNING: Cannot find File:", zfile
           cycle
        end if
        rfile = line
        rfile(nchar-1:nchar) = '.r'
        inquire(file=rfile,exist=lfile)
        if (.not. lfile) then
           write(*,*)"WARNING: Cannot find File:", rfile
           cycle
        end if
        tfile = line
        tfile(nchar-1:nchar) = '.t'
        inquire(file=tfile,exist=lfile)
        if (.not. lfile) then
           write(*,*)"WARNING: Cannot find File:", tfile
           cycle
        end if

        ntrc = ntrc + 1
        if (ntrc > ntrc_max) then
           write(0,*)"ERROR: increase maximum trace number and " // &
                & "recompile program!"
           stop
        end if
        
        fnames(ntrc,iz) = zfile
        fnames(ntrc,ir) = rfile
        fnames(ntrc,it) = tfile
        
     end if
  end do
  close(io_unit)
  write(*,*)
  write(*,*)"--- read from list ---"
  write(*,*)(trim(fnames(itrc,iz)), itrc = 1, ntrc)
  write(*,*)(trim(fnames(itrc,ir)), itrc = 1, ntrc)
  write(*,*)(trim(fnames(itrc,it)), itrc = 1, ntrc)
  write(*,*)"Number of traces : ", ntrc
  write(*,*)
  return 
end subroutine read_list
   
!-----------------------------------------------------------------------

subroutine get_line(io_unit,line,ierr)
  implicit none 
  integer, intent(in) :: io_unit
  character(*), intent(out) :: line
  integer, intent(out) :: ierr
  integer :: i, nchar
  
  lines: do 
     read(io_unit,'(a)',iostat=ierr)line
     if (ierr /= 0) then
        return
     end if
     nchar = len_trim(line)
     chars: do i = 1, nchar
        if (line(i:i) == " ") then
           cycle chars
        else if (line(i:i) == "#") then
           cycle lines
        else
           exit chars
        end if
     end do chars
     exit lines
  end do lines
  
  return 
end subroutine get_line

!-----------------------------------------------------------------------

subroutine get_param(list, a_gus, t0, tlen, pcnt, out_dir, tpre, &
     & stack_flag, freq_flag, sp_flag)
  use const
  implicit none 
  integer :: narg, iarg
  real(8), intent(out) :: a_gus, t0, tlen, pcnt, tpre
  character(clen_max), intent(out) :: list, out_dir
  logical, intent(out) :: stack_flag, freq_flag, sp_flag
  character(clen_max) :: opt, val

  ! defuault parameter
  a_gus   = 1.d0
  tlen    = 50.d0
  t0      = 200.d0
  pcnt    = 0.05d0
  tpre    = 3.d0
  out_dir = "./"
  list    = ""
  stack_flag = .false.
  sp_flag = .false.
  
  ! read from arguments
  narg = iargc()
  do iarg = 1, narg
     call getarg(iarg,opt)
     if (opt(1:2) == "a=") then
        val = opt(3:clen_max)
        read(val,*)a_gus
     else if (opt(1:2) == "t=") then
        val = opt(3:clen_max)
        read(val,*)t0
     else if (opt(1:2) == "l=") then
        val = opt(3:clen_max)
        read(val,*)tlen
     else if (opt(1:2) == "p=") then
        val = opt(3:clen_max)
        read(val,*)pcnt
     else if (opt(1:2) == "o=") then
        out_dir = opt(3:clen_max)
     else if (opt(1:2) == "s=") then
        val = opt(3:clen_max)
        read(val,*)stack_flag
     else if (opt(1:2) == "f=") then
        val = opt(3:clen_max)
        read(val,*)freq_flag
     else if (opt(1:2) == "n=") then
        val = opt(3:clen_max)
        read(val,*)tpre
     else if (opt(1:3) == "sp=") then
        val = opt(4:clen_max)
        read(val,*)sp_flag
     else
        list = opt
     end if
  end do
  
  ! Display help message
  if (list == "") then
     write(0,*) "USAGE: spec_deconv [Input SAC file list] " // &
          & "(a=<Gaussian parameter> " // &
          & "t=<Onset time> l=<Time window length> " // &
          & "p=<fraction of cos taper> o=<Output directory> "// &
          & "n=<negative time length in RF> "// &
          & "s=<Stack flag> f=<Frequency flag> sp=<sp flag>)" 
     stop
  end if
  
  ! display 
  write(*,*)
  write(*,*)"--- read from arument (or default value) ---"
  write(*,*)"Gaussian parameter     : ", a_gus
  write(*,*)"Onset time (s)         : ", t0
  write(*,*)"Time window length (s) : ", tlen
  write(*,*)"Fraction of cos taper  : ", pcnt
  write(*,*)"Output directory       : ", trim(out_dir)
  write(*,*)"Input SAC file list    : ", trim(list)
  write(*,*)"Negative time len. (s) : ", tpre
  write(*,*)"Stack flag             : ", stack_flag
  write(*,*)"Frequency flag         : ", freq_flag
  write(*,*)"Sp flag                : ", sp_flag
  write(*,*)
  return 
end subroutine get_param
