!Created=Wed 13 Dec 2017 01:56:55 PM STD
!Last Modified=Mon 24 Sep 2018 11:34:19 AM DST

      write(inputfile_k,'(a,i4.4)')trim(inputfile)//'_',k
      open(11,file=trim(inputfile_k)//".dat",&
          iostat=ierr,status="old",access="stream")

      ! read the content of the k-th realization
      read(11,iostat=ierr) Nc
      write(*,*)"Nc=",Nc
      if (k_proc.eq.0)then
          allocate(mu_avg(0:Nc-1,0:Nc-1),mu2_avg(0:Nc-1,0:Nc-1),mu_tilde(0:Nc-1,0:Nc-1))
      endif
      read(11,iostat=ierr) norm_a,norm_b
      read(11,iostat=ierr) mu_avg,mu2_avg
      close(11,iostat=ierr)
      if (ierr .ne. 0) then 
          badfiles = badfiles+1
          close(11)
          cycle
      endif
      write(*,*)"Done reading"

      Ntilde = 2*Nc 


      include "get_kernel.f90"




        write(*,*) norm_a,norm_b


