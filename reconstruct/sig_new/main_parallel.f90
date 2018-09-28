! Created=Wed 13 Dec 2017 01:47:49 PM STD
! Last Modified=Mon 24 Sep 2018 11:32:20 AM DST

      program main

          implicit none
          include "mpif.h"
          include "header.h"

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

          call getarg(1,inputfile)
          call getarg(2,arg_tmp)
          read(arg_tmp,*) RLZmin
          call getarg(3,arg_tmp)
          read(arg_tmp,*) RLZmax
          call getarg(4,arg_tmp)
          read(arg_tmp,*) Noutput
          write(*,*)RLZmin,RLZmax,Noutput
          badfiles=0
        
          include "get_eps_grid.f90"
        
          do k_proc = 0,(RLZmax-RLZmin)/num_procs+1
          k = k_proc*num_procs+my_id+RLZmin
          if (k.gt.RLZmax) then
              cycle
          endif
          write(*,*)k

          !first read, and prepare Egrid
          include "read_input.f90"
          write(*,*) k, "read done, kernel prepared"
          include "get_mu_tilde.f90"
          include "get_gamma_mu.f90"


          open(13,file=trim(inputfile_k)//"out.dat",status="replace",&
              form="unformatted",access="stream")
          write(13)RLZmax+1-badfiles-RLZmin,&
              Noutput,eps_grid,gamma_mu,norm_a,norm_b
          close(13)
          open(22,file=trim(inputfile_k)//"out.txt")
          write(22,*)RLZmax+1-badfiles-RLZmin
          write(22,*)'------Noutput'
          write(22,*)Noutput
          write(22,*)'------eps_grid,gamma_mu'
          do i=1,Noutput
          write(22,*)eps_grid(i),gamma_mu(i)
          enddo
          write(22,*)'------norms'
          write(22,*)norm_a,norm_b

          close(22)
          
          enddo

          deallocate(mu_avg,mu2_avg,mu_tilde,eps_grid,gamma_mu,gJ,hm)
        call MPI_FINALIZE(ierr)

      contains 
          include "Chebyshev.f90"
          include "GammaMN.f90"
          include "StatMech.f90"
      End program main



