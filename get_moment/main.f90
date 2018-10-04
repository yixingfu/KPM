! Created=Tue 12 Dec 2017 02:59:28 PM STD
! Last Modified=Thu 04 Oct 2018 05:45:53 PM DST
      program main
          use lapack95
          use f95_precision
          implicit none
          include "mkl.fi"
          include "mpif.h"
          include "header.h"
          my_id = 0
          call getarg(2,arg_tmp)
          read(arg_tmp,*) seq_rep
          ! initialize
          call MPI_INIT(ierr)
          call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
          call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
          ! read inputs
          include "read_input.f90"

          ! prepare (including the realization dependent)
          include "prepare.f90"

          allocate(Jxrp(N+1),JxA(JNNZ),Jxcol(JNNZ))
          allocate(Jyrp(N+1),JyA(JNNZ),Jycol(JNNZ))
          do seq_i=0,seq_rep-1

          rlz_id = REALIZATION0+my_id
          write(outputfile_final,&
              '(a,i4.4)')trim(outputfile)//"_",rlz_id

          ! make H (and normalize)
          if (BHZ) then
              include "make_h_BHZ.f90"
          else if (PIFLUX) then
              include "make_h_PI_SM.f90"
          else if (GRAPHENE) then
              include "make_h_GRAPHENE.f90"
          else if (LRH1D) then
              include "make_h_LRH1D.f90"
          else if (IQHE_SQ) then
              include "make_h_IQHE_SQ.f90"
          else
              include "make_h.f90"
          endif 
          ! check
          !call printCSR(N,NNZ,A,col,rp)


          if (task.eq.OPTCOND) then
                ! make J
                include "make_j.f90"
          endif

          ! find moment and save result
          if ((task.eq.RHO) .or. (task.eq.RHODER))then
              include "get_mu.f90"
                write(*,*) rlz_id, "done mu calculation"

              ! save results
              ! in the binary file we only need to record data needed
              ! for reconstruction. Other information should be kept
              ! in log file
              if (ExactSpectrum .or. ExactStates) then
                  write(*,*) "Exact Diag"
              else

                  open(15,file=trim(outputfile_final)//".dat",&
                      status="replace",&
                      form="unformatted",access="stream")
                  write(15) Nc
                  write(15) norm_a,norm_b
                  write(15) mu_avg,mu2_avg
                  close(15)
                  open(16,file=trim(outputfile_final)//".log",&
                      status="replace")
                  write(16,*)"D=",D,",L=",L,",Nc=",Nc,",Rep=",Rep
                  write(16,*)"W=",W,",norm_a=",norm_a,",norm_b=",norm_b
                  close(16)
              endif

          else if (task.eq.OPTCOND) then
              if (CondTensor) then
                  write(*,*) "include hall conductivity"
                  include "get_mu2d_general.f90"
              else
                  include "get_mu2d.f90"
              endif

              ! save result
              rlz_id = REALIZATION0+my_id
              write(outputfile_final,&
                  '(a,i4.4)')trim(outputfile)//"_",rlz_id
              open(18,file=trim(outputfile_final)//".dat",&
                  status="replace",&
                  form="unformatted",access="stream")
              write(18) Nc
              write(18) norm_a,norm_b
              write(18) real(mu2d_avg),imag(mu2d_avg)
              write(18) real(mu2d2_avg),imag(mu2d2_avg)
              close(18)
              open(19,file=trim(outputfile_final)//".log",&
                  status="replace")
              write(19,*)"D=",D,",L=",L,",Nc=",Nc,",Rep=",Rep
              write(19,*)"W=",W,",norm_a=",norm_a,",norm_b=",norm_b
              close(19)




          endif


          deallocate(A,rp,col)
          if (task .eq. RHO) then
              deallocate(mu_avg,mu2_avg)
          else if (task.eq.OPTCOND) then
              deallocate(mu2d_avg,mu2d2_avg)
          endif

          my_id = my_id+num_procs
          enddo

          if (ExactSpectrum) then
              call MPI_REDUCE(EigValTot, EigValTotALL,EIGVALCOUNT, &
                  MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
              if (my_id .eq. (seq_rep*num_procs)) then
                  write(*,*) "SAVING EIGEN VALUES..."
                  open(62, file=trim(outputfile_final)//".eigval",&
                      status="replace",form="unformatted",&
                      access="stream",action="write")
                  write(62) EigValTotALL
                  close(62)


                  open(63, file=trim(outputfile_final)//".eigvalr",&
                      status="replace",action="write")
                  write(63,*) EigValTotALL
                  close(63)
              endif
          endif

          deallocate(Jxrp,JxA,Jxcol)
          deallocate(Jyrp,JyA,Jycol)
          call MPI_FINALIZE(ierr)

      contains
          include "Sparse.f90"
          include "index_convert.f90"
          include "random.f90"
          include "fibonacci.f90"
          include "Lanczos.f90"
          include "Rescale.f90"
          include "open_bc.f90"
          include "safe_mod.f90"
          include "dual.f90"
          include "operators.f90"

      end program main
