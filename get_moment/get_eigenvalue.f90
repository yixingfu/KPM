! peeled from make_h.f90=Thu 10 May 2018 02:58:02 PM DST
! Last Modified=Tue 23 Oct 2018 06:57:35 PM DST
      ! This includes looking for exact eigenvalues (spectrum), exact
      ! states if requested, and IPR of real/momentum space if
      ! requested. Each depend on previous.
      if (ExactSpectrum .or. ExactStates .or. ExactIPR) then
          allocate(H_dense(N,N))
          call Sparse2Dense(N,NNZ,A,rp,col,H_dense)

          ! Find states?
          if (ExactStates .or. ExactIPR) then
              JOBZ = 'V'
          else
              JOBZ = 'N'
          endif


          !Query
          allocate(work(2))
          allocate(rwork(3*N-2))
          allocate(EigVal(N))
          call zheev(JOBZ,'U',int(N),H_dense,int(N),EigVal,work,&
              -1,rwork,info)
          if (info.ne.0) then
              write(*,*)"something wrong, exact diag query"
          endif
          lwork = work(1)
          write(*,*) "query result: ",lwork
          deallocate(work)
          !End of query


          ! Working on ED
          allocate(work(lwork))
          call zheev(JOBZ,'U',int(N),H_dense,int(N),EigVal,work,&
              lwork,rwork,info)
          if (info.ne.0) then
              write(*,*)"something wrong, exact diag "
          endif

          ! compare with new Lanzcos
!        call LanczosLowest(N,NNZ,A,rp,col,&
!                      1000*EIGVALCOUNT,EIGVALCOUNT,EigValLanc)
          ! empirical: 1.5*eigvalcount should be enough
          ! come back later. Not resolving the lowest ones
          deallocate(work,rwork)
          ! End of ED

          if (ExactStates) then
!        allocate(psi_test1(N),psi_test2(N))
!        do i=1,N
!        psi_test1 = H_dense(:,i)
!        call CSRmultVc16(N,NNZ,A,rp,col,psi_test1,psi_test2)
!       temp_val = sum(psi_test1/psi_test2)/N
!        write(*,*) temp_val,&
!                sqrt(sum(abs(psi_test1/psi_test2-temp_val)**2)/N)
!        enddo
              write(*,*) "saving exact states testing"
              open(62,file=trim(outputfile_final)//".eigvec",&
                  status="replace",access="stream",action="write")
              write(62) N
              write(62) dreal(H_dense),dimag(H_dense)
              close(62)
          endif

          if (ExactIPR) then
                include "get_IPR.f90"
          endif


          if (ExactSpectrum) then
              STARTPOINT=ceiling((N-EIGVALCOUNT)/real(2))
              ENDPOINT = STARTPOINT+EIGVALCOUNT-1
              EigValTot = EigValTot + EigVal(STARTPOINT:ENDPOINT)
              EigValLancTot = EigValLancTot + EigValLanc
!              do i_tmp=1,EIGVALCOUNT
!              write(*,*) EigValTot(i_tmp),'vs',EigValLancTot(i_tmp),&
!                  '@',my_id
!              enddo
!                open(62,file=trim(outputfile_final)//".eigval",&
!                        status="replace",access="stream",action="write")
!                write(62) N
!                write(62) EigVal
!                close(62)
          endif
          deallocate(EigVal)

          deallocate(H_dense)
      endif        
