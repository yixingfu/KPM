!

! N,NNZ,A,col,rp: sparse matrix information

       subroutine arpack_wrap_ends(N,NNZ,A,rp,col,&
                OUTPUTRESULT,OUTPUT_COUNT,&
                tol,maxitr) 

        integer*8:: N,NNZ
        complex*16,dimension(NNZ)::A
        integer*8,dimension(NNZ)::col
        integer*8,dimension(N+1)::rp

        integer::OUTPUT_COUNT
        complex*16,dimension(OUTPUT_COUNT)::OUTPUTRESULT
        complex*16,dimension(:,:),allocatable::v
        complex*16,dimension(:),allocatable::workl,workd,workev,d,resid
        integer::lworkl,nx,ncv,nev,nn
        character(10)::bmat,which
        real*8::tol
        integer::ido,info,ishfts,maxitr,mode
        integer,dimension(11)::iparam
        integer,dimension(14)::ipntr
        integer::ldv,ierr
        real*8,dimension(:),allocatable::rwork
        Logical::rvec
        Logical,dimension(:),allocatable::select
        complex*16::sigma
        integer::countr
        integer::i
        
           nx = N
           nn = nx*nx
           ldv = nn
           nev = OUTPUT_COUNT
           ncv = floor(2.5*OUTPUT_COUNT)
        allocate(workd(3*nn),rwork(ncv),select(ncv),d(ncv))
        d=0
        countr = 1
          bmat = 'I'
          which = 'LM'
          lworkl=3*ncv**2+5*ncv
        allocate(workl(lworkl))
        allocate(resid(nn))
        allocate(v(ldv,ncv))
        allocate(workev(3*ncv))
                  
!          tol=0.00001
          ido = 0
                info = 0  

      ishfts = 1
!      maxitr = 200
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
10   continue 
        countr = countr+1
        write(*,*)countr
           call znaupd ( ido, bmat, nn, which, nev, tol, resid, ncv,&
           v, ldv, iparam, ipntr, workd, workl, lworkl,&
           rwork,info )
!        write(*,*)"step1",ido
        if (ido .eq. -1 .or. ido .eq. 1) then
                call CSRmultVc16(N,NNZ,A,rp,col,&
                 workd(ipntr(1)), workd(ipntr(2)))
        goto 10
        endif
!        write(*,*)"step2"
                
        if ( info .lt. 0 ) then
                write(*,*)"error1", info
        endif
        rvec = .true.
             call zneupd (rvec, 'A', select, d, v, ldv, sigma,&
             workev, bmat, nn, which, nev, tol, resid, ncv,&
             v, ldv, iparam, ipntr, workd, workl, lworkl,&
             rwork, ierr)

        write(*,*)"error?",ierr
        OUTPUTRESULT = d


        return


        end subroutine arpack_wrap_ends

                
       subroutine arpack_wrap_center(N,NNZ,A,rp,col,&
                OUTPUTRESULT,OUTPUT_COUNT,&
                tol,maxitr) 

        integer*8:: N,NNZ
        complex*16,dimension(NNZ)::A
        integer*8,dimension(NNZ)::col
        integer*8,dimension(N+1)::rp

        integer::OUTPUT_COUNT
        complex*16,dimension(OUTPUT_COUNT)::OUTPUTRESULT
        complex*16,dimension(:,:),allocatable::v
        complex*16,dimension(:),allocatable::workl,workd,workev,d,resid
        integer::lworkl,nx,ncv,nev,nn
        character(10)::bmat,which
        real*8::tol
        integer::ido,info,ishfts,maxitr,mode
        integer,dimension(11)::iparam
        integer,dimension(14)::ipntr
        integer::ldv,ierr
        real*8,dimension(:),allocatable::rwork
        Logical::rvec
        Logical,dimension(:),allocatable::select
        complex*16::sigma
        integer::countr
        integer::i
        
           nx = N
           nn = nx*nx
           ldv = nn
           nev = OUTPUT_COUNT
           ncv = floor(2.5*OUTPUT_COUNT)
        allocate(workd(3*nn),rwork(ncv),select(ncv),d(ncv))
        d=0
        countr = 1
          bmat = 'I'
          which = 'LM'
          lworkl=3*ncv**2+5*ncv
        allocate(workl(lworkl))
        allocate(resid(nn))
        allocate(v(ldv,ncv))
        allocate(workev(3*ncv))
                  
!          tol=0.00001
          ido = 0
                info = 0  

      ishfts = 1
!      maxitr = 200
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
12   continue 
           call znaupd ( ido, bmat, nn, which, nev, tol, resid, ncv,&
           v, ldv, iparam, ipntr, workd, workl, lworkl,&
           rwork,info )
!        write(*,*)"step1",ido
        if (ido .eq. -1 .or. ido .eq. 1) then
                call CSRmultVc16_square_offset(N,NNZ,&
                  A,rp,col,10000d0,&
                 workd(ipntr(1)),workd(ipntr(2)))
        goto 12
        endif
!        write(*,*)"step2"
                
        if ( info .lt. 0 ) then
                write(*,*)"error1", info
        endif
        rvec = .true.
             call zneupd (rvec, 'A', select, d, v, ldv, sigma,&
             workev, bmat, nn, which, nev, tol, resid, ncv,&
             v, ldv, iparam, ipntr, workd, workl, lworkl,&
             rwork, ierr)

        write(*,*)"error?",ierr
        OUTPUTRESULT = d


        return


        end subroutine arpack_wrap_center

                
