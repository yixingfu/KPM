! Created=Thu 20 Sep 2018 05:01:31 PM DST
! Last Modified=Thu 20 Sep 2018 05:52:23 PM DST
      ! Makes x and y operator
        allocate(XA(XYNNZ),YA(XYNNZ))
        allocate(XYcol(XYNNZ),XYrp(N+1))

      XYrp = 0
      XYcol = 0
      XA = 0
      YA = 0

      rp_ind = 1
      
      ! go through all indices up to N (not NNZ)
      if (my_id.eq.0) then
          write(*,*)"2D with spin"
      endif
      do j=1,L !y
      do i=1,L !x
      do s=0,1 !spin
        XYrp(rp_ind) = rp_ind
        XYcol(rp_ind) = rp_ind

        XA(rp_ind) = i 
        YA(rp_ind) = j
        ! this scaling is arbitrary but should be consistent for x and y


        rp_ind = rp_ind + 1
      enddo
      enddo
      enddo
      XYrp(rp_ind) = XYNNZ+1


