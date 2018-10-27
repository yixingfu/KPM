! Created=Fri 26 Oct 2018 01:01:00 AM DST
! Last Modified=Fri 26 Oct 2018 11:40:56 AM DST

      ! This part calculates level statistics. 
      ! it prepares 2 arrays, one for aggregating each bin and the
      ! other keep track of the number of data points added. Bin is
      ! divided by energy (determined by LSEmax and LSBinCount) and data to
      ! be aggregated is max(gap_{i+1},gap_i)/min(gap_{i+1},gap_i)

        ! initiate the arrays and consts on first run
        if (seq_i.eq.0) then
            write(*,*) "Level stats"
            allocate(GapRatio_sum(BinCount),GapRatio_num(BinCount))
            GapRatio_sum = 0d0
            GapRatio_num = 0
        endif

        ! 
        temp_denom = (2*BinCount*LSEmax)
        do i=2,N-1
          curr_E = EigVal(i)
          bin_ind = ceiling((curr_E - LSEmax)/temp_denom
          prev_E = EigVal(i-1)
          next_E = EigVal(i+1)
          
          prev_delta = abs(curr_E-prev_E)
          next_delta = abs(curr_E-next_E)
          temp_gapratio = min(prev_delta,next_delta)&
              /max(prev_delta,next_delta)

          GapRatio_sum(bin_ind) = GapRatio_sum(bin_ind) + temp_gapratio
          GapRatio_num(bin_ind) = GapRatio_num(bin_ind) + 1

        enddo



