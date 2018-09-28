	! Input file
	character(200)::inputfile, matrixfile

	! All inputs
	integer::D,L,Nc,Rep
	real*8::W
	integer::task,RandType
	character(200)::outputfile
	character(200)::outputfile_final
	integer::fiboN,fiboM,commC! fiboN is for L, fiboM is for F_{n-m}/L
	integer,dimension(2)::FiboBasis
	real*8::setQ
	integer::Qn
	integer::REALIZATION0
	character(200)::arg_tmp
	real*8,dimension(3)::inputPhase
	real*8::BHZ_M
	integer::BHZ_SPIN
	integer::OPEN_BC_x,OPEN_BC_y,OPEN_BC_z
        integer::MODEL_TYPE
        integer,parameter::TYPE_SLNN=0,TYPE_BHZ=1,TYPE_PI_SM=2,TYPE_GRAPHENE=3
        integer,parameter::TYPE_LRH1D=4
        integer::PIECE,temp_i,temp_start,total_count,set_count
        real*8::temp_i_real
	integer::HONEYCOMB_BASIS
	integer,parameter::HC_LAT=0,HC_XY=1
	integer,parameter::HC_set=2,HC_set_theta=3
	real*8::HC_Jx,HC_Jy
        real*8,dimension(2)::HC_d1,HC_d2,HC_d3
	real*8::HC_d3_abs
        real*8,dimension(2)::HC_a1,HC_a2,HC_a
        real*8,dimension(2)::HC_b,HC_r,HC_r_,HC_r_ref
	integer::controller,SET_HEXAGON,SET_PARA
	integer::FLAKE_SHAPE
        real*8::LRH_sigma
	integer,parameter::HEXAGON=2,PARAGRAM=1
	real*8::HC_theta_in
	real*8::a1,a2,b1,b2,HC_denom

	
        



	! Options
	Logical::QP, BHZ, LIMIT_CORRELATION
	Logical::slowOPTCOND
	Logical::fixedTwist
	Logical::Inherit,SaveAll
	Logical::ExactSpectrum,ExactStates
	Logical::RandPhase
	Logical::PIFLUX,GRAPHENE,LRH1D
	Logical::Scramble, CondTensor

	! Matrix dimensions
	integer*8::N, NNZ, JNNZ,XYNNZ


	! Twist 
	
	real*8,dimension(3)::Twist,OrigTwist
	real*8,dimension(:),allocatable::TwistAll
	complex*16,dimension(3)::expTwist

	! Disorder
	real*8::WQP,Wrnd,P,Q,R,TQP,Trnd,t0,t
	real*8::k_i
	integer*8::iQ
	complex*16::t_tmp
	real*8,dimension(:),allocatable::eps
	integer::eps_ind

	! temp index (r for corresponding to rhs, l for lhs)
	integer*8::ind_r

	
	! parameters
	real*8,parameter::pi=3.1415926535897932384626433832795d0
	integer,parameter::OPTCOND=2,RHO=1,RHODER=3
	integer,parameter::RandHOP=2,RandPOT=1
	complex*16,parameter::III=dcmplx(0d0,1d0)
	
	
	integer::PerSite
	! H matrix (CSR)
	complex*16,dimension(:),allocatable::A
	integer*8,dimension(:),allocatable::col,rp
	integer*8::rp_ind,col_ind
	integer::s,s_
	! H matrix (dense)
	complex*16,dimension(:,:),allocatable::H_dense
	complex*16,dimension(:),allocatable::work
	real*8,dimension(:),allocatable::rwork,EigVal
	integer,parameter::EIGVALCOUNT=3000
	integer::STARTPOINT,ENDPOINT,i_tmp
	real*8,dimension(EIGVALCOUNT)::EigValTot,EigValTotALL
	real*8,dimension(EIGVALCOUNT)::EigValLanc,EigValLancTot
	integer::lwork,info
	character(1)::JOBZ

	! Normalization
	real*8::norm_a,norm_b,Emin,Emax,Set_Norm_a,Set_Norm_b

	! J matrix (CSR)
        integer,parameter::DIR_X=1,DIR_Y=2
        integer::Dir_a,Dir_b
	complex*16,dimension(:),allocatable::JA,JxA,JyA,JaA,JbA
	integer*8,dimension(:),allocatable::Jcol,Jrp
	integer*8,dimension(:),allocatable::Jxcol,Jxrp,Jycol,Jyrp


	! local variables
	integer::i,j,k,i_,j_,k_
	integer*8::i8,j8,k8

	! 2D QP: phase
	real*8,dimension(3)::phase
        real*8,dimension(:,:),allocatable::phase_vals,phase_all

	! moment mu
	real*8,dimension(:),allocatable::mu_tot,mu2_tot,mu,psi0R
	real*8,dimension(:),allocatable::mu_avg,mu2_avg
	complex*16,dimension(:),allocatable::psi0,psi1
	complex*16,dimension(:),allocatable::psi_tmp
	complex*16,dimension(:),allocatable::psi_j,psi_j_p,psi_j_pp
	!saveload
	integer::prevNc,j0

	! 2d moment mu2d
	complex*16,dimension(:,:),allocatable::mu2d_tot,mu2d2_tot
	complex*16,dimension(:,:),allocatable::mu2d_avg,mu2d2_avg
	complex*16,dimension(:,:),allocatable::mu2d

	! fast (high memory consumption)
	complex*16,dimension(:,:),allocatable::psi_all_out
        complex*16,dimension(:),allocatable::JaPsi,JbTmJaPsi
        complex*16,dimension(:),allocatable::TmJaPsi,TmpJaPsi,TmppJaPsi
        complex*16,dimension(:),allocatable::TnPsi,TnpPsi,TnppPsi
	integer::cond_m,cond_n

	! old
        complex*16,dimension(:),allocatable::psi_out,&
	psi0_out,psi1_out,psi_tmp_out,psi_in,&
	psi_p_in,psi_pp_in,psi_p_out,psi_pp_out
	

!!	! 2d moment reduced memory
!!	integer::proj_N
!!	complex*16,dimension(:,:),allocatable::proj_psi_all
!!	real*8,dimension(:,:),allocatable::proj_psi_all_R
!!	complex*16,dimension(:),allocatable::psi_r_rdc,psi_l_rdc
	! 2d moment save all
	
	

	! parameters not in the form of parameter
	complex*16,dimension(0:1,0:1)::pauli_x,pauli_y,pauli_z
	complex*16,dimension(0:1,0:1)::xf,xb,yf,yb,zf,zb,zz
	complex*16,dimension(0:1,0:1)::txf,txb,tyf,tyb,tzf,tzb
	complex*16,dimension(0:1,0:1)::Jtxf,Jtxb
        complex*16,dimension(0:1)::eiAx,eiAy
        complex*16::U_fwd,U_bwd
        complex*16::tx,tx_,ty,ty_
        real*8::ix,iy,jx,jy,xx,yy,ww,xx_,yy_,ABx,ABy
        real*8::real_ix,real_iy,real_jx,real_jy,real_ABx,real_ABy

	complex*16,dimension(0:1,1:3)::texp_theta,Jtexp_theta ! graphene
        integer,dimension(0:1,1:3)::NNx,NNy
        integer::NNi
        



	! MPI
	integer::status,ierr,num_procs,my_id,rlz_id
	integer::seq_rep,seq_i
	integer,parameter::twist_tag=2002
	
	! RHODER
	integer::iT
	real*8::tempTsum
	real*8,dimension(:),allocatable::orig_mu
	
	! randoms
	real*8::randomtest1,randomtest2
	real*8::ran2,psirand
	integer::idum
	integer::ipsi

	
	! 0-up; 1-down
	pauli_x = 0d0
	pauli_x(0,1) = 1d0
	pauli_x(1,0) = 1d0
	pauli_y = 0d0
	pauli_y(0,1) = -III
	pauli_y(1,0) = III
	pauli_z = 0d0
	pauli_z(0,0) = 1d0
	pauli_z(1,1) = - 1d0

	! ordinary
	xf= 0.5d0*III*pauli_x
	xb=-0.5d0*III*pauli_x

	yf= 0.5d0*III*pauli_y
	yb=-0.5d0*III*pauli_y

	zf= 0.5d0*III*pauli_z
	zb=-0.5d0*III*pauli_z
	
	zz= 0d0
