!Created=Wed 13 Dec 2017 01:48:29 PM STD
!Last Modified=Sun 23 Sep 2018 11:09:08 PM DST
	character(200)::inputfile,inputfile_k
	character(200)::outputfile_k
	character(200)::arg_tmp
	integer::RLZmax,RLZmin
	integer::Nc,Noutput,ForceNc
	real*8::norm_a,norm_b
	complex*16,dimension(:,:),allocatable::mu_avg,mu2_avg,mu_tilde
	real*8,dimension(:,:),allocatable::mu_avg_real,mu_avg_imag
	real*8,dimension(:,:),allocatable::mu2_avg_real,mu2_avg_imag
	
	integer::Ntilde
	real*8,dimension(:),allocatable::eps_grid,Gamma_mu
        real*8::eps,Gamma_mu_tmp
        integer::eps_ind

	integer::my_id, num_procs,k_proc


	
	! Jackson Kernel
	real*8,dimension(:),allocatable::gJ,hm
	real*8::a_



	! parameters
	real*8,parameter::pi=3.1415926535897932384626433832795d0

	! index & stats
	integer*8::i8,j8,k8
	integer::i,j,k,m,n
	integer::badfiles
	integer::ierr

