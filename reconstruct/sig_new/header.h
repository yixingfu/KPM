!Created=Wed 13 Dec 2017 01:48:29 PM STD
!Last Modified=Sun 23 Sep 2018 11:09:08 PM DST
	character(200)::inputfile,inputfile_k
	character(200)::outputfile_k
	character(200)::arg_tmp
	integer::RLZmax,RLZmin
	integer::Nc,Noutput
	real*8::norm_a,norm_b
	real*8,dimension(:,:),allocatable::mu_avg,mu2_avg,mu_tilde
	
	integer::Ntilde
	real*8,dimension(:),allocatable::eps_grid,gamma_mu
        real*8::eps,Gamma_mu_tmp
        integer::eps_ind


	
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

