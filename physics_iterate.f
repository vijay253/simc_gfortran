 	real*8 function peepi(vertex,main)

* Purpose:
* This function determines the kinematics in the PHOTON-NUCLEON center of mass
* frame, calculates some of the kinematical variables (s,t, and CM quantities
* in the 'main' structure), and returns the pion cross section.
*
*   output:
*	peepi		!d5sigma/dEe'dOmegae'Omegapi	(microbarn/MeV/sr^2)

	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex

* NOTE: when we refer to the center of mass system, it always refers to the
* photon-NUCLEON center of mass, not the photon-NUCLEUS!  The model gives
* the cross section in the photon-nucleon center of mass frame.

	real*8 sigma_eepi		!final cross section (returned as peepi)
	real*8 k_eq			!equivalent photon energy.
	real*8 sig219,sig,fac
	real*8 wfactor
	real*8 sigt,sigl,siglt,sigtt	!components of dsigma/dt
c	real*8 efer			!energy of target particle
	real*8 epsi			!epsilon of virtual photon
	real*8 gtpr			!gamma_t prime.
	real*8 tcos,tsin		!cos/sin of theta between ppi and q
	real*8 tfcos,tfsin		!cos/sin of theta between pfermi and q
	real*8 ppix,ppiy,ppiz		!pion momentum in lab.
	real*8 zero

	real*8 s,s_gev,t,t_gev,Q2_g			!t,s, Q2 in (GeV/c)**2
	real*8 mtar_gev

	real*8 new_x_x,new_x_y,new_x_z
	real*8 new_y_x,new_y_y,new_y_z
	real*8 new_z_x,new_z_y,new_z_z
	real*8 dummy,p_new_x,p_new_y
	real*8 davesig,phipq
	real*8 square_root,dp_dcos_num,dp_dcos_den,dp_dcos
	real*8 dp_dphi_num,dp_dphi_den,dp_dphi
	real*8 dt_dcos_lab,dt_dphi_lab,psign
	real*8 dpxdphi,dpydphi,dpxdcos,dpydcos,dpzdcos,dpzdphi
	real*8 dpxnewdphi,dpynewdphi,dpxnewdcos,dpynewdcos
	real*8 dpznewdphi,dpznewdcos
	real*8 dphicmdphi,dphicmdcos

	real*8 pbeam,beam_newx,beam_newy,beam_newz
c	real*8 pbeamcmx,pbeamcmy,pbeamcmz,ebeamcm,pbeamcm
	real*8 ppicm_newx,ppicm_newy,ppicm_newz

	real*8 dEcmdcos,dEcmdphi,dcoscmdcos,dcoscmdphi
c	real*8 tmp
	real*8 q2_set,tav,ftav,ft

	logical first

c Variables calculated in transformation to gamma-NUCLEON center of mass.
        real*8 gstar,bstar,bstarx,bstary,bstarz		!beta of boost to C.M.
        real*8 nustar,qstar,qstarx,qstary,qstarz	!q in C.M.
        real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz		!p_hadron in C.M.
        real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
        real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz	!p_fermi in C.M.
        real*8 thetacm,phicm,phiqn,jacobian,jac_old
        real*8 pm2_tmp
	real*8 s_fer,t_old,thetacm_fer,phicm_fer,davejac_fer

*=====================================================================
c       Fit parameters.
	integer npar,ipar
	parameter (npar=12)	!number of fit parameters for H, pi+ and D, pi-
	real*8 fitpar(npar),par,par_er
	save fitpar
	logical first_call
	save first_call

	data first_call/.true./
*=====================================================================

* Initialize some stuff.
	Q2_g = vertex%Q2/1.d6
c NOTE: phipq calculation in event.f reverted to original version.
	phipq= main%phi_pq
	mtar_gev = targ%Mtar_struck/1000.

* calculate energy of initial (struck) nucleon, using the same assumptions that
* go into calculating the pion angle/momentum (in event.f).  For A>1, the struck
* nucleon is off shell, the 2nd nucleon (always a neutron) is on shell, and has
* p = -p_fermi, and any additional nucleons are at rest.
* NOTE pfer, efer appear to be in MeV at this point.
	efer = sqrt(pfer**2+targ%Mtar_struck**2)
	if(doing_deutpi.or.doing_hepi) then
	  efer = targ%M-sqrt(mn**2+pfer**2)
	  if(doing_hepi)efer=efer-mp
c	  mtar_offshell = sqrt(efer**2-pfer**2)
	endif

* calculate some kinematical variables
* 'f' and 'fer' indicate fermi momenta. 'star' or 'cm' indicate CM system

	tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y
     1         +vertex%up%z*vertex%uq%z
	if(tcos-1..gt.0..and.tcos-1..lt.1.d-8)tcos=1.0
	tsin=sqrt(1.-tcos**2)

	tfcos = pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z
	if(tfcos-1..gt.0..and.tfcos-1..lt.1.d-8)tfcos=1.0
	tfsin=sqrt(1.-tfcos**2)

	epsi = 1./(1.+2*(1.+vertex%nu**2/vertex%Q2)*tan(vertex%e%theta/2.)**2)

	s_fer = (vertex%nu+efer)**2-(vertex%q+pfer*tfcos)**2-(pfer*tfsin)**2
	ntup%W_fer=sqrt(s_fer)/1.d3      

c GH:
c calculate with target nucleon at rest (as in experimental replay)
	s = (vertex%nu+targ%Mtar_struck)**2-vertex%q**2
	main%wcm = sqrt(s)
	s_gev = s/1.d6

	t_old = vertex%Q2-Mh2+2.*vertex%nu*vertex%p%E-2.*vertex%p%P*vertex%q*tcos
	ntup%t_old=t_old/1.d6

c GH:
c calculate t the same way as in experimental replay
	pm2_tmp = (vertex%q*vertex%uq%x-vertex%p%P*vertex%up%x)**2 +
     1            (vertex%q*vertex%uq%y-vertex%p%P*vertex%up%y)**2 +
     2            (vertex%q*vertex%uq%z-vertex%p%P*vertex%up%z)**2
        t = -(mn-mp)**2 +2*targ%Mtar_struck*
     1          ( sqrt(targ%Mrec_struck**2+pm2_tmp)-targ%Mrec_struck )

	t_gev = t/1.d6			!CONVERT TO (GeV/c)**2
	main%t = t

*******************************************************************************
c OLD VERSION WHERE TARGET NUCLEON HAS FERMI MOMENTUM
* Calculate velocity of PHOTON-NUCLEON C.M. system in the lab frame. Use beta
* and gamma of the cm system (bstar and gstar) to transform particles into
* c.m. frame.  Define z along the direction of q, and x to be along the
* direction of the pion momentum perpendicular to q.

* DJG: Get pfer components in the lab "q" system
	dummy=sqrt((vertex%uq%x**2+vertex%uq%y**2)*
     1             (vertex%uq%x**2+vertex%uq%y**2+vertex%uq%z**2))
	new_x_x = -(-vertex%uq%x)*vertex%uq%z/dummy
	new_x_y = -(-vertex%uq%y)*vertex%uq%z/dummy
	new_x_z = ((-vertex%uq%x)**2 + (-vertex%uq%y)**2)/dummy

	dummy   = sqrt(vertex%uq%x**2 + vertex%uq%y**2)
	new_y_x =  (-vertex%uq%y)/dummy
	new_y_y = -(-vertex%uq%x)/dummy
	new_y_z =  0.0

	p_new_x = pfer*((-pferx)*new_x_x + (-pfery)*new_x_y + pferz*new_x_z)
	p_new_y = pfer*((-pferx)*new_y_x + (-pfery)*new_y_y + pferz*new_y_z)

	if(p_new_x.eq.0.)then
	  phiqn=0.
	else
	  phiqn = atan2(p_new_y,p_new_x)
	endif
	if(phiqn.lt.0.)phiqn = phiqn+2.*pi

* get beam in "q" system.

	pbeam = sqrt(vertex%Ein**2-me**2)
	beam_newx = pbeam*new_x_z
	beam_newy = pbeam*new_y_z
	beam_newz = pbeam*vertex%uq%z

	bstar=sqrt((vertex%q+pfer*tfcos)**2+(pfer*tfsin)**2)/(efer+vertex%nu)
	gstar=1./sqrt(1. - bstar**2)

	bstarz = (vertex%q+pfer*tfcos)/(efer+vertex%nu)
	bstarx = p_new_x/(efer+vertex%nu)
	bstary = p_new_y/(efer+vertex%nu)

* DJG: Boost virtual photon to CM.

	zero =0.d0
	call loren(gstar,bstarx,bstary,bstarz,vertex%nu,
     >		zero,zero,vertex%q,nustar,qstarx,qstary,qstarz,qstar)

* DJG: Boost pion to CM.
	
	ppiz = vertex%p%P*tcos
	ppix = vertex%p%P*tsin*cos(phipq)
	ppiy = vertex%p%P*tsin*sin(phipq)
	call loren(gstar,bstarx,bstary,bstarz,vertex%p%E,
     >		ppix,ppiy,ppiz,epicm,ppicmx,ppicmy,ppicmz,ppicm)
	thetacm = acos((ppicmx*qstarx+ppicmy*qstary+ppicmz*qstarz)/ppicm/qstar)
	main%pcm = ppicm

* DJG Boost the beam to CM.

	call loren(gstar,bstarx,bstary,bstarz,vertex%Ein,beam_newx,
     >		beam_newy,beam_newz,ebeamcm,pbeamcmx,pbeamcmy,pbeamcmz,pbeamcm)

* Thetacm is defined as angle between ppicm and qstar.
* To get phicm, we need out of plane angle relative to scattering plane
* (plane defined by pbeamcm and qcm).  For stationary target, this plane
* does not change.  In general, the new coordinate system is defined such
* that the new y direction is given by (qcm x pbeamcm) and the new x
* is given by (qcm x pbeamcm) x qcm.

	dummy = sqrt((qstary*pbeamcmz-qstarz*pbeamcmy)**2+
     >		(qstarz*pbeamcmx-qstarx*pbeamcmz)**2
     >		+(qstarx*pbeamcmy-qstary*pbeamcmx)**2)
	new_y_x = (qstary*pbeamcmz-qstarz*pbeamcmy)/dummy
	new_y_y = (qstarz*pbeamcmx-qstarx*pbeamcmz)/dummy
	new_y_z = (qstarx*pbeamcmy-qstary*pbeamcmx)/dummy

	dummy = sqrt((new_y_y*qstarz-new_y_z*qstary)**2
     >		+(new_y_z*qstarx-new_y_x*qstarz)**2
     >		+(new_y_x*qstary-new_y_y*qstarx)**2)
	new_x_x = (new_y_y*qstarz-new_y_z*qstary)/dummy
	new_x_y = (new_y_z*qstarx-new_y_x*qstarz)/dummy
	new_x_z = (new_y_x*qstary-new_y_y*qstarx)/dummy

	new_z_x = qstarx/qstar
	new_z_y = qstary/qstar
	new_z_z = qstarz/qstar

	ppicm_newx = ppicmx*new_x_x + ppicmy*new_x_y + ppicmz*new_x_z
	ppicm_newy = ppicmx*new_y_x + ppicmy*new_y_y + ppicmz*new_y_z
	ppicm_newz = ppicmx*new_z_x + ppicmy*new_z_y + ppicmz*new_z_z

	phicm = atan2(ppicm_newy,ppicm_newx)
	if(phicm.lt.0.) phicm = 2.*3.141592654+phicm

	thetacm_fer = thetacm
	ntup%thetacm_fer=thetacm_fer
	phicm_fer = phicm
	ntup%phicm_fer=phicm_fer

*******************************************************************************
c NEW VERSION WHERE TARGET NUCLEON IS AT REST (AS IN EXPERIMENTAL REPLAY)
* Calculate velocity of PHOTON-NUCLEON C.M. system in the lab frame. Use beta
* and gamma of the cm system (bstar and gstar) to transform particles into
* c.m. frame.  Define z along the direction of q, and x to be along the
* direction of the pion momentum perpendicular to q.

	dummy=sqrt((vertex%uq%x**2+vertex%uq%y**2) * 
     1             (vertex%uq%x**2+vertex%uq%y**2+vertex%uq%z**2))
	new_x_x = -(-vertex%uq%x)*vertex%uq%z/dummy
	new_x_y = -(-vertex%uq%y)*vertex%uq%z/dummy
	new_x_z = ((-vertex%uq%x)**2 + (-vertex%uq%y)**2)/dummy

	dummy   = sqrt(vertex%uq%x**2 + vertex%uq%y**2)
	new_y_x =  (-vertex%uq%y)/dummy
	new_y_y = -(-vertex%uq%x)/dummy
	new_y_z =  0.0

* get beam in q+nucleon CM system.

	pbeam = sqrt(vertex%Ein**2-me**2)
	beam_newx = pbeam*new_x_z
	beam_newy = pbeam*new_y_z
	beam_newz = pbeam*vertex%uq%z

	bstar=vertex%q/(vertex%nu+targ%Mtar_struck)
	gstar=1./sqrt(1. - bstar**2)

	zero =0.d0
	bstarx = zero
	bstary = zero
	bstarz = vertex%q/(targ%Mtar_struck+vertex%nu)

* DJG: Boost virtual photon to q+nucleon CM.

	call loren(gstar,bstarx,bstary,bstarz,vertex%nu,
     >		zero,zero,vertex%q,nustar,qstarx,qstary,qstarz,qstar)

* DJG: Boost pion to q+nucleon CM.
	
	ppiz = vertex%p%P*tcos
	ppix = vertex%p%P*tsin*cos(phipq)
	ppiy = vertex%p%P*tsin*sin(phipq)
	call loren(gstar,bstarx,bstary,bstarz,vertex%p%E,
     >		ppix,ppiy,ppiz,epicm,ppicmx,ppicmy,ppicmz,ppicm)
	thetacm = acos((ppicmx*qstarx+ppicmy*qstary+ppicmz*qstarz)/ppicm/qstar)
	main%pcm = ppicm

* DJG Boost the beam to q+nucleon CM.

	call loren(gstar,bstarx,bstary,bstarz,vertex%Ein,beam_newx,
     >		beam_newy,beam_newz,ebeamcm,pbeamcmx,pbeamcmy,pbeamcmz,pbeamcm)

* Thetacm is defined as angle between ppicm and qstar.
* To get phicm, we need out of plane angle relative to scattering plane
* (plane defined by pbeamcm and qcm).  For stationary target, this plane
* does not change.  In general, the new coordinate system is defined such
* that the new y direction is given by (qcm x pbeamcm) and the new x
* is given by (qcm x pbeamcm) x qcm.

	dummy = sqrt((qstary*pbeamcmz-qstarz*pbeamcmy)**2+
     >		(qstarz*pbeamcmx-qstarx*pbeamcmz)**2
     >		+(qstarx*pbeamcmy-qstary*pbeamcmx)**2)
	new_y_x = (qstary*pbeamcmz-qstarz*pbeamcmy)/dummy
	new_y_y = (qstarz*pbeamcmx-qstarx*pbeamcmz)/dummy
	new_y_z = (qstarx*pbeamcmy-qstary*pbeamcmx)/dummy

	dummy = sqrt((new_y_y*qstarz-new_y_z*qstary)**2
     >		+(new_y_z*qstarx-new_y_x*qstarz)**2
     >		+(new_y_x*qstary-new_y_y*qstarx)**2)
	new_x_x = (new_y_y*qstarz-new_y_z*qstary)/dummy
	new_x_y = (new_y_z*qstarx-new_y_x*qstarz)/dummy
	new_x_z = (new_y_x*qstary-new_y_y*qstarx)/dummy

	new_z_x = qstarx/qstar
	new_z_y = qstary/qstar
	new_z_z = qstarz/qstar

	ppicm_newx = ppicmx*new_x_x + ppicmy*new_x_y + ppicmz*new_x_z
	ppicm_newy = ppicmx*new_y_x + ppicmy*new_y_y + ppicmz*new_y_z
	ppicm_newz = ppicmx*new_z_x + ppicmy*new_z_y + ppicmz*new_z_z

	phicm = atan2(ppicm_newy,ppicm_newx)
	if(phicm.lt.0.) phicm = 2.*3.141592654+phicm

	main%thetacm = thetacm
	main%phicm = phicm

c	write(6,*)'  '
c 	write(6,*)' pfer ',pfer
c	write(6,*)' t ',t,t_old
c	write(6,*)' s ',s,s_fer
c	write(6,*)' thetacm ',thetacm*180./3.14159,thetacm_fer*180./3.14159
c	write(6,*)' phicm ',phicm*180./3.14159,phicm_fer*180./3.14159,phipq*180./3.14159

*******************************************************************************
* Read fit parameters when first called.

	   if(first_call) then

	      first_call=.false.
	
	      if(doing_pion) then
		 open(88,file='par.pl',status='old')
	      else
		 open(88,file='par.mn',status='old')
	      end if

	      do while(.true.)
		 read(88,*,end=99) par,par_er,ipar
		 fitpar(ipar)=par
	      end do
 99	      close(88)

	      print*,doing_pion
	      do ipar=1,npar
		 print*,fitpar(ipar),ipar
	      end do

	   end if		!first_call.

* Models for sigL, sigT, sigLT, sigTT  for Deuterium.

***
* Parameterization revised for IT26, 12.11.09
c	   q2_set=2.45
	   q2_set=0.375
	   tav=(0.0735+0.028*log(q2_set))*q2_set
	   ftav=(abs(t_gev)-tav)/tav
	   ft=t_gev/(abs(t_gev)+0.139570**2)**2

	   sigl=(fitpar(1)+fitpar(2)*log(Q2_g))
     1           *exp((fitpar(3)+fitpar(4)*log(Q2_g))*(abs(t_gev)-0.2))
c	   print*,"sigl", sigl
	   sigt=fitpar(5)+fitpar(6)*log(Q2_g)
     1           +(fitpar(7)+fitpar(8)*log(Q2_g))*ftav
c	   print*,"sigt", sigt
	   siglt=((fitpar(9)*exp(fitpar(10)*abs(t_gev))
     1           +fitpar(11)/abs(t_gev))*sin(thetacm))
c	   print*,"siglt", siglt
	   sigtt=-((fitpar(12)/Q2_g*exp(-Q2_g))*ft*sin(thetacm)**2)
c	   print*,"sigtt", sigtt
	   tav=(-0.178+0.315*log(Q2_g))*Q2_g
	   
	   sig219=(sigt+main%epsilon*sigl+main%epsilon*cos(2.*phicm)*sigtt
     >		+sqrt(2.0*main%epsilon*(1.+main%epsilon))*cos(phicm)*siglt)/1.d0
	  
c	   print*,"sig", sig219

c now convert to different W
c W dependence given by 1/(W^2-M^2)^2
c factor 15.333 is value of (w**2-ami**2)**2 at W=2.19

c	  wfactor=15.333/(s-(targ.Mtar_pion/1000.)**2)**2
*	  wfactor=8.539/(s-(targ.Mtar_pion/1000.)**2)**2

c	  wfactor=1.D0/(s-(targ.Mtar_pion/1000.)**2)**2
	  wfactor=1.D0/(s_gev-mtar_gev**2)**2
	  sig=sig219*wfactor
	  sigl=sigl*wfactor
	  sigt=sigt*wfactor
	  sigtt=sigtt*wfactor
	  siglt=siglt*wfactor

C--->Debug
c 	  write(*,*) 's =',s
c 	  write(*,*) 'wfactor =',wfactor
c 	  write(*,*) 'sig =',sig
c 	  write(*,*) 'sigL  =',sigL
c 	  write(*,*) 'sigT  =',sigT
c 	  write(*,*) 'sigLT =',sigLT
c 	  write(*,*) 'sigTT =',sigTT
c 	  write(*,*) '-----------------------------------------------------'
C--->Debug

	  sig=sig/2./pi/1.d+06   !dsig/dtdphicm in microbarns/MeV**2/rad

	  ntup%dsigdt = sig
C--->Debug
c	  write(*,*) 'sig =',sig
c	  write(*,*) '====================================================='
C--->Debug

*******************************************************************************
* sigma_eepi is two-fold C.M. cross section: d2sigma/dt/dphi_cm [ub/MeV**2/rad]
* Convert from dt dphi_cm --> dOmega_lab using 'jacobian' [ub/sr]
* Convert to 5-fold by multiplying by flux factor, gtpr [1/MeV]
* to give d5sigma/dOmega_pi/dOmega_e/dE_e [ub/Mev/sr].
*******************************************************************************
c OLD VERSION WHERE TARGET NUCLEON HAS FERMI MOMENTUM
* DJG We need a full blown Jacobian with target not at rest:
* DJG | dt/dcos(theta) dphicm/dphi - dt/dphi dphicm/dcos(theta) |
* DJG: Calculate dt/d(cos(theta)) and dt/dphi for the general case.
c
c	psign = cos(phiqn)*cos(phipq)+sin(phiqn)*sin(phipq)
c
c	square_root = vertex%q + pfer*tfcos - vertex%p%P*tcos
c	dp_dcos_num = vertex%p%P + (vertex%p%P**2*tcos -
c     >		psign*pfer*vertex%p%P*tfsin*tcos/tsin)/square_root
c	dp_dcos_den = ( (vertex%nu+efer-vertex%p%E)*vertex%p%P/vertex%p%E +
c     >		vertex%p%P*tsin**2-psign*pfer*tfsin*tsin )/square_root - tcos
c	dp_dcos = dp_dcos_num/dp_dcos_den
c
c	dp_dphi_num = pfer*vertex%p%P*tsin*tfsin*(cos(phiqn)*sin(phipq)-
c     >		sin(phiqn)*cos(phipq))/square_root
c	dp_dphi_den = tcos + (pfer*tsin*tfsin*psign - vertex%p%P*tsin**2
c     >		- (vertex%nu+efer-vertex%p%E)*vertex%p%P/vertex%p%E)/square_root
c	dp_dphi = dp_dphi_num/dp_dphi_den
c
c	dt_dcos_lab = 2.*(vertex%q*vertex%p%P +
c     >		(vertex%q*tcos-vertex%nu*vertex%p%P/vertex%p%E)*dp_dcos)
c
c	dt_dphi_lab = 2.*(vertex%q*tcos-vertex%nu*vertex%p%P/vertex%p%E)*dp_dphi
c
c* DJG: Now calculate dphicm/dphi and dphicm/d(cos(theta)) in the most
c* DJG: excruciating way possible.
c
c        dpxdphi = vertex%p%P*tsin*(-sin(phipq)+(gstar-1.)*bstarx/bstar**2*
c     >		(bstary*cos(phipq)-bstarx*sin(phipq)) ) +
c     >		( (ppicmx+gstar*bstarx*vertex%p%E)/vertex%p%P -
c     >		gstar*bstarx*vertex%p%P/vertex%p%E)*dp_dphi
c        dpydphi = vertex%p%P*tsin*(cos(phipq)+(gstar-1.)*bstary/bstar**2*
c     >		(bstary*cos(phipq)-bstarx*sin(phipq)) ) +
c     >		( (ppicmy+gstar*bstary*vertex%p%E)/vertex%p%P -
c     >		gstar*bstary*vertex%p%P/vertex%p%E)*dp_dphi
c	dpzdphi =  vertex%p%P*(gstar-1.)/bstar**2*bstarz*tsin*
c     >		(bstary*cos(phipq)-bstarx*sin(phipq)) +
c     > 		((ppicmz+gstar*bstarz*vertex%p%E)/vertex%p%P-
c     >		gstar*bstarz*vertex%p%P/vertex%p%E)*dp_dphi
c
c	dpxdcos = -vertex%p%P*tcos/tsin*(cos(phipq)+(gstar-1.)*bstarx/bstar**2*
c     >		(bstarx*cos(phipq)+bstary*sin(phipq)-bstarz*tsin/tcos)) +
c     >		( (ppicmx+gstar*bstarx*vertex%p%E)/vertex%p%P -
c     >		gstar*bstarx*vertex%p%P/vertex%p%E)*dp_dcos
c	dpydcos = -vertex%p%P*tcos/tsin*(sin(phipq)+(gstar-1.)*bstary/bstar**2*
c     >		(bstarx*cos(phipq)+bstary*sin(phipq)-bstarz*tsin/tcos)) +
c     >		( (ppicmy+gstar*bstary*vertex%p%E)/vertex%p%P -
c     >		gstar*bstary*vertex%p%P/vertex%p%E)*dp_dcos
c	dpzdcos = vertex%p%P*(1.-(gstar-1.)/bstar**2*bstarz*tcos/tsin*
c     >		(bstarx*cos(phipq)+bstary*sin(phipq)-tsin/tcos*bstarz))
c     >		+((ppicmz+gstar*bstarz*vertex%p%E)/vertex%p%P-gstar*bstarz*
c     >		vertex%p%P/vertex%p%E)*dp_dcos
c
c	dpxnewdphi = dpxdphi*new_x_x+dpydphi*new_x_y+dpzdphi*new_x_z
c	dpynewdphi = dpxdphi*new_y_x+dpydphi*new_y_y+dpzdphi*new_y_z
c	dpznewdphi = dpxdphi*new_z_x+dpydphi*new_z_y+dpzdphi*new_z_z
c
c	dphicmdphi = (dpynewdphi*ppicm_newx - ppicm_newy*dpxnewdphi)/
c     >			(ppicm_newx**2+ppicm_newy**2)
c
c	dpxnewdcos = dpxdcos*new_x_x+dpydcos*new_x_y+dpzdcos*new_x_z
c	dpynewdcos = dpxdcos*new_y_x+dpydcos*new_y_y+dpzdcos*new_y_z
c	dpznewdcos = dpxdcos*new_z_x+dpydcos*new_z_y+dpzdcos*new_z_z
c
c	dphicmdcos = (dpynewdcos*ppicm_newx - ppicm_newy*dpxnewdcos)
c     >			/(ppicm_newx**2+ppicm_newy**2)
c
c	dcoscmdcos = dpznewdcos/ppicm-ppicm_newz*epicm*dEcmdcos/ppicm**1.5
c	dcoscmdphi = dpznewdphi/ppicm-ppicm_newz*epicm*dEcmdphi/ppicm**1.5
c
c	jacobian = abs(dt_dcos_lab*dphicmdphi-dt_dphi_lab*dphicmdcos)
c	davejac_fer = jacobian
c
c	main%johnjac = 2*(efer-2*pferz*pfer*vertex%p%E/vertex%p%P*tcos)*
c     >		(vertex%q+pferz*pfer)*vertex%p%P /
c     >		( efer+vertex%nu-(vertex%q+pferz*pfer)
c     1           *vertex%p%E/vertex%p%P*tcos )- 2*vertex%p%P*pfer
c
*******************************************************************************
c NEW VERSION WHERE TARGET NUCLEON IS AT REST (AS IN EXPERIMENTAL REPLAY)

c	square_root = vertex%q - vertex%p%P*tcos
c	dp_dcos_num = vertex%p%P + vertex%p%P**2*tcos/square_root
c	dp_dcos_den = ( (vertex%nu+targ%Mtar_struck-vertex%p%E)*vertex%p%P/vertex%p%E +
c     >		vertex%p%P*tsin**2)/square_root - tcos
c	dp_dcos = dp_dcos_num/dp_dcos_den
c
c	dt_dcos_lab = 2.*(vertex%q*vertex%p%P +
c     >		(vertex%q*tcos-vertex%nu*vertex%p%P/vertex%p%E)*dp_dcos)

	dt_dcos_lab = 2.*(vertex%q*vertex%p%P*targ%Mtar_struck) / 
     1            ( targ%Mtar_struck+vertex%nu
     2              -vertex%q*tcos*(vertex%p%E/vertex%p%P) )

	jacobian = abs(dt_dcos_lab)
	main%davejac = jacobian

c	write(6,*)' jac ',davejac_fer,jacobian

*******************************************************************************
* Get photon flux factor (two options, see comments below).
*
* DJG,2000: Replace targ.Mtar_struck in denominator of gammaflux with more 
* general efer-pfer*tfcos, for pfer =0 this reverts to old form
*	k_eq = (s-targ%Mtar_struck**2)/2./(efer-pfer*tfcos)
*
* JRA,2001: Go back to original version - more consistent with phase space used
* in the subroutine (according to DJG - see gaskell_model.ps)
*	k_eq = (main%wcm**2-targ%Mtar_struck**2)/2./targ%Mtar_struck

* Note that there is an additional factor 'fac' included with gtpr.   This
* takes into account pieces in the flux factor that are neglected (=1) in
* colinear collisions.  The flux factor is |v_1-v_2| * 2E_1 * 2E_2.
* For a stationary target, v_2=0 and so velocity term is v_1=1 (electron
* beam), and E_2=M_2.  For collinear boost, the flux factor can be expressed
* in a way that is lorenz invariant, and so can be used for lab or C.M.
* For a NON-COLLINEAR boost, there are two changes.  First, the |v| term
* becomes 1 - (z component of pfer)/efer.  Second, E_2 isn't just the mass,
* it becomes E_fermi, so we have to remove targ.Mtar_struck (which is used
* for E_2 by default) and replace it with efer.  Since the flux factor 
* comes in the denominator, we replace the usual flux factor (gtpr) with
* gtpr*fac, where fac = 1/ ( (1-pfer_z/efer)* (efer/mtar_struck) ).
*
*       fac = 1./(1.-pferz*pfer/efer) * targ%Mtar_struck/efer
*	gtpr = alpha/2./(pi**2)*vertex%e%E/vertex%Ein*k_eq/vertex%Q2/(1.-epsi)
*
*	davesig = gtpr*fac*sig*jacobian
*
*******************************************************************************
* DJG: Replace targ.Mtar_pion in denominator of gammaflux with more general
* DJG: efer-pfer*tfcos, for pfer =0 this reverts to old form

c	gtpr = alpha/2./(pi**2)*vertex%e%E/vertex%Ein*(s_gev-mtar_gev**2)/2./
c     >		(efer-pfer*tfcos)/Q2_g/(1.-epsi)

	gtpr = alpha/2./(pi**2)*vertex%e%E/vertex%Ein*(s_gev-mtar_gev**2)/2./
     >		(targ%Mtar_struck)/Q2_g/(1.-epsi)

	davesig = gtpr*sig*jacobian
*******************************************************************************

	sigma_eepi = davesig
	peepi = sigma_eepi

	ntup%sigcm = sigma_eepi		!sig_cm
c	write(6,*)' 1 ',jacobian,thetacm,phicm,ppicm
c	write(6,*)'   ',efer,pfer
c       write(6,*)'   ',gtpr,sig,sigma_eepi
c	write(6,*)'   ',sigl,sigt,siglt,sigtt,wfactor

202	format(/11X,f5.1/)
203	format(11X,f5.0)
204	format(6(/9X,7f8.3))

	return
	end
