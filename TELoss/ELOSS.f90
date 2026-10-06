!======================================================================
!== Energy loss in material
!==
!== File: ELOSS.f90
!== Date: 08/12/2000
!======================================================================
!  ELOSS.f90 
!
!  FUNCTIONS/SUBROUTINES exported from ELOSS.dll:
!	ELOSS      - subroutine 
!
subroutine ELOSS(NEL,ZEL,AEL,WEL,DEN,ZP,AP,NE,ETAB,RE,ZW,AW)

  ! Expose subroutine ELOSS to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::ELOSS

	IMPLICIT REAL*8(A-H,O-Z)
	INTERFACE
	SUBROUTINE rde(E, R, R1, R2, IRE)
		REAL*8   E, R, R1, R2
		INTEGER  IRE
	END SUBROUTINE rde
	END INTERFACE

  ! Variables

  ! Input:
	INTEGER NEL, NE
	REAL*8  ZEL(NEL),AEL(NEL),WEL(NEL), DEN,ZP,AP, ETAB(NE)

  ! Output:
	REAL*8  RE(NE),ZW,AW

  ! Local
	COMMON adj(5),as(5),z(5),f(5),um(5),oz,ro,azm,aj,zion,aion,&
	       zmed,amed,jc

  ! Body of ELOSS

	zion = ZP
	aion = AP
	jc   = NEL
	ro   = DEN
	
	!== two parameters with question
	ire  = 1
	oz   = 1.
	DO i = 1, NEL
		as(i) = AEL(i)
		z(i)  = ZEL(i)
		um(i) = WEL(i)
		IF(z(i) .LE. 13.) THEN
			adj(i) = 12.*z(i) + 7.
		ELSE
			adj(i) = 9.76*z(i) + 58.8 / z(i)**0.19
		END IF
	END DO

	uma=0.
	zam=0.
	umed=0.
	zme=0.
	aj=0.
	DO i = 1, jc
		uma=uma+um(i)*as(i)
		umed=umed+um(i)
		zme=zme+um(i)*z(i)
	END DO

	amed=uma/umed
	zmed=zme/umed
	
	DO i = 1, jc
		f(i)=um(i)*as(i)/uma
	END DO

	DO i = 1, jc
		zam=zam+f(i)*z(i)/as(i)
	END DO

	azm=1./zam
	aj=0.
	
	DO i = 1, jc
		aj=aj+f(i)*z(i)*log(adj(i))/as(i)
	END DO
	
	aj=aj*azm
	aj=exp(aj)
	ZW = zmed
	AW = amed

	DO j = 1, NE
		!== result [mg/cm^2]
		CALL rde(ETAB(j), RE(j), rel1, rel, ire)
	END DO
end subroutine ELOSS

SUBROUTINE rde(e,range,rel1,rel,ix)
	IMPLICIT REAL*8(A-H,O-Z)

	INTERFACE 
		REAL*8 FUNCTION c(x)
			REAL*8 x
		END FUNCTION c

		REAL*8 FUNCTION c1(x)
			REAL*8 x
		END FUNCTION c1
	END INTERFACE



! calculates range and de/dx for compounds
	dimension ala(3),ala1(3)
	dimension a(3,3),a1(4,4),a2(4,4),b(2,3),cc(5)
	common	adj(5),as(5),z(5),f(5),um(5),oz,ro,azm,aj,zion,aion,&
			zmed,amed,jc
	dimension alaj1(4),altau1(4)

	data a/	-.75265, .073736, .040556, 2.5398,-.312,.018664, &
			-.24598, .11548, -.0099661/

	data a1/-8.0155,     0.36916,    -1.4307e-02,  3.4718e-03,	&
			 1.8371,    -1.452e-02,  -3.0142e-02,  2.3603e-03,	&
			 0.045233,  -9.5873e-04,  7.1303e-03, -6.8538e-04,	&
			-5.9898e-03,-5.2315e-04, -3.3802e-04,  3.9405e-05/

	data a2/-8.725,		 0.8309,	 -0.13396,	    0.012625,	&
			 1.8797,	 0.11139,	 -0.064808,		0.0054043,	&
			 0.74192,	-0.528805,	  0.1264232,   -0.00934142,	&
			 0.752005,  -0.5558937,   0.12843119,  -0.009306336/

	data b/	0.422377e-06,   3.858019e-09, 0.0304043e-06, -0.1667989e-09, &
           -0.00038106e-06, 0.00157955e-09/

	if(e.le.0.002)then
		range=0.
		rel1=0.
		rel=0.
		return
	endif

! this a flag
	ibis=-1

	en	=	e/aion
	tau	=	en*1.008

	altau	=	log(tau)
	alaj	=	log(aj)
	alaj1(1)=	1.
	altau1(1)=	1.
	DO kk = 2,4
		alaj1(kk)=alaj**(kk-1)
		altau1(kk)=altau**(kk-1)
	END DO
	t	=	tau/938.256
	tt	=	2.*t + t*t
	beta=	sqrt(tt)/(1.+t)


	s1=0.
	DO i=1,4
		DO j=1,4
			s1	=	s1 + a2(j,i)*alaj1(j)*altau1(i)
		END DO
	END DO

	ala(1) = azm*exp(s1)
	
	s2=0.
	DO i=2,4
		DO j=1,4
			s2	=	s2 + (i-1)*a2(j,i)*alaj1(j)*altau1(i-1)
		END DO
	END DO
	ala1(1)=ala(1)*s2/tau

	s1=0.
	DO i=1,4
		DO j=1,4
			s1	=	s1 + a1(j,i)*alaj1(j)*altau1(i)
		END DO
	END DO

	ala(3) = azm*exp(s1)

	s2=0.
	DO i=2,4
		DO j=1,4
			s2	=	s2 + (i-1)*a1(j,i)*alaj1(j)*altau1(i-1)
		END DO
	END DO
	ala1(3)=ala(3)*s2/tau

	s1=0.
	DO i=1,3
		DO j=1,3
			s1	=	s1	+	a(j,i)*alaj1(j)*altau1(i)
		END DO
	END DO
	ala(2)=azm*exp(s1)/1000.

	s2=0.
	DO i=2,3
		DO j=1,3
			s2	=	s2 + (i-1)*a(j,i)*alaj1(j)*altau1(i-1)
		END DO
	END DO
	ala1(2)=ala(2)*s2/tau

25	continue
	coefa=3.
	coefb=1.

	alaa=(ala(1)*(tanh(coefa*(.98 - en))+1.)/2.			&
		+ ala(2)*(tanh(coefa*(en - .98))+1.)/2.)		&
		*        (tanh(coefb*(8.0 - en))+1.)/2.			&
		+ ala(3)*(tanh(coefb*(en  - 8.))+1.)/2.


	alim1=0.
	alim2=0.
	if(coefa*(en-.98).lt.85)then
		alim1=1.008/cosh(coefa*(.98-en))/cosh(coefa*(.98-en))
	endif
	if(coefb*(en-8.).lt.85)then
		alim2=1.008/cosh(coefb*(8.-en))/cosh(coefb*(8.-en))
	endif

	alaa1=(ala1(1)*(tanh(coefa*(.98-en))+1.)/2.+		&
		ala1(2)*(tanh(coefa*(en-.98))+1.)/2.)*			&
		(tanh(coefb*(8.-en))+1.)/2.+					&
		ala1(3)*(tanh(coefb*(en-8.))+1.)/2.+			&
		coefa/2.*(tanh(coefb*(8.-en))+1.)/2.*			&
		(ala(2)*alim1-ala(1)*alim1)+					&
		coefb/2.*(ala(3)*alim2-							&
		(ala(1)*(tanh(coefa*(.98-en))+1.)/2.+			&
		ala(2)*(tanh(coefa*(en-.98))+1.)/2.)*alim2)


	hi=137.*beta/zion
	bz=(31.8+3.86*(aj**.625))*azm*.000001
	bz=bz*(zion**2.666667)*c(hi)
	bz1=(4.357+.5288*(aj**.625))*azm*.001
	bz1=bz1*(zion**1.666667)*c1(hi)
	bep=beta*beta
	rel1=zion*zion/(alaa1+bz1*((1.-bep)**1.5)/931.141/beta)
	rel1=rel1/1000.
	range=(alaa+bz)*aion/1.008/zion/zion
	range=range*1000.

! Atention!! this version do not work correctly for ix.ne.1
	if(ix.ne.1)return
	ank=.153574*ro/azm
	z23=zion**.666667
	abet=beta*125./z23
	zef=zion*(1.-exp(-abet))
	zef2=zef*zef
	om=1022000.*bep/(1.-bep)
	cbet=0.
	DO k=1,jc
		cc(k)=0.
		DO i=1,2
			DO j=1,3
				cc(k)=cc(k)+b(i,j)*((1./bep-1.)**j)*((adj(k)**(i+1)))
			END DO
		END DO
		cbet=cbet+f(k)*cc(k)/as(k)
	END DO
	cbet=cbet*azm
	del1=ank*zef2*(log(om/oz)-bep)/bep
	del1=del1/ro
	del2=2.*ank*zef2*(log(oz/aj)-cbet)/bep
	del2=del2/ro
	rel2=rel1-del1
	rel3=rel1-rel2+del2
	if(del1)8,8,9
8	rel=rel1
	goto 12
9	if(del1+del2-rel2)10,10,11
10	rel=rel3
	goto 12
11	if(del1.lt.rel1)then
	 rel=rel2
	 else
	 rel=rel1
	endif
12	return
	END


	REAL*8 FUNCTION c(x) 
	REAL*8 x
	IF(x .LE. 0.2) THEN
		c = -0.00006 + (0.05252 + 0.12857*x)*x
	ELSE IF(x .LE. 2.) THEN
		c = -0.00185 + (0.07355 + (0.07172 - 0.02723*x)*x)*x
	ELSE IF(x .LE. 3.) THEN
		c = -0.0793  + (0.3323  - (0.1234  - 0.0153*x)*x)*x
	ELSE
		c = 0.22
	END IF
	END

	REAL*8 FUNCTION c1(x)
	REAL*8 x

	IF(x .LE. 0.2) THEN
		c1 = 0.05252+.25714*x
	ELSE IF(x .LE. 2.0) THEN
		c1 = 0.07355 + (0.14344 - 0.08169*x)*x
	ELSE IF(x .LE. 3.0) THEN
		c1 = 0.3323  - (0.2468  - 0.0459*x)*x
	ELSE
		c1 = 0.
	END IF
	
	END


