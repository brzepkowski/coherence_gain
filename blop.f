	    program jeden
	    implicit none
	    real wynik,temp1,temp2,tau,wyn3,dtau,w
	    integer n
	    external w

	    open(1,file='blop0.dat',status='new')

	
	    tau=0.
        dtau=0.005
    
        temp1=0.000001
	    temp2=1.

	    do n=1,1000
	    wynik=0.
	    wyn3=0.
	    wynik=w(tau,temp1)
	    wyn3=w(tau,temp2)
	    write(1,*) tau,wynik,wyn3
	    tau=tau+dtau
	    enddo


	    close(1)

	    end

c	delta rho

    	real function w(tau,temp)
    	implicit none
    	real tau,temp,dampw2,fazaw2
    	external dampw2,fazaw2

    	w=0.5*exp(dampw2(temp))*abs(sin(fazaw2(tau)))

    	end

c	damping

    	real function dampw2(temp)
    	implicit none
    	real dampw1,omega,wynik
    	real zero,temp,tempdatad2
    	integer inf
	    external dampw1
	    parameter(zero=0.,inf=1)
	    common /datad2/ tempdatad2


    	tempdatad2=temp

    	call calkai(dampw1,zero,inf,wynik)
	    dampw2=wynik
    	end

	    real function dampw1(omega)
	    implicit none
	    real omega,dampw,tempdatad2
	    external dampw
	    common /datad2/ tempdatad2

    	dampw1=dampw(omega,tempdatad2)
    	end
	

    	real function dampw(omega,temp)
	    implicit none
    	real omega,g,temp,nb
    	external g,nb

    	dampw=-g(omega)*(2.*nb(omega,temp)+1.)

    	end


c	faza

    	real function fazaw2(czas)
    	implicit none
    	real fazaw1,czas,czasdata2,omega,wynik
    	real zero
    	integer inf
	    external fazaw1
	    parameter(zero=0.,inf=1)
	    common /data2/ czasdata2

    	czasdata2=czas

    	call calkai(fazaw1,zero,inf,wynik)
    	fazaw2=wynik
	    end

	    real function fazaw1(omega)
	    implicit none
	    real omega,czasdata2,fazaw
	    external fazaw
	    common /data2/ czasdata2

    	fazaw1=fazaw(omega,czasdata2)
    	end
    	

    	real function fazaw(omega,czas)
    	implicit none
    	real omega,czas,g
    	external g

    	fazaw=g(omega)*sin(omega*czas)
    
    	end





c	omega razy eksponenta; przecalkowanie po theta

	    real function g(omega)
	    implicit none
	    real omega,stala,g02,poprawka
	    external g02
	    parameter(stala=0.016335232,poprawka=1.653061224)

	    g=poprawka*stala*g02(omega)

	    end

    	real function g02(omega)
	    implicit none
    	real omega,g02b,omegadata1,a,b,wynik
    	external g02b
    	common /data1/ omegadata1

    	a=0.
    	b=1.

    	omegadata1=omega

    	call calka(g02b,a,b,wynik)
	    g02=wynik

	    end


	    real function g02b(sintheta)
	    implicit none
	    real omegadata1,sintheta,g02a
	    external g02a
	    common /data1/ omegadata1

	    g02b=g02a(omegadata1,sintheta)
	    end

	    real function g02a(omega,sintheta)
	    implicit none 
	    real omega,sintheta,lp,lz,cl
	    parameter(lp=5.,lz=1.,cl=5.1)
	
	    g02a=omega*
     +	exp(-0.5*lp*lp*omega*omega*(1.-sintheta*sintheta)/(cl*cl)
     +	-0.5*lz*lz*omega*omega*sintheta*sintheta/(cl*cl))

c	omega razy eksponenta
	
    	end






c	nb

    	real function nb(w,temp)
	    implicit none
	    real war,bt,igrek,w,temp

!	bt=h/k
    	bt=7.6386533675
    	igrek=bt*w/temp
    	war=abs(igrek)

		if (igrek.gt.20000.) then
		nb=0.
		else
			if (igrek.eq.0.) then
			nb=0.
			else
				if (war.lt.0.0000001) then
				nb=temp/(bt*w)
				else
				nb=1./(exp(igrek)-1.)
				endif
			endif
		endif
    	end


	    subroutine calka(f,a,b,wynik)
	    implicit none

	    real a,abserr,b,epsabs,epsrel,f,res,work,wynik
      integer ier,iwork,key,lenw,limit,neval,last
	    parameter(epsabs=0.1,epsrel=0.1)
	    parameter(key=2,limit=10000)
	    parameter(lenw=4*limit+1)
      dimension iwork(limit),work(lenw)
      external f
	
	    call qag(f,a,b,epsabs,epsrel,key,res,abserr,neval,ier,
     *    limit,lenw,last,iwork,work)

	    wynik=res
	    end


	    subroutine calkai(f,bound,inf,wynik)
	    implicit none

	    real abserr,epsabs,epsrel,f,res,work,wynik,bound
      integer ier,iwork,key,lenw,limit,neval,inf,last
	    parameter(epsabs=0.000001,epsrel=0.000001)
	    parameter(key=2,limit=10000)
	    parameter(lenw=4*limit+1)
      dimension iwork(limit),work(lenw)
      external f
	
	    call qagi(f,bound,inf,epsabs,epsrel,res,abserr,neval,
     *   ier,limit,lenw,last,iwork,work)

	    wynik=res
    	end
