      subroutine qag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     *    limit,lenw,last,iwork,work)
c***begin prologue  qag
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result)le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real version
c
c            f      - real
c                     function subprogam defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                       7 - 15 points if key.lt.2,
c                      10 - 21 points if key = 2,
c                      15 - 31 points if key = 3,
c                      20 - 41 points if key = 4,
c                      25 - 51 points if key = 5,
c                      30 - 61 points if key.gt.5.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c                      error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yield no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulaties.
c                             if the position of a local difficulty can
c                             be determined (i.e.singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set
c                             to zero.
c                             except when lenw is invalid, iwork(1),
c                             work(limit*2+1) and work(limit*3+1) are
c                             set to zero, work(1) is set to a and
c                             work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end with
c                    ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers to the error
c                    estimates over the subintervals, such that
c                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
c                    form a decreasing sequence with k = last if
c                    last.le.(limit/2+2), and k = limit+1-last otherwise
c
c            work  - real
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left end
c                    points of the subintervals in the partition of
c                     (a,b),
c                    work(limit+1), ..., work(limit+last) contain the
c                     right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last) contain
c                     the error estimates.
c
c***references  (none)
c***routines called  qage,xerror
c***end prologue  qag
c
      real a,abserr,b,epsabs,epsrel,f,result,work
      integer ier,iwork,key,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of lenw.
c
c***first executable statement  qag
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for qage.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call qage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,
     *  ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
c      if(ier.ne.0) call xerror(26habnormal return from  qag ,
c     *  26,ier,lvl)
	if(ier.eq.0) then
	goto 122
	else 
	write(*,*) ier
	endif
122   return
      end
