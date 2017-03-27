       parameter(jmax=180,jmin=0)
       integer bs,bn,bsc,bnc
       integer jmin1,jmax0,jmax00,bn0,jdim,bdim,bs0,bs1,bn00
       integer dt,mm,nstep,time,tmax,ncal
       integer jbs,jbn0,j1,k1,j2,k2
c      apal : palaeotime (propagated forward in this routine)
       real apal,abid,fpal, ppmgz2
       real aa,da,g,lambda,pi,rho,dphi,dlat,r,xdel
       real alpha,beta,gamma
       real h(jmin:jmax),ht(jmin:jmax),b(jmin:jmax),b0(jmin:jmax)
       real acc(jmin:jmax),accum(1:18,5:7)
       real accwri(1:18,5:7),ablwri(1:18,5:7)
       real v(jmin:jmax),vv(jmin:jmax),bbh(jmin:jmax)
       real cosv(jmin:jmax),cosp(jmin:jmax)
       real sinv(jmin:jmax),dsinp(jmin:jmax),del(jmin:jmax)
       real aht(jmin:jmax),bht(jmin:jmax),cht(jmin:jmax),
     &        dht(jmin:jmax)
       real ab(jmin:jmax),bb(jmin:jmax),cb(jmin:jmax),db(jmin:jmax)
       real wb(jmin:jmax),wht(jmin:jmax)
       real dx(jmin:jmax),xdx(jmin:jmax),xx(jmin:jmax),xfix(jmin:jmax)
       real lat(jmin:jmax),vol,volum
       real szp01(1:18,3:7),ajsu0(1:18,5:7)
c***********************************************************************
        open(1,status='old',file='bacc.d')
        open(2,status='unknown',access='append',file='vol.dat')
        open(5,status='old',file='cal5.dat')
        open(15,status='unknown',access='append',file='vol5.dat')
        open(6,status='old',file='cal6.dat')
        open(16,status='unknown',access='append',file='vol6.dat')
        open(8,status='old',file='cal7.dat')
        open(18,status='unknown',access='append',file='vol7.dat')
        open(10,status='unknown',file='cal.dat')
        open(3,status='unknown',file='calo.dat')
c***********************************************************************
c constantes
       jmin1=jmin+1
       jmax0=jmax-1
       jmax00=jmax0-1
       jdim=jmax0-jmin
       g=9.81
       rho=900.
       pi=3.14159265
       r=6371000.
c v---- perturbation stochastique (seed from clock time)
       am = 0
       CALL SYSTEM_CLOCK(ix)
       ix = mod(ix,7919)    


       open (unit=9,file="calctr",status="old")
       read (9,*) tmax    ! 1000 in JGR1
       read (9,*) sd      ! 0.00 in JGR1
       close (9)

c ^---- perturbation stochastique

c ... lecture de l'annee
       open(unit=9,status='unknown',file='bcmin')
       read(9,3104) inpu
c
       read(9,3102) ipal
      read (9,*) apal, ppmgz2
       close(9)
c ... apal : annee du cycle glacial (ky)
c
c _STD da= 45000000. 
       da=100000000. 
       dt=1
       dlat=0.5
       beta=0.5
       calv=-20.
       dphi=dlat*pi/180.
c***********************************************************************
c constantes dependantes de la latitude
       do 200 j=jmin,jmax0
         cosp(j)=cos((j-jmin)*dphi)
         cosv(j)=cos(dphi/2.+(j-jmin)*dphi)
         sinv(j)=sin(dphi/2.+(j-jmin)*dphi)
         lat(j)=j+dlat/2.-j*dlat
 200     continue
       cosp(jmax)=cosp(jmax0)
       sinv(jmax)=0.
       cosv(jmax)=1.
       lat(jmax)=jmax+dlat/2.
       do 210 j=jmin,jmax0
         dsinp(j)=sinv(j+1)-sinv(j)
 210     continue
       dsinp(jmax)=dsinp(jmax0)
c***********************************************************************
c initialisations
c

c ... hubert v
         read(1,1008)apal,(accwri(j,5  ),j=1,18),(ablwri(j,5   ),j=1,18)
         do 302 j=1,18
         accum(j,5)      = accwri(j,5   ) + ablwri(j,5)
 302     continue
c ... hubert ^
c
       do 300 ncal=5,7
c
c ... hubert v
         read(1,1008)apal,(accwri(j,ncal),j=1,18),
     c                    (ablwri(j,ncal),j=1,18)
         do 303 j=1,18
         accum(j,ncal)   = accwri(j,ncal) + ablwri(j,ncal)
 303     continue
c ... hubert ^
c
 300     continue
       volum=0.
       do 301 j=1,18
         szp01(j,3)=0.0
         szp01(j,4)=0.0
 301     continue

c***********************************************************************
c application a chacune des calottes: 5=groen.,6=n.am.,7=euras.
       do 1 ncal=7,5,-1
c***********************************************************************
c constantes
       if(ncal.eq.5)then
         bn=170
         bs=120
         read(5,2011)
         do 350 j=jmin,jmax
           read(5,1011)lat(j),h(j),ht(j),b0(j),b(j),v(j),xfix(j),
     &                 xx(j),acc(j)
 350       continue
c _STD   aa=1.5e-15
c        aa=0.8e-15
c ...    cfr. Hyde et Peltier, 1985, 1987, JAS
c _t01   aa=0.5e-15
         aa=0.2e-15
         dt=1
         mm=3
         lambda=12.4
         else
         if(ncal.eq.6)then
           bn=160
           bs=60
           read(6,2011)
           do 360 j=jmin,jmax
             read(6,1011)lat(j),h(j),ht(j),b0(j),b(j),v(j),xfix(j),
     &                   xx(j),acc(j)
 360         continue
c _STD     aa=1.5e-15
c          aa=0.8e-15
c ...      cfr. Hyde et Peltier, 1985, 1987, JAS
           aa=0.6e-15
           dt=1
           mm=3
           lambda=3.5
           else
             bn=150
             bs=80
             read(8,2011)
             do 370 j=jmin,jmax
               read(8,1011)lat(j),h(j),ht(j),b0(j),b(j),v(j),xfix(j),
     &                     xx(j),acc(j)
 370           continue
c _STD       aa=1.5e-15
c            aa=0.8e-15
c ...        cfr. Hyde et Peltier, 1985, 1987, JAS
             aa=0.6e-15
             dt=1
             mm=3
             lambda=5.5
           endif
         endif            
       bn0=bn-1
       bn00=bn0-1
       bs0=bs-1
       bs1=bs+1
       bdim=bn-bs
       alpha=aa*((rho*g)**mm)*2./(mm+2.)
       gamma=dt*alpha*(lambda**2)/4./(r**(mm-1))
       jbs=ifix(float(bs/10))+1
       jbn0=ifix(float((bn0)/10))+1
       do 380 j=jmin,jmax
         aht(j)=0.0
         bht(j)=0.0
         cht(j)=0.0
         dht(j)=0.0
         dx(j)=0.0
         bbh(j)=0.0
         db(j)=0.0
 380     continue   
c***********************************************************************
c parametres dependant du pas de temps
       xdel=dt*da/r/r/dphi
       do 220 j=jmin,jmax0
         del(j)=xdel/dsinp(j)
         ab(j)=-del(j)*beta*cosv(j)
         bb(j)=1.+del(j)*beta*(cosv(j)+cosv(j+1))
         cb(j)=-del(j)*beta*cosv(j+1)
 220     continue
       bb(jmin)=bb(jmin)+ab(jmin)
       ab(jmin)=0.0
       bb(jmax0)=bb(jmax0)+cb(jmax0)
       cb(jmax0)=0.0
       do 230 j=jmin,jmax
         vv(j)=dt*cosv(j)/2./r/dsinp(j)
 230     continue
c***********************************************************************
c integration sur 1000 ans
       nstep=0
       time=0
c***********************************************************************
c boucle princapale
 2     nstep=nstep+1
       time=time+dt
c***********************************************************************
c calcul des accumulations
       acc(bs  )=(accum(jbs-1,ncal)*4.+accum(jbs,ncal)*6.)/10.*.33/.9
       acc(bs+1)=(accum(jbs-1,ncal)*3.+accum(jbs,ncal)*7.)/10.*.33/.9
       acc(bs+2)=(accum(jbs-1,ncal)*2.+accum(jbs,ncal)*8.)/10.*.33/.9
       acc(bs+3)=(accum(jbs-1,ncal)   +accum(jbs,ncal)*9.)/10.*.33/.9
       do 402 j=jbs,jbn0-1,1
         jj=(j-1)*10+4
         do 403 k=0,9
           acc(jj+k)=(accum(j,ncal)*(10-k)+accum(j+1,ncal)*k)/10.*.33/.9
 403       continue
 402       continue
       acc(bn0-5)=(accum(jbn0,ncal)*.33/.9)
       acc(bn0-4)=(accum(jbn0,ncal)*9.+accum(jbn0+1,ncal)*1.)/10.*.33/.9
       acc(bn0-3)=(accum(jbn0,ncal)*8.+accum(jbn0+1,ncal)*2.)/10.*.33/.9
       acc(bn0-2)=(accum(jbn0,ncal)*7.+accum(jbn0+1,ncal)*3.)/10.*.33/.9
       acc(bn0-1)=(accum(jbn0,ncal)*6.+accum(jbn0+1,ncal)*4.)/10.*.33/.9
       acc(bn0  )=(accum(jbn0,ncal)*5.+accum(jbn0+1,ncal)*5.)/10.*.33/.9


       do 404 j=jmin,bs0
         acc(j)=0.0
 404     continue
       do 405 j=bn,jmax
         acc(j)=0.0 
 405     continue
c***********************************************************************
c calcul de v
       do 500 j=bs,bn
         v(j)=-alpha/(2**(mm+1))/((dphi*r)**mm)
     &        *((ht(j)+ht(j-1))**(mm+1))
     &        *((abs(h(j)-h(j-1))**(mm-1)))*(h(j)-h(j-1))
 500     continue
c***********************************************************************
c calcul du nouveau ht
       do 600 j=bs,bn0
         aht(j)=-beta*v(j)*vv(j)
         bht(j)=1.+beta*(v(j+1)*vv(j+1)-v(j)*vv(j))
         cht(j)=beta*v(j+1)*vv(j+1)
         dht(j)=ht(j-1)*(1.-beta)*v(j)*vv(j)
     &         +ht(j)*(1.+(1.-beta)*(v(j)*vv(j)-v(j+1)*vv(j+1)))
     &         -ht(j+1)*(1.-beta)*v(j+1)*vv(j+1)
 600     continue
       call tridia(aht,bht,cht,dht,wht,bs,bn0,jmin,jmax)
       do 610 j=bs,bn0
         ht(j)=dht(j)
         if(ht(j).lt.1.e-4)ht(j)=0.0
 610     continue
       do 620 j=bs1,bn00
         dx(j)=gamma*(ht(j)**(mm-1))
     &        *(abs(h(j+1)-h(j-1))/2./dphi)**(mm-1)
 620     continue
       dx(bn0)=gamma*(ht(bn0)**(mm-1))
     &        *(abs(h(bn0)-h(bn00))/dphi)**(mm-1)
       dx(bs)=gamma*(ht(bs)**(mm-1))
     &        *(abs(h(bs1)-h(bs))/dphi)**(mm-1)
       do 630 j=bs,bn0
c
c _HUB .v.
         htmax = 4. *h(j) /3.
         if (ht(j).ge.htmax) then
          cor = (4.   /3.  )**3
         else
          cor = (ht(j)/h(j))**3
         end if
         dx(j)=dx(j)*cor
         ht(j)=ht(j)-dx(j)
c ...    Contribution negative du terme de decharge laterale
c _HUB .^.
c
c _STD   ht(j)=ht(j)+dx(j)
         if((-acc(j)*dt).gt.ht(j))acc(j)=-ht(j)/dt
         ht(j)=ht(j)+acc(j)*dt
chubert        if (ht(j).lt..5) ht(j)=0.
 630     continue
c
c do while
c       j=bs
c 6311  continue
c       if (ht(j).gt.0.) go to 6310
c       j=j+1
c       go to 6311
c 6310  continue
c       bsc=j
c end do while
c
c do while
c       j=bn0
c 6321  continue
c       if (ht(j).gt.0.) go to 6320
c       j=j-1
c       go to 6321
c 6320  continue
c       bnc=j
c end do while
c***********************************************************************
c "calving surge" au sud de la calotte
c       j = bsc + 1
c  640  if ((0.9*ht(j-1)).lt.(-b(j-1)).and.b(j).lt.0.0) then
c       ht(j) = ht(j) + calv * dt
c       j = j + 1
c       go to 640
c       endif
c***********************************************************************
c "calving surge" au nord de la calotte
c       j = bnc - 1
c  650  if ((0.9*ht(j+1)).lt.(-b(j+1)).and.b(j).lt.0.0) then
c       ht(j) = ht(j) + calv * dt
c       j = j - 1
c       go to 650
c       endif
c***********************************************************************
c calcul du nouveau bedrock
       do 700 j=jmin,jmax0
         bbh(j)=(1.-beta)*b(j)-b0(j)+ht(j)/3.
 700     continue
       db(jmin)=b(jmin)+del(jmin)*cosv(jmin1)*(bbh(jmin1)-bbh(jmin))
       do 710 j=jmin1,jmax00
         db(j)=b(j)+del(j)*(cosv(j+1)*(bbh(j+1)-bbh(j))
     &                      -cosv(j)*(bbh(j)-bbh(j-1)))   
 710     continue
       db(jmax0)=b(jmax0)+del(jmax0)*cosv(jmax0)
     &                   *(bbh(jmax0)-bbh(jmax00))
       call tridia(ab,bb,cb,db,wb,jmin,jmax0,jmin,jmax)
       do 720 j=jmin,jmax0
           b(j)=db(j)
 720       continue     
c***********************************************************************
c calcul de l'altitude
       do 800 j=jmin,jmax
         h(j)=ht(j)+b(j)
 800     continue 
c***********************************************************************
       if(time.ge.tmax)goto 99
       goto 2
 99    continue
c************************************************************************
c calcul de l'extension en longitude et du volume de glace
       do 900 j=jmin,bs0
         xx(j)=0.0
 900     continue
       vol=0.0
       do 901 j=bs,bn0
         if(ht(j).gt.0.0)then
           xx(j)=2.*h(j)*h(j)/lambda
           vol=vol+2./3.*xx(j)*ht(j)*r*dphi
           xx(j)=xx(j)/2./pi/r/cosp(j)
           if(xx(j).gt.xfix(j))then
             vol=vol-4./3.*sqrt(lambda)*(1.-b(j)/h(j))
     &           *((xx(j)-xfix(j))*pi*r*cosp(j))**1.5*r*dphi
             xx(j)=xfix(j)
             endif
           else
           xx(j)=0.0
           endif
 901     continue
       do 902 j=bn,jmax
         xx(j)=0.0
 902     continue
c output
       if(ncal.eq.5)then
         close(5)
         open(5,status='unknown',file='cal5.dat')
         write(5,2012)ncal,apal
         do 15 j=jmin,jmax
           write(5,1012)lat(j),h(j),ht(j),b0(j),b(j),v(j),xfix(j),
     &                  xx(j),acc(j)
 15        continue
         close(15)
         open(15,status='old',access='append',file='vol5.dat')
         write(15,1009)apal,vol/1.e15
         volum=volum+vol
         do 5 j=1,18
           ajsu0(j,5)=0.0
           szp01(j,5)=0.0
           do 25 k=0,9
             l=(j-1)*10+k
             ajsu0(j,5)=ajsu0(j,5)+xx(l)/10.
             szp01(j,5)=szp01(j,5)+h(l)*xx(l)
25           continue
c
c ... hubert v
cifte
         if (.not.ajsu0(j,5).ne.0.) go to 250
cthen
             szp01(j,5)=szp01(j,5)/(10.*ajsu0(j,5))
         go to 251
celse
250      continue
             l=(j-1)*10+4
             szp01(j,5)=h(l)
251      continue
cendif
 5           continue
c ... hubert ^
c
         write(10,1010)ncal,(szp01(j,5),j=1,18) 
         write(10,1010)ncal,(ajsu0(j,5),j=1,18) 
         else
         if(ncal.eq.6)then
           close(6)
           open(6,status='unknown',file='cal6.dat')
           write(6,2012)ncal,apal
           do 16 j=jmin,jmax
             write(6,1012)lat(j),h(j),ht(j),b0(j),b(j),v(j),xfix(j),
     &                    xx(j),acc(j)
 16          continue
           close(16)
           open(16,status='old',access='append',file='vol6.dat')
           write(16,1009) apal,vol/1.e15
           volum=volum+vol
           do 6 j=1,18
             ajsu0(j,6)=0.0
             szp01(j,6)=0.0
             do 26 k=0,9
               l=(j-1)*10+k
               ajsu0(j,6)=ajsu0(j,6)+xx(l)/10.
               szp01(j,6)=szp01(j,6)+h(l)*xx(l)
26             continue
c
c ... hubert v
cifte
           if (.not.ajsu0(j,6).ne.0.) go to 260
cthen
               szp01(j,6)=szp01(j,6)/(10.*ajsu0(j,6))
           go to 261
celse
260        continue
               l=(j-1)*10+4
               szp01(j,6)=h(l)
261        continue
cendif
 6             continue
c ... hubert ^
c
           write(10,1010)ncal,(szp01(j,6),j=1,18) 
           write(10,1010)ncal,(ajsu0(j,6),j=1,18) 
           else
             close(8)
             open(8,status='unknown',file='cal7.dat')
             write(8,2012)ncal,apal
             do 17 j=jmin,jmax
               write(8,1012)lat(j),h(j),ht(j),b0(j),b(j),v(j),xfix(j),
     &                      xx(j),acc(j)
 17            continue
             close(18)
             open(18,status='old',access='append',file='vol7.dat')
             write(18,* ) apal,vol/1.e15, ppmgz2
             volum=volum+vol
             do 7 j=1,18
               ajsu0(j,7)=0.0
               szp01(j,7)=0.0
               do 27 k=0,9
                 l=(j-1)*10+k
                 ajsu0(j,7)=ajsu0(j,7)+xx(l)/10.
                 szp01(j,7)=szp01(j,7)+h(l)*xx(l)
27               continue
c
c ... hubert v
cifte
             if (.not.ajsu0(j,7).ne.0.) go to 270
cthen
                 szp01(j,7)=szp01(j,7)/(10.*ajsu0(j,7))
             go to 271
celse
270          continue
                 l=(j-1)*10+4
                 szp01(j,7)=h(l)
271          continue
cendif
 7               continue
c ... hubert ^
c
             write(10,1010)ncal,(szp01(j,7),j=1,18) 
             write(10,1010)ncal,(ajsu0(j,7),j=1,18) 
           endif
         endif  
c***********************************************************************
c valeurs transmises au bcm
       do 3 j=1,18 
         do 4 k=0,9
           l=(j-1)*10+k
           szp01(j,3)=szp01(j,3)+h(l)*xfix(l)
 4         continue
         szp01(j,4)=szp01(j,3)
 3       continue
 1     continue
       do 8 k=4,3,-1
         write(10,1010)k,(szp01(j,k),j=1,18)
 8       continue
       close(2)
        open(2,status='old',access='append',file='vol.dat')
c***********************************************************************

       open(unit=9,status='unknown',file='bcmin')
       read(9,3104) inpu
 104   format(l12)
c ... utilisation ou non etat d'une simulation precedente
c
       read(9,3102) ipal
 102   format(i12)
      read (9,*) apal, ppmgz2
c ... apal : annee du cycle glacial (ky)
c
      read(9,103) fpal
 103   format(f9.3)
c ... apal : annee de fin du run (utile pour calot uniquement) (ky)
       rewind 9
c
      write(9,3104) inpu
c ... utilisation ou non etat d'une simulation precedente
c
      write(9,3102) ipal
c ... ipal : no du run
c
c
C update time (for output)
         apal = apal + tmax / 1.e3 
         if (apal.ge.fpal) then
          open (unit=1,file="stop",status="unknown")
          close(1)
         endif

      write(9,*) apal, ppmgz2
c ... apal : annee du cycle glacial (ky)
      write(9,3103) fpal
c ... fpal : fin du cycle glacial (ky)
c
      close(unit=9)
      close(unit=10)
c
       write(2,1009) apal,volum/1.e15

      open (10, status='unknown',file='volco2.dat')
      write (10,*) volum/1.e15, ppmgz2
      close(10)

c1010 format(i5,9e12.4,/5x,9e12.4)
 1010 format(i5,2x,9e12.4,/,7x,9e12.4)

 1008   format(f9.3,2x,9e12.4,/,7x,9e12.4,
     .        4x,2(/,7x,9e12.4),4x)
 1009 format(f9.3,2f14.6)
 1011 format(1x,f5.2,8(1x,e10.4))
 1012 format(1x,f5.2,8(1x,e10.4))
 2011 format(1x,/,1x,/,1x)
 2012 format(' calotte nr',i2,' annee',f9.3,
     &//,' lat   altitude   epaisseur  bed.'
     &'init.  bedrock    vitesse   % cont.    % glace   accumulation')
 3102   format(i12)
 3103   format(f9.3)
 3104   format(l12)
       stop
       end
c???????????????????????????????????????????????????????????????????????
c{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
	subroutine tridia(a,b,c,d,w,kmin,kmax,km,kp)
        real a(km:kp),b(km:kp),c(km:kp),d(km:kp),w(km:kp)
        ksup=kmax
        if(b(ksup).eq.0.0)then
          if(a(ksup).eq.0.0)then
            write(3,*)'la ',ksup,'ieme equation est id. nulle'
            return
            endif
          if(c(ksup-1).eq.0.0)then
            write(3,*)'le syst. est ind. pour la ',ksup,'ieme valeur'
            return
            endif
          d(ksup)=d(ksup)/a(ksup)
          d(ksup-1)=d(ksup-1)-d(ksup)*b(ksup-1)
          d(ksup-2)=d(ksup-2)-d(ksup)*c(ksup-2)
          b(ksup-1)=c(ksup-1)
          c(ksup-2)=0.
          ksup=ksup-1
          endif
        w(ksup)=1./b(ksup)
        do 20 k=ksup-1,kmin,-1
          w(k)=1./(-c(k)*w(k+1)*a(k+1)+b(k))
20       continue
        d(ksup)=w(ksup)*d(ksup)
        do 120 k=ksup-1,kmin,-1
          d(k)=w(k)*(d(k)-c(k)*d(k+1))
 120      continue
        do 180 k=kmin+1,ksup
          d(k)=-w(k)*a(k)*d(k-1)+d(k)
 180      continue
        if(ksup.ne.kmax)then
          dd=d(ksup)
          d(ksup)=d(ksup+1)
          d(ksup+1)=dd
          endif
        return
        end
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c la sous routine tridiv tient compte de la presence de nombreux 0 dans
c la matrice tridiagonale pour simplifier la methode du compagnon.
c wk1 : matrice a inverser perdue en fin de routine de dim jdim.
c wk2 : matrice resultat.
c lmin,lmax : indices des matrices ; lmin:ind.minimum; lmax:ind.maximum
c
      subroutine tridiv(wk1,wk2,lmin,lmax)
c
      real wk1(lmin:lmax,lmin:lmax),wk2(lmin:lmax,lmin:lmax)
c
      lmin1=lmin+1
      lmax0=lmax-1
c
c calcul des nouveaux elements de la diagonale superieure suite aux
c operations effectuees pour obtenir 0 sur la diagonale inferieure et
c 1 sur la diagonale centrale. calcul de l'element diagonal central 
c suite a l'operation effectuee pour obtenir 0 sur la diagonale 
c inferieure.
c 
      do 50 l=lmin,lmax0
        wk1(l,l+1)=wk1(l,l+1)/wk1(l,l)
        wk1(l+1,l+1)=wk1(l+1,l+1)-wk1(l+1,l)*wk1(l,l+1)
 50     continue
c
c modification de la matrice compagnon unite en une matrice triangulaire
c inferieure suite au operations effectuees sur la matrice de depart
c pour obtenir des 0 sur la diagonale inf. et des 1 sur la diagonale
c centrale.
c
      wk2(lmin,lmin)=1./wk1(lmin,lmin)
      do 30 l=lmin1,lmax
        wk2(l,l)=1./wk1(l,l)
        lm1=l-1
        do 30 k=lm1,lmin,-1
          wk2(l,k)=-wk2(l,k+1)*wk1(k+1,k)/wk1(k,k)
 30       continue
c
c modification de la matrice compagnon suite au operations effectuees
c sur la matrice de depart pour obtenir des 0 sur la diagonale sup.
c
      do 60 l=lmax0,lmin,-1
        do 60 k=lmin,lmax
          wk2(l,k)=wk2(l,k)-wk1(l,l+1)*wk2(l+1,k)
 60       continue
      return
      end
