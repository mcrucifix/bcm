c ---------------------------------------------------------------------
c   climat present  cycle saisonnier  run 20 ans (0h40 cpu sur Dec3100)
c ---------------------------------------------------------------------
c         simplification routine infra
c   calcul des flux i.r. a 4 niveaux seulement
c   calcul des flux solaires par simps = true
c   degradation du nombre de niveaux dans bcmr2.for (subroutine iniras)
c   couches de 100 mb au lieu de couches de 75 mb
c   sauf pres de la surface  situation inchangee  ( 3 couches )
c        ------   ciel normal  -----------
c ---------------------------------------------------------------------
c   pour adapter le programme au systeme u n i x 
c                 a partir de la version v m s, faire les operations :
c                 1) substitute/       open/c mvax open/whole
c                 2) substitute/c unix/      /whole
c                 3) ne pas oublier la correction dans le block data :
c                    substitute/       data iwr1/c mvax data iwr1
c ---------------------------------------------------------------------
c cray .if [ $concord .eq. 1 ]
c cray nzyx  bcm=====
c cray .endif
      program bcm
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
c +++ bcm00.for +++
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical       lgm18k,taiga,otter,summer,co2var
c
      double precision  cpwth2,rhwth2,ocath2
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      double precision  hssn2,hspsn2
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/varl02/lgm18k,taiga,otter,summer,co2var
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,jour,nts
      common/time2/delt,tht,it
      common/pal/apal,fpal,ipal,ipaleo
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      common/th1/cpath1,rath1,akth1
      common/th2/cpwth2,rhwth2,ocath2
c
c --- variables dynamiques
c
      common/p01/szp01(np,nilw),spp01(np,nilw)
      common/p02/szp02(np,nilw)
      common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
      common/dy2/t2ady2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .                 pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
      common/dy4/qdy4(np,nhd),dqdy4(np,nhd)
      common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
      common/dy7/zsdy7(np) ,ztdy7(np),
     .                 usdy7(npp),utdy7(npp),udy7(npp,nhd)
c
c --- variables surface
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .                ajosu0(np),ajsu0(np,nilw)
      common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .                ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .             hiic2(np),hipic2(np),awic2(np),awpic2(np)
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .                 hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      common/la1/tcmla1(np),tcdla1(np),csmla1(np),csdla1(np)
      common/la2/tcmla2(np),tcdla2(np)
      common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
      common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .           vopsn3(np,nilw), vosn3(np,nilw),
     .           hspsn3(np,nilw), hssn3(np,nilw),
     .           ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
      common/ca1/ajca1(np,nilw)
      common /so10/ perh,ecc,so
c pour exp. co2
c
c --- variables rayonnement
c
       common/fac/facsc,facco2,faccl
       common/gz2/ppmgz2
       real palco2, co2inter, ppmgz2, ppmgz2t, volum
c
c **********************************************************************
c
c --- ouverture des fichiers d'output
c
       open(unit=iwr6,status='unknown',file='bcmlog')
       rewind iwr6
       open(unit=iwr2,status='old',file='bcmfn')
       rewind iwr2
c mVax open(unit=iwr1,status='new',recl=150,file='bcmout.d')
       open(unit=iwr1,status='unknown',file='bout.d')
       rewind iwr1
c mVax open(unit=3,status='new',recl=150,file='bcmio.d')
       open(unit=3,status='unknown',file='bcmo.d')
       rewind 3
c mVax open(unit=4,status='new',recl=150,file='bcmls.d')
       open(unit=4,status='unknown',file='bcml.d')
       rewind 4
c mVax open(unit=iwr5,status='new',recl=150,file='bcmc.d')
       open(unit=iwr5,status='unknown',file='bcmc.d')
       rewind iwr5
c mVax open(unit=iwr7,status='new',recl=150,file='bcmr.d')
       open(unit=iwr7,status='unknown',file='bcmr.d')
       rewind iwr7
c mVax open(unit=8,status='new',recl=150,file='bcmm.d')
       open(unit=8,status='unknown',file='bcmm.d')
       rewind 8
c mVax open(unit=9,status='new',recl=150,file='bcma.d')
       open(unit=9,status='unknown',file='bcma.d')
       rewind 9
c
c **********************************************************************
c
c --- lecture de quelques parametres concernant le cycle glacial
c
c mVax open(unit=10,status='old',recl=150,file='bcmin.dat')
       open(unit=10,status='old',file='bcmin')
       rewind 10
       read(10,104) inpu
 104   format(l12)
c ... utilisation ou non etat d'une simulation precedente
c
       read(10,102) ipal
 102   format(i12)
c ... ipal : no du run
c
       read(10,103) apal0
       apal = apal0
 103   format(f8.3)
       read(10,103) fpal0
       fpal = fpal0
c ... apal : annee du cycle glacial (ky)
c
c pour exp. co2
c     apali=apal
c pour exp. co2
      close(unit=10)
      ipaleo=ipal-1
c
c **********************************************************************
c
c --- choix des parametres de la simulation
c
c mVax open(unit=18,status='old',recl=150,file='bcmctr.dat')
       open(unit=18,status='old',file='bcmctr')
       rewind 18
c
       read(18,18) nts
 18    format(i8)
c ...  nts    = number of time steps run in this program
c     (nts    = 120  =>  1 an)
c
       read(18,18) iprint
c      iprint = 0
       if(nts.eq.1)  iprint = 1
c ...  iprint = 1 => program print every time step
c
c ----------------------------------------------------------------------
c
      read(18,180) season
 180  format(7x,l1)
c     season =.true.
c ... cycle saisonnier ou non
c
      read(18,180) ocen
c     ocen   =.true.
c ... transports de chaleur oceanique ou non
c
      read(18,180) ekm
c     ekm    =.true.
c ... ekm  : decide ou non de la parametrisation du courant d'ekman
c
      read(18,180) tsclim
c     tsclim =.false.
c ... forcage de l'atmosphere via la climatologie de t4 ou non
c
      read(18,180) taclim
c     taclim =.false.
c ... forcage de la surface   via la climatologie de t2 et u4 ou non
c
      read(18,180) simpir
c     simpir =.false.
c ... retroaction temperature - extension glace continentale ou non
c
      read(18,180) turbu
c     turbu  =.true.
c ... parametrisation du flux tourbillonnaire macroechelle ou non
c
      read(18,180) lgm18k
c     lgm18k  =.false.
c ... simulation dernier maximum glaciaire avec donnees de CLIMAP ou non
c
      read(18,180) taiga
c     taiga  =.false.
c ... prise en compte taiga feedback ou non pour le climat present
c
      read(18,180) otter
c     otter  =.false.
c ... prise en compte taiga feedback d'Otterman ou non
c
      read(18,180) summer
c                summer =.false.
      if (otter) summer =.false.
c ... taiga feedback base sur la temperature estivale ou non
c
      read(18,180) co2var
c     co2var =.false.
c ... prise en compte variation CO2(Vostok) ou non
c
      read(18,181) facco2
 181  format(f9.0)
c     facco2 = 1.
c ... Simulation facco2 X CO2 ou non
c
      read(18,181) facsc
c     facsc  = 0.
c ... Simulation  (1 + facsc/100) X 1368Wm-2 
c
      close(unit=18)
c
c --- ecriture des parametres de la simulation
c
      write(iwr1,20)nts,season,tsclim,taclim,ocen,ekm,
     .                  turbu,facsc,facco2,simpir,inpu,
     .                  apal,lgm18k,co2var,taiga,otter,summer
      write(iwr6,20)nts,season,tsclim,taclim,ocen,ekm,
     .                  turbu,facsc,facco2,simpir,inpu,
     .                  apal,lgm18k,co2var,taiga,otter,summer
   20 format (
     ./' nts   ',7x,'= ',i9,
     ./' season',7x,'= ',l9,
     ./' tsclim',7x,'= ',l9,18x,' taclim',7x,'= ',l9,
     ./' ocen  ',7x,'= ',l9,18x,' ekm   ',7x,'= ',l9,
     ./' turbu ',7x,'= ',l9,
     ./' facsc ',7x,'= ',f13.3,
     ./' facco2',7x,'= ',f13.3,
     ./' simpir',7x,'= ',l9,
     ./' inpu  ',7x,'= ',l9,
     ./' apal  ',7x,'= ',f13.3,
     ./' lgm18k',7x,'= ',l9,
     ./' co2var',7x,'= ',l9,
     ./' taiga ',7x,'= ',l9,
     ./' otter ',7x,'= ',l9,
     ./' summer',7x,'= ',l9)
c
c **********************************************************************
c
c --- constantes dynamiques et thermodynamiques
c
      call iniphy
c
c --- discretisation meridienne
c
      ca = 6.371e6
c ... ca : rayon de la terre (m)
c
      pi = 3.1415926567
      dg = 2.*pi/360.
c ... dg converts degrees to radians
c
      cp = pi/(real(np)*2.)
c ... cp = angle subtended by one interval
c
      dels = ca*cp
c ... dels : largeur d'une bande de latitude
      ds2  = dels**2
c
      gs  =  90./np
c ... gs  :  largeur d'une bande de latitude (degres)
      gs2 =  gs /2.
c
      jp  = np
c ... jp  =  nombre de bandes de latitude
      jpp = jp+1
c
c --- initialisation des variables auxiliaires spatiales
c
      call inigr
c
c --- discretisation temporelle
c
      delt = 72.*3600.
c ... delt = pas de temps
      tht  = delt / 2.
c
      td   = 86400./delt
c ... td   = nombre de pas de temps par jour
      td1  = 1./td
c
       tsm = 360.*24.*3600./(delt*12.)
c ...  tsm = nb de pas de temps par mois
      itsm = int(tsm)
      ints = itsm/2
c
      itint = 1
c ... itint = time steps between interpolation
c
c --- initialisation des variables temporelles
c
      it   = 0
c ... it   = no of time step
c   - attention : ne pas modifier l'assignation it=0 *
c                 certaines variables sont initialisees sous la conditio
      mo    = 1
      mon   = 0
      mm    =12
      monpr = 1
c ... monpr = the first month being printed
c
      ian   = 1
c ... ian   = no du cycle annuel en cours
c
c --- recalibrage du nombre d'iterations de telle maniere
c             a avoir un nombre entier de cycles annuels
c
      if (season.and.it.gt.itsm*12) then
       nts=nts/(itsm*12)
       nts=nts*(itsm*12)
      end if
c
      write(iwr1,21)dels,td,itint
      write(iwr6,21)dels,td,itint
   21 format (
     .//' grid spacing                     = ',e12.5,
     .//' ts per day                       = ',f12.3,
     .//' time steps between interpolation = ',i12)
c
c pour exp. co2
c     write(iwr1,22)dels,td,itint,ipal,apal
c     write(iwr6,22)dels,td,itint,ipal,apal
c  22 format (
c    .//' grid spacing                     = ',e12.5,
c    .//' ts per day                       = ',f12.3,
c    .//' time steps between interpolation = ',i12,
c    .//' annee du run pour insolation (ky)= ',i12,
c    .//' annee du debut pour scenario co2 = ',f12.3)
c pour exp. co2
c
c --- initialisation des capacites calorifiques en surface
c
      call inics
c
c **********************************************************************
c
c --- initialisation par lecture de fichier d'input
c
      if (inpu) then
c
c --- initialisation du climat
c
c mVax open(unit=11,status='old',recl=150,file='bcmsim.d')
       open(unit=11,status='old',file='bsim.d')
       rewind 11
       read(11,1)
     . t2dy2,t4su0,t4ssu0,t4rsu0,tmsu1,t4dsu1,tcdla1,tcdla2,ajsu0,
     . ajasu1,
     . twoc2,sstoc2,hoc2,heqoc2,ajosu0,
     . twic2,twiic2,tsiic2,awpic2,awic2,hipic2,hiic2,
     . ajca1,hspsn2,hssn2,
     . szp01,pstdy2,qdy4,akody5,akdy5,zsdy7,ztdy7,assn3
    1  format(/,3(3(1x,e12.5),/))
       close(unit=11)

c --- initialisation de la calotte
c
c mVax open(unit=12,status='old',recl=150,file='cal.dat')
       open(unit=12,status='old',file='cal.dat')
       rewind 12
       do 3 iilw=nilw,5,-1
       read(12,2)ical,(szp01(j,iilw),j=1,jp)
       read(12,2)ical,(ajsu0(j,iilw),j=1,jp)
 2     format(i5,2x,9e12.4,/7x,9e12.4)
       do 30 j=1,jp
       szp01(j,iilw) = szp01(j,iilw) * 2. / 3.
       if (ajasu1(j,iilw).gt.0.) then
        if (ajsu0(j,iilw).lt.0.0001) then 
         ajsu0(j,iilw) = 0.0001
c ...    Traitement des montagnes :
c        if (iilw.eq.5) then
c          if (szp01(j,iilw).lt.1200.) 
c    .     szp01(j,iilw) = 0.5 * (szp01(j,iilw) +1200.)
c ...     Groenland
c        else
c         if (j.gt.11.or.iilw.eq.6) then
c          if (szp01(j,iilw).lt. 600.) 
c    .     szp01(j,iilw) = 0.5 * (szp01(j,iilw) + 600.)
c         else
c          if (szp01(j,iilw).lt. 200.) 
c    .     szp01(j,iilw) = 0.5 * (szp01(j,iilw) + 200.)
c         end if
c ...     Scandinavie / Cotes Est de l'Amerique du Nord 
c        end if
        end if   
       end if
       szp02(j,iilw) = szp01(j,iilw)
c ...  rappel : le traitement des autres variables caracterisant le 
c               sous-systeme iilw se fait dans la routine inila2.
 30    continue
 3     continue
      end if

       open (10, file='volco2.dat', status='old')
       read (10,*) volum, ppmgz2
       close(10)

      write(*,*) 'volum =',volum
      if (co2var) then
         ppmgz2=palco2(apal)
c         call celest(double(apal))
c        ppmgz2=co2inter (ppmgz2, volum, perh, so, ecc)

        write (*,*) "co2 =",ppmgz2
       else
        ppmgz2 = 210
       end if
c

c
c **********************************************************************
c
c --- beginning new time step
c
  501  day = real(it)/td
      iday =  day
c ... iday = total no of days calculated
c
c pour exp. co2
c     apal=real(it)/120. + apali
c pour exp. co2
c
c --- parametres determinant les output
c
c
      mts  = mod(it+ints,itsm)
c ... mts  = time steps since middle of month
c
      momt = mod(mts,itint)
c ... momt = count of time steps for interpolation
c
c --- bilan annuel
c
      if (season) then
       if (amod(day,360.).lt.td1) call enbil
      else
       if (mts           .eq.0  ) call enbil
      end if
c
      jour = iday - (ian-1) * 360
      if (momt.gt.0) go to 50
      if (mts .ne.0) go to 50
          mo  = mo  + 1
          mon = mon + 1
      if (mo.eq.13) mo = 1
          mm  = mo  - 1
   50 if (mo.eq.1)  mm = 12
c
      if (iprint.eq.1) go to 41
      if (mts.ne.0)    go to 40
      if (monpr.gt.mm) go to 40
   41 continue
      write(iwr1,42)season,it,ian,mm,jour
   42 format (///1x,120('+'),
     . /,'  season =',2x,l1,5x,'iter.no',i6,
     . /,'  annee  =',i3,5x,'mois  =',i6,5x,'jour  =',i5,/1x,120('+'))
   40 continue
c
c **********************************************************************
c
c --- thermodynamique et modeles de surface
c
      call htng
c pour exp. co2
c     write(iwr6,24)it,apal,ppmgz2
c  24 format (
c    .//' numero de l''iteration    it    = ',i12,
c    .//' temps en annee    pour   it    = ',f12.3,
c    .//' concentration co2 pour   it    = ',f12.3)
c pour exp. co2
c
c **********************************************************************
c
c --- dynamique atmospherique
c
      call dyn
c
c **********************************************************************
c
      call prtflu
c
      it =   it+1
      ipaleo = ipal
      if (it.le.nts) go to 501
c
c **********************************************************************
c
c --- ecriture de quelques parametres concernant le cycle glacial
c
      inpu = .true.
      ipal = ipal + 1
c      apal = ipal   ! apal est maintenant gere depuis cal.v4
c paleo
c
c pour exp. co2
c     inpu = .false.
c     apal = apali
c pour exp. co2
c
       open(unit=10,status='unknown',file='bcmin')
       rewind 10
c
      write(10,104) inpu
c ... utilisation ou non etat d'une simulation precedente
c
      write(10,102) ipal
c ... ipal : no du run
c
      write(10,*) apal, ppmgz2
c ... apal : annee du cycle glacial (ky)
c
      write(10,103) fpal
c ... apal : annee de fin du run (utile pour calot uniquement) (ky)
c
      close(unit=10)
c
c pour exp. co2
      write(iwr1,23)              ipal,apal,ppmgz2
      write(iwr6,23)              ipal,apal,ppmgz2
   23 format (
     .//' annee du run pour insolation (ky)= ',i12,
     .//' annee de fin pour scenario co2 = ',f12.3,
     .//' concentration co2 en fin de run= ',f12.3)
c pour exp. co2
c **********************************************************************
c
c --- fermeture des fichiers d'output
c
       close(unit=iwr1)
       close(unit=iwr2)
       close(unit=3)
       close(unit=4)
       close(unit=iwr5)
       close(unit=iwr7)
       close(unit=8)
       close(unit=9)
c
       close(unit=iwr6)
c
      stop
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  block===
c cray .endif
      block data
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
       data iwr1/ 1/,iwr2/ 2/,iwr5/ 5/,iwr6/ 6/,iwr7/ 7/
c mvax data iwr1/21/,iwr2/22/,iwr5/25/,iwr6/26/,iwr7/27/
c
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  phy=====
c cray .endif
      subroutine iniphy
c
c +++ bcm01.for +++
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/ir0/stfir0
      common/th0/cawth0
      common/th1/cpath1,rath1,akth1
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
      double precision  a1th5,a2th5,c2th5,c4th5
      double precision m1th5,m2th5,m3th5,m4th5,m5th5,k1th5,k2th5,k3th5
      common/th5/m1th5,m2th5,m3th5,m4th5,m5th5,a1th5,a2th5,c2th5,c4th5,
     . k1th5,k2th5,k3th5
c
       stfir0 = 5.67e-8
c ...          constante de stefan-boltzman                     (w.m-2)
c      emwir0 = 0.97
       emwir0 = 1.
c ...          emw = longwave emissivity of sea water
c      emiir0 = 0.97
       emiir0 = 1.
c ...          emi = longwave emissivity of ice
c
       cawth0 =     4.1855
c ...          constante de conversion (cal) -> (j)
c
       cpath1 =  1004.
c ...          chaleur specifique de l'air a pression constante (j.kg-1)
        rath1 =   287.05
c ...          constante des gaz parfaits pour l'air sec        (j.kg-1)
        akth1 = rath1/cpath1
c
       cpwth2 = 4000.
c ...          capacite calorifique de l'eau de mer             (j.kg-1)
       rhwth2 = 1025.
c ...          masse volumique de l'eau de mer                  (kg.m-3)
        ocath2 = cpwth2 * rhwth2
c
        tfth3 = 271.2
c ...          tf=freezing point of sea water (k)
       tfith3 = 273.05
c ...          tfi=melting point of surface ice (k)
       hfith3 = 302.e+06
c ...          hfi=heat of fusion of ice (j.m-3)
       cdith3 =   2.04
c ...          cdi=conductivity of ice (w.m-1.k-1)
       c88th3 =   0.
c      c88th3 =   0.88
c ...          c88=density of ice
c
       tfsth4 = 273.15
c ...          tfs=melting point        of snow                 (k)
       hfsth4 = 3.3e+08
c ...          hfs=heat of fusion       of snow                 (j.m-3)
       hvsth4 = 2.5e+09
c ...          hfs=heat of vaporisation of snow                 (j.m-3)
       rhsth4 = hfsth4 / hvsth4
c ...          rapport
       cdsth4 =    .3
c ...          cds=conductivity         of snow                 (w.m-1)
       c33th4 =    .33
c ...          c33=density              of snow
       chfth4 = c33th4 * hfsth4 / (360.d00*86400.d00)
c
       m1th5=0.45d0
       m2th5=2.6d0
       m3th5=1.9d0
       m4th5=2.3d0
       m5th5=0.6d0
       a1th5=0.6d0
       a2th5=0.3d0
       c2th5=((3.d0-2.d0*m5th5)*(m2th5+m3th5)-m5th5*m3th5)/3.d0
       c4th5=2.d0*m4th5/(m1th5*m1th5)
       k1th5=(3.d0*c2th5)/(1.d0-m5th5)
       k2th5=m4th5/(2.d0-2.d0*m5th5)
       k3th5=m3th5+m2th5/(1.d0-m5th5)
c ...          constantes pour le modele oml
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  htng====
c cray .endif
      subroutine htng
c
c +++ bcm02s.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nir=50,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical fulwri
      logical hydroa
      logical simps
c ... variable pour le choix routine solaire
c
      double precision  chhv0,cshv0,htshv0
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c pour exp. co2
      common/pal/apal,fpal,ipal,ipaleo
c pour exp. co2
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
c
c --- constantes physiques
c
      common/th1/cpath1,rath1,akth1
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables dynamiques
c
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy2/t2hdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
c
c --- variables vapeur d'eau
c
       common/cl1/clcl1(np),clhcl1,ztcl1(np),zbcl1(np),
     .                             ptcl1(np),pbcl1(np)
       common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
       common/ql2/qlql2(np),fsql2(np)
c
c --- variables rayonnement
c
       common/gz2/ppmgz2
       common/ra1/ttra1(nir),hhra1(nir),ozra1(nir),cldra1(nir),
     .            pmbra1(nir),zkmra1(nir)
       common/rf0/rgsrf0(np,nilw)
       common/so1/sotso1(np),cmuso1(np),bhso1(np,6)
       common/so2/gtso2(np),
     .            atso2(np),ab0so2(np),abcso2(np),albso2(np)
       common/so3/saso3(np),sgso3(np)
       dimension alb(np),solag(np)
       common/ir0/stfir0
       common/ir1/fimir1(nir),fidir1(nir),finir1(nir)
       common/ir2/fupir2(np),fdnir2(np),fnir2(np),divir2(np)
       common/ir3/fupir3(np),fdnir3(np),fnir3(np)
       common/ir4/c3ir4(np,nilw)
        common/as0/tas0(np,nilw,3),oas0(np,nilw,3),gas0(np,nilw,3)
        common/as1/tas1(np,nilw,3,12),oas1(np,nilw,3,12),
     .             gas1(np,nilw,3,12)
        common/as2/tas2(np,nilw,3),oas2(np,nilw,3),gas2(np,nilw,3)
c
c --- variables interface
c
       common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
       common/hv1/htshv1(6,np),shshv1(np)
       common/hv2/htshv2(6,np,nilw),zhshv2(6,np),shshv2(np),shahv2(np),
     .  albhv2(np),sahv2(np),fuphv2(np),fdnhv2(np),hslhv2(np),t2hv2(np)
       common/hv3/ialhv3(np)
c
c --- variables surface
c
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       dimension            t4s   (np,nilw)
       common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .           ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
       common/su3/hswsu3(np,nilw),cksu3(np,nilw)
       common/su4/t4su4(np),t4asu4(np,nilw)
       common/sn6/hfasn6(np,nilw),  nsn6(np,nilw),nnsn6(np,nilw)
c
c --- variables bilan general
c
       common/bi0/surbi0(np,nilw),atmbi0(np),topbi0(np)
       common/bi1/surbi1(nilw)   ,atmbi1    ,topbi1    ,radbi1(np)
       common/ham/t4ham,t2ham,albham,fupham,fdnham,shham,alhham
       common/dt4/dt4ep,dt4t2
       dimension sens(np),evap (np),hshl(np)
       dimension          evapa(np)
c
c **********************************************************************
c
c --- print at middle of the month
c
      if (iprint.eq.1) go to 11
      if (mts.ne.0) go to 10
      if (monpr.gt.mm) go to  10
   11 continue
      write(iwr1,12)
   12 format(//' ---- htng   ----')
   10 continue
c
c ----------------------------------------------------------------------
c
c --- initialisation des constantes locales
c
      fulwri = .false.
c ... fulwri : parametre controlant le volume des output
c
      simps  = .true .
c ---  choix routine radiative solaire simps =true => routine simplifiee
c ---  choix routine radiative i.r.    simpir=true => routine simplifiee
c
      hydroa = .false.
c ... hydroa = .true. => cycle hydrologique d'apres ohring & adler, 1978
c
      gcl = 240.
c     gcl = 238.5
c ... gcl : coefficient de proportionalite
c           nebulosite - precipitations (w.m-2) (ohring et adler, 1978)
c
c **********************************************************************
c
c --- initialisation
c      des climatologies propres aux flux de chaleur verticaux et
c      des extensions des divers secteurs en surface
c
      if (it.eq.0) then
                      call inira1
                      call iniso(facco2)
                      call iniir
                      call inihyd
       if (   tsclim) call initsz
       if (.not.inpu) call iniaja
      end if
c
c **********************************************************************
c pour exp. co2
c co2  ppmgz2= 341.396+1.539*apal*exp(0.009173*apal)
c co2  ppmgz2= 330.000
c pour exp. co2
c **********************************************************************
c
c --- insolation au sommet atmosphere
c
      call insol
c
c **********************************************************************
c
c --- temperature observee en surface
c
      if (tsclim) call tsolzv
c
c **********************************************************************
c
c --- mise a jour climatologie utilisee pour flux de chaleur verticaux
c
      if (it.ne.0) then
c
c **********************************************************************
c
c
       if (season) call clmrad
                   call reflio
                   call reflls(3,4)
       do 510 ic=5,nilw
                   call reflca(ic,ic)
 510   continue
c
c ----------------------------------------------------------------------
c
c --- mise a zero des variables flux vericaux moyennes zonales
c
       do 30 j=1,jp
       do 31 n=1,6
       htshv1(n,j)=0.
   31  continue
       do 32 n=1,3
       htady2(n,j)=0.
   32  continue
       do 33 n=5,6
       htady2(n,j)=0.
   33  continue
       do 34 iilw=1,nilw
                                 c3ir4(j,iilw) = 0.
c      if (ajsu0(j,iilw).eq.0.) hswsu3(j,iilw) = 0.
   34  continue
       alb   (j)  = 0.
       saso3 (j)  = 0.
       sgso3 (j)  = 0.
       fupir3(j)  = 0.
       fdnir3(j)  = 0.
       fnir3(j)   = 0.
       t4su4(j)   = 0.
   30  continue
c
c ----------------------------------------------------------------------
c
c --- flux verticaux de chaleur : decoupage
c     glace marine - ocean - continent - glace continentale
c
       do 101 j=1,jp
       t4su4(j) = 0.
c
       do 100 iilw=1,nilw
       call evarad(j,iilw)
c ...  calcul de l'humidite relative en surface
c
       if (ajsu0(j,iilw).gt.0.) then
c
c -------------------------------------------------------------------
c *126* call iniras(dt0,j,iilw,nl4)
        call iniras(dt0,j,iilw,ntrop,nl4)
c -------------------------------------------------------------------
c
        if (cmuso1(j).ne.0.) then
            nl4m=nl4-1
c -------------------------------------------------------------------
c
         if (simps) then
          call sola3g (ysol,j,iilw,nl4m)
         else
c *126*   call solari (ysol,j,iilw,nl4m)
          call solari (ysol,j,iilw,ntrop,nl4m)
         end if
c
c -------------------------------------------------------------------
          alb  (j)    = alb (j) +(1.-atso2(j)-gtso2(j))*ajsu0(j,iilw)
        else
          alb  (j)    = 0.
        end if
c
          chhv0(1,j,iilw) = gtso2(j)*sotso1(j)
          cshv0(1,j,iilw) = 0.
          htady2(1,j) = htady2(1,j)+ atso2(j)*sotso1(j)*ajsu0(j,iilw)
          saso3(j)    = saso3(j)   + atso2(j)*sotso1(j)*ajsu0(j,iilw)
          sgso3(j)    = sgso3(j)   + gtso2(j)*sotso1(j)*ajsu0(j,iilw)
c
        if ((.not.tsclim).or.iilw.le.2) then
         if (simpir) then
          call infras(dt0,j,iilw)
         else
          call infra(dt0,nl4)
         end if
        end if
c
        if (mod(j-1,4).eq.0.and.fulwri) then
         write(iwr1,72)j
   72    format(/,i4,2x,'pmb',5x,'zkm',5x,'tt',10x,'hh',
     .         10x,'oz',5x,'cld',4x,'firm',4x,'fird',4x,'firn')
         write(iwr1,73)(pmbra1(i),zkmra1(i),ttra1(i),hhra1(i),ozra1(i),
     .                cldra1(i),fimir1(i),fidir1(i),finir1(i),i=1,nl4)
   73    format(f9.2,f8.3,f8.2,2x,e10.3,2x,e10.3,f8.3,3f8.2)
        end if
c
        divir2(j)=finir1(nl4)-finir1(1)
c _cess divir2(j)=finir1(ntrop)-finir1(1)
        fupir2(j)=fimir1(nl4)
         fnir2(j)=finir1(1)
        fdnir2(j)=fidir1(1)
        chhv0(2,j,iilw) = fdnir2(j) + stfir0*3.*t4ssu0(j,iilw)**4
        cshv0(2,j,iilw) =             stfir0*4.*t4ssu0(j,iilw)**3
        htady2(2,j) = htady2(2,j) + divir2(j)*ajsu0(j,iilw)
        fupir3(j)   = fupir3(j)   + fupir2(j)*ajsu0(j,iilw)
        fdnir3(j)   = fdnir3(j)   + fdnir2(j)*ajsu0(j,iilw)
         fnir3(j)   =  fnir3(j)   +  fnir2(j)*ajsu0(j,iilw)
c
        call c3rad(j,iilw)
        call hlat(u4dy8(j),u4dy8(j+1),rhzev1(j,iilw),t4ssu0(j,iilw),
     .            hse,hsf,j,iilw)
c
        call hsens(t4ssu0(j,iilw),hsb,hsc,hse,hsf,j,iilw)
        chhv0(3,j,iilw) = hsb*t2dy2(j)-hsc
        cshv0(3,j,iilw) = hsb
c
        chhv0(4,j,iilw)=hswsu3(j,iilw)*(hsf+hse*chhv0(3,j,iilw))
        cshv0(4,j,iilw)=hswsu3(j,iilw)* hse*hsb
c
       else
        do 102 n=1,6
            chhv0(n,j,iilw) = 0.
            cshv0(n,j,iilw) = 0.
  102   continue
       end if
c
       t4su4(j)=t4su4(j) + ajsu0(j,iilw) * t4asu4(j,iilw)
  100  continue
  101  continue
c
      end if
c
c **********************************************************************
c
c --- determination des nouvelles temperatures en surface
c
c  -- a) glace marine et ocean
c
      call tocean
c
c  -- b) continent libre ou recouvert de neige
c
      call tland(3,4)
c
c  -- c) calotte groenland
c
      call tcalot(5)
c
c  -- c) calotte laurentide
c
      call tcalot(6)
c
c  -- e) calotte finno-scandinave
c
      call tcalot(7)
c
c --- incrementation bilan des temperatures de chaque secteur
c
      do 115 j=1,jp
      do 115 iilw=1,nilw
      t4msu1(j,iilw) = t4msu1(j,iilw) + t4ssu0(j,iilw)
      if (t4ssu0(j,iilw).ne.0.) bjmsu1(j,iilw) = bjmsu1(j,iilw) + 1.
      ajmsu1(j,iilw) = ajmsu1(j,iilw) +  ajsu0(j,iilw)
  115 continue
c
c ----------------------------------------------------------------------
c
c --- moyennes zonales des temperatures
c
      do 120 j=1,jp
      t4su0(j) = 0.
      do 120 iilw=1,nilw
      t4su0(j) = t4su0(j) + ajsu0(j,iilw)*t4ssu0(j,iilw)
  120 continue
      do 122 j=1,jp
      tmsu1(j) = tmsu1(j) + t4su0(j)
  122 continue
c
c ----------------------------------------------------------------------
c
c --- bilans moyens zonaux
c
      do 130 j=1,jp
          shshv1(j) =           htshv1(1,j)
      do 130 n=2,6
          shshv1(j) = shshv1(j)+htshv1(n,j)
  130 continue
      do 131 j=1,jp
       sens(j)=htshv1(3,j)
       evap(j)=htshv1(4,j)
       hshl(j)=sens(j)+evap(j)
  131 continue
c
c ----------------------------------------------------------------------
c
c --- incrementation des bilans moyens annuels
c
      do 132 j=1,jp
      do 133 iilw=1,nilw
      do 134 n=1,6
      htshv2(n,j,iilw) = htshv2(n,j,iilw) + htshv0(n,j,iilw)
  134 continue
  133 continue
      do 135 n=1,6
      zhshv2(n,j     ) = zhshv2(n,j     ) + htshv1(n,j)
  135 continue
      shshv2(  j     ) = shshv2(  j     ) + shady2(  j)
      shahv2(  j     ) = shahv2(  j     ) + shady2(  j)
c
      if (saso3(j).gt.0.) then
       albhv2(  j     ) = albhv2(  j     ) + alb   (  j)
       ialhv3(  j     ) = ialhv3(  j     ) + 1
      endif
       sahv2(  j     ) =  sahv2(  j     ) +  saso3(  j)
      fuphv2(  j     ) = fuphv2(  j     ) + fupir3(  j)
      fdnhv2(  j     ) = fdnhv2(  j     ) + fdnir3(  j)
      hslhv2(  j     ) = hslhv2(  j     ) + hshl  (  j)
       t2hv2(  j     ) =  t2hv2(  j     ) +  t2dy2(  j)
  132 continue
c
c ----------------------------------------------------------------------
c
c --- bilan au sommet de l'atmosphere
c
       do 150 j=1,jp
       solag (j) = sgso3(j) + saso3(j)
       topbi0(j) = sgso3(j) + saso3(j) - fupir3(j)
       radbi1(j) = radbi1(j) + topbi0(j)
  150  continue
       call av(sotso1,sohm ,cp)
       call av(solag ,saghm,cp)
       call av(topbi0,topbhm,cp)
       aghm      = 1. - saghm/sohm
       albham    = albham + aghm
       topbi1    = topbi1 + topbhm
c
c ----------------------------------------------------------------------
c
c --- moyennes hemispheriques des flux de chaleur verticaux
c
      call av(t4su0,t4hsu1,cp)
      t4ham = t4ham + t4hsu1
      dt4ep = dt4ep + t4su0(1) - t4su0(jp)
      dt4t2 = dt4t2 + t4hsu1 - t2hdy2
      t2ham = t2ham + t2hdy2
c
      call av(fupir3,fuphm,cp)
      fupham=fupham+fuphm
      call av(fdnir3,fdnhm,cp)
      fdnham=fdnham+fdnhm
c
      call av(sens,senshm,cp)
      shham = shham + senshm
      call av(evap,alhhm,cp)
      alhham= alhham + alhhm
c
c ----------------------------------------------------------------------
c
c --- chauffage diabatique atmospherique
c
c --- calcul des precipitations
c
      if (hydroa) then
       do 142 j = 1,jp
       evapa(j) = -htshv1(4,j)+gcl*(clcl1(j)-clhcl1)
       qlql2(j) =  evapa(j) * 86400.e3 /hvsth4
c ...  qlql2    : precipitations, d'apres ohring et adler, 1978.
c
 142   continue
      else
       call precip(alhhm,evapa)
c ...  cycle hydrologique d'apres les observations de jaeger, 1976
       call hadley
c ...  calcul du transport de chaleur par la cellule de hadley
c
      end if
c
       do 140 j=1,jp
          shady2(j)   = 0.
          htady2(2,j) = -htady2(2,j)
          htady2(3,j) = -htshv1(3,j)
          htady2(4,j) =  evapa(j)
c
      do 141 n=1,6
          shady2(j)   =  shady2(j) + htady2(n,j)
  141 continue
          atmbi0(j) = atmbi0(j) +shady2(j)
          htdy2(j)  = q2dy0 *crfdy0 *gdy0 *1.e-5 *shady2(j) /cpath1
c ...     le facteur gdy0 * 1.e-5 permet le passage des w.m-2 aux w.kg-1
  140 continue
c
      call av(shady2,htahm,cp)
      call av(shshv1,htshm,cp)
      bitot=topbhm-htahm-htshm
c
c ----------------------------------------------------------------------
c
c --- print at middle of the month
c
      if (iprint.eq.1) go to 41
      if (mts.ne.0)    go to 40
      if (monpr.gt.mm) go to 40
   41 continue
c ------- pour exp. co2  thierry     ------ debut ----------------------
co2      write(iwr7,50) apal,ppmgz2
co2      write(iwr7,51) (t4su0(j),j=1,jp)
co2      write(iwr7,52) (t4ssu0(j,2),j=1,jp)
co2      write(iwr7,53) (t4ssu0(j,3),j=1,jp)
co2      write(iwr7,54) (t4ssu0(j,4),j=1,jp)
50    format(5x,'annee =',f6.2,3x,'co2 =',f8.3)
51    format (1x,'t4su0 (j)   =',8f8.3,/,10f8.3)
52    format (1x,'t4ssu0(j,2) =',8f8.3,/,10f8.3)
53    format (1x,'t4ssu0(j,3) =',8f8.3,/,10f8.3)
54    format (1x,'t4ssu0(j,4) =',8f8.3,/,10f8.3)
c ------- pour exp. co2  thierry     ------ fin   ----------------------
c
      write(iwr1,42)t4hsu1
   42 format (
     .//,' moyennes des temperatures et flux de chaleur verticaux :',
     .//,' temperature          a 1000mb                      : ',f12.2)
c
      if (it   .eq.0 ) go to 40
c
      write(iwr1,43)t2hdy2,clhcl1,aghm,fuphm,fdnhm,senshm,alhhm,
     .            htshm,htahm,topbhm,bitot
   43 format (
     .  ' temperature          a  500mb                      : ',f12.2,
     ./,' nebulosite moyenne hemispherique                   : ',f12.4,
     ./,' albedo planetaire                                  : ',f12.4,
     ./,' flux i.r. montant    a    0mb                      : ',f12.2,
     ./,' flux i.r. descendant a 1000mb                      : ',f12.2,
     ./,' flux chaleur sensible (- = vers l`atmosphere)      : ',f12.2,
     ./,' flux chaleur latente  (- = vers l`atmosphere)      : ',f12.2,
     ./,' bilan surface         (+ = absorbe par la surface) : ',f12.3,
     ./,' bilan atmosphere      (+ = absorbe par atmosphere) : ',f12.3,
     ./,' bilan sommet atmosphere    (+ = rentrant )         : ',f12.3,
     ./,' bilan energetique total    (doit etre nul)         : ',f12.3)
c
      write(iwr1,44)
   44 format(/,3x,'saso3',3x,'sgso3',5x,'alb',2x,'fupir3',3x,'fnir3',
     . 2x,'fdnir3',4x,'sens',4x,'evap',3x,'hs+hl',4x,'shts',
     . 6x,'ql',3x,'hta 4',4x,'shta','  t4 zon')
      write(iwr1,45)
     .(saso3(j),sgso3(j),alb(j),fupir3(j),fnir3(j),fdnir3(j),
     . sens(j),evap(j),hshl(j),shshv1(j),qlql2(j),htady2(4,j),shady2(j),
     . t4su0(j),j=1,jp)
   45 format((2f8.2,f8.4,11f8.2))
c
      do 48 iilw=1,nilw
      do 49 j=1,jp
                               t4s(j,iilw)=t4ssu0(j,iilw)
      if (ajsu0(j,iilw).eq.0.) t4s(j,iilw)=0.
   49 continue
      write(iwr1,46)iilw,(n,n=1,6)
   46 format(/' secteur',i2,6x,'flux de chaleur verticaux',
     . ' et horizontaux (w.m-2)',16x,'w',3x,'h (%)  c3 (i.r.)',3x,
     . 'albedo  temper.   fraction',/,10x,6i10)
      write(iwr1,47)(deggr1(j),(htshv0(n,j,iilw),n=1,6),hswsu3(j,iilw),
     . rhzev1(j,iilw),c3ir4(j,iilw),rgsrf0(j,iilw),t4s(j,iilw),
     . ajsu0(j,iilw),j=1,jp)
   47 format((7f10.2,f11.3,f8.3,f11.2,f9.3,f10.2,f10.5))
   48 continue
c
   40 continue
c
      if (mts.eq.0)    then
       write(iwr7,471)((rgsrf0(j,iilw),j=1,jp),iilw,mm,iilw=1,nilw)
  471  format((9f8.3,/,9f8.3,'  rg',i2,'   month',i3))
       write(iwr7,472)((hfasn6(j,iilw),j=1,jp),iilw,mm,iilw=4,nilw)
  472  format((9f8.3,/,9f8.3,' hf*',i2,'   month',i3))
       write(8,450)(sotso1(j),j=1,jp)
  450  format(9f8.2,/,9f8.2,'  sotso1')
       write(8,451)(saso3(j),j=1,jp)
  451  format(9f8.2,/,9f8.2,'  saso3')
       write(8,452)(sgso3(j),j=1,jp)
  452  format(9f8.2,/,9f8.2,'  sgso3')
       write(8,453)(alb(j),j=1,jp)
  453  format(9f8.4,/,9f8.4,'  alb')
       write(8,454)(fupir3(j),j=1,jp)
  454  format(9f8.2,/,9f8.2,'  fupir3')
       write(8,455)(fnir3(j),j=1,jp)
  455  format(9f8.2,/,9f8.2,'  fnir3')
       write(8,456)(fdnir3(j),j=1,jp)
  456  format(9f8.2,/,9f8.2,'  fdnir3')
       write(8,458)(htshv1(2,j),j=1,jp)
  458  format(9f8.2,/,9f8.2,'  h2z')
       write(8,459)(sens(j),j=1,jp)
  459  format(9f8.2,/,9f8.2,'  sens')
       write(8,460)(evap(j),j=1,jp)
  460  format(9f8.2,/,9f8.2,'  evap')
       write(8,461)(qlql2(j),j=1,jp)
  461  format(9f8.2,/,9f8.2,'  ql')
       write(8,462)(t4su0(j),j=1,jp),mm
  462  format(9f8.2,/,9f8.2,'  t1000mb   month',i3)
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  prtf====
c cray .endif
      subroutine prtflu
c
c +++ bcm03.for +++
c
      parameter(np=18,npp=19)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr2/argr2(np),atgr2
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
c
       common/so3/saso3(np),sgso3(np)
       common/ir3/fupir3(np),fdnir3(np),fnir3(np)
c
       common/oc6/otoc6(npp),otmoc6(npp),flooc6(npp)
c
       common/dy6/atdy6(npp),atmdy6(npp),flady6(npp),azdy6,akzdy6
c
       dimension rbt(np),rbg(np),flt(npp)
c
       do 200 j=1,jp
       rbt(j)=saso3(j) +sgso3(j) -fupir3(j)
  200  rbg(j)=sgso3(j) -fnir3(j)
       call av(rbt,rbtav,cp)
       call av(rbg,rbgav,cp)
c
       do 220 j=1,jpp
       flt  (j)=atdy6(j)+flooc6(j)
  220  continue
      if (iprint.eq.1) go to 196
      if (mts.ne.0   ) go to 199
      if (monpr.gt.mm) go to 196
  196 continue
       write(iwr1,230) it,ian,mm,iday
  230 format(//'   -- prtflu --  ',
     . 3x,'it =',i6,3x,i3,'e annee',3x,i2,'e mois',3x,i4,'e jour',
     . //'     j        rbt(w/m2)   rbg(w/m2)     fla(pw)',
     . 8x,'at(w)',8x,'ot(w)',7x,'flt(w)')
  240  continue
       write(iwr1,250)(j,rbt(j),rbg(j),flady6(j),atdy6(j),
     .                   flooc6(j),flt(j),j=1,jp)
  250  format(i12,6e12.4)
       write(iwr1,251) rbtav,rbgav
  251  format(/' mean : ',4x,2e12.4)
  199  continue
       return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tszi====
c cray .endif
      subroutine initsz
c
c +++ bcm04.for +++
c
      parameter (np=18)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
c
      common/su6/t4ssu6(12,np)
c
      dimension t4s0(12,18),x(18),y(np)
c
      data ((t4s0(m,j),m=1,12),j=1,9)/
     . 299.20,299.70,299.90,300.10,299.90,299.55,
     . 299.00,298.95,299.15,299.20,299.30,299.25,
     . 298.90,299.35,299.80,300.25,300.30,300.00,
     . 299.65,299.50,299.60,299.60,299.45,299.05,
     . 297.90,298.35,299.25,300.05,300.75,300.55,
     . 300.25,300.05,300.20,300.00,299.30,298.30,
     . 295.80,296.30,297.65,298.90,300.30,300.70,
     . 300.60,300.60,300.55,299.75,298.25,296.60,
     . 292.80,293.30,295.20,296.95,298.85,300.20,
     . 300.55,300.70,300.30,298.60,296.35,294.00,
     . 288.70,289.25,291.45,293.80,296.20,298.50,
     . 299.65,299.65,298.85,296.35,293.15,290.30,
     . 283.80,284.45,286.65,289.60,292.45,295.70,
     . 297.80,297.90,296.30,292.85,288.70,285.60,
     . 279.10,280.00,282.35,286.05,289.30,293.00,
     . 295.90,295.95,293.60,289.20,284.40,280.90,
     . 273.80,275.05,277.80,282.55,286.50,290.25,
     . 293.35,293.40,290.50,285.55,279.95,275.90/
      data ((t4s0(m,j),m=1,12),j=10,18)/
     . 268.20,269.60,273.05,278.90,283.70,287.55,
     . 290.40,290.40,286.95,281.65,275.35,270.50,
     . 264.10,265.35,269.40,276.00,281.45,285.65,
     . 288.45,287.95,284.25,278.50,271.70,266.35,
     . 260.15,261.25,265.85,273.00,279.10,284.10,
     . 287.05,286.15,282.05,275.65,268.00,262.40,
     . 254.15,255.65,260.55,268.75,276.55,282.60,
     . 285.95,284.45,279.45,271.65,262.45,256.40,
     . 249.95,250.50,254.55,262.55,271.80,278.90,
     . 282.80,281.20,276.05,267.35,257.70,252.00,
     . 248.85,248.15,250.25,256.80,266.55,274.25,
     . 278.05,276.70,272.30,264.10,255.15,250.80,
     . 246.50,246.15,246.45,252.85,263.80,271.85,
     . 275.35,273.90,269.10,260.30,252.35,248.20,
     . 241.95,242.95,243.25,249.80,261.95,271.25,
     . 274.70,273.15,265.75,256.35,249.20,244.85,
     . 238.00,241.00,241.90,248.10,260.75,271.05,
     . 274.20,272.90,263.20,253.80,246.30,243.00/
c ------------------------------------------------------------
c  source : s.g.warren & s.h.schneider, p.1379, 1979 (gts.obs)
c  modification : interpolation 90-85, 85-80 ...
c ------------------------------------------------------------
c
      do 1 i=1,12
c
      do 2 j=1,18
      x(j)=t4s0(i,j)
    2 continue
c
      call pol(x,y,18,np)
c
      do 5 j=1,np
      t4ssu6(i,j)=y(j)
    5 continue
    1 continue
c
      write(1,60)(i,i=1,12)
   60 format(//'   -- initsz -- ',//,
     .   6x,' temperature moyenne zonale',/,7x,26('-'),
     . /,6x,'j',4x,12i8)
      write(1,61)(j,(t4ssu6(i,j),i=1,12),j=1,np)
   61 format((i7,4x,12f8.2))
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tszc====
c cray .endif
      subroutine tsolzv
c
c +++ bcm05.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables surfaces
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/su3/hswsu3(np,nilw),cksu3(np,nilw)
c
      common/su5/t4ssu5(np)
      common/su6/t4ssu6(12,np)
c
c **********************************************************************
c
c --- parametres temporels de l'interpolation
c
      ita=int(12.*tsm)
      itsm=int(tsm)
      itm=mod(it,ita)
      itj=mod(itm,itsm)
      i  =1+itm/itsm
      ip1=mod(i+1,12)
      if(ip1.eq.0) ip1=12
c
c --- interpolation des temperatures
c
      do 1 j   =1,jp
      t4ssu5(j)=t4ssu6(i,j)+itj*(t4ssu6(ip1,j)-t4ssu6(i,j))/tsm
c
c **********************************************************************
c
      do 2 iilw=1,nilw
      if (ajsu0(j,iilw).gt.0.) then
       t4ssu0(j,iilw) = t4ssu5(j)
      else
       t4ssu0(j,iilw) = 0.
      end if
 2    continue
c
c **********************************************************************
c
    1 continue
c
c **********************************************************************
c
c --- print
c
      if (iprint.eq.1) go to 601
      if (mts.ne.0) go to 600
      if (monpr.gt.mm) go to 600
  601 continue
      write(1,602)(t4ssu5(j),j=1,jp)
  602 format(//'   -- tsolzv --',//,' surface temperature :',2(/,9f8.3))
  600 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tazi====
c cray .endif
      subroutine initaz
c
c +++ bcm06.for +++
c
      parameter (np=18,npp=19)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
c
      common/dyc/t2dyc(12,np),u4dyc(12,npp)
c
      dimension t4a0(12,18),u4a0(12,19),x(18),y(np),x1(19),y1(npp)
c
      data ((t4a0(m,j),m=1,12),j=1,9)/
     . 267.60,267.80,267.85,267.90,267.90,267.45,
     . 267.15,267.00,267.35,267.50,267.60,267.60,
     . 267.55,267.85,267.80,267.75,267.80,267.40,
     . 267.10,267.10,267.45,267.60,267.70,267.70,
     . 267.05,267.30,267.25,267.25,267.50,267.40,
     . 267.15,267.25,267.55,267.55,267.45,267.30,
     . 265.70,265.80,265.90,266.15,266.85,267.40,
     . 267.30,267.45,267.50,267.20,266.60,266.10,
     . 263.20,263.20,263.50,264.35,265.70,267.10,
     . 267.40,267.55,267.30,266.25,264.90,263.95,
     . 259.55,259.55,260.15,261.95,264.00,266.25,
     . 267.35,267.40,266.70,264.55,262.20,260.70,
     . 255.45,255.30,256.40,259.00,261.75,264.60,
     . 266.65,266.70,265.25,262.15,258.90,256.75,
     . 251.25,250.95,252.55,255.60,259.00,262.30,
     . 265.00,265.15,262.90,259.25,255.40,252.65,
     . 247.40,247.05,249.00,252.35,256.25,259.85,
     . 262.85,262.95,260.10,256.25,252.05,248.95/
      data ((t4a0(m,j),m=1,12),j=10,16)/
     . 244.20,243.80,245.90,249.45,253.60,257.50,
     . 260.65,260.50,257.20,253.25,248.90,245.75,
     . 241.60,241.30,243.30,246.75,251.00,255.40,
     . 258.70,258.15,254.50,250.40,245.95,243.10,
     . 239.50,239.40,241.10,244.25,248.75,253.65,
     . 257.05,256.15,252.15,247.80,243.30,241.00,
     . 237.70,237.80,239.10,242.00,246.80,252.15,
     . 255.50,254.40,250.10,245.40,240.95,239.15,
     . 236.05,236.30,237.20,239.90,245.00,250.70,
     . 254.10,252.75,248.30,243.15,238.90,237.50,
     . 234.65,234.95,235.45,238.05,243.45,249.30,
     . 252.85,251.25,246.80,241.05,237.30,236.10,
     . 233.50,233.75,233.85,236.50,242.10,248.05,
     . 251.70,249.90,245.55,239.05,236.05,234.80/
      data u4a0 /
     .  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  
     . -2.5,-2.4,-2.2,-2.0,-1.4,-1.1,-1.0,-0.6,-0.8,-1.2,-1.7,-2.2,
     . -3.7,-3.7,-3.3,-2.9,-2.2,-1.4,-1.1,-0.6,-0.9,-1.7,-2.8,-3.6,
     . -4.0,-4.1,-3.7,-3.4,-2.9,-2.3,-2.1,-1.5,-1.8,-2.4,-3.6,-4.2,
     . -3.1,-3.1,-3.2,-3.1,-2.9,-2.7,-2.6,-2.2,-2.4,-2.8,-3.5,-3.5,
     . -1.4,-1.3,-1.8,-2.1,-2.1,-2.0,-2.1,-2.1,-2.2,-2.2,-2.5,-2.0,
     .  0.7, 0.8, 0.2,-0.4,-0.6,-0.7,-1.0,-1.3,-1.5,-0.9,-0.9, 0.2,
     .  2.3, 2.4, 1.9, 1.4, 1.0, 0.8, 0.2,-0.1,-0.1, 0.5, 0.9, 2.0,
     .  2.7, 2.8, 2.4, 2.1, 1.6, 1.5, 1.0, 0.9, 1.1, 1.3, 2.1, 2.7,
     .  2.5, 2.4, 2.0, 2.0, 1.5, 1.7, 1.5, 1.6, 2.0, 2.1, 2.7, 2.9,
     .  2.0, 1.7, 1.3, 1.3, 1.1, 1.3, 1.4, 1.8, 2.1, 2.3, 2.6, 2.6,
     .  1.4, 1.3, 0.8, 0.7, 0.6, 0.8, 1.0, 1.5, 1.8, 2.0, 2.1, 1.8,
     .  1.0, 0.8, 0.6, 0.3, 0.1, 0.3, 0.4, 0.9, 1.2, 1.1, 1.3, 1.0,
     .  0.6, 0.5, 0.4, 0.1,-0.2,-0.1, 0.1, 0.2, 0.5, 0.6, 0.8, 0.4,
     .  0.2, 0.0,-0.2,-0.3,-0.4,-0.5,-1.0,-0.4, 0.0, 0.0, 0.5, 0.1,
     . -0.1,-0.2,-0.6,-0.6,-0.6,-0.6,-0.7,-0.4,-0.1,-0.3, 0.5,-0.1,
     . -0.3,-0.3,-1.0,-0.9,-0.7,-0.7,-0.3,-0.4,-0.2,-0.5, 0.5,-0.2,
     . -0.1,-0.2,-0.5,-0.4,-0.3,-0.3,-0.1,-0.2,-0.1,-0.2, 0.2,-0.1,
     .  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
c ...  source : oort and rasmusson,1971 + interpolation a 75 et 85 deg.n
c
      do 1 i=1,12
c
      t4a0(i,17) = t4a0(i,16) - .75 *(t4a0(i,15)-t4a0(i,16)) 
      t4a0(i,18) = t4a0(i,17) - .30 *(t4a0(i,15)-t4a0(i,16)) 
c
      do 2 j=1,18
      x  (j)=t4a0(i,j)
      x1 (j)=u4a0(i,j)
    2 continue
      x1(19)=u4a0(i,19)
c
      call    pol(x ,y ,18,np)
      call intpol(x1,y1,19,npp)
c
      do 5 j=1,np
      t2dyc(i,j)  =y (j)
      u4dyc(i,j)  =y1(j)
    5 continue
      u4dyc(i,npp)=y1(npp)
    1 continue
c
      write(1,60)(i,i=1,12)
   60 format(//'   -- initaz -- ',//,
     .   6x,' temperature moyenne zonale',/,7x,26('-'),
     . /,6x,'j',4x,12i8)
      write(1,61)(j,(t2dyc(i,j),i=1,12),j=1,np)
   61 format((i7,4x,12f8.2))
      write(1,62)(i,i=1,12)
   62 format(                     //,
     .   6x,' vent moyen zonal',/,7x,16('-'),
     . /,6x,'j',4x,12i8)
      write(1,61)(j,(u4dyc(i,j),i=1,12),j=1,npp)
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tazc====
c cray .endif
      subroutine tatmzv
c
c +++ bcm07.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables dynamiques
c
      common/dy2/t2ady2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .                 pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
      common/dy8/u4dy8(npp),u4mdy8(npp),
     .          th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
      common/dyc/t2dyc(12,np),u4dyc(12,npp)
c
c **********************************************************************
c
c --- parametres temporels de l'interpolation
c
      ita=int(12.*tsm)
      itsm=int(tsm)
      itm=mod(it,ita)
      itj=mod(itm,itsm)
      i  =1+itm/itsm
      ip1=mod(i+1,12)
      if(ip1.eq.0) ip1=12
c
c --- interpolation des temperatures
c
      do 1 j   =1,jp
      t2dy2(j)=t2dyc(i,j)+itj*(t2dyc(ip1,j)-t2dyc(i,j))/tsm
      u4dy8(j)=u4dyc(i,j)+itj*(u4dyc(ip1,j)-u4dyc(i,j))/tsm
    1 continue
      u4dy8(jpp)=0.
c
c **********************************************************************
c
c --- print
c
      if (iprint.eq.1) go to 601
      if (mts.ne.0) go to 600
      if (monpr.gt.mm) go to 600
  601 continue
      write(1,602)(t2dy2(j),j=1,jp)
  602 format(//'   -- tatmzv --',//,' atmosp. temperature :',2(/,9f8.3))
  600 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  dy0=====
c cray .endif
      subroutine inidy0
c
c +++ bcma0.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      double precision  cpwth2,rhwth2,ocath2
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
c
c --- constantes physiques
c
       common/th1/cpath1,rath1,akth1
c
c --- variables dynamiques
c
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/p00/p4p00,p3p00,p2p00,p1p00
       common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
       common/dy3/copdy3(np),comdy3(np),coedy3(np),idy3
       common/dy7/zsdy7(np) ,ztdy7(np),
     .            usdy7(npp),utdy7(npp),udy7(npp,nhd)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
       common/th2/cpwth2,rhwth2,ocath2
       common/oc7/ceoc7
c
c **********************************************************************
c
c --- initialisation des constantes physiques propres a la dynamique
c
       comdy0 = 1.454e-4
c ...           2 * vitesse angulaire de rotation de la terre    (s-1)
       sigdy0 = 2.80e-6
c ...           parametre de stabilite statique de l'atmosphere  (m4.s2.
         gdy0 =     9.81
c ...           acceleration de gravite                          (m.s-2)
        f0dy0 = 1.e-4
c ...           parametre de coriolis a 43 deg.n                 (s-1)
         phi0 = asin(f0dy0 /comdy0)
        betdy0 = comdy0 * cos(phi0) / ca
c ...           parametre beta pour 43.deg.n
        p2p00 = .5e5
c ...           pression du niveau 2                             (pa)
        q2dy0 = 2.*(f0dy0**2) / (sigdy0*(p2p00**2))
       crfdy0 = rath1 / (2.*f0dy0)
c
        cedy0 = 4.e-6
c ...           parametre de frottement en surface               (  )
c
        ceoc7 = sqrt(cedy0*p2p00/(rhwth2*gdy0))
        write(iwr6,10) ceoc7
 10     format(/,' ceoc7 = ',e12.4,/,1x)
c
        p4p00 = 1.00e5
        p3p00 =  .75e5
        p1p00 =  .25e5
c ...            pression autres niveaux du modele quasigeostrophique (p
c
c **********************************************************************
c
c --- initialisation des variables auxiliaires de la dynamique
c
      idy3 = 1
c ... idy3 : parametre du schema numerique de resolution de la dynamique
c
      write(iwr1,28)
   28 format (//'   -- inidy0 --',//
     .5x,'j',9x,'deg',11x,'fj',10x,'den',7x,'comdy3',7x,'copdy3',
     .7x,'coedy3',8x,'degc',10x,'fjj  ')
      do 30 j=1,jpp
          fjjdy1(j) = comdy0*sin(cp*real(j-1))
c ...     fjjdy1    : parametre de coriolis a chaque frontiere de bande
          betdy1(j) = comdy0*cos(cp*real(j-1)) /ca
c ...     betdy1    : parametre beta        a chaque frontiere de bande
      if (j.eq.jpp) go to 30
          fjdy1(j)  = comdy0*sin(cp*(real(j)-.5))
c ...     fjdy1     : parametre de coriolis au milieu de bande
          den       = ds2*cos(cp*(real(j)-0.5))
          comdy3(j)  = cos(cp*(j-1.0))/den
          copdy3(j)  = cos(cp*(j    ))/den
          coedy3(j)  = copdy3(j)+comdy3(j)+q2dy0
      write(iwr1,29)j,deggr1(j),fjdy1(j),den,comdy3(j),copdy3(j),
     .             coedy3(j),dgcgr1(j),fjjdy1(j)
   29 format(i6,f12.2,5(1x,e12.5),f12.2,1x,e12.5)
   30 continue
      write(iwr1,31)jpp,dgcgr1(jpp),fjjdy1(jpp)
   31 format(i6,77x,f12.2,1x,e12.5)
c
c **********************************************************************
c
c --- initialisation du vent aux frontieres laterales
c
      do 25 j=1,jpp,jp
          usdy7(j)   = 0.
          utdy7(j)   = 0.
           udy7(j,1) = 0.
           udy7(j,2) = 0.
          u4dy8(j)   = 0.
   25 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  dy1=====
c cray .endif
       subroutine inidy1(ic)
c
c +++ bcma1.for +++
c
       parameter(np=18,npp=19,nilw=7)
c _dp implicit double precision (a-h,o-z)
c
       common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
       common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
       common/p01/szp01(np,nilw),spp01(np,nilw)
       common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
c
c --- initialisation des pressions continentales
c
       do 71 j=1,jp
       call ttpp(spp01(j,ic),t2dy2(j),z2a,tt,szp01(j,ic))
   71  continue
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ttpp====
c cray .endif
      subroutine ttpp(pp,t2,z2,tt,zz)
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/th1/cpath1,rath1,akth1
      common/p00/p4p00,p3p00,p2p00,p1p00
      common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
c
c --- constantes
c
       p4k = p4p00 ** akth1
       p2k = p2p00 ** akth1
       p21k= p2p00 **(1.-akth1)
c
c --- calcul
c
       s0  = sigdy0 * p2p00
       sr1k= s0 /((1.-akth1) *rath1)
c
c  -- altitude 500 mb
       z2 =-(rath1/gdy0) * (sr1k*(p4p00-p2p00)
     .     -(t2/p2p00 +sr1k) *p21k *(p4k-p2k) /akth1)
c
c  -- pression niveau courant
       cc = (t2/p2p00 + sr1k) * p21k / akth1
       aa = (gdy0/rath1) * (zz-z2) + sr1k*p2p00 - cc*p2k
       bb = sr1k
       it = 0
       t0 = 288.
       gam= .0065
       tt = t0-gam*zz
       pp = 1.e5 * (tt/t0)**(gdy0/(rath1*gam))
c dowhile
   21  continue
       it = it + 1
c      write(iwr6,22)it,pp
   22  format(i4,2x,e12.5)
       ppo= pp
       pp = ( aa + cc*(1.-akth1)*(pp**akth1) )
     .     / ( bb - cc*akth1*(pp**(akth1-1.)) )
       if (abs(pp-ppo).lt.1.) go to 20
       go to 21
   20  continue
c
c  -- temperature niveau courant
       pp2k  = (pp/p2p00)**akth1
       tt =  t2 * pp2k - sr1k * ( pp - p2p00 * pp2k)
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  dy2=====
c cray .endif
      subroutine inidy2
c
c +++ bcma2.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9,nilw=7)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical corako
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
c
c --- constantes physiques
c
      common/th1/cpath1,rath1,akth1
c
c --- variables dynamiques
c
       common/p00/p4p00,p3p00,p2p00,p1p00
       common/p01/szp01(np,nilw),spp01(np,nilw)
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
       common/dy2/t2ady2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/dy4/qdy4(np,nhd),dqdy4(np,nhd)
       dimension  qdyo(np)    ,dqdyo(npp)   ,zeta(np),vq(np)
       dimension  qdy (np)    ,dqdy (npp)
       dimension  aux1(np)    ,aux2 (np)
       common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
       common/dy7/zsdy7(np) ,ztdy7(np),
     .            usdy7(npp),utdy7(npp),udy7(npp,nhd)
c
      dimension t2a0(9),t2s0(9),t2a(np),u10(10),u30(10),uz(npp)
c     dimension ak1(10),ak2(10)
      dimension ak10(19),ak20(19),akk1(npp),akk2(npp)
c
c **********************************************************************
c
c --- observations moyenne annuelle
c
c      data t2a0 / 269.63, 267.60, 263.97,
c     .            259.21, 252.54, 248.24,
c     .            245.04, 242.05, 239.34 /
c ... ancienne initialisation (mv30)
c
      data t2a0 /-5.6,-6.0, -8.6,-13.9,-20.0,-25.1,-29.0,-32.1,-35.1/
c ... obs. oort (1983) moyenne annuelle
      data t2s0 /-5.6,-6.6,-11.6,-19.8,-27.5,-32.7,-36.3,-39.1,-41.4/
c ... obs. oort (1983) moyenne du mois de janvier
c
c oa  data  ak1            /.55e6,.55e6,.55e6,.50e6,.35e6,1.2e6,2.5e6,
c oa 1    2.7e6,1.5e6,0.0/
c std data  ak1            /.30e6,.30e6,.30e6,.26e6,.18e6,1.2e6,2.5e6,
c std1    2.7e6,1.5e6,0.0/
c p00 data ak10 /0.5e6,0.5e6,0.5e6,0.5e6,0.5e6,0.5e6,0.5e6,
c p00.                 0.4e6,0.3e6,0.6e6,1.2e6,1.8e6,2.3e6,
c p00.                 2.6e6,2.7e6,2.1e6,1.4e6,0.7e6,0.0  /
      data ak10 /0.5e6,0.5e6,0.5e6,0.5e6,0.5e6,0.5e6,0.5e6,
     .                 0.5e6,0.6e6,0.8e6,1.2e6,1.8e6,2.3e6,
     .                 2.6e6,2.7e6,2.1e6,1.4e6,0.7e6,0.0  /
c oa  data  ak2            /1.3e6,1.3e6,1.3e6,1.3e6,1.6e6,6.8e6,5.8e6,
c oa 2    3.6e6,1.5e6,0.0/
c std data  ak2            /1.5e6,1.5e6,1.5e6,1.5e6,1.8e6,6.8e6,5.8e6,
c std2    3.6e6,1.5e6,0.0/
      data ak20 /1.3e6,1.3e6,1.3e6,1.3e6,1.3e6,1.3e6,1.3e6,
     .                 1.3e6,2.6e6,3.5e6,6.5e6,7.4e6,5.8e6,
     .                 4.5e6,3.7e6,2.9e6,2.0e6,1.0e6,0.0  /
c
      data  u10 /-3.8,  .7,13.6,24.8,24.0,16.7,10.5, 7.3, 4.2, 0./
      data  u30 /-3.4,-3.9, 0. , 4.8, 6.6, 6.5, 4.3, 2.7, 1.5, 0./
      corako = .false.
c
      write(iwr1,6)
    6 format(//'   -- inidy2 --')
c
c **********************************************************************
c
c --- initialisation des variables dependantes en cas de non lecture de
c     d'input
c
      call pol(t2a0,t2a,9,np)
      if (season) then
       call pol(t2s0,t2dy2,9,np)
      else
       do 102 j=1,jp
       t2dy2(j)=t2a(j)
 102   continue
      end if
c
      do 120 j=1,jp
      t2a  (j) = t2a  (j) + 273.15
      t2dy2(j) = t2dy2(j) + 273.15
          zsdy7 (j) = 0.
          ztdy7 (j) = 0.
          pstdy2(j) = crfdy0*t2dy2(j)
          qdy4(j,1) = fjdy1(j)-q2dy0*crfdy0*t2dy2(j)
          qdy4(j,2) = fjdy1(j)+q2dy0*crfdy0*t2dy2(j)
  120 continue
c
c --- coefficient de diffusion macroechelle
c
c std  call intpol(ak1,akk1,10,npp)
c std  call intpol(ak2,akk2,10,npp)
       call intpol(ak10,akk1,19,npp)
       call intpol(ak20,akk2,19,npp)
          akody5(jpp,1)=akk1(jpp)
          akody5(jpp,2)=akk2(jpp)
      do 10 j=1,jp
          akody5(j  ,1)=akk1(j)
          akody5(j  ,2)=akk2(j)
   10 continue
c
c --- correction du coefficient de diffusion macroechelle en fonction de
c
      if (corako) then
       jmin  = 1 + int(40./gs)
       do 21 k=1,nhd
       jminc = jmin + np *(2-k)/9
c
c begin case
       go to (210,211) k
 210   continue
       call intpol(u10,uz,10,npp)
       akp = akth1 * ((p1p00/p2p00)**(akth1-2.))
       go to 212
 211   continue
       call intpol(u30,uz,10,npp)
       akp = akth1 * ((p3p00/p2p00)**(akth1-2.))
 212   continue
c end case
c
       do 22 j=jminc,jp
       zeta(j)= - (c10gr1(j+1)*uz(j+1) - c10gr1(j)*uz(j))
     .        /   ( dels * cos(cp*(j-.5)) )
       qdyo(j)= fjdy1(j) + zeta(j)
     .        - q2dy0 *(sigdy0/2.e-6) * crfdy0 * akp * t2a(j)
       qdy (j)= fjdy1(j) + zeta(j)
     .        - q2dy0 *                 crfdy0 * akp * t2a(j)
 22    continue
       do 23 j=1,jp
       dqdyo(j) = (qdyo(j) -qdyo(j-1)) / dels
       dqdy (j) = (qdy (j) -qdy (j-1)) / dels
 23    continue
       do 24 j=2,jp
       aux1(j) = (dqdyo(j-1) + 2.*dqdyo(j) + dqdyo(j+1) ) / 4.
       aux2(j) = (dqdy (j-1) + 2.*dqdy (j) + dqdy (j+1) ) / 4.
 24    continue
       do 25 j=2,jp
       dqdyo(j) = aux1(j)
       dqdy (j) = aux2(j)
 25    continue
       do 26 j=jminc,jp
       vq (j)   = - akody5(j,k) * dqdyo(j)
       akody5(j,k) =    - vq(j) / dqdy (j)
 26    continue
       akody5(1,k) = akody5(2,k)
       akody5(1,k) = akody5(2,k)
       write(iwr1,60)k
 60    format(//'   -- inidy2 -- k=',i1)
       write(iwr1,61)
 61    format(/5x,'j',10x,'fj',8x,'zeta',7x,'dqdyo',8x,'dqdy',
     .                10x,'vq',10x,'ak')
       write(iwr1,62)
     . (j,fjdy1(j),zeta(j),dqdyo(j),dqdy(j),vq(j),akody5(j,k),j=1,jp)
 62    format(i6,6e12.4)
 21    continue
      end if
c
      do 30 k=1,nhd
      do 30 j=1,jpp
      akdy5(j,k)=akody5(j,k)
 30   continue
c
      call shfdyn
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  dyn=====
c cray .endif
      subroutine dyn
c
c +++ bcma3.for +++
c
      
      parameter(np=18,npp=19,nhd=2,nho=9,nilw=7)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
c
c --- variables dynamiques
c
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
       common/dy2/t2hdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/z00/z2z00(np),z3z00(np),t3z00(np)
       common/p00/p4p00,p3p00,p2p00,p1p00
       common/p01/szp01(np,nilw),spp01(np,nilw)
       common/dy3/copdy3(np),comdy3(np),coedy3(np),idy3
       common/dy4/qdy4(np,nhd),dqdy4(np,nhd)
       common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
       dimension a(np,nhd),b(np,nhd),c(np,nhd),fq(np,nhd)
       dimension q1(np),q2(np),dq1(np),dq2(np),qt(np),dp(np)
       dimension a1(np),a2(np),b1(np),b2(np),c1(np),c2(np)
       common/dy7/zsdy7(np) ,ztdy7(np),
     .            usdy7(npp),utdy7(npp),udy7(npp,nhd)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
       dimension  z4dy (np)
       common/dy9/u1jdy9,jjdy9,itjdy9
c
       data c2a/1.2e-6/
       jp1 = jp-1
c
c **********************************************************************
c
c --- print
c
      if (iprint.eq.1) go to 86
      if (mts.ne.0) go to 87
      if (monpr.gt.mm) go to 87
   86 continue
      write(iwr1,84) it,ian,mm,iday,sigdy0
   84 format(//,' ---- dyn    ----',
     . 3x,'it =',i6,3x,i3,'e annee',3x,i2,'e mois',3x,i4,'e jour',
     . /,6x,'st.st.:',e12.4)
   87 continue
c
c --- initialisation
c
      if (it.eq.0) then
                      call inidy0
       if (.not.inpu) call inidy2
       if (   taclim) call initaz
       do 502 iilw=1,nilw
                      call inidy1(iilw)
  502  continue
        write(iwr1,506)(iilw,iilw=1,4),
     .  (deggr1(j),(szp01(j,iilw),spp01(j,iilw),iilw=1,4),j=1,jp)
  506   format(/' altitude-pression de surface',/' secteur',4i18,
     .      18(/,f8.1,4(' (',f7.0,',',f7.0,') ')))
        write(iwr1,507)(iilw,iilw=5,7),
     .  (deggr1(j),(szp01(j,iilw),spp01(j,iilw),iilw=5,7),j=1,jp)
  507   format(/' secteur',3i18,
     .      18(/,f8.1,3(' (',f7.0,',',f7.0,') ')))
c
c **********************************************************************
c
      else
c
c --- calcul de la dynamique
c
c --- coefficient de la matrice tridiagonale
c
       do 34 k=1,2
       do 34 j=1,jp
           a(j,k)=tht*akdy5(j+1,k)*copdy3(j)
           c(j,k)=tht*akdy5(j  ,k)*comdy3(j)
           b(j,k)=a(j,k)+c(j,k)+1.
 34    continue
c
c --- prognostic equation
c
       do 75 j=1,jp
       z4dy (j) = - (tsdy8(j+1)*c10gr1(j+1) - tsdy8(j)*c10gr1(j))
     .      / (ca * (           s10gr1(j+1) -          s10gr1(j)))
 75    continue
       call av(z4dy,z4hm,cp)
       do 76 j=1,jp
       z4dy (j) = z4dy (j) - z4hm
 76    continue
c
       do 77 k=1,2
       do 77 j=1,jp
           th  = htdy2(j)
           t3  = c2a*ztdy7(j)
           z4  =(    zsdy7(j)-2.    *ztdy7(j))
c          z4  =(    zsdy7(j)-1.5   *ztdy7(j))
c          z4  =(    zsdy7(j)-1.5   *ztdy7(j)) *2./3.
c          th4 = -cedy0*z4dy(j)
           th4 = -cedy0*z4
           if (k.eq.1) fq(j,k)=(-th-t3    ) *delt
           if (k.eq.2) fq(j,k)=( th+t3+th4) *delt
       if (j.eq.1.or.j.eq.jp) go to 78
       dqdy4(j,k) =
     .    tht *copdy3(j) *akdy5(j+1,k) *(qdy4(j+1,k)-qdy4(j,k))
     .  - tht *comdy3(j) *akdy5(j  ,k) *(qdy4(j  ,k)-qdy4(j-1,k))
     .  + fq(j,k) + qdy4(j,k)
 78    if (j.eq.1)  dqdy4(j,k) =
     .    tht *copdy3(j) *akdy5(j+1,k) *(qdy4(j+1,k)-qdy4(j,k))
     .  + fq(j,k) + qdy4(j,k)
       if (j.eq.jp) dqdy4(j,k) =
     .  - tht *comdy3(j) *akdy5(j,k)   *(qdy4(j,k)  -qdy4(j-1,k))
     .  + fq(j,k) + qdy4(j,k)
 77    continue
       do 85 j=1,jp
           a1 (j)=a (j,1)
           a2 (j)=a (j,2)
           b1 (j)=b (j,1)
           b2 (j)=b (j,2)
           c1 (j)=c (j,1)
           c2 (j)=c (j,2)
           dq1(j)=dqdy4(j,1)
           dq2(j)=dqdy4(j,2)
           q1(j) = qdy4(j,1)
           q2(j) = qdy4(j,2)
 85    continue
       call pq(a1,b1,c1,dq1,q1,idy3)
       call pq(a2,b2,c2,dq2,q2,idy3)
       do 90 j=1,jp
           qdy4(j,1) = q1(j)
 90        qdy4(j,2) = q2(j)
c **********************************************************************
c
c --- diagnostic equations
c
       do 82 j=1,jp
           zsdy7(j) = 0.5*(qdy4(j,1)+qdy4(j,2))-fjdy1(j)
           qt(j)    = 0.5*(qdy4(j,1)-qdy4(j,2))
           dp(j)    = -qt(j)
 82    continue
       logequ = 1
       call pq(copdy3,coedy3,comdy3,dp,pstdy2,logequ)
       do 92 j=1,jp
           ztdy7(j) = q2dy0*pstdy2(j)+qt(j)
           t2ody2(j)= t2dy2(j)
           t2dy2(j) = pstdy2(j)/crfdy0
        if (j.gt.1)
     .    dt2dy2(j) = t2dy2(j)-t2dy2(j-1)
 92    continue
          dt2dy2(1) = 0.
          dt2dy2(jpp)=0.
       do 99 j=1,jp1
           j1 = j+1
           usdy7(j1)  = (usdy7(j)*c10gr1(j) -
     .            ca*zsdy7(j) *(s10gr1(j+1)-s10gr1(j)) ) / c10gr1(j+1)
           utdy7(j1)  = (utdy7(j)*c10gr1(j) -
     .            ca*ztdy7(j) *(s10gr1(j+1)-s10gr1(j)) ) / c10gr1(j+1)
           udy7(j1,1) = usdy7(j1)+utdy7(j1)
           udy7(j1,2) = usdy7(j1)-utdy7(j1)
 99    continue
c
c --- ajustement du coefficient de diffusion macroechelle
c
       if (taclim) then
        call tatmzv
       else
        call turbul
        call usol
       end if
c
c --- jet
c
       do 12 j=1,jp
c
       if (udy7(j,1).gt.u1jdy9) then
          u1jdy9 = udy7(j,1)
           jjdy9 = j
          itjdy9 = iday
       end if
   12 continue
c
      end if
c
c --- caracteristiques niveau 750 mb
c
      do 100 j=1,jp
      call ttzz(p3p00,t2dy2(j),z2z00(j),t3z00(j),z3z00(j))
  100 continue
c
      call av(t2dy2,t2hdy2,cp)
      call prtdyn
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  pq======
c cray .endif
      subroutine pq(a,b,c,d,r,logequ)
c
c +++ bcma4.for +++
c
c     resolution d'un systeme tridiagonal
c     -a(j) *x(j+1) + b(j) *x(j) + c(j) *x(j-1) = d(j)
c     pour lequel b(j) = a(j) + c(j) + 1  
c      (schema crank-nicholson applique a 
c       l'equation de diffusion ecrite en coordonees spheriques)
c       conditions limites : 
c        a l'equateur : dx / dy = 0. 
c          => b(1) = a(1) + 1 c-a-d b(1) := b(1) - c(1)  => logequ = 1
c        ou encore    :  x(0) fixe
c          => b(1) = a(1) + c(1) + 1                     => logequ = 0
c        au pole      : dx / dy = 0.
c          => b(bp)= c(jp)+ 1 (rem.: a(jp) = 0. par construction)
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
      dimension a(np),b(np),c(np),d(np),r(np),e(np),f(np)
          e(1) = a(1)/(b(1)- logequ*c(1))
      do 10 j=2,jp
   10     e(j) = a(j)/(b(j)-        c(j)*e(j-1))
          f(1) = d(1)/(b(1)- logequ*c(1))
      do 20 j=2,jp
   20     f(j) = (d(j)+c(j)*f(j-1))/(b(j)-c(j)*e(j-1))
          r  (jp) = f(jp)/(1.-e(jp))
          jpm = jp-1
      do 30 j=1,jpm
          m = jp-j
          r  (m) = e(m)* r (m+1)+f(m)
   30 continue
          return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  usol====
c cray .endif
      subroutine usol
c
c +++ bcma5.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
c
c --- variables dynamiques
c
       common/p00/p4p00,p3p00,p2p00,p1p00
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
       common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
       common/dy7/zsdy7(np) ,ztdy7(np),
     .            usdy7(npp),utdy7(npp),udy7(npp,nhd)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
c
      u4dy8(1)=0.
c
      do 1 j=1,jp-1
      z4 =(zsdy7(j)-2. *ztdy7(j))
c     z4 =(zsdy7(j)-1.5*ztdy7(j))
c     z4 =(zsdy7(j)-1.5*ztdy7(j)) *2./3.
      u4dy8(j+1) =
     .   ( u4dy8(j)*c10gr1(j) - ca *(s10gr1(j+1)-s10gr1(j)) *z4 )
     . /            c10gr1(j+1)
c     u4dy8(j)=(1.5*u(j,2)- .5*u(j,1)) *2./3.
      u4mdy8(j) =u4mdy8(j)  + u4dy8(j)
    1 continue
c
      u4mdy8(jp)=u4mdy8(jp) + u4dy8(jp)
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  turb====
c cray .endif
      subroutine turbul
c
c +++  bcma6.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9,nd=3)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,snofdb,inpu,turbu
      logical bransc
c
      common/varl00/season,ocen,ekm,tsclim,taclim,snofdb,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
      common/gr2/argr2(np),atgr2
c
c --- constantes physiques
c
      common/th1/cpath1,rath1,akth1
c
c --- variables rayonnement
c
      common/ir3/fupir3(np),fdnir3(np),fnir3(np)
c
c --- variables flux de chaleur verticaux
c
      common/hv1/htshv1(6,np),shshv1(np)
c
c --- variables dynamiques
c
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
       common/dy2/t2ady2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/z00/z2z00(np),z3z00(np),t3z00(np)
       common/p00/p4p00,p3p00,p2p00,p1p00
       common/dy3/copdy3(np),comdy3(np),coedy3(np),idy3
       common/dy4/qdy4(np,nhd),dqdy4(np,nhd)
       common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
       common/dy6/atdy6(npp),atmdy6(npp),flady6(npp),azdy6,akzdy6
       common/ql3/dqyql3(npp,nhd),omgql3(np)
       common/dy7/zsdy7(np) ,ztdy7(np),
     .            usdy7(npp),utdy7(npp),udy7(npp,nhd)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
c
c --- variables locales reutilisees
c
       common/a60/dqya60(nd)
c
c --- variables locales
c
       dimension zk(np),t1(np),th(np),rho(np),dth2(npp)
       dimension akhb(npp,nhd),gam(npp,nhd),akh(npp)
       dimension aka (npp,nhd)
       dimension sumi1(npp),sumi(np)
       dimension jwsmx(nhd)
c
       dimension w(np)
       dimension dt2dt(np),h2(np),h2d(np)
       dimension pstd(np),angmot(npp)
       dimension akdqu(np),akdq(np)
       dimension utmb2(np)
       dimension htw(np)
       dimension zint(np)
c
c --- variables surface
c
       common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .           ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
c
c **********************************************************************
c
      bransc = .true.
c
      jp1      = jp - 1
      jmih     = 1 + int( 5./gs)
      jmin     = 1 + int(30./gs)
c
      jwsmx(1) = 1 + int(70./gs)
      jwsmx(2) = 1 + int(55./gs)
c
c **********************************************************************
c
c --- variation meridienne moyenne
c      de la vorticite potentielle quasigeostrophique
c
      do 200 k=1,2
      do 212 j=1,jpp
c
      if (j.gt.1.and.j.lt.jpp) then
       dqyql3(j  ,k)=( qdy4(j,k)   -qdy4(j-1,k)  ) / cp
      end if
  212 continue
       dqyql3(1  ,k)=0.
       dqyql3(jpp,k)=0.
c
  200 continue
c
c --- spatial moving average (elimination parasite numerique)
c
      do 202 k=1,nhd
      sumi1(1)= dqyql3(1,k)
      do 203 j=2,jp
      sumi1(j)=(dqyql3(j-1,k)+2.*dqyql3(j,k)+dqyql3(j+1,k))/4.
 203  continue
      sumi1(jpp)=dqyql3(jpp,k)
 202  continue
c
c ... dowhile - determination de la zone dqyql3(j,2) > 0
      jminh = 0
 221  continue
      jminh = 1 + jminh
      if(.not.       jminh   .lt.jpp) go to 220
      if(.not.dqyql3(jminh,2).ge.0. ) go to 220
c ... jminh - 1 : indice derniere frontiere de bande avec dqy2(j,2) > 0
      goto 221
 220  continue
      if (jminh.gt.jp-2) jminh = jp-2
c ... end dowhile
c
      do 222 id=1,nd-1
      dqya60(id)= dqya60(id+1)
 222  continue
      jminp = jminh - 1
      dqya60(nd)=(qdy4(jp1,2)-qdy4(jminp,2)) / ((jp1-jminp)*cp)
      dqy2t     = dqya60(1)
      do 223 id=2,nd-1
      dqy2t     = dqy2t + 2*dqya60(id)
 223  continue
      dqy2t     =(dqy2t +   dqya60(nd)) / (2*(nd-1))
c
c **********************************************************************
c
c --- log temperature potentielle niveau 500mb
c
      do 10 j=1,jp
      th2dy2(j) = akth1*alog(p4p00/p2p00) + alog(t2dy2(j))
 10   continue
      dth2t     = th2dy2(jmin)    -th2dy2(jp)
c
c --- temperature potentielle niveau 500mb
c
      do 11 j=1,jp
      th2dy2(j) = exp(th2dy2(j))
 11   continue
      call av(th2dy2,th2av,cp)
c
c --- calcul du gradient meridien de temperature potentielle a 500mb
c
      dth2(  1)= 0.
      do 12  j=1,jp1
      dth2(j+1)= (th2dy2(j+1)-th2dy2(j)) / cp
 12   continue
      dth2(jpp)= 0.
c
c **********************************************************************
c
c --- parametrisation du transport de la chaleur d'apres branscome, 1983
c
      hatm     = rath1 *t4hsu1 /gdy0
c
      do 50 k=1,2
      do 51 j=1,jp
c
c begin case
      go to (511,512) k
c
  511 continue
      call ttzz(p1p00,t2dy2(j),z2z00(j),t1(j),zk(j))
      th (j)   = t1   (j) * ((p4p00/p1p00)**akth1)
      rho(j)   = p1p00 / (rath1 *t1   (j))
      go to 515
c
  512 continue
      zk (j)   = z3z00(j)
      th (j)   = t3z00(j) * ((p4p00/p3p00)**akth1)
      rho(j)   = p3p00 / (rath1 *t3z00(j))
c
  515 continue
c end case
c
   51 continue
c
      do 52 j=1,jp1
      an2      = sigdy0 * ((gdy0 *rho(j))**2)
      an       = sqrt(an2)
c ... an2      : carre de la frequence de brunt-vaisala
c
      dthdy    =            (th(j+1) -th(j)) / dels
      dthdz    = an2 * .5 * (th(j+1) +th(j)) / gdy0
      agam     = - (fjjdy1(j+1) *dthdy) / (betdy0 *hatm *dthdz)
c
      if (agam.lt.0.02) then
       agam     = 0.02
      end if
c
      gam (j+1,k) = 1. / agam
      dd          = hatm / (1. +gam(j+1,k))
      akhb(j+1,k) =(exp(-(zk(j)/dd))) * (abs(dthdy)) * gdy0 * an
     .            * (dd**2) / (th(j) *(fjjdy1(j+1)**2))
   52 continue
c
   50 continue
c
c **********************************************************************
c
c --- choix de la parametrisation : initialisation prealable de aka
c
      jmatch = jminp + 1
      djmax  = jpp - jmin
      djmin  = jpp - jmatch
      if (djmin.gt.0.) rdj = djmax / djmin
c
      if (bransc) then
c
c --- parametrisation de Branscome pour tenir compte de l'effet beta
c
       do 102 j=2,jmatch
       aka(j,2) =(akhb(j-1,2)+2.*akhb(j,2)+akhb(j+1,2))/4.
c ...  moving average (elimination parasite numerique)
c
       aka(j,1) = aka(j,2) * akody5(1,1) / akody5(1,2)
 102   continue
c
      else
c
c --- parametrisation basee sur les observations de Wiin-nielsen et Sela
c
       do 112 k=1,nhd
       do 113 j=2,jmatch
       aka(j,k) = akody5(j,k)
c ...  calibration de wiin-nielsen et sela (m.w.r.,1971)
c
 113   continue
 112   continue
c
      end if
       do 104 k=1,nhd
       aka(1,k) = aka(2,k)
c ...  condition limite equatoriale sur kq
c
 104   continue
c
      if (jmatch.lt.jpp) then 
       lrak = 1
       do 1041 j=jmatch+1,jpp
       do 1042 k=1,nhd
       dj = rdj *(j -jmatch)
        j1 = jmin + int(dj)
       if (j1.ge.jpp) then
        akj = 0.
       else
        j2 = j1 + 1
        akj= akody5(j1,k) +(dj -j1 +jmin) *(akody5(j2,k) -akody5(j1,k))
       end if
c       
c      aka(j,k) = akj * dth2t / .12
c      aka(j,k) = akj *(dth2t**3) * 579.
c      aka(j,k) =-akj * dqy2t / .06d-4
       aka(j,k) =-akj * dqy2t / .12d-4
c ...  calibration observations de wiin-nielsen et sela (m.w.r.,1971)
c
 1042  continue
 1041  continue
      else
       lrak = 0
       do 1043 j=1,jpp
       do 1043 k=1,nhd
       aka(j,k) = 0.
 1043  continue
      end if
c
      if (turbu.and.lrak.eq.1) call wheat(aka,theat,amu1,jmin,jmatch)
c     if (turbu.and.lrak.eq.1) call bheat(aka,theat,amu1)
c
c **********************************************************************
c
c --- flux moyen hemispherique de chaleur (green)
c
      do 141 j=1,jpp
c
      sumi1(j)= c10gr1(j) * dth2(j) * dth2t * 6.5e7
     .        * s10gr1(j) * c10gr1(j) * c10gr1(j)
c ... parametrisation de wu & white ( 6.5e7 = 2.5e7 * sqrt(3.) *3. /2.)
c
  141 continue
      do 142 j=1,jp
      sumi(j)=.5*(sumi1(j)+sumi1(j+1))
  142 continue
      call av(sumi,tgreen,cp)
      tgreen = - tgreen * q2dy0 * crfdy0 / (2.**akth1)
c
c **********************************************************************
c
c --- moment angulaire total de l'atmosphere
c
      do 340 j=1,jpp
      sumi1(j)=0.
 340  continue
      do 341 k=1,2
      do 342 j=1,jpp
      sumi1(j)= sumi1(j) + c10gr1(j)*aka(j,k)*dqyql3(j,k)
 342  continue
 341  continue
      do 343 j=1,jp
      sumi(j)=.5*(sumi1(j)+sumi1(j+1))
 343  continue
      call av(sumi,amb,cp)
c
c **********************************************************************
c
c --- ecriture des resultats 
c
      if (iprint.eq.1 ) go to 600
      if (mts   .ne.0 ) go to 601
      if (monpr .gt.mm) go to 601
  600 continue
      write(iwr1,602)it
  602 format(//'   -- turbul --',5x,'it =',i6)
      write(iwr1,603)tgreen,dth2t,theat,t2ady2,amu1,amb
  603 format(/5x,'t (green) = ',e12.5,6x,'d theta   = ',e12.5,
     .       /5x,'t (inter) = ',e12.5,6x,'t 500mb   = ',e12.5,
     .       /5x,'amu1      = ',f12.5,6x,'ang.m.bal.= ',e12.5)
      write(iwr1,604)
  604 format(/'  j',9x,'dq1',9x,'dq3',8x,'dth2',8x,'aka1',8x,'aka2',
     . '  gam1  gam2',4x,'akh1 (b)',4x,'akh2 (b)',9x,'deg',9x,'th2')
      write(iwr1,605)(j,
     . (dqyql3(j,k),k=1,2),dth2(j),(aka(j,k),k=1,2),(gam(j,k),k=1,nhd),
     . (akhb(j,k),k=1,nhd),dgcgr1(j),th2dy2(j),j=1,jp)
  605 format((i3,5e12.4,2f6.2,4e12.5))
  601 continue
c
c **********************************************************************
c
      do 402 j=1,jpp
      do 404 k=1,nhd
      akdy5(j,k)=aka(j,k)
 404  continue
c
      if (j.gt.1.and.j.lt.jpp)
     .dpsi= 2.*(pstdy2(j) -pstdy2(j-1)) / (-p2p00*dels)
      tsidy8(j) = ( f0dy0**2) * (akdy5(j,2)-akdy5(j,1)) *dpsi
     .          / (sigdy0*p2p00*cedy0)
       tsdy8(j) = -(akdy5(j,1) *dqyql3(j,1) +akdy5(j,2) *dqyql3(j,2))
     .          /  (ca *cedy0)
      tsody8(j) = tsdy8(j) - tsidy8(j)
 402  continue
c
c **********************************************************************
c
c --- calcul du flux de chaleur tourbillonnaire atmospherique instantann
c
      atdy6(jpp)= 0.
c
      do 60 j=jp,1,-1
      akh   (j)= .5*(akdy5(j,1) + akdy5(j,2))
c     flady6(j)=   1.e-15 *cpath1 *2. *pi *ca *c10gr1(j)
c    . * (p4p00 /gdy0) * akh(j) * utdy7(j) /crfdy0
      flady6(j)= - 1.e-15 *cpath1 *2. *pi *ca *c10gr1(j)
     . * (p4p00 /gdy0) * akh(j) *dt2dy2(j) /dels
c
      dt2dt(j) = (t2dy2(j) - t2ody2(j))/delt
       atdy6(j) = atdy6(j+1) - argr2(j)
     .  * (shady2(j) - htady2(6,j) - dt2dt(j) *cpath1 *p4p00 /gdy0) 
      atmdy6(j) = atmdy6(j) + flady6(j)
   60 continue
      akh (jpp) = 0.
c
c --- calcul de la vitesse verticale omega
c
      do 70 j=1,jp
      h2   (j)=gdy0 *1.e-5 *shady2(j)
c ... le facteur gdy0 * 1.e-5 permet le passage des w.m-2 aux w.kg-1
      ha   = crfdy0 * ( dt2dt(j) - h2(j)/cpath1 )
                   adv = 0.
      if (j.lt.jp) adv =
     .       copdy3(j) *akh(j+1) *(dt2dy2(j+1)+t2ody2(j+1)-t2ody2(j))
     .     * crfdy0 /2.
      if (j.gt.1 ) adv = adv
     .     - comdy3(j) *akh(j)   *(dt2dy2(j  )+t2ody2(j)-t2ody2(j-1))
     .     * crfdy0 /2.
      omgql3(j)=(q2dy0*p4p00*crfdy0/rath1) * (ha-adv)
      w     (j)=-1.e3*omgql3(j)*rath1*(t2dy2(j)+t2ody2(j))/(p4p00*gdy0)
   70 continue
c
c --- bilan energetique de l'atmosphere
c
c --- energie potentielle disponible
c
      call av(pstdy2,pstav,cp)
      do 81 j=1,jp
      pstd(j)=pstdy2(j)-pstav
      zint(j)=(p4p00/2.)*(q2dy0/gdy0)*(pstd(j)**2)
   81 continue
      call av(zint,azdy6,cp)
c
c --- energie cinetique du flux moyenne verticalement
c
      do 82 j=1,jp
      zint(j)=(p4p00/(8.*gdy0))*((usdy7(j)+usdy7(j+1))**2)
   82 continue
      call av(zint,akzs,cp)
c
c --- energie cinetique du cisaillement vertical du flux
c
      do 83 j=1,jp
      utmb2(j)=((utdy7(j)+utdy7(j+1))/2.)**2
      zint (j)=(p4p00/(2.*gdy0))*utmb2(j)
   83 continue
      call av(zint,akzt,cp)
      akzo   = akzdy6
      akzdy6 = akzs + akzt
c
c --- conversion de l'energie potentielle disponible en energie cinetiqu
c
      do 84 j=1,jp
      zint(j)=-2.*(f0dy0/gdy0)*pstdy2(j)*omgql3(j)
   84 continue
      call av(zint,cazkz,cp)
c
c --- conversion de l'energie potentielle disponible moyennee zonalement
c                  en energie potentielle disponible des tourbillons
c     (wiin-nielsen et fuenzalida,1975,tellus,p205)
c
      do 85 j=1,jp
      zint(j)=q2dy0*(p4p00/gdy0)*utmb2(j)*(akh(j)+akh(j+1))/2.
   85 continue
      call av(zint,cazae,cp)
c
c --- generation de l'energie potentielle disponible
c     (wiin-nielsen et fuenzalida,1975,tellus,p205)
c
      call av(h2,h2av,cp)
      do 86 j=1,jp
      h2d (j)=h2(j)-h2av
      zint(j)=4.*crfdy0*q2dy0*p4p00*h2d(j)*pstd(j)/(gdy0*cpath1)
   86 continue
      call av(zint,gaz,cp)
c
c --- conversion de l'energie cinetique des tourbillons
c                  en energie cinetique moyennee zonalement
c
      do 87 j=1,jp
      akdq (j)=0.
      akdqu(j)=0.
      do 88 k=1,2
      avla = akdy5(j  ,k)*dqyql3(j  ,k)
      avlb = akdy5(j+1,k)*dqyql3(j+1,k)
      avlo =(avla *c10gr1(j) + avlb *c10gr1(j+1))/2.
      akdq (j) = akdq (j) + avlo
      akdqu(j) = akdqu(j) + avla*udy7(j,k) + avlb*udy7(j+1,k)
   88 continue
      akdqu(j) = akdqu(j) / 2.
      zint(j) =-(p4p00/(2.*gdy0))*akdqu(j)/ca
   87 continue
      call av(zint,ckekz,cp)
      ckekz=ckekz+cazae
c
c --- dissipation de l'energie cinetique
c
      dkz = cazkz + ckekz - (akzdy6-akzo)/delt
c
c --- transport du moment angulaire total
c
      angmot(jpp)=0.
      do 89 jj=2,jpp
      j=jpp+1-jj
      angmot(j)=angmot(j+1) - akdq(j)*pi*(ca**2)*(p4p00/gdy0)*1.e-18
     .         *(sin(cp*j)-sin(cp*(j-1)))
c  -- 1.e-18 : conversion (kg.m2.s-2) ==> (1.e15 ton.m2.s-2)
   89 continue
c ----------------------------------------------------------------------
c
c --- impression des resultats
c
c
      if (iprint.eq.1 ) go to 621
      if (mts   .ne.0 ) go to 620
      if (monpr .gt.mm) go to 620
  621 continue
      write(iwr1,22)it,ian,mm,iday,azdy6,gaz,akzdy6,akzs,akzt,dkz,
     . cazkz,cazae,ckekz
   22 format(//,'   -- caracteristiques dynamiques --',
     .3x,'it =',i6,3x,i3,'e annee',3x,i2,'e mois',3x,i4,'e jour',
     .//,14x,'az   = ',e12.5,3x,'gaz  = ',e12.5,3x,
     . /,14x,'akz  = ',e12.5,3x,'akzs = ',e12.5,3x,'akzt = ',e12.5,
     .    3x,'dkz  = ',e12.5,
     . /,14x,'cazkz= ',e12.5,3x,'cazae= ',e12.5,3x,'ckekz= ',e12.5)
c
  620 continue
c
      if (mts.eq.0) then
c
       write(9,91)(t2dy2(j),j=1,jp)
   91  format(9(f8.3,2x),3x,'t2')
       write(9,92)(udy7(j,1),j=1,jp)
   92  format(9(f8.3,2x),3x,'u1')
       write(9,93)(udy7(j,2),j=1,jp)
   93  format(9(f8.3,2x),3x,'u3')
       write(9,94)(u4dy8(j),j=1,jp)
   94  format(9(f8.3,2x),3x,'u4')
       write(9,95)(w(j),j=1,jp)
   95  format(9(f8.3,2x),'w')
c
       write(9,96)((akdy5(j,k),j=1,jp),k=1,nhd)
   96  format(9(e12.5,1x),'ak')
       write(9,97)((dqyql3(j,k),j=1,jp),k=1,nhd)
   97  format(9(e12.5,1x),'dqy')
c
       write(9,98)(1.e-15*atdy6(j),j=1,jp)
   98  format(9(f8.3,2x),'atm')
c
       write(9,99)(angmot(j),j=1,jp)
   99  format(9(e12.5,1x),'angmot')
c
       write(9,190)(shady2(j),j=1,jp)
  190  format(9(e12.5,1x),'shady2')
c
       do 191 j=1,jp
       htw(j) = htady2(4,j)
  191  continue
       write(9,192)(htw(j),j=1,jp)
  192  format(9(f8.3,2x),'htw')
c
c  -- conversion de j.m-2 en 1.e3 kj.m-2
       azdy6 =azdy6 *1.e-6
       akzdy6=akzdy6*1.e-6
c
       write(9,193)azdy6,akzdy6,cazkz,cazae,gaz,ckekz,dkz,it,mm
  193  format(
     .  7(e12.5,1x),/,1x,'az',11x,'akz',10x,'cazkz',8x,'cazae',8x,
     .  'gaz',10x,'ckekz',8x,'dkz',11x,'it =',i5,3x,'mm =',i3)
      end if
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  wheat===
c cray .endif
      subroutine wheat(aka,theat,rakatp,jmin,jminp)
c
c  -- contrainte integrale sur la vorticite, 
c      respectee en modifiant kq dans les regions tropicales
c
      parameter(np=18,npp=19,nilw=7,nhd=2)
c _dp implicit double precision (a-h,o-z)
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
      common/ql3/dqyql3(npp,nhd),omgql3(np)
      dimension sumi1 (npp),sumi3 (npp),sumi(np)
      dimension theats(nhd),theatn(nhd)
      dimension aka(npp,nhd)
c
      do 1  k=1,nhd
      do 10 j=1,jminp
      sumi1(j)=c10gr1(j)*dqyql3(j,k)*aka(j,k)
      sumi3(j)=0.
 10   continue
      do 11 j=jminp+1,jpp
      sumi1(j)=0.
      sumi3(j)=c10gr1(j)*dqyql3(j,k)*aka(j,k)
 11   continue
      do 12 j=1,jp
      sumi (j)=sumi1(j)+sumi1(j+1)
 12   continue
      call av(sumi,theats(k),cp)
      do 13 j=1,jp
      sumi (j)=sumi3(j)+sumi3(j+1)
 13   continue
      call av(sumi,theatn(k),cp)
 1    continue
c
      if (theatn(1).le.0.or.theatn(2).ge.0.) then
       laka0     = 0
      else
       rakatp   = -.8 * theatn(2) / theatn(1)
       do 15 j=jminp+1,jpp
       aka(j,1) = rakatp   * aka(j,1)
 15    continue
       aakatr   = -  (rakatp *theatn(1) + theatn(2))
     .             / (        theats(1) + theats(2))
c ...  aakatr   : amplitude Kq regions tropicales
c
       if (aakatr.lt.0.) then
        laka0   = 0
       else
        laka0   = 1
       end if
      end if
      if (laka0.eq.0) then
       do 20 k=1,nhd
       do 20 j=1,jpp
       aka(j,k) = 0.
 20    continue
       theat    = 0.
      else
       do 21 k=1,nhd
       do 21 j=1,jminp
       aka(j,k) = aakatr * aka(j,k)
 21    continue
       theat    = -2.*(aakatr*theats(2) + theatn(2))
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  bheat===
c cray .endif
      subroutine bheat(aka,theat,amu1)
c
c  -- contrainte integrale sur la vorticite, 
c      respectee en modifiant kq au niveau 250mb
c
      parameter(np=18,npp=19,nilw=7,nhd=2)
c _dp implicit double precision (a-h,o-z)
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
      common/ql3/dqyql3(npp,nhd),omgql3(np)
      dimension sumi1 (npp),sumi(np)
      dimension aka(npp,nhd)
c
c **********************************************************************
c
c  -- flux moyen hemispherique de chaleur
c
      do 120 j=1,jpp
      sumi1(j)= c10gr1(j)*dqyql3(j,2)*aka(j,2)
 120  continue
      do 121 j=1,jp
      sumi(j) = sumi1(j)+sumi1(j+1)
 121  continue
      call av(sumi,si3,cp)
      theat   =-si3
c
c --- pas de transfert meridien atmospherique
c
      if (theat.lt.0.) then
       do 122 j=1,jpp
       aka(j,2)=0.
 122   continue
       theat   =0.
       amu1    =0.
      else
c
c --- application de la contrainte sur kq a 250mb
c
       do 124 j=1,jpp
       sumi1(j)= c10gr1(j)*aka(j,1)*dqyql3(j,1)
 124   continue
c
       do 125 j=1,jp
       sumi(j) = sumi1(j)+sumi1(j+1)
 125   continue
       call av(sumi,si1,cp)
       amu1    = -si3 / si1
      end if
c
      do 130 j=1,jpp
      aka(j,1)= amu1 * aka(j,1)
 130  continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ttzz====
c cray .endif
      subroutine ttzz(pp,t2,z2,tt,zz)
c
c +++ bcma7.for +++
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
      common/th1/cpath1,rath1,akth1
c
c --- constantes
c
       data p2/5.e4/,p4/1.e5/
c
       p2k = p2 ** akth1
       p21k= p2 **(1.-akth1)
c
c --- calcul
c
       s0  = sigdy0 * p2
       sr1k= s0 /((1.-akth1) *rath1)
c
c  -- altitude 500 mb
       z2 =-(rath1/gdy0) *
     .  (sr1k*(p4-p2) -(t2/p2 +sr1k) *p21k *(p4**akth1-p2k) /akth1)
c
c  -- temperature niveau courant
       pp2k  = (pp/p2)**akth1
       tt =  t2 * pp2k - sr1k * ( pp - p2 * pp2k)
c
c  -- altitude niveau courant
       zz = (rath1/gdy0) *
     .  (sr1k*(pp-p2) -(t2/p2 +sr1k) *p21k *(pp**akth1-p2k) /akth1)
     .  + z2
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  prtdyn==
c cray .endif
      subroutine prtdyn
c
c +++ bcma8.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables dynamiques
c
       common/dy2/t2hdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/dy4/qdy4(np,nhd),dqdy4(np,nhd)
       common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
       common/dy7/zsdy7(np) ,ztdy7(np),
     .            usdy7(npp),utdy7(npp),udy7(npp,nhd)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
c
c --- variables bilan
c
       common/ham/t4ham,t2ham,albham,fupham,fdnham,shham,alhham
c
c ----------------------------------------------------------------------
c
      if (iprint.eq.1) go to 101
      if (mts.ne.0   ) go to 100
      if (monpr.gt.mm) go to 100
  101 continue
      write(iwr1,97)it,ian,mm,iday
   97 format(//'   -- prtdyn --',
     . 3x,'it =',i6,3x,i3,'e annee',3x,i2,'e mois',3x,i4,'e jour')
      write(iwr1,30)
   30 format(/,6x,'j',11x,'t2',10x,'ht',10x,
     . 'q1',10x,'q2',10x,'dq1',9x,'dq2',9x,'z*',10x,'zt',9x,'pst')
      write(iwr1,31)(j,t2dy2(j),htdy2(j),(qdy4(j,k),k=1,nhd),
     . (dqdy4(j,k),k=1,nhd),zsdy7(j),ztdy7(j),pstdy2(j),j=1,jp)
   31 format((i7,4x,9(e12.4)))
c
      write(iwr1,32)
   32 format(/3x,'j',9x,'ak1',9x,'ak2',10x,'u*',10x,'ut',10x,
     . 'u1',10x,'u3',10x,'u4',10x,'us',5x,'us btro',5x,'us bcli')
      write(iwr1,33)(j,(akdy5(j,k),k=1,nhd),usdy7(j),utdy7(j),
     .              (udy7(j,k),k=1,nhd),u4dy8(j),
     .               tsdy8(j),tsody8(j),tsidy8(j),j=1,jpp)
   33 format((i4,10(e12.4)))
      write(iwr1,109)t2hdy2
  109 format (/,' temperature moyenne a 500 mb = ',e12.5)
  100 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  shfdyn==
c cray .endif
      subroutine shfdyn
c
c +++ bcma9.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
c
      dimension jminc(nhd),jminv(nhd)
      dimension aint(npp),bint(npp)
c
      jmin=1+int(30./gs)
c
c  -- largeur de la zone kq petit
c
      jminh=7
c
      write(iwr1,101)jminh
  101 format(//'   -- shfdyn --  largeur kq petit : ',i5)
c
      do 10 k=1,nhd
      jminc(k)=jmin + np*(k - 2)/9
      jminv(k)=jminh+ np*(k - 2)/9
c
      do 18 j=1,jpp
      akody5(j,k)=0.
   18 continue
c
      do 19 j=jminc(k),jpp
      jj      =j-jminc(k)+1
      aint(jj)=akdy5(j,k)
   19 continue
c
      do 20 j=1,jminv(k)
      akody5(j,k)=akdy5(1,k)
   20 continue
c
      ja=jpp-jminc(k)
      jb=jpp-jminv(k)
      call intp(aint,bint,ja,jb)
c
      do 30 j=jminv(k),jp
      jj      =j-jminv(k)+1
      akody5(j,k)=bint(jj)
   30 continue
c
      write(iwr1,60)k,(akody5(j,k),j=1,jpp)
   60 format(/' ak',i1,' :',3(/7(1x,e10.3)))
c
   10 continue
c
      return
      end
c ######################################################################
c cray .if [ $concord .eq. 1 ]
c cray nzyx  intp====
c cray .endif
      subroutine intp(a,b,naa,nbb)
      parameter(npp=19)
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      dimension a(npp),b(npp)
c
      na  = naa-1
      nb  = nbb-1
      xt  = 1.
      dxa = xt /     real(na)
      dxb = xt /     real(nb)
      b(  1)=a(  1)
c
      if (nb.ge.na) then
       ja = 0
       jb = 0
c
c do until
 110   continue
       ja = ja + 1
       xa = dxa*real(ja-1)
       xa1= dxa*real(ja)
       da = a(ja+1)-a(ja)
c
c   do until
 1100  continue
       jb = jb + 1
       xb = dxb*real(jb-1)
       dx = xb-xa
       b(jb) = a(ja) + dx * da / dxa
       if (.not.(xb.ge.xa1.or.jb.ge.nbb-1) ) go to 1100
c
       if (.not.(ja.ge.naa-1)               ) go to 110
c
c   do while
 1110   continue
        if (jb.ge.nbb) go to 1111
        jb = jb + 1
        xb = dxb*real(jb-1)
        dx = xb-xa
        b(jb) = a(ja) + dx * da / dxa
        go to 1110
 1111   continue
c
      else
       ja = 1
       jb = 1
c
c
c do until
  100  continue
       jb = jb + 1
       xb = dxb*real(jb-1)
c
c  do  until
 1000  continue
       ja = ja + 1
       xa = dxa*real(ja-1)
       xa1= dxa*real(ja)
       if (.not.((xb.ge.xa.and.xb.lt.xa1).or.(ja.ge.naa-1))) go to 1000
c
       if (jb.ge.nbb) go to 1001
       da = a(ja+1)-a(ja)
       dx = xb-xa
       b(jb) = a(ja) + dx * da / dxa
 1001  continue
c
       if (.not.(jb.ge.nbb-1)) go to 100
c
       b(nbb)=a(naa)
c
      end if
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  av======
c cray .endif
      subroutine av(x,xav,cp)
      parameter (np=18)
c _dp implicit double precision (a-h,o-z)
c
c +++ bcmav.for +++
c
      double precision  y,ycos,yav,dcp
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      dimension  x(np),y(np)
c
          dcp  = cp
          yav  = 0.d00
      do 5 j = 1,np
          y(j) = x(j)
          ycos = y(j)*(dsin(dcp*dble(j))-dsin(dcp*dble(j-1)))
          yav  = yav+ycos
    5 continue
          xav  = yav
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  enbil===
c cray .endif
      subroutine enbil
c
c +++ bcmb0.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical       lgm18k,taiga,otter,summer,co2var
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      double precision  hssn2,hspsn2
c
       character*1 ch(nilw)
       character*4 label
       character*8 fnam
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/varl02/lgm18k,taiga,otter,summer,co2var
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/pal/apal,fpal,ipal,ipaleo
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
      common/gr2/argr2(np),atgr2
c
c --- constantes physiques
c
       common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables radiatives
c
       common/fac/facsc,facco2,faccl
       common/gz2/ppmgz2
c
c --- variables dynamiques
c
       common/p01/szp01(np,nilw),spp01(np,nilw)
       common/p02/szp02(np,nilw)
       common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/dy4/qdy4(np,nhd),dqdy4(np,nhd)
       common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
       common/dy6/atdy6(npp),atmdy6(npp),flady6(npp),azdy6,akzdy6
       dimension             atm   (npp)
       common/dy7/zsdy7(np) ,ztdy7(np),
     .            usdy7(npp),utdy7(npp),udy7(npp,nhd)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
       common/dy9/u1jdy9,jjdy9,itjdy9
c
c --- variables surface
c
       common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .        hiic2(np),hipic2(np),awic2(np),awpic2(np)
       common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .                 hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
       common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
c
       common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .            vopsn3(np,nilw), vosn3(np,nilw),
     .            hspsn3(np,nilw), hssn3(np,nilw),
     .            ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
       dimension  dsnnet(np,nilw),dsnajs(np,nilw),dsnvol(np,nilw)
       dimension  densn (   nilw)
       common/sn4/accsn4(np,4:nilw,5),ablsn4(np,4:nilw,5)
       dimension  accwri(np,4:nilw)  ,ablwri(np,4:nilw)
       common/ca1/ajca1(np,nilw)
c
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .           ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
       common/su2/hocsu2(np),hicsu2(np),awsu2(np),haisu2(np)
       dimension ajn(np),ajmav(nilw),t4n(np),t4mav(nilw)
c
c --- variables ocean
c
       common/oc6/otoc6(npp),otmoc6(npp),flooc6(npp)
c
c --- variables continent
c
       common/la1/tcmla1(np),tcdla1(np),csmla1(np),csdla1(np)
       common/la2/tcmla2(np),tcdla2(np)
c
c --- variables bilan general
c
       common/bi0/surbi0(np,nilw),atmbi0(np),topbi0(np)
       common/bi1/surbi1(nilw)   ,atmbi1    ,topbi1    ,radbi1(np)
       common/ham/t4ham,t2ham,albham,fupham,fdnham,shham,alhham
       common/dt4/dt4ep,dt4t2
       common/hv2/htshv2(6,np,nilw),zhshv2(6,np),shshv2(np),shahv2(np),
     .  albhv2(np),sahv2(np),fuphv2(np),fdnhv2(np),hslhv2(np),t2hv2(np)
       common/hv3/ialhv3(np)
       dimension aux(np)
c
       data ch /'i','w','l','*','g','l','f'/
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      if (it.ne.0) then
       read(iwr2,2)label
    2  format(a4)
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
       write(iwr1,10)
   10  format (////1h ,120('*'),/1h ,10('bilans',6x),/1h ,120('*'))
c
c ----------------------------------------------------------------------
c
c   -- determination du diviseur de moyenne
       if (season) then
                           bb=tsm*12.
        if(ian.eq.1)       bb=bb-1.

       else
                           bb=tsm/2.
        if (it.gt.int(bb)) bb=tsm
       end if
c
c --- nouvelles altitudes
c
c     do 100 j=1,jp
c     do 100 iilw=4,nilw
c     szp01(j,iilw) = szp02(j,iilw) + .367 *hssn2(j,iilw) 
c     call ttpp(spp01(j,iilw),t2dy2(j),z2a,tt,szp01(j,iilw))
c100  continue
c
c --- bilan au sommet de l'atmosphere
c
       topbi1 = topbi1 / bb
c
c --- bilan de l'energie emmagasinnee par les differents secteurs (w.m-2
c
       do 132 j=1,jp
       do 133 iilw=1,nilw
       do 134 n=1,6
       htshv2(n,j,iilw) = htshv2(n,j,iilw) / bb
 134   continue
 133   continue
       do 135 n=1,6
       zhshv2(n,j     ) = zhshv2(n,j     ) / bb
 135   continue
       shshv2(  j     ) = shshv2(  j     ) / bb
       shahv2(  j     ) = shahv2(  j     ) / bb
       if (ialhv3(j).gt.0)
     . albhv2(  j     ) = albhv2(  j     ) / ialhv3(j)
        sahv2(  j     ) =  sahv2(  j     ) / bb
       fuphv2(  j     ) = fuphv2(  j     ) / bb
       fdnhv2(  j     ) = fdnhv2(  j     ) / bb
       hslhv2(  j     ) = hslhv2(  j     ) / bb
        t2hv2(  j     ) =  t2hv2(  j     ) / bb
  132  continue
       call av(shahv2,shaham,cp)
       call av( sahv2, saham,cp)
       call av(hslhv2,hslham,cp)
c
       surbih       = 0.
       do 200 iilw=1,nilw
       do 201 j=1,jp
       aux(j)=surbi0(j,iilw)
  201  continue
       call av(aux,surbi1(iilw),cp)
       surbi1(iilw) = surbi1(iilw) / bb
       surbih       = surbih       + surbi1(iilw)
  200  continue
c
c --- moyenne annuelle dans le systeme ocean - glace marine
c
       do 205 j=1,jp
       hocsu2(j) = hocsu2(j) / bb
       if (awsu2(j).ne.0.) hicsu2(j) = hicsu2(j) / awsu2(j)
        awsu2(j) =  awsu2(j) / bb
       haisu2(j) = haisu2(j) * ajosu0(j) * argr2(j) / bb
  205  continue
c
c --- moyenne annuelle des transports de chaleur atmospherique et oceani
c
       atm (jpp)= 0.
       do 210 ji=1,jp
       j = jpp - ji
       otmoc6(j)= otmoc6(j)/bb
       radbi1(j)= radbi1(j)/bb
       atm   (j)= atm   (j+1) - radbi1(j) *argr2(j)
       atmdy6(j)= atmdy6(j)/bb
c
       u4mdy8(j)= u4mdy8(j)/bb
c
  210  continue
       do 211 j=1,jp
         atm   (j)= atm   (j)   - otmoc6(j)
  211  continue
c
c --- volume global de neige
c
       do 233 iilw=4,nilw
        densn(  iilw) =  0.
       do 232 j=1,jp
        vosn3(j,iilw) =  vosn3(j,iilw) * argr2(j) / bb
        ajsn3(j,iilw) =  ajsn3(j,iilw) * argr2(j) / bb
       dsnnet(j,iilw) = accsn3(j,iilw) - ablsn3(j,iilw)
       dsnajs(j,iilw) =  ajsn3(j,iilw) - ajpsn3(j,iilw)
       dsnvol(j,iilw) =  vosn3(j,iilw) - vopsn3(j,iilw)
c
        densn(  iilw) =  densn(  iilw) + dsnvol(j,iilw)
c
       vopsn3(j,iilw) =  vosn3(j,iilw)
       ajpsn3(j,iilw) =  ajsn3(j,iilw)
        do 234 i=1,4
        accsn4(j,iilw,i)= accsn4(j,iilw,i+1)
        ablsn4(j,iilw,i)= ablsn4(j,iilw,i+1)
 234    continue
        accsn4(j,iilw,5)= accsn3(j,iilw)
        ablsn4(j,iilw,5)= ablsn3(j,iilw)
        accwri(j,iilw)  = 0.
        ablwri(j,iilw)  = 0.
        do 235 i=1,5
        accwri(j,iilw)   = accwri(j,iilw) + accsn4(j,iilw,i)
        ablwri(j,iilw)   = ablwri(j,iilw) + ablsn4(j,iilw,i)
 235    continue
                      div = ian
        if (ian.gt.5) div = 5.
        accwri(j,iilw)   = accwri(j,iilw) / div
        ablwri(j,iilw)   = ablwri(j,iilw) / div
 232   continue
        densn(   iilw)   = chfth4 * densn(  iilw) / atgr2
 233   continue
c
c --- moyennes annuelles des temperatures de surface continent sans calo
c
       do 222 j=1,jp
          tmsu1 (j)= tmsu1(j)  / bb
          if (csmla1(j).ne.0.) tcdla1(j)= tcmla1(j) / csmla1(j)
          if (summer) then
                               tcdla2(j)=tcmla2(j) 
c                              tcdla2(j)=tcmla2(j)  / tsm
                      else
                               tcdla2(j)=tcmla2(j)  / bb
                      end if
       do 223 iilw=1,nilw
c
c
       if (bjmsu1(j,iilw).ne.0.) then
          t4msu1(j,iilw)=t4msu1(j,iilw) / bjmsu1(j,iilw)
          t4dsu1(j,iilw)=t4msu1(j,iilw)
       else
          t4msu1(j,iilw)=0.
          t4dsu1(j,iilw)=tmsu1(j)
       endif
c
  223  continue
  222  continue
c
       if ((taiga.or.inpu).and.(.not.otter).and.mod(ian,20).eq.0) then
       is = 4
          if (summer) then
           tmax = 290.
           tran =  12.
          else
           tmax = 274.
           tran =  14.
          end if
          call taifdb(tmax,tran,is)
          call av(tcdla1,tcdham,cp)
        end if
c
c ------   pour exp. co2   thierry   -----  debut  ---------------------
          call av(haisu2,voliav,cp)
c ------   pour exp. co2   thierry   -----  fin    ---------------------
       do 225 iilw=1,nilw
       do 226 j=1,jp
          ajn(j)=ajmsu1(j,iilw)
          t4n(j)=t4msu1(j,iilw) * ajmsu1(j,iilw)
  226  continue
          call av(ajn,ajmav(iilw),cp)
          call av(t4n,t4mav(iilw),cp)
c
       if (ajmav(iilw).gt.0.) then
        t4mav(iilw) = t4mav(iilw) / ajmav(iilw)
       else
        t4mav(iilw) = 0.
       end if
c
  225  continue
c
c --- calcul des moyennes annuelles globales
c
          t4ham  = t4ham  / bb
          t2ham  = t2ham  / bb
          dt4ep  = dt4ep  / bb
          dt4t2  = t4ham - t2ham
c
       call av(atmbi0,atmbi1,cp)
          atmbi1 = atmbi1 / bb
          albham = albham / bb
          fupham = fupham / bb
          fdnham = fdnham / bb
          shham  = shham  / bb
          alhham = alhham / bb
c
c --- ecriture des resultats
c
       if (mod(ian-1,20).eq.0.and.season) then
        write(iwr6,117)
 117    format(4x,'moyennes hemispheriques (w.m-2)',27x,
     .   '  accum./ablat.*(m.an-1)      ts g  t4 Hem. alb.pla.',
     .   '   h ice   oc.ml.d',/,
     .   '  an  sommet atm  net atm -surf.   -groen.   -land *',
     .   '  -groen.*    * 65-70n groe 60-65n    60-65n',18x,
     .   '  80-85n     5-10n')
       end if
       write(iwr6,116)ian,topbi1,atmbi1,surbih,surbi1(5),
     .  densn(4),densn(5),
     .  0.367*accsn3(14,4),0.367*ablsn3(14,4),
     .  0.367*accsn3(13,5),0.367*ablsn3(13,5),
     .  t4ssu0(13,5),t4ham,albham,hiic2(17),hoc2(2)
c ...  facteur 0.367 : rapport de densite neige tassee / glace
 116   format(i4,e12.4,f9.5,f7.3,3f10.6,
     .           4f6.3,f10.5,f9.4,f9.5,f8.5,f10.5)
c
       write(iwr1,60)ian,facsc,facco2,faccl
   60  format(/i4,'e annee : ',/,
     . '  (',f4.2,' x sol.const.)',/,
     . '  (',f4.2,' x co2       )',/,
     . '  (',f4.2,' x cloudiness)')
       write(iwr1,61)topbi1
   61  format(/'  bilans relatifs a l`atmosphere : ',/,2x,32('*'),
     .  //,3x,'sommet de l`atmosphere :',f12.4,' (w.m-2)',4x,
     .        'dans et a la base de l`atmosphere')
       write(iwr1,62)
   62  format(4x,'latitude',6x,'albedo',3x,'f.sol.dwn',4x,'f.i.r.up',
     .        9x,'hta',3x,'f.i.r.dwn',2x,'h. s. & l.',4x,'t 500 mb')
       write(iwr1,63)(deggr1(j),albhv2(j),sahv2(j),fuphv2(j),shahv2(j),
     .             fdnhv2(j),hslhv2(j),t2hv2(j),j=1,jp)
   63  format((f12.1,f12.4,6f12.2))
       write(iwr1,64)albham,saham,fupham,shaham,fdnham,hslham,t2ham
   64  format(' hem.an.av. :',f11.4,6f12.2)
c
       write(iwr1,65)shham,alhham,dt4ep,dt4t2
   65  format(72x,f12.2,' (sensible)',/,72x,f12.2,' (latent)',
     . /,72x,'dtsol eq-po ',f12.2,
     . /,72x,'dt 1000-500 ',f12.2)
c
c
       write(iwr1,66)(surbi1(iilw),iilw=1,nilw),surbih
   66  format(/,2x,'bilans pour chacun des secteurs :',/2x,32('*'),
     . //,3x,'input glace marine  :',f12.4,' (w.m-2)',
     .  /,3x,'input ocean         :',f12.4,' (w.m-2)',
     .  /,3x,'input continent   1 :',f12.4,' (w.m-2)',
     .  /,3x,'input champ neige 1 :',f12.4,' (w.m-2)',
     .  /,3x,'input calotte groen.:',f12.4,' (w.m-2)',
     .  /,3x,'input calotte laur. :',f12.4,' (w.m-2)',
     .  /,3x,'input calotte finsc.:',f12.4,' (w.m-2)',
     .  /,3x,'input hemisphere    :',f12.4,' (w.m-2)')
       write(iwr1,67)
   67  format(/,2x,'flux de chaleur a l`interface atmosphere - secteur')
       write(iwr1,68)(n,n=1,6),(n,n=1,6)
   68  format(/' latitude  si',i6,5i10,'  oc',i6,5i10)
       write(iwr1,69)(deggr1(j),((htshv2(n,j,iilw),n=1,6),iilw=1,2),
     .                j=1,jp)
   69  format((f9.2,12f10.2))
       write(iwr1,70)(n,n=1,6),(n,n=1,6)
   70  format(/' latitude  ln 1',i4,5i10,'  sn 1',i4,5i10)
       write(iwr1,69)(deggr1(j),((htshv2(n,j,iilw),n=1,6),iilw=3,4),
     .                j=1,jp)
       write(iwr1,71)(n,n=1,6),(n,n=1,6)
   71  format(/' latitude  ca g',i4,5i10,'  ca l',i4,5i10)
       write(iwr1,69)(deggr1(j),((htshv2(n,j,iilw),n=1,6),iilw=5,6),
     .                j=1,jp)
       write(iwr1,72)(n,n=1,6)
   72  format(/' latitude  ca f',i4,5i10)
       write(iwr1,76)(deggr1(j),(htshv2(n,j,7),n=1,6),j=1,jp)
   76  format((f9.2,6f10.2))
c
       write(iwr1,73)
   73  format(/2x,'temperatures moyennes annuelles',/,2x,31('*'),
     . //,8x,'tice',8x,'tocn',4x,'t land 1',4x,'t snow 1',5x,'t cal g',
     .    5x,'t cal l',5x,'t cal f',4x,'t cont 1',8x,'tzon')
c
       write(iwr1,74)((t4msu1(j,iilw),iilw=1,nilw),tcdla2(j),tmsu1(j),
     .                     j=1,jp)
   74  format((9f12.2))
       write(iwr1,75)(t4mav(iilw),iilw=1,nilw),tcdham,t4ham
   75  format(  '               moyennes annuelles hemispheriques',
     .      /, 9f12.2)
c
c --- output accumulations
c
c mVax open(unit=13,status='new',file='bcmacc.d')
       fnam(1:4) ='bacc'
       fnam(5:8) = label
       open(unit=13,status='unknown',file=fnam)
       rewind 13
c
       if (season) then
        write(iwr1,92)ian
   92   format(/2x,'caracteristiques modele glace marine-ocean',
     .  12x,'modele de neige',33x,'annee',i5,/2x,42('*'),12x,15('*'),
     .  /2x,'lat.',9x,'hoc',10x,'aw',3x,'vol.glace',4x,
     .  '  acc(m.y-1)  abl(m.y-1)  net(m.y-1)     aj* min',
     .  '   d j* (m2)   d v* (m3)   vol* (m3)')
        write(iwr1,93)(deggr1(j),hocsu2(j),awsu2(j),haisu2(j),
     .   accsn3(j,4),ablsn3(j,4),dsnnet(j,4),assn3(j,4),
     .   dsnajs(j,4),dsnvol(j,4), vosn3(j,4),j=1,jp)
   93   format((f6.1,f12.2,f12.4,e12.4,4x,7e12.4))
        write(13,94)apal,ch(4),(accwri(j,   4),j=1,jp),
     .                   ((-1.)*ablwri(j,   4),j=1,jp)
   94   format(f8.3,1x,a1,9e12.4,/,7x,9e12.4,
     .        ' acc',2(/,7x,9e12.4),' abl')
c
        do 96 iilw=5,nilw
        write(iwr1,97) iilw,ian
   97   format(/56x,'caracteristiques calotte no',i2,31x,'annee',i5,
     .         /56x,29('*'),/2x,'lat.',36x,
     .  '  acc(m.y-1)  abl(m.y-1)  net(m.y-1)     aj* min',
     .  '               d v* (m3)   vol* (m3)')
        do 99 j=1,jp
c
        if (ajsu0(j,iilw).gt.0.) then
         write(iwr1,98)deggr1(j),accsn3(j,iilw),ablsn3(j,iilw),
     .    dsnnet(j,iilw),assn3(j,iilw),dsnvol(j,iilw),vosn3(j,iilw)
   98    format((f6.1,36x,4e12.4,12x,2e12.4))
        else
         accwri(j,iilw) = accwri(j,4)
         ablwri(j,iilw) = ablwri(j,4)
        end if
c
   99   continue
        write(13,94)apal,ch(iilw),(accwri(j,iilw),j=1,jp),
     .                      ((-1.)*ablwri(j,iilw),j=1,jp)
   96   continue
c
        write(13,940)(otmoc6(j),j=1,jp)
 940    format('    to ',9e12.4,/,7x,9e12.4,' to')
        write(13,941)(atm   (j),j=1,jp)
 941    format('    ta ',9e12.4,/,7x,9e12.4,' ta')
        write(13,942)(albhv2(j),j=1,jp),albham,ppmgz2,topbi1
 942    format('   alb ',9e12.4,/,7x,10e12.4,/,
     .         '   CO2 ', f12.2,81x,'bil.top atm. : ',e12.4)
c
        close(unit=13)
c
       end if
c
       write(iwr1,80)(otmoc6(j),j=1,jpp)
   80  format(/14x,'o.f.',10x,': ',/,(10(1x,e11.4)))
       write(iwr1,81)(atmdy6(j),j=1,jpp)
   81  format(/14x,'a.f.(eddies)',2x,': ',/,(10(1x,e11.4)))
       write(iwr1,82)(atm   (j),j=1,jpp)
   82  format(/14x,'a.f.',10x,': ',/,(10(1x,e11.4)))
c
       write(iwr1,83)u1jdy9,jjdy9,itjdy9,(u4mdy8(j),j=1,jpp)
   83  format(/14x,'jet stream maximum : ',e12.5,4x,'j = ',i2,4x,
     .  'jour : ',i7,//14x,'u1000mb',7x,': ',/,(10(1x,e11.4)))
c
c ----------  pour exp. co2 de thierry   ------  debut   ---------------
c co2      write(iwr7,50) apal,ppmgz2,ajmav(1),voliav
 50        format (3x,'moyennes annuelles',5x,'annee =',f6.2,5x,
     .             'co2 =',f8.3,/,15x,
     .             'ice surf. =',e12.4,3x,'ice vol. =',e12.4)
c co2      write(iwr7,51) (t4msu1(j,2),j=1,jp)
c co2      write(iwr7,52) (t4msu1(j,4),j=1,jp)
c co2      write(iwr7,53) (tmsu1(j),j=1,jp)
c co2      write(iwr7,54) (t4msu1(j,3),j=1,jp)
c co2      write(iwr7,55) (ablsn3(j,5),j=1,jp)
c co2      write(iwr7,56) (accsn3(j,5),j=1,jp)
c co2      write(iwr7,57) (hssn2(j,5),j=1,jp)
 51        format (1x,'t4msu1(j,2) =',8f8.3,/,10f8.3)
 52        format (1x,'t4msu1(j,4) =',8f8.3,/,10f8.3)
 53        format (1x,' tmsu1(j)   =',8f8.3,/,10f8.3)
 54        format (1x,'t4msu1(j,3) =',8f8.3,/,10f8.3)
 55        format (1x,'ablsn3(j,5) =',5e11.3,/,7e11.3,/,6e11.3)
 56        format (1x,'accsn3(j,5) =',5e11.3,/,7e11.3,/,6e11.3)
 57        format (1x,' hssn2(j,5) =',5e11.3,/,7e11.3,/,6e11.3)
c ----------  pour exp. co2 de thierry   ------  fin     ---------------
c
c **********************************************************************
c
c --- fermeture des fichiers d'output
c
       if (it.lt.nts) then
        close(unit=iwr1)
        close(unit=3)
        close(unit=4)
        close(unit=iwr5)
        close(unit=iwr7)
        close(unit=8)
        close(unit=9)
       end if
c
c **********************************************************************
c
c --- sauvetage de la simulation
c
c mVax open(unit=20,status='new',file='bcmsim.d')
       fnam(1:4) ='bsim'
       fnam(5:8) = label
       open(unit=20,status='unknown',file=fnam)
       rewind 20
       write(20,1)
     . t2dy2,t4su0,t4ssu0,t4rsu0,tmsu1,t4dsu1,tcdla1,tcdla2,ajsu0,
     . ajasu1,twoc2,sstoc2,hoc2,heqoc2,ajosu0,
     . twic2,twiic2,tsiic2,awpic2,awic2,hipic2,hiic2,
     . ajca1,hspsn2,hssn2,
     . szp01,pstdy2,qdy4,akody5,akdy5,zsdy7,ztdy7 ,ajasu1
    1  format(/,3(3(1x,e12.5),/))
       close(unit=20)
c
c **********************************************************************
c
c --- ouverture de nouveaux fichiers d'output
c
       if (it.lt.nts) then
c mVax open(unit=iwr1,status='new',file='bcmout.d')
       fnam(1:4) ='bout'
       fnam(5:8) = label
       open(unit=iwr1,status='unknown',file=fnam)
       rewind iwr1
c mVax open(unit=3,status='new',file='bcmio.d')
       fnam(1:4) ='bcmo'
       fnam(5:8) = label
       open(unit=3,status='unknown',file=fnam)
       rewind 3
c mVax open(unit=4,status='new',file='bcmls.d')
       fnam(1:4) ='bcml'
       fnam(5:8) = label
       open(unit=4,status='unknown',file=fnam)
       rewind 4
c mVax open(unit=iwr5,status='new',file='bcmc.d')
       fnam(1:4) ='bcmc'
       fnam(5:8) = label
       open(unit=iwr5,status='unknown',file=fnam)
       rewind iwr5
c mVax open(unit=iwr7,status='new',recl=150,file='bcmr.d')
       fnam(1:4) ='bcmr'
       fnam(5:8) = label
       open(unit=iwr7,status='unknown',file=fnam)
       rewind iwr7
c mVax open(unit=8,status='new',file='bcmm.d')
       fnam(1:4) ='bcmm'
       fnam(5:8) = label
       open(unit=8,status='unknown',file=fnam)
       rewind 8
c mVax open(unit=9,status='new',file='bcma.d')
       fnam(1:4) ='bcma'
       fnam(5:8) = label
       open(unit=9,status='unknown',file=fnam)
       rewind 9
       end if
c
c **********************************************************************
c
c --- ecriture des parametres de la simulation
c
c
      write(iwr1,20)nts,season,tsclim,taclim,ocen,ekm,
     .                  turbu,facsc,facco2,simpir,inpu,
     .                  apal,lgm18k,co2var,taiga,otter,summer
   20 format (
     ./' nts   ',7x,'= ',i9,
     ./' season',7x,'= ',l9,
     ./' tsclim',7x,'= ',l9,18x,' taclim',7x,'= ',l9,
     ./' ocen  ',7x,'= ',l9,18x,' ekm   ',7x,'= ',l9,
     ./' turbu ',7x,'= ',l9,
     ./' facsc ',7x,'= ',f13.3,
     ./' facco2',7x,'= ',f13.3,
     ./' simpir',7x,'= ',l9,
     ./' inpu  ',7x,'= ',l9,
     ./' apal  ',7x,'= ',f13.3,
     ./' lgm18k',7x,'= ',l9,
     ./' co2var',7x,'= ',l9,
     ./' taiga ',7x,'= ',l9,
     ./' otter ',7x,'= ',l9,
     ./' summer',7x,'= ',l9)
       delth=delt/3600.
       write(iwr1,21)dels,delth
  21   format (
     . /' grid spacing  = ',e12.5,' (h)',
     . /' time spacing  = ',f12.3,' (m)')
c
       if (season) ian = ian + 1
c
      else
       do 11 iilw=4,nilw
       do 11 i=1,5
       do 11 j=1,jp
       accsn4(j,iilw,i) = 0.
       ablsn4(j,iilw,i) = 0.
 11    continue
      end if
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      call inibil
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  inibi===
c cray .endif
      subroutine inibil
c
c +++ bcmb1.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      double precision  hssn2,hspsn2
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables dynamique
c
       common/dy6/atdy6(npp),atmdy6(npp),flady6(npp),azdy6,akzdy6
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
       common/dy9/u1jdy9,jjdy9,itjdy9
c
c --- variables surface
c
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .           ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
       common/su2/hocsu2(np),hicsu2(np),awsu2(np),haisu2(np)
       common/oc6/otoc6(npp),otmoc6(npp),flooc6(npp)
       common/la1/tcmla1(np),tcdla1(np),csmla1(np),csdla1(np)
       common/la2/tcmla2(np),tcdla2(np)
       common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
       common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .            vopsn3(np,nilw), vosn3(np,nilw),
     .            hspsn3(np,nilw), hssn3(np,nilw),
     .            ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
c
c --- variables bilan general
c
       common/hv2/htshv2(6,np,nilw),zhshv2(6,np),shshv2(np),shahv2(np),
     .  albhv2(np),sahv2(np),fuphv2(np),fdnhv2(np),hslhv2(np),t2hv2(np)
       common/hv3/ialhv3(np)
       common/ham/t4ham,t2ham,albham,fupham,fdnham,shham,alhham
       common/dt4/dt4ep,dt4t2
       common/bi0/surbi0(np,nilw),atmbi0(np),topbi0(np)
       common/bi1/surbi1(nilw)   ,atmbi1    ,topbi1    ,radbi1(np)
c
c **********************************************************************
c
c --- initialisation des variables bilan moyen hemispherique
c
          t4ham          = 0.
          t2ham          = 0.
          dt4ep          = 0.
          albham         = 0.
          fupham         = 0.
          fdnham         = 0.
          shham          = 0.
          alhham         = 0.
c
c --- initialisation des variables flux de chaleur et vents
c
          otmoc6 (jpp)   = 0.
          atmdy6 (jpp)   = 0.
c
          u4mdy8(jpp)    = 0.
          u1jdy9         = 0.
c
          topbi1         = 0.
c
      do 10 j=1,jp
          otmoc6(j)      = 0.
          atmdy6(j)      = 0.
          u4mdy8(j)      = 0.
          atmbi0(j)      = 0.
          radbi1(j)      = 0.
      do 11 iilw=1,nilw
          surbi0(j,iilw) = 0.
          ajmsu1(j,iilw) = 0.
          bjmsu1(j,iilw) = 0.
c
c --- initialisation des variables champ de neige
c
          ablsn3(j,iilw) = 0.
          accsn3(j,iilw) = 0.
           vosn3(j,iilw) = 0.
           ajsn3(j,iilw) = 0.
           assn3(j,iilw) = ajsu0(j,iilw)
c
c --- initialisation des variables temperature
c
          t4msu1(j,iilw) = 0.
   11 continue
c
           tmsu1(j)      = 0.
          tcmla1(j)      = 0.
          csmla1(j)      = 0.
          tcmla2(j)      = 0.
c
c --- initialisation des variables ocean - glace marine
c
          hocsu2(j)      = 0.
          hicsu2(j)      = 0.
           awsu2(j)      = 0.
          haisu2(j)      = 0.
   10 continue
c
c --- initialisation des bilans moyens annuels
c
      do 132 j=1,jp
      do 133 iilw=1,nilw
      do 134 n=1,6
      htshv2(n,j,iilw) = 0.
  134 continue
  133 continue
      do 135 n=1,6
      zhshv2(n,j     ) = 0.
  135 continue
      shshv2(  j     ) = 0.
      shahv2(  j     ) = 0.
      albhv2(  j     ) = 0.
      ialhv3(  j     ) = 0
       sahv2(  j     ) = 0.
      fuphv2(  j     ) = 0.
      fdnhv2(  j     ) = 0.
      hslhv2(  j     ) = 0.
       t2hv2(  j     ) = 0.
  132 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ca1=====
c cray .endif
      subroutine inica1(ic)
c
c +++ bcmc1.for +++
c
      parameter(np=18,npp=19,nilw=7)
c _dp implicit double precision (a-h,o-z)
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/p01/szp01(np,nilw),spp01(np,nilw)
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/ca1/ajca1(np,nilw)
      common/ca2/caca2(np,nilw),daca2(np,nilw)
c
       data dmax  /500.e3/
c ...       dmax = 500 km : rayon de continentalite
            dma2 = dmax*dmax
       data dx   /  10.e3/
c ...       dx   =  10 km : increment pour calcul de continentalite
c
       data h1   /   2.  /
c ...       h1   =   2 km : seuil pour le Desert - Altitude Effect
            alo2 =   1.d0 / dlog(2.d0)
c
c --- facteur de continentalite de la calotte
c
      write(iwr6,736)ic
 736  format(/,' Parametres de la calotte',i2,' (--- inica1 ---)',
     .       /,'  j    jc H.crest D.crest  lambda D.cont.',
     .         '      x1      x2','     x10     x20',
     .         '      Sl',7x,'C     SCH      DA')
c
       do 10 j=1,jp
       if (ajsu0(j,ic).ne.0.) then
        ajg  = ajsu0(j,ic)
        ahg  = szp01(j,ic) * 1.5
c ...   correction : altitude moyenne => altitude max
c
        dpa  = 2.*pi*ca*cos(cp*(j-.5))
        djg  = ajg/2.
        dcrest = dpa * djg
        dcont  = dpa * ajca1(j,ic)
        al1  = ahg*ahg/dcrest
c ...   al1  : parametre lambda de la calotte dans bcm (m)
c
        zz1  = 0.
  731   continue
        zz2  = ahg
c
        if (zz2.gt.0.d0) then
         slcal = al1 / (zz2 + zz1)
         if (ajsu0(j-1,ic).le.0.) slcal = slcal * 4./3.
         if (ajsu0(j+1,ic).le.0.) slcal = slcal * 4./3.
         if (slcal.gt..07d0) slcal = .07d0
        else
         slcal = 0.d0
        end if
c ...   pente moyenne Ouest-Est 
c
        x1             = zz1         * zz1         / al1
        x2             = zz2         * zz2         / al1
        xt             = x2 * 2.
        if (xt.ge.dcont) then
          x10          = x1
        else
         if (xt.ge.dcont-dmax) then
          x10          = amax1(0.,dcont-xt)
          x10          = amin1(x10,x2)
         else
          x10          = x2
         end if
        end if
          x20          = x10   + x2
c
c  do while
         xx             = x1
         xc             = 0.
         cont           = 0.
c  begin dowhile
 7301   continue
         if (xx.lt.dmax) then
          cont          = cont +2.*dx
     .     /(1.+(acos(xx/dmax) -xx*sqrt(dma2-xx*xx)/dma2) /pi)
         else
          cont          = cont +2.*dx
         end if
         xx             = xx + dx
         xc             = xc + dx
        if (xx.gt.x2 ) go to 7300
        go to 7301
 7300   continue
c  end dowhile
c
c
c  do while
         xx             = x10
c  begin dowhile
 7311   continue
         if (xx.lt.dmax) then
          cont          = cont +2.*dx
     .     /(1.+(acos(xx/dmax) -xx*sqrt(dma2-xx*xx)/dma2) /pi)
         else
          cont          = cont +2.*dx
         end if
         xx             = xx + dx
         xc             = xc + dx
        if (xx.gt.x20) go to 7310
        go to 7311
 7310   continue
c  end dowhile
c
         cocal  = cont / xc
c ...    continentalite moyenne sur la tranche de denivellation
c
         if (ajsu0(j,ic).gt.0.0001) then
          caca2(j,ic) =1.5*(1. +14.*slcal -1.7e-4*szp01(j,ic)) /cocal
c ...     parametrisation d'Oerlemans, J.of Climatology, 1982.
c              facteur 1.5 pour convertir moyenne zonale 
c                          des precipitations en moyenne oceanique.
c
         else
          caca2(j,ic) = 1.0
         end if
         h2           = 1.d-3 * ahg
         if (h2.lt.h1) then
          daca2(j,ic) = 1.d0
         else
          daca2(j,ic) =
     .    (h1*h1 +(2*alo2) *(h1+alo2 -(2**(h1-h2)) *(h2+alo2))) /(h2*h2)
         end if
c
        write(iwr6,737)j,ajg,ahg,1.d-3*dcrest,al1,1.e-3*dcont,
     .  1.e-3*x1,1.e-3*x2,1.e-3*x10,1.e-3*x20,
     .  100.*slcal,cocal,caca2(j,ic),daca2(j,ic)
 737    format(i3,f6.3,2f8.1,f8.2,5f8.1,f8.2,3f8.3)
c
       end if
 10   continue
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ca2=====
c cray .endif
      subroutine inica2(ic)
c
c +++ bcmc2.for +++
c
      parameter(np=18,npp=19,nilw=7)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      double precision  hssn2,hspsn2
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables altitude et pression surface
c
       common/p01/szp01(np,nilw),spp01(np,nilw)
c
c --- variables surface
c
       common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
       common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .            vopsn3(np,nilw), vosn3(np,nilw),
     .            hspsn3(np,nilw), hssn3(np,nilw),
     .            ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
       common/sn5/rsusn5(np,nilw),rslsn5(np,nilw),rssn5(np,nilw)
       common/sn6/hfasn6(np,nilw),  nsn6(np,nilw),nnsn6(np,nilw)
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .           ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
       common/la3/t4sla3(np,nilw)
c
      dimension t4a0(9),t4s0(9),t4a(np),t4s(np)
c
c     data t4a0  /300.1,297.7,293.5,288.3,281.8,276.2,270.1,262.4,255.1/
c ... ancienne initialisation (mv30)
c
      data t4a0/26.3,26.3,23.2,15.9, 8.4,  2.2, -5.5,-12.7,-18.0/
      data t4s0/26.0,24.0,18.0, 8.2,-2.4,-10.6,-22.6,-24.8,-33.9/
c ... source : warren & schneider (1979) moyenne annuelle et de janvier
c
      if (.not.inpu) then
       call pol(t4a0,t4a,9,np)
       call pol(t4s0,t4s,9,np)
c
       if (season) then
        do 112 j=1,jp
          t4ssu0(j,ic) = t4s(j) + tfsth4
          t4sla3(j,ic) = t4s(j) + tfsth4
          t4rsu0(j,ic) = t4s(j) + tfsth4
 112    continue
       else
        do 102 j=1,jp
          t4ssu0(j,ic) = t4a(j) + tfsth4
          t4sla3(j,ic) = t4a(j) + tfsth4
          t4rsu0(j,ic) = t4a(j) + tfsth4
 102    continue
       end if
c
       do 12 j=1,jp
          t4dsu1(j,ic) = t4a(j) + tfsth4
 12    continue
      else
       do 103 j=1,jp
       t4s(j) = 0.
       ajc    = 0.
       do 104 iilw=4,nilw
       t4s(j) = t4s(j) + ajsu0(j,iilw) * t4ssu0(j,iilw)
       if (t4ssu0(j,iilw).ne.0.) ajc = ajc + ajsu0(j,iilw)
 104   continue
c
       if (ajc.gt.0.) then
        t4s(j) = t4s(j) / ajc
       else
        t4s(j) = t4su0(j)
        if (t4s(j).gt.tfsth4) t4s(j) = tfsth4
       end if
c
       if ( ajsu0(j,ic).gt.0..and.t4ssu0(j,ic).eq.0.) then
        t4ssu0(j,ic) = t4s(j)
        t4sla3(j,ic) = t4s(j)
       end if
c
 103   continue
c
      end if
c
      do 13 j=1,jp
c
       if (ajsu0(j,ic).gt.0.) then
        hssn2(j,ic) =   0.5
       else
        hssn2(j,ic) =   0.
       end if
c
        hssn3(j,ic) = hssn2(j,ic)
c ...      epaisseur de la neige (m)
c
       rsusn5(j,ic) = .850
c petz rsusn5(j,ic) = .850
       rslsn5(j,ic) = .7
c
         nsn6(j,ic) = 0
        nnsn6(j,ic) = 0
       hfasn6(j,ic) = 0.
c
   13 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tcal====
c cray .endif
      subroutine tcalot(ic)
c
c +++ bcmc3.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical fulwri
      double precision  chhv0,cshv0,htshv0
      double precision  hssn2,hspsn2
      double precision  r2p,r3p,r2n,r3n,flux,bil
      dimension r2p(np),r3p(np)
      dimension r2n(np),r3n(np),flux(np),bil(np)
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables flux de chaleur
c
       common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
       common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
       common/hv1/htshv1(6,np),shshv1(np)
c
c --- variables surface
c
       common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
       common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .            vopsn3(np,nilw), vosn3(np,nilw),
     .            hspsn3(np,nilw), hssn3(np,nilw),
     .            ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .           ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
       common/la0/cslla0(np),csila0(np),csla0(np),tcla0(np)
       common/la1/tcmla1(np),tcdla1(np),csmla1(np),csdla1(np)
       common/la3/t4sla3(np,nilw)
       common/bi0/surbi0(np,nilw),atmbi0(np),topbi0(np)
       dimension fsnos(np)
c
c **********************************************************************
c
c --- initialisation
c
      fulwri = .false.
c
      if (it.eq.0) then
       call inica1(ic)
       call inica2(ic)
c
c **********************************************************************
c
      else
c
       if (season) then
c
c --- energy before computation
c
        do 1 j=1,jp
        r2p(j) = t4rsu0(j,ic) * ajsu0(j,ic) / csila0(j)
        r3p(j) =  hssn2(j,ic) * ajsu0(j,ic) * hfsth4 * c33th4
        do 2 n=1,6
        htshv0(n,j,ic) = 0.
 2      continue
 1      continue
c
        pr = 1.
c ...   pr = 1. => schema numerique pronostique
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
       else
c
        pr = 0.
c ...   pr = 0. => schema numerique diagnostique
c
       end if
c
c **********************************************************************
c
c --- modele de sol (deardorff, 1978)
c
c --- continent recouvert de neige et champ de neige
c
       do 212 j= 1,jp
       dtsnow =  0.
c
       if (ajsu0(j,ic).ne.0.) then
        scs = 0.
        sch = 0.
        cshv0(6,j,ic) = cshv0(4,j,ic)*rhsth4
        chhv0(6,j,ic) = chhv0(4,j,ic)*rhsth4
c ...   pour tenir compte de la partie chaleur de fonte retiree 
c       du manteau neigeux lors de la sublimation
        do 213 n=1,6
        scs = scs + cshv0(n,j,ic)
        sch = sch + chhv0(n,j,ic)
 213    continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if (season) then
         t4sc          = t4sla3(j,ic)
c ...    temperature du champ de neige
c                    a la fin du pas de temps precedent
c
         cdsh          =  cdsth4  *.6   /  hssn3(j,ic)
         schr          =  sch *    cdsh / (scs+cdsh)
         scsr          = cdsh *(1.-cdsh / (scs+cdsh))
         t4ro          = t4rsu0(j,ic)
         call tsol(t4rsu0(j,ic),csila0(j),scsr,schr,t4dsu1(j,ic),pr,
     .             alpha,beta)
         t4rn          = alpha *t4rsu0(j,ic) + beta *t4ro
         t4sla3(j,ic)  = (sch +cdsh *t4rn) / (scs +cdsh)
c
         if (t4sla3(j,ic).gt.tfsth4) then
          t4rsu0(j,ic)  = t4ro
          schr          = cdsh * tfsth4
          scsr          = cdsh
          call tsol(t4rsu0(j,ic),csila0(j),scsr,schr,t4dsu1(j,ic),pr,
     .              alpha,beta)
          t4rn          = alpha *t4rsu0(j,ic) + beta *t4ro
          f1            = (sch +cdsh *t4rn) - (scs +cdsh) *tfsth4
          dtsnow        = t4sla3(j,ic) - t4sc
          if (f1.lt.0.d0) t4sla3(j,ic) = tfsth4
          t4hts         = tfsth4
         else
          f1            = 0.
          t4hts         = t4sla3(j,ic)
         end if
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        else
         t4ro          = t4sla3(j,ic)
         call tsol(t4sla3(j,ic),csila0(j),scs,sch,tcdla1(j),pr,
     .             alpha,beta)
         t4hts         = alpha *t4sla3(j,ic) + beta *t4ro
        end if
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if (tsclim) t4hts = t4ssu0(j,ic)
c
        do 214 n=1,6
        htshv0(n,j,ic)= chhv0(n,j,ic) - cshv0(n,j,ic) *t4hts
 214    continue
c
       else
        f1           = 0.
        dtsnow       = 0.
        t4sla3(j,ic) = 0.
        t4rsu0(j,ic) = t4sla3(j,ic)
       endif
c
c **********************************************************************
c
       if (season) then
c
c  -- variation des fractions continent-neige --
c
        ajcont=   ajsu0(j,ic)
        call ajland(dtsnow,f1,csila0(j),ajcont,hfal,hmel,htasn,fsnos,
     .              j,ic,ic)
c
        if (ajcont.gt.0.) then
         accsn3(j,ic) = accsn3(j,ic) + hfal
         ablsn3(j,ic) = ablsn3(j,ic) + hmel
c
         vosn3(j,ic) =  vosn3(j,ic) + ajsu0(j,ic) * hssn2(j,ic)
         ajsn3(j,ic) =  ajsn3(j,ic) + ajsu0(j,ic)
c
         if (hmel.gt.0.d0) then
          hspsn3(j,ic)=  0.01d0
         else
          hspsn3(j,ic)=  hssn3(j,ic) + hfal
         end if
c
         hssn3(j,ic) = hspsn3(j,ic) + .1d0*( 0.01d0 -hspsn3(j,ic))
     .                                    *(10.00d0 +t4sc -tfsth4)
         if (hssn3(j,ic).lt.0.01d0      ) hssn3(j,ic) = 0.01d0
         if (hssn3(j,ic).gt.hspsn3(j,ic)) hssn3(j,ic) = hspsn3(j,ic)
        end if
c
       end if
c
  212  continue
c
c **********************************************************************
c
       if (season) then
c
c --- energy input
c
        do 3 j=1,jp
        flux(j) = 0.
        do 4 n=1,6
        flux(j) = flux(j) + ajsu0(j,ic)*htshv0(n,j,ic)
 4      continue
        flux(j) = flux(j) * delt
 3      continue
c
c --- energy after  computation
c
        do 7 j=1,jp
        r2n(j) = t4rsu0(j,ic) * ajsu0(j,ic) / csila0(j)
        r3n(j) =  hssn2(j,ic) * ajsu0(j,ic) * hfsth4 * c33th4
 7      continue
c
c --- energy balance
c
        do 5 j=1,jp
        bil(j) = (r2n(j) - r2p(j) - r3n(j) + r3p(j) - flux(j)) / delt
c       if (.not.tsclim.and.abs(bil(j)).gt..05)
c    .  write(iwr6,6)it,ic,j,bil(j)
 6      format('  -- tcalo no',i2,' --  * error *  iter',i4,2x,'j =',i2,
     .          2x,'bilan =',e12.4)
 5      continue
c
c --- print
c
        if (.not.fulwri) go to 60
        if (iprint.eq.1) go to 61
        if (mts.ne.0)    go to 60
        if (monpr.gt.mm) go to 60
 61     continue
        write(iwr1,62)
 62     format(/'  -- tcalot --')
        write(iwr1,63)
 63     format(5x,'j',9x,'r2p',9x,'r3p',8x,'flux',
     .                9x,'r2n',9x,'r3n',9x,'bil')
        write(iwr1,64)(j,r2p(j),r3p(j),flux(j),r2n(j),r3n(j),bil(j),
     .                 j=1,jp)
 64     format((i6,6e12.4))
 60     continue
c
       end if
c
c **********************************************************************
c
c --- moyennes
c
       do 400 j=1,jp
c
       do 402 n=1,6
       htshv1(n,j)  = htshv1(n,j)  + ajsu0(j,ic)*htshv0(n,j,ic)
       surbi0(j,ic) = surbi0(j,ic) + ajsu0(j,ic)*htshv0(n,j,ic)
 402   continue
 400   continue
c
       if (.not.tsclim) then
        do 410 j=1,jp
        t4ssu0(j,ic) = t4sla3(j,ic)
 410    continue
       end if
c
      end if
c
c **********************************************************************
c
c --- print sur fichier bcmcal.d
c
      if (.not.mts.eq.0) go to 600
c
      write(iwr5,641)(t4sla3(j,ic),j=1,jp)
  641 format((9f8.2,4x,'tsno'))
      write(iwr5,642)(ajsu0(j,ic),j=1,jp)
  642 format((9f8.3,4x,'aj c'))
      write(iwr5,648)(hssn3(j,ic),j=1,jp)
  648 format((9f8.3,4x,'hs c'))
      write(iwr5,643)(htshv0(1,j,ic),j=1,jp)
  643 format((9f8.2,4x,'h1 c'))
      write(iwr5,644)(htshv0(2,j,ic),j=1,jp)
  644 format((9f8.2,4x,'h2 c'))
      write(iwr5,645)(htshv0(3,j,ic),j=1,jp)
  645 format((9f8.2,4x,'h3 c'))
      write(iwr5,646)(htshv0(4,j,ic),j=1,jp)
  646 format((9f8.2,4x,'h4 c'))
      write(iwr5,647)ian,mm,ic
  647 format(' cycle annuel no',i3,5x,'mois',i3,5x,'calotte no',i3)
c
  600 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  rfcal===
c cray .endif
      subroutine reflca(il,is)
c
c +++ bcmcr.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      double precision  chhv0,cshv0,htshv0
      double precision  hssn2,hspsn2
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables rayonnement
c
      common/so11/zso11(np),czso11(np)
c
c --- variables flux de chaleur verticaux
c
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables vapeur d'eau
c
      common/ql2/qlql2(np),fsql2(np)
c
c --- variables surface
c
      common/rf0/rgsrf0(np,nilw)
      common/rf1/rgwrf1(np,12),rgarf1(np,nilw)
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/su3/hswsu3(np,nilw),cksu3(np,nilw)
      common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
      common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .           vopsn3(np,nilw), vosn3(np,nilw),
     .           hspsn3(np,nilw), hssn3(np,nilw),
     .           ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
      common/sn5/rsusn5(np,nilw),rslsn5(np,nilw),rssn5(np,nilw)
      common/sn6/hfasn6(np,nilw),  nsn6(np,nilw),nnsn6(np,nilw)
c
c **********************************************************************
c
c --- initialisation
c
      do 20 j=1,jp
c
c --- secteurs champ de neige sur continent et calottes
c
      if (ajsu0(j,is).ne.0.) then
c
c --- w factor (sa76)
c
           hswsu3(j,is) =   .4  + .06  *(t4ssu0(j,is)-263.05)
       if (hswsu3(j,is).gt.1. )   hswsu3(j,is)=1.
       if (hswsu3(j,is).lt.0.4)   hswsu3(j,is)=0.4
c
c --- heat diffusivity (sa76)
c
            cksu3(j,is) = 16.e2
c           cksu3(j,is) =  6.e2 + 1.e2 *(t4ssu0(j,is)-263.05)
c      if ( cksu3(j,is).gt.16.e2) cksu3(j,is)=16.e2
c      if ( cksu3(j,is).lt. 6.e2) cksu3(j,is)= 6.e2
c
c --- albedo
c
       if (season) then
c
        if (t4ssu0(j,is).lt.tfsth4) then
         a0 = rsusn5(j,is)
         ai = rsusn5(j,is) - .20
c        ai = rsusn5(j,is) - .35
c
         if (hfasn6(j,is).gt.0.001) then
          nd = .5 + 1./hfasn6(j,is)
          nm = min(nd,nnsn6(j,is))
         else
          nm =        nnsn6(j,is)
         end if
         an = a0
c
         if (nm.ge.1) then
          do 322 n=1,nm
          an = an - (10.**(.78-.069*n  )) /100.
 322      continue
         end if
         an = amax1(an,ai)
        else
         a0  = rslsn5(j,is)
c
         if (hfasn6(j,is).gt.0.001) then
          nd = .5 + 1./hfasn6(j,is)
          nm = min(nd,nsn6(j,is))
         else
          nm =        nsn6(j,is)
         end if
         an  =      a0
c
         if (nm.ge.1) then
          do 302 n=1,nm
          an = an - (10.**(1.05-.07*n  )) /100.
 302      continue
         end if
        end if
c
        amin = .40
        amin = amin1(amin,rslsn5(j,is))
        an   = amax1(an  ,amin        )
c
        rgsrf0(j,is) = an
c
        znsa = zso11(j)
c ...   znsa : noon solar elevation angle
c              (danard et al., mwr, 1984, p.1162)
c
        if (znsa.lt.10.) znsa = 10.
        if (znsa.gt.40.) znsa = 40.
        das  = (10. - (znsa - 10.) / 3.) / 100.
        rgsrf0(j,is) = rgsrf0(j,is) + das
        if (rgsrf0(j,is).gt..9999) rgsrf0(j,is)=.9999
c
        if (il.ne.is) then
         dsnow = hssn2(j,is)
         dsnow = amin1(.1,dsnow) 
         rgsrf0(j,is) = rgsrf0(j,il) 
     .                +(rgsrf0(j,is) - rgsrf0(j,il)) * dsnow / .1
c ...    raccord avec albedo sol sans neige et sature d'eau
c              dans regions polaires =.16
c
        end if
       end if
      else
       hswsu3(j,is)=1.
        cksu3(j,is)=0.
       rgsrf0(j,is)=0.
      end if
c
   20 continue
c
c **********************************************************************
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  dav=====
c cray .endif
      subroutine dav(x,xav,cp)
c
c +++ bcmdav.for +++
c
      parameter (np=18)
c _dp implicit double precision (a-h,o-z)
      double precision  x,xcos,xav,dcp
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      dimension x(np)
          dcp  = cp
          xav  = 0.d00
      do 5 j = 1,np
          xcos = x(j)*(dsin(dcp*dble(j))-dsin(dcp*dble(j-1)))
          xav  = xav+xcos
    5 continue
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  esat====
c cray .endif
      function esat(tt)
c
c +++ bcmes.for +++
c
c ...          esat : tension de vapeur saturante (mb)
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      if (tt.ge.273.) then
       esl  = 9.4051-2354./tt
       esat = 10.**esl
      else
       esl  = 10.553-2667./tt
       esat = 10.**esl
      end if
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  inigr===
c cray .endif
      subroutine inigr
c
c +++ bcmg0.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
      common/gr2/argr2(np),atgr2
c
c **********************************************************************
c
c --- initialisation des variables auxiliaires de la grille
c
          pi2r =  ca * pi * 2.
          p2r2 =  ca * pi2r
          atgr2=  p2r2
c
      write(iwr1,20)
   20 format (//' ---- inigr ----',//,
     . 5x,'j',9x,'deg',8x,'degc',6x,'s10gr1',6x,'c10gr1',8x,'area')
      do 32 jj=1,jpp
                 j  = jpp + 1 - jj
          dgcgr1(j) = (j-1)*gs
            prad    = cp*real(j-1)
          c10gr1(j) = cos(prad)
          s10gr1(j) = sin(prad)
c
      if (j.eq.jpp) then
       write(iwr1,21)jpp,dgcgr1(jpp),s10gr1(jpp),c10gr1(jpp)
 21    format(i6,12x,f12.1,2f12.4)
      else
          deggr1(j) = j*gs - gs2
          argr2 (j) = p2r2*(s10gr1(j+1)-s10gr1(j))
       write(iwr1,22)j,deggr1(j),dgcgr1(j),s10gr1(j),c10gr1(j),argr2(j)
 22    format(i6,2f12.1,2f12.4,1x,e11.4)
      end if
 32   continue
c
c **********************************************************************
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  hyd=====
c cray .endif
      subroutine inihyd
c
c +++ bcmh1.for +++
c
c *** modifier egalement les arguments d'entree de intpol
      parameter(np=18,npp2=37)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c --- variables precipitations
c
       common/ql1/qlql1(12,np),esql1(12,np),qlaql1(np),esaql1(np)
c
       dimension qlv0(24,12),qld0(24),qld1(46),qld2(npp2)
       dimension t4s0(12,9),t4a0(9),x(9),y(np)
       dimension qlv0a(144),qlv0b(144)
       equivalence(qlv0(1,1),qlv0a(1))
       equivalence(qlv0(1,7),qlv0b(1))
c
       data qlv0a/
     .4.514,4.629,3.642,2.615,1.782,1.523,1.170,1.434,1.747,2.085,2.483,
     .2.387,2.292,2.184,2.021,1.674,1.173, .782, .727, .714, .564, .257,
     . .176, .167,
     .4.752,4.779,3.315,2.309,1.619,1.408,1.176,1.540,1.901,2.250,2.591,
     .2.531,2.343,2.220,2.057,1.686,1.175, .825, .757, .720, .582, .606,
     . .637, .641,
     .5.016,4.938,3.261,2.157,1.405,1.145,1.006,1.293,1.624,1.984,2.402,
     .2.148,2.003,1.944,1.811,1.381,1.022, .691, .540, .467, .467, .563,
     . .541, .538,
     .5.419,6.213,4.309,2.832,1.744,1.335,1.076,1.519,1.789,1.991,2.176,
     .2.145,1.986,1.871,1.738,1.488,1.087, .844, .605, .461, .521, .135,
     . .029, .017,
     .4.822,6.540,6.037,4.530,2.966,1.961,1.796,1.858,1.875,1.861,1.819,
     .2.128,2.085,1.982,1.839,1.592,1.191, .903, .633, .480, .563, .143,
     . .067, .059,
     .3.792,6.331,6.863,5.990,4.747,3.530,2.643,2.529,2.326,2.085,1.870,
     .2.151,2.327,2.311,2.231,2.100,1.794,1.258, .802, .488, .425, .181,
     . .089, .078 /
      data qlv0b /
     .3.264,4.980,6.184,5.861,4.937,3.543,3.215,2.771,2.613,2.414,1.833,
     .1.849,2.267,2.384,2.410,2.410,1.987,1.587,1.116, .719, .503, .221,
     . .103, .089,
     .3.268,4.860,6.041,6.155,5.636,4.095,3.273,2.804,2.509,2.200,1.737,
     .1.967,2.151,2.261,2.343,2.386,2.121,1.811,1.347, .898, .541, .311,
     . .124, .099,
     .3.393,4.752,6.157,6.199,5.615,4.295,2.848,2.361,2.178,2.038,1.830,
     .2.067,2.112,2.163,2.265,2.463,2.000,1.563,1.183, .865, .640, .633,
     . .485, .464,
     .3.572,5.178,5.894,5.176,4.107,3.172,1.832,1.648,1.667,1.732,1.838,
     .2.082,2.132,2.136,2.151,2.188,1.685,1.283,1.004, .774, .577, .574,
     . .481, .469,
     .4.118,4.872,5.212,4.332,3.182,2.274,1.307,1.374,1.594,1.854,2.164,
     .2.361,2.460,2.347,2.228,2.204,1.381,1.026, .895, .808, .716, .939,
     . .977, .981,
     .4.293,4.717,4.514,3.481,2.384,1.820,1.328,1.383,1.638,1.969,2.380,
     .2.335,2.360,2.223,2.027,1.734,1.238, .860, .807, .782, .619, .591,
     . .565, .561/
c
c ... qlv0 : moyenne zonale precipitations (mm.day-1) (jaeger, 1976)
c            moyenne mois par mois
c
       data qld0/
     .4.181,5.241,5.128,4.314,3.353,2.514,1.894,1.878,1.955,2.037,2.090,
     .2.176,2.209,2.169,2.093,1.945,1.490,1.121, .869, .681, .560, .428,
     . .353, .344/
c
c ... qld0 : moyenne zonale precipitations (mm.day-1) (jaeger, 1976)
c            moyenne annuelle
c
       data t4s0/
     . 26.0, 26.5, 26.7, 27.0, 26.8, 26.5,
     . 26.1, 26.0, 26.1, 26.2, 26.2, 26.1,
c ...  data  5 deg n.
     . 24.0, 24.5, 25.6, 26.6, 27.7, 27.6,
     . 27.3, 27.1, 27.3, 27.0, 25.9, 24.6,
c ...  data 15 deg n.
     . 18.0, 18.5, 20.7, 22.7, 24.8, 26.6,
     . 27.2, 27.3, 26.8, 24.7, 22.1, 19.4,
c ...  data 25 deg n.
     .  8.2,  8.9, 11.1, 14.3, 17.3, 21.0,
     . 23.5, 23.8, 21.7, 17.7, 13.2, 10.0,
c ...  data 35 deg n.
     .- 2.4,- 1.0,  2.0,  7.3, 11.7, 15.5,
     . 18.4, 18.7, 15.5, 10.4,  4.3,  0.0,
c ...  data 45 deg n.
     .-10.6,- 9.5,- 5.3,  1.5,  7.2, 11.7,
     . 14.5, 13.8, 10.1,  4.1,- 3.0,- 8.3,
c ...  data 55 deg n.
     .-22.6,-20.7,-15.9, -7.0,  2.1,  8.7,
     . 12.3, 10.4,  4.9, -3.9,-14.1,-20.3,
c ...  data 65 deg n.
     .-24.8,-25.4,-24.5,-18.5,- 8.4,- 0.6,
     .  2.8,  1.4,- 2.6,-10.4,-19.2,-22.7,
c ...  data 75 deg n.
     .-33.9,-31.8,-30.9,-24.6,-12.1, -1.8,
     .  1.5, -0.1, -9.3,-18.3,-25.5,-29.4/
c ...  data 85 deg n.
c
       data t4a0/
     . 26.3, 26.3, 23.2, 15.9,  8.4,  2.2, -5.5,-12.7,-18.0/
                   mt=1
       if (season) mt=12
       do 10 mm=1,mt
c
       if (season) then
        do 22 j=1,23
        qld1(2*j-1)=.5*(qlv0(j,mm)+qlv0(j+1,mm))
        qld1(2*j  )=               qlv0(j+1,mm)
 22     continue
       else
        do 23 j=1,23
        qld1(2*j-1)=.5*(qld0(j)    +qld0(j+1))
        qld1(2*j  )=                qld0(j+1)
 23     continue
       end if
       call intpol(qld1,qld2,46,37)
c
       if (season) then
        do 32 j=1,np
        qlql1(mm,j)=qld2(2*j)
 32     continue
       else
        do 33 j=1,np
        qlaql1(j)=qld2(2*j)
 33     continue
       end if
c
       if (season) then
        do 42 j=1,9
        x(j) = t4s0(mm,j) + 273.15
 42     continue
        call pol(x,y,9,np)
        do 43 j=1,np
        esql1(mm,j) = esat(y(j))
 43     continue
       else
        call pol(t4a0,y,9,np)
        do 44 j=1,np
        y     (   j) =      y(j)  + 273.15
        esaql1(   j) = esat(y(j))
 44     continue
       end if
 10    continue
c
      write(iwr1,62)(i,i=1,12)
   62 format(
     . /,6x,' precipitations moyennes zonales',/,7x,31('-'),
     . /,6x,'j',4x,12i8)
      write(iwr1,63)(j,(qlql1(i,j),i=1,12),j=1,np)
   63 format((i7,4x,12f8.3))
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  precip==
c cray .endif
      subroutine precip(alhhm,evapa)
c
c +++ bcmh3.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables vapeur d'eau
c
      common/ql1/qlql1(12,np),esql1(12,np),qlaql1(np),esaql1(np)
      common/ql2/qlql2(np   ),fsql2(np)
      dimension qlevp(np),evapa(np)
c
c --- variables surface
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/ca2/caca2(np,nilw),daca2(np,nilw)
c
c **********************************************************************
c
      if (season) then
c
c --- parametres temporels de l'interpolation
c
       ita=int(12.*tsm)
       itsm=int(tsm)
       itm=mod(it,ita)
       itj=mod(itm,itsm)
       i  =1+itm/itsm
       ip1=mod(i+1,12)
       if(ip1.eq.0) ip1=12
c
c --- interpolation des precipitations
c
       do 12  j=1,jp
       qlql2(j)= qlql1(i,j) + itj *(qlql1(ip1,j) -qlql1(i,j)) /tsm
c _pal esdat   = esql1(i,j) + itj *(esql1(ip1,j) -esql1(i,j)) /tsm
c _pal qlql2(j)= qlql2(j) * esat(t4su0(j)) / esdat
c
       qlevp(j)= 0.
       ajcal   = 0.
       do 120 ic = 5,nilw
       qlevp(j)= qlevp(j) + qlql2(j) *caca2(j,ic) *ajsu0(j,ic)
       ajcal   = ajcal    +                        ajsu0(j,ic)
 120   continue
       qlevp(j)= qlevp(j) + qlql2(j) *(1.-ajcal)
 12    continue
      else
       do 13  j=1,jp
       qlql2(j)= qlaql1(j)
       qlevp(j)= qlaql1(j)
c _pal esdat   = esaql1(j)
c _pal qlql2(j)= qlql2(j) * esat(t4su0(j)) / esdat
 13    continue
      end if
c
c --- ajustement sur le total d'evaporation
c
      qls = - alhhm * 86400.e3 / hvsth4
      call av(qlevp,qlav,cp)
      qlr = qls / qlav
      if (it.eq.0) qlr=1.
      do 14 j=1,jp
      qlql2(j) = qlql2(j) * qlr
      qlevp(j) = qlevp(j) * qlr
 14   continue
c
c --- liberation moyenne de chaleur latente
c
      do 210 j=1,jp
      evapa(j) = (hvsth4 /86400.e3) * qlevp(j)
 210  continue
c
c **********************************************************************
c
c --- print at middle of the month
c
      if (iprint.eq.1) go to 61
      if (mts.ne.0)    go to 60
      if (monpr.gt.mm) go to 60
   61 continue
      write(iwr1,62)qlr,(qlql2(j),j=1,jp),(evapa(j),j=1,jp)
   62 format(/'   -- precip --  ql corr. = ',f12.3,2(/,9f8.3),
     .       /'      evapa  :',2(/,9f8.3))
   60 continue
c
c **********************************************************************
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  hadley==
c cray .endif
      subroutine hadley
c
c +++ bcmh4 +++
c
c ... parametrisation du transport de chaleur du a la cellule de hadley
c     ref. : peng, chou & arking, j.g.r., 1987, pp.5519-5520
c            held & hou,          j.a.s., 1980, pp. 515- 533
c
      parameter(np=18,npp=np+1,nilw=7,nhd=2)
      parameter(nd=7,nd2=(nd-1)/2,nda=120,nda2=1+nda/2)
      parameter(ndi1=nda2-nd2+1,ndi2=nda2+nd2)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical fulwri
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,jour,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
      common/gr2/argr2(np),atgr2
c
c --- variables radiatives
c
      common/ir0/stfir0
      common/ir3/fupir3(np),fdnir3(np),fnir3(np)
c
c --- variables thermodynamique
c
      common/th1/cpath1,rath1,akth1
      common/ql2/qlql2(np   ),fsql2(np)
      common/ql3/dqyql3(npp,nhd),omgql3(np)
c
c --- variables dynamique
c
      common/p00/p4p00,p3p00,p2p00,p1p00
      common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
      common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
      common/dy2/t2ady2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .                 pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
      common/dy3/copdy3(np),comdy3(np),coedy3(np),idy3
      dimension th2mav(0:np),vth (npp),hta6(np)
      dimension tpq   (  np),akth(npp)
      dimension apq(np),bpq(np),cpq(np),dpq(np)
      common/h40/ph0h40(nd ),phih40(nd ),pmnh40,pmxh40
      common/h41/vthh41(nda),th2h41(nda),th0h41,phih41
c
c --- variables surface
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/hv1/htshv1(6,np),shshv1(np)
c
      data a1 /8.7e-4/
      data a2 /2.1e-2/
      data cb /2.0e-4/
      data ch /1.5e+7/
      data hh /15.e+3/
      data ct /1.4e+1/
c
      fulwri = .false.
      pi2 = pi * 2.d0
      if (it.gt.1) then
c
c ---  champ temperature a 500mb : lissage meridien
c
       th2mav(1) =                         th2dy2(1)
       do 13 j=2,jp-1
       th2mav(j) = .25 * (th2dy2(j-1) + 2.*th2dy2(j) + th2dy2(j+1))
 13    continue
       th2mav(jp)=                         th2dy2(jp)
c
c ---  temperature de rappel
c
       do 130 ida = nda-1,1,-1
       th2h41(ida+1) = th2h41(ida)
 130   continue
       th2h41(    1) = th2mav(1)
       th2mav(    0) = th2h41(1)
       do 131 ida = 2,nda
       th2mav(    0) = th2mav(0) + th2h41(ida)
 131   continue
       th2mav(    0) = th2mav(0) / nda
c
c ---  localisation du maximum de temperature potentielle
c
        thmax = 0.
       do 140 j=1,jp
       if (th2mav(j).gt.thmax) then
        thmax = th2mav(j)
         jmax = j
       end if
 140   continue
c
c --- latitude de la zone de convergence intertropicale (itcz)
c
       phi0 = cp
     .      + a1 * ( (htshv1(1,2) +htshv1(2,2)) 
     .              -(htshv1(1,1) +htshv1(2,1)) ) / cp 
c    .      + a2 * (   t4su0(  2) - t4su0(  1)  ) / cp
c ...  phio : latitude itcz (parametrisation de peng & al.)
c
c      phi0 = (jmax -.5) * cp
c ...  phi0 : latitude itcz (parametrisation gall)
c
       do 101 id = nd-1,1,-1
       ph0h40(id+1) = ph0h40(id)
 101   continue
       ph0h40(1   ) = phi0
c
       phi0 =        ph0h40(1)
       do 102 id = 2,nd
       phi0 = phi0 + ph0h40(id)
 102   continue
       phi0 = phi0 / real(nd)
c
       phi0d= phi0 * gs / cp
c
c --- calcul du gradient moyen equateur-pole de temperature potentielle
c
       j = 1
 111   continue
       if (       j .ge.jp   ) go to 110
       if (deggr1(j).gt.phi0d) go to 110
       j = j + 1
       go to 111
 110   continue
       if (j.eq.1) then
        phimin = deggr1(j)
        thitcz = th2dy2(j)
       else
        phimin = phi0d
        thitcz = th2dy2(j) + (phi0d     - deggr1(j-1))
     .                      *(th2dy2(j) - th2dy2(j-1)) / gs
       end if
       dth2  = (thitcz - th2dy2(jp)) * (deggr1(1) - deggr1(jp))
     .                               / (phi0d     - deggr1(jp)) / 300.
c
c ---  latitude de la branche nord de la cellule de hadley
c
       phiip = sqrt(cb * hh * dth2)
c ...  phiip : latitude branche nord de la cellule de hadley
c              parametrisation de peng & al.
c
       j = int(70./gs)
c ...          70 deg.n : limite nord zone barocline (green, 1970)
c
 121   continue
       if (dqyql3(j,2).gt.0.) go to 120
       j = j - 1
       if (j          .le.2 ) go to 120
       go to 121
 120   continue
       if (j          .gt.2) then
        ddqy  = dqyql3(j+1,2) -dqyql3(j,2)
        if (ddqy.ne.0.) then
         phiid = dgcgr1(j) - dqyql3(j,2) * gs /ddqy
        else
         phiid = phih41
        end if
       end if
        phih41= phiid
        phii  = phiid * pi / 180.d0
c ...   phii : latitude branche nord de la cellule de hadley
c              link sur bordure sud de la zone barocline
c
       phii = amax1(phii,phiip)
       phii = amin1(phii,pmxh40)
c
       do 122 id = nd-1,1,-1
       phih40(id+1) = phih40(id)
 122   continue
       phih40(1   ) = phii
c
       phii =        phih40(1)
       do 123 id = 2,nd
       phii = phii + phih40(id)
 123   continue
       phii = phii / real(nd)
       phiid= phii * gs / cp
c
       dphi = phii - phi0
c
       if (iprint.eq.1) go to 301
       if (mts.ne.0)    go to 300
       if (monpr.gt.mm) go to 300
 301   continue
       write(iwr1,30)it,phi0d,phiid,dth2
       if (fulwri) write(iwr6,30)it,phi0d,phiid,dth2
 30    format(/,' caracteristiques cellule de hadley :',i5,
     .        /,' latitude  itcz   : ',f10.2,' deg.n',
     .        /,' extension hadley : ',f10.2,' deg.n',
     .        /,' dh               = ',f12.4)
 300   continue
c
c --- calcul du coefficient de diffusion de chaleur meridien
c
        j = 0
 151   continue
        j = j + 1
        phi = cp * real(j-1)
       if (j  .ge.jp  ) go to 150
       if (phi.gt.phii) go to 150
        yi  = phi         / phii
        yis =(phi - phi0) / dphi
        akth(j) =  ch * phiip*phiip * (1.d0 -yis*yis) * (1.d0 -yi*yi)
       go to 151
 150   continue
        jh  = j
        do 152 j=jh,jpp
        akth(j) = 0.d0
 152    continue
c
c --- calcul contribution de hadley a l'evolution de la temperature
c     utilisation du schema crank-nicholson
c
       do 153 j=1,jp
       tpq(j) = th2mav(j)
       apq(j) = tht * akth(j+1) * copdy3(j)
       cpq(j) = tht * akth(j  ) * comdy3(j)
       bpq(j) = apq(j) + cpq(j) + 1.d0
 153   continue
       do 154 j=2,jp-1
       dpq(j) = apq(j)*tpq(j+1) -(bpq(j)-2.d0)*tpq(j) + cpq(j)*tpq(j-1)
 154   continue
       dpq(1) = apq(1)*tpq(  2) -(bpq(1)-2.d0)*tpq(1) 
     .        + cpq(1)*th2mav(0)*2.d0
       dpq(jp)= cpq(jp)*tpq(jp-1)-(bpq(jp)-2.d0)*tpq(jp)
       logequ = 0
       call pq(apq,bpq,cpq,dpq,tpq,logequ)
       do 155 j=1,jp
       th2mav(j) = .5 *(th2mav(j) + tpq(j))
 155   continue
c
       do 156 j=1,jp
        vth( j) =- akth(j) * ((th2mav(j) -th2mav(j-1)) /cp)
     .          *  pi2 * c10gr1(j) * cpath1 * p4p00 / gdy0
c ...   vth     : flux de chaleur du a hadley [w]
c
 156   continue
c
       do 157 ida =nda-1,1,-1
       vthh41(ida+1) = vthh41(ida)
 157   continue
       vthh41(    1) = vth   (1)
       vthav         = vthh41(nda)
       do 158 ida =1,nda-1
       vthav         = vthh41(ida) + vthav
 158   continue
       vthav         = vthav / real(nda)
c
       vth(1 )       = vth(1) - vthav
c
c ---  calcul du terme de chauffage diabatique mmc
c
       do 20 j=1,jh-1
       htady2(6,j) =(vth(j) - vth(j+1)) / argr2(j)
       hta6  (  j) =(p4p00 /gdy0) * cpath1 * sigdy0 * p2p00 * omgql3(j)
     .                            * (f0dy0/fjdy1(j) - 1.) / rath1
 20    continue
c
       if (iprint.eq.1) go to 311
       if (mts.ne.0)    go to 310
       if (monpr.gt.mm) go to 310
 311   continue
       write(iwr1,31)
       if (fulwri) write(iwr6,31)
 31    format(/,' latitude   flux chaleur(w)',
     .' latitude    temp. 500mb(k)   smooth 500mb(k)',
     .'  chauffage(w.m-2)   chauffage (om.)')
       write(iwr1,33) vthav,th2mav(0)
       if (fulwri) write(iwr6,33) vthav,th2mav(0)
 33    format(9x,e18.5,27x,f18.2)
       write(iwr1,32)(dgcgr1(j),vth(j),deggr1(j),t2dy2(j),th2mav(j),
     .                htady2(6,j),hta6(j),j=1,jh)
       if (fulwri) write(iwr6,32)(dgcgr1(j),vth(j),deggr1(j),
     .              t2dy2(j),th2mav(j),htady2(6,j),hta6(j),j=1,jh)
 32    format((f9.1,e18.5,f9.1,4f18.2))
 310   continue
      else
       do 40 id = 1,nd
       ph0h40(id)   = 0.d0
       phih40(id)   = 0.d0
 40    continue
c
       pmnh40       = cp * 30./gs
       pmxh40       = cp * 50./gs
       th0h41       = 266.*(2.**akth1)
       do 41 ida = 1,nda
       th2h41(ida)  = th0h41
 41    continue
      end if
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  intpol==
c cray .endif
      subroutine intpol(a,b,naa,nbb)
c
c +++ bcmipl.for +++
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      dimension a(naa),b(nbb)
c
      na  = naa-1
      nb  = nbb-1
      xt  = 1.
      dxa = xt /     real(na)
      dxb = xt /     real(nb)
      b(  1)=a(  1)
c
      if (nb.ge.na) then
       ja = 0
       jb = 0
c
c
c do until
 110   continue
       ja = ja + 1
       xa = dxa*real(ja-1)
       xa1= dxa*real(ja)
       da = a(ja+1)-a(ja)
c
c  do until
 1100  continue
       jb = jb + 1
       xb = dxb*real(jb-1)
       dx = xb-xa
       b(jb) = a(ja) + dx * da / dxa
       if (.not.(xb.ge.xa1.or.jb.ge.nbb-1) ) go to 1100
c
       if (.not.(ja.ge.naa-1)               ) go to 110
c
c  do while
 1110  continue
       if (jb.ge.nbb) go to 1111
       jb = jb + 1
       xb = dxb*real(jb-1)
       dx = xb-xa
       b(jb) = a(ja) + dx * da / dxa
       go to 1110
 1111  continue
c
c
      else
       ja = 1
       jb = 1
c
c
c do until
  100  continue
       jb = jb + 1
       xb = dxb*real(jb-1)
c
c  do until
 1000  continue
       ja = ja + 1
       xa = dxa*real(ja-1)
       xa1= dxa*real(ja)
       if (.not.((xb.ge.xa.and.xb.lt.xa1).or.(ja.ge.naa-1))) go to 1000
c
       if (jb.ge.nbb) go to 1001
       da = a(ja+1)-a(ja)
       dx = xb-xa
       b(jb) = a(ja) + dx * da / dxa
 1001  continue
c
       if (.not.(jb.ge.nbb-1)) go to 100
c
       b(nbb)=a(naa)
c
c
      end if
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  iniir===
c cray .endif
      subroutine iniir
c
c +++ bcmir0.for +++
c
c *** infrared radiation data ***
c
c _dp implicit double precision (a-h,o-z)
      logical lpri
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common /padel/ ga(16,6),gb(16,6),npad(16)
      common /ir5/ rt1(2),wg1(2),ng1
      common /ir6/ at(10),bt(6),o1h,o2h
      common /ir7/ xp(2,6,4),tstand,nint
c
c *** for infra routine, morcrette,these,1984 ***
c
      lpri=.false.
c mVax open(unit=14,file='bcmir0.dat',status='old')
       open(unit=14,file='bcmir0',status='old')
       rewind 14
      read (14,100) ng1,nint,tstand,o1h,o2h
      read (14,110) (rt1(i),i=1,2)
      read (14,120) (wg1(i),i=1,2)
      read (14,130) (npad(i),i=1,16)
      if (lpri) write(iwr6,100) ng1,nint,tstand,o1h,o2h
      if (lpri) write(iwr6,110) (rt1(i),i=1,2)
      if (lpri) write(iwr6,120) (wg1(i),i=1,2)
      if (lpri) write(iwr6,130) (npad(i),i=1,16)
 100  format (1x,i1,1x,i1,1x,f5.1,1x,f6.1,1x,f6.4)
 110  format (2(1x,f12.9))
 120  format (2(1x,f3.1))
 130  format (16(1x,i1))
      read (14,140) (at(i),i=1,6)
      read (14,140) (at(i),i=7,10)
      if (lpri) write(iwr6,140) (at(i),i=1,6)
      if (lpri) write(iwr6,140) (at(i),i=7,10)
 140  format (6(e10.3))
      do 10 k=1,4
      do 10 i=1,2
      read (14,150) (xp(i,j,k),j=1,3)
      read (14,150) (xp(i,j,k),j=4,6)
      if (lpri) write(iwr6,151) (xp(i,j,k),j=1,3)
      if (lpri) write(iwr6,151) (xp(i,j,k),j=4,6)
   10 continue
 150  format (3(e14.7))
 151  format (3(1x,e14.7))
      do 20 i=1,6
      read (14,160) (ga(i,j),j=1,3)
      read (14,160) (ga(i,j),j=4,6)
      if (lpri) write(iwr6,160) (ga(i,j),j=1,3)
      if (lpri) write(iwr6,160) (ga(i,j),j=4,6)
  20  continue
      do 30 i=7,8
      read (14,170) (ga(i,j),j=1,4)
      if (lpri) write(iwr6,170) (ga(i,j),j=1,4)
  30  continue
      do 40 i=9,12
      read (14,160) (ga(i,j),j=1,3)
      read (14,160) (ga(i,j),j=4,6)
      if (lpri) write(iwr6,160) (ga(i,j),j=1,3)
      if (lpri) write(iwr6,160) (ga(i,j),j=4,6)
  40  continue
      do 50 i=13,16
      read (14,170) (ga(i,j),j=1,4)
      if (lpri) write(iwr6,170) (ga(i,j),j=1,4)
  50  continue
      do 60 i=1,6
      read (14,160) (gb(i,j),j=1,3)
      read (14,160) (gb(i,j),j=4,6)
      if (lpri) write(iwr6,160) (gb(i,j),j=1,3)
      if (lpri) write(iwr6,160) (gb(i,j),j=4,6)
  60  continue
      do 70 i=7,8
      read (14,170) (gb(i,j),j=1,4)
      if (lpri) write(iwr6,170) (gb(i,j),j=1,4)
  70  continue
      do 80 i=9,12
      read (14,160) (gb(i,j),j=1,3)
      read (14,160) (gb(i,j),j=4,6)
      if (lpri) write(iwr6,160) (gb(i,j),j=1,3)
      if (lpri) write(iwr6,160) (gb(i,j),j=4,6)
  80  continue
      do 90 i=13,16
      read (14,170) (gb(i,j),j=1,4)
      if (lpri) write(iwr6,170) (gb(i,j),j=1,4)
  90  continue
 160  format (3(e13.7))
 170  format (4(e13.7))
      close(unit=14)
      write(iwr6,2000)
2000  format (///,'    end of the ir radiative data lecture ',//)
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  infra===
c cray .endif
      subroutine infra (dt0,nl)
c
c +++ bcmir3.for +++
c
      parameter(nir=50,nir3=151)
c _dp implicit double precision (a-h,o-z)
      logical gauss
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c *** input ***
c
      common/ir0/stfir0
      common/gz1/ih1,ih2,ih3,ihc1,ihc2,ihc3,ico2,io3
      common/gz2/ppmgz2
c
      common/ra1/ttra1(nir),hhra1(nir),ozra1(nir),cldra1(nir),
     .           pmbra1(nir),zkmra1(nir)
c
      common /padel/ ga(16,6),gb(16,6),npad(16)
      common /ir5/ rt1(2),wg1(2),ng1
      common /ir6/ at(10),bt(6),o1h,o2h
      common /ir7/ xp(2,6,4),tstand,nint
c -------
        common/di1/it1,it2,it3,ic1,ic2,ic3,inu
c -------
c
c --- variables locales input
c
      dimension z(nir),b(4,nir),bint(nir),f(2,nir),tave(nir)
     1 ,indcl(nir),cntrb(nir,nir),fup(nir,nir),fdn(nir,nir),
     2 dz(nir),bsol(4),u(8,nir3),v(10,nir3),w(4,nir3),xt(nir3),btop(4)
     3 ,uu(8),tt(8)
      equivalence(uu(1),uh),(uu(2),uc),(uu(3),uo),(uu(4),uw),(uu(5),uho)
     1,(uu(6),uf),(uu(7),uw63),(uu(8),uw15)
      equivalence(tt(1),th),(tt(2),tc),(tt(3),to),(tt(4),tw),(tt(5),t15)
     1,(tt(6),tf),(tt(7),tw63),(tt(8),tw15)
c
c *** output ***
c
       common/ir1/fimir1(nir),fidir1(nir),finir1(nir)
c
c ----------------------------------------------------------------------
      gauss = .true .
c ----------------------------------------------------------------------
c
      iw2 = 2
c ... parametres d'ecriture
c
      em0 = 1.
c ... emissivite surface fixee provisoirement a 1.
c
      nc=nl-1
      ng1p1=ng1+1
      nglp1=nc*ng1p1+1
      ngl=nglp1-1
c ccc  if (iw2.eq.1) write(iwr6,1023) (cldra1(i),i=1,nc)
1023  format (/,10x,'cld(i)',/,(5x,e12.5))
      do 3001 i=1,nglp1
      do 3001 j=1,nint
      w(j,i)=0.0
3001  continue
c-- ozone modifie
      do 104 i=1,nl
      bint(i)=0.00
      do 104 j=1,nl
      fup(j,i)=0.0
      fdn(j,i)=0.0
      finir1(i)=0.0
104   continue
      bsolin=0.0
c-- ppm en ppm   mco2= 44 divise par mair= 28.97
      co2= 1.519e-06*ppmgz2
      mx=0
      ncloud=0
      do 106 i=1,nc
      tave(i)=(ttra1(i)+ttra1(i+1))/2.
      if(cldra1(i).eq.0.  )go to 106
      ncloud=ncloud+1
      indcl(ncloud)=i
106   continue
      if(ncloud.ne.0)mx=indcl(ncloud)
      psol=pmbra1(1)
      do 2 i=1,nl
      z(i)=-alog(pmbra1(i)/psol)
      q=hhra1(i)
      w(1,i)=q
2     continue
      if (iw2.eq.1) write(iwr6,1000) (w(1,i),ozra1(i),ttra1(i),i=1,nl)
c ccc  if (iw2.eq.1) write(iwr6,1001) (tave(i),i=1,nc)
1001  format (/,10x,'tave',/,(5x,e14.6))
1000  format (/,10x,'  w  ',10x,'  o  ',10x,'  t  ',/,
     *(4x,e12.4,5x,e12.4,5x,e12.4))
      do 3 i=1,nl
      ip1=i+1
      j=(i-1)*ng1p1+1
      u(1,j)=ozra1(i)
      u(7,j)=pmbra1(i)/psol
      if(i.eq.nl) go to 3
      zp=(z(i)+z(ip1))/2.0
      zm=(z(ip1)-z(i))/2.0
      dz(i)=zm
      zm2=2.0*zm
      do 29 ig1=1,ng1
      j=j+1
      xz=rt1(ig1)*zm+zp
      u(7,j)=exp(-xz)
      xt(j)=(ttra1(ip1)-ttra1(i))/zm2
      xr=(xz-z(i))/zm2
      u(1,j)=ozra1(i)+xr*(ozra1(ip1)-ozra1(i))
29    continue
 3    continue
c-- 101325/9.80665/10
      coef=1033.2275
      do 49 i=1,ngl
      ip1=i+1
c-- 22400/48
      xo=(u(1,i)-u(1,ip1))
      xxo=xo/466.968
      upm=(u(7,i)+u(7,ip1))/2.0
      u(3,i)=(u(7,i)-u(7,ip1))*coef
      u(5,i)=upm
      u(1,i)=xxo
      u(2,i)=xxo*upm
      xxc=co2*u(3,i)
      u(4,i)=xxc*upm
49    continue
      do 555 i=1,nc
      ip1=i+1
      xq=(w(1,i)+w(1,ip1))/2.0
      j=(i-1)*ng1p1+1
      jpn=j+ng1
      do 554 k=j,jpn
      xxu=xq*u(3,k)
      u6=xxu*u(5,k)
c-- 28.9647 / 18.0153
      fppw= 1.6078 *xq/(1.0+0.608*xq)
      u(6,k)=u6
      u(7,k)=u6*fppw
c ccc  if (iw2.eq.1) write(iwr6,888)u(2,k),u(4,k),u(6,k),u(7,k)
554   continue
555   continue
      do 11 i=1,10
      v(i,nglp1)=0.0
11    continue
c-- introduction de l'effet de temperature sur les quantites d'absorbant
      do 121 i=1,nc
      j=(i-1)*ng1p1+1
      jpn=j+ng1
      l=nl-i
      tavic=tave(l)
      do 120 k=j,jpn
      ic=nglp1-k
      icp1=ic+1
      factt=exp(6.08*(296.0/tavic   -1.0))
      fact63=1.
      fact15=factt
      v(7,ic)=v(7,icp1)+u(7,ic)*factt*1.66
      v(2,ic)=v(2,icp1)+u(7,ic)*fact63*1.66
      v(3,ic)=v(3,icp1)+u(7,ic)*fact15*1.66
      tx=tavic-250.0
      tx2=tx*tx
      psih=at(1)*tx+at(2)*tx2
      epsih=exp(psih)
      ps15=at(3)*tx+at(4)*tx2
      eps15=exp(ps15)
      psic=at(5)*tx+at(6)*tx2
      epsic=exp(psic)
      pshw=at(7)*tx+at(8)*tx2
      epshw=exp(pshw)
      phio=at(9)*tx+at(10)*tx2
      ephio=exp(phio)
      v(1,ic)=v(1,icp1)+u(1,ic)*ephio*1.66
      v(4,ic)=v(4,icp1)+u(4,ic)*epsic*1.66
      v(6,ic)=v(6,icp1)+u(6,ic)*epsih*1.66
      v(9,ic)=v(9,icp1)+u(6,ic)*eps15*1.66
      v(10,ic)=v(10,icp1)+u(6,ic)*epshw*1.66
c ccc  if (iw2.ne.1) goto 1020
c ccc  write(iwr6,889)ic,v(1,ic),v(4,ic),v(6,ic),v(7,ic),v(9,ic),v(10,ic)
c ccc1020  continue
120   continue
121   continue
888   format (10e12.6)
c-- calcul des fonctions de planck b(t) et db/dt (polynome degre 5 )
      do 89 int=1,nint
      do 8 i=1,nl
      ti=(ttra1(i)-tstand)/tstand
      res=xp(1,6,int)
      do 86 ll=1,5
       l=6-ll
      res=res*ti+xp(1,l,int)
86    continue
      b(int,i)=res
      bint(i)=bint(i)+res
8     continue
      bsol(int)=b(int,1)
      if (dt0.eq.0.0) go to 9
      ti=(ttra1(1)+dt0-tstand)/tstand
      res=xp(1,6,int)
      do 87 ll=1,5
      l=6-ll
      res=res*ti+xp(1,l,int)
87    continue
      bsol(int)=res
9     btop(int)=b(int,nl)
      bsolin=bsolin+bsol(int)
      do 10 i=1,nc
      j=(i-1)*ng1p1+1
      ti=(tave(i)-tstand)/tstand
      res=xp(2,6,int)
      do 88 ll=1,5
      l=6-ll
      res=res*ti+xp(2,l,int)
88    continue
      do 10 ig=1,ng1
      jpig=j+ig
      w(int,jpig)=res
10    continue
c cc  if (iw2.eq.1) write(iwr6,888)(b(int,l),l=1,nl)
c cc  if (iw2.eq.1) write(iwr6,888)(w(int,l),l=1,nglp1)
89    continue
c cc  if (iw2.eq.1) write(iwr6,888)(bsol(int),int=1,nint)
c cc  if (iw2.eq.1) write(iwr6,888)(btop(int),int=1,nint)
      ind1=0
      ind2=1
c -- integration verticale
c -- se fait par trapeze sauf pour les couches adjacentes au niveau de
c    calcul ou on utilise une quadrature de gauss a ng1 points
c ini  do 21 i=1,nl
c new------------
c     do 21 ii=1,4
c     i=1
c     if (ii.eq.2) i=it1
c     if (ii.eq.3) i=it2
c     if (ii.eq.4) i=it3
      do 21 ii=1,2
      i=1
      if (ii.eq.2) i=it3
c new------------
      in=(i-1)*ng1p1+1
      ip1=i+1
      xmont=0.0
      xdesc=0.0
      xadjd=0.0
      xadjm=0.0
c -- flux descendants
      if(i.eq.nl) go to 15
      do 14 j=i,nc
      jx=(j-1)*ng1p1+1
      xdg=0.
c ----------------------------------------------------------------------
      if (.not.gauss) goto 131
c ----------------------------------------------------------------------
      if (j.ne.i) go to 131
      do 13 ig=1,ng1
      id=jx+ig
      uh=v(6,in)-v(6,id)
      uho=v(9,in)-v(9,id)
      uf=v(10,in)-v(10,id)
      uc=v(4,in)-v(4,id)
      uo=v(1,in)-v(1,id)
      uw=v(7,in)-v(7,id)
      uw63=v(2,in)-v(2,id)
      uw15=v(3,in)-v(3,id)
      call hornel( uu, tt, ind1)
c-- 140 / 490
      wtr=w(1,id)*th*tw63+w(2,id)*t15*tc*tw15+tw*tf*(w(3,id)+w(4,id)*to)
      xdg=xdg+wtr*xt(id)*wg1(ig)
13    continue
      dzxdg=dz(j)*xdg
      xadjd=dzxdg
      go to 139
131   if (j.ne.ip1) go to 133
      id1=jx
      uh=v(6,in)-v(6,id1)
      uho=v(9,in)-v(9,id1)
      uf=v(10,in)-v(10,id1)
      uc=v(4,in)-v(4,id1)
      uo=v(1,in)-v(1,id1)
      uw=v(7,in)-v(7,id1)
      uw63=v(2,in)-v(2,id1)
      uw15=v(3,in)-v(3,id1)
      call hornel( uu, tt, ind1)
      th1=th
      tc1=tc
      to1=to
      tw1=tw
      t151=t15
      tf1=tf
      tw631=tw63
      tw151=tw15
133   id2=jx+ng1p1
      uh=v(6,in)-v(6,id2)
      uho=v(9,in)-v(9,id2)
      uf=v(10,in)-v(10,id2)
      uc=v(4,in)-v(4,id2)
      uo=v(1,in)-v(1,id2)
      uw=v(7,in)-v(7,id2)
      uw63=v(2,in)-v(2,id2)
      uw15=v(3,in)-v(3,id2)
      call hornel( uu, tt, ind1)
      th2=th
      tc2=tc
      to2=to
      tw2=tw
      t152=t15
      tf2=tf
      tw632=tw63
      tw152=tw15
      th=(th1+th2)*0.5
      tc=(tc1+tc2)*0.5
      to=(to1+to2)*0.5
      tw=(tw1+tw2)*0.5
      t15=(t151+t152)*0.5
      tf=(tf1+tf2)*0.5
      tw63=(tw631+tw632)*0.5
      tw15=(tw151+tw152)*0.5
      k2=id2-1
c ccc  if (iw2.ne.1) goto 1021
c ccc  write(iwr6,890)i,th,tc,to,tw,t15,tf,tw63,tw15,(w(l,k2),l=1,3),
c ccc *xt(k2),dz(j)
c ccc1021  continue
890   format(1x,i2,8f7.4,5e12.6)
      ww=w(1,k2)*th*tw63+w(2,k2)*tc*t15*tw15+tw*tf*(w(3,k2)+w(4,k2)*to)
      xdg=ww*xt(k2)*2.0
      dzxdg=dz(j)*xdg
      xdesc=xdesc+dzxdg
      th1=th2
      tc1=tc2
      to1=to2
      tw1=tw2
      t151=t152
      tf1=tf2
      tw631=tw632
      tw151=tw152
139   jp1=j+1
      cntrb(i,jp1)=dzxdg
14    continue
c-- flux montants
15    if(i.eq.1) go to 18
      im1=i-1
      im2=i-2
      do 17 l=1,im1
      j=i-l
      jx=(j-1)*ng1p1+1
      xmg=0.
c ----------------------------------------------------------------------
      if (.not.gauss) goto 161
c ----------------------------------------------------------------------
      if(j.ne.im1) go to 161
      do 16 ig=1,ng1
      im=jx+ig
      uh=v(6,im)-v(6,in)
      uho=v(9,im)-v(9,in)
      uf=v(10,im)-v(10,in)
      uc=v(4,im)-v(4,in)
      uo=v(1,im)-v(1,in)
      uw=v(7,im)-v(7,in)
      uw63=v(2,im)-v(2,in)
      uw15=v(3,im)-v(3,in)
      call hornel( uu, tt, ind1)
      wtr=w(1,im)*th*tw63+w(2,im)*t15*tc*tw15+tw*tf*(w(3,im)+w(4,im)*to)
      xmg=xmg+wtr*xt(im)*wg1(ig)
16    continue
      dzxmg=dz(j)*xmg
      xadjm=dzxmg
      go to 169
161   if(j.ne.im2) go to 163
      iu1=jx+ng1p1
      uh=v(6,iu1)-v(6,in)
      uho=v(9,iu1)-v(9,in)
      uf=v(10,iu1)-v(10,in)
      uc=v(4,iu1)-v(4,in)
      uo=v(1,iu1)-v(1,in)
      uw=v(7,iu1)-v(7,in)
      uw63=v(2,iu1)-v(2,in)
      uw15=v(3,iu1)-v(3,in)
      call hornel( uu, tt, ind1)
      th1=th
      tc1=tc
      to1=to
      tw1=tw
      t151=t15
      tf1=tf
      tw631=tw63
      tw151=tw15
163   iu2=jx
      uh=v(6,iu2)-v(6,in)
      uho=v(9,iu2)-v(9,in)
      uf=v(10,iu2)-v(10,in)
      uc=v(4,iu2)-v(4,in)
      uo=v(1,iu2)-v(1,in)
      uw=v(7,iu2)-v(7,in)
      uw63=v(2,iu2)-v(2,in)
      uw15=v(3,iu2)-v(3,in)
      call hornel( uu, tt, ind1)
      th2=th
      tc2=tc
      to2=to
      tw2=tw
      t152=t15
      tf2=tf
      tw632=tw63
      tw152=tw15
      th=(th1+th2)*0.5
      tc=(tc1+tc2)*0.5
      to=(to1+to2)*0.5
      tw=(tw1+tw2)*0.5
      t15=(t151+t152)*0.5
      tf=(tf1+tf2)*0.5
      tw63=(tw631+tw632)*0.5
      tw15=(tw151+tw152)*0.5
      k2=jx+1
      ww=w(1,k2)*th*tw63+w(2,k2)*tc*t15*tw15+tw*tf*(w(3,k2)+w(4,k2)*to)
      xmg=ww*xt(k2)*2.0
      dzxmg=dz(j)*xmg
      xmont=xmont+dzxmg
      th1=th2
      tc1=tc2
      to1=to2
      tw1=tw2
      t151=t152
      tf1=tf2
      tw631=tw632
      tw151=tw152
169   cntrb(i,j)=dzxmg
17    continue
18    cnsol=0.0
c **  contribution du sommet de l'atmosphere  ***
      cntop=btop(1)+btop(2)+btop(3)+btop(4)
      if(i.eq.nl) go to 19
      uh=v(6,in)
      uho=v(9,in)
      uf=v(10,in)
      uc=v(4,in)
      uo=v(1,in)
      uw=v(7,in)
      uw63=v(2,in)
      uw15=v(3,in)
      call hornel( uu, tt, ind2)
      cntop=btop(1)*th*tw63+btop(2)*tc*t15*tw15+tw*tf*(btop(3)+btop(4)*
     *to)
19    fd=cntop-bint(i)-xdesc-xadjd
      if(i.eq.1) fd1=-fd
c **   contribution de la surface   ***
      bgnd=bsolin*em0+(1.0-em0)*fd1-bint(1)
      if(bgnd.eq.0.) go to 20
      uh=v(6,1)-v(6,in)
      uho=v(9,1)-v(9,in)
      uf=v(10,1)-v(10,in)
      uc=v(4,1)-v(4,in)
      uo=v(1,1)-v(1,in)
      uw=v(7,1)-v(7,in)
      uw63=v(2,1)-v(2,in)
      uw15=v(3,1)-v(3,in)
      call hornel( uu, tt, ind2)
      cnsol=bsol(1)*th*tw63+bsol(2)*tc*t15*tw15+tw*tf*(bsol(3)+bsol(4)*
     *to)
      cnsol=cnsol*bgnd/bsolin
20    continue
      fm=cnsol+bint(i)-xmont-xadjm
      fup(1,i)=fm
      fdn(1,i)=fd
c ccc  if (iw2.ne.1) goto 1022
c ccc  write(iwr6,889)i,cntop,xadjd,xdesc,fd,bint(i),fm,cnsol,xadjm,xmont
c ccc1022  continue
21    continue
      if(ncloud.ne.0) go to 220
      do 219 i=1,nl
      f(1,i)=fup(1,i)
      f(2,i)=fdn(1,i)
      fimir1(i)=f(1,i)
      fidir1(i)=-f(2,i)
      finir1(i)=f(1,i)+f(2,i)
889   format(i3,10e12.6)
219   continue
      go to 221
220   call cloudy (bsolin,cntrb,f,fup,fdn,bint,em0,indcl,mx
     *,nl,ncloud)
      do 3000 i=1,nl
      fimir1(i)=f(1,i)
      fidir1(i)=-f(2,i)
      finir1(i)=f(1,i)+f(2,i)
3000  continue
221   do 222 i=1,nc
      ip1=i+1
222   continue
      outgir=f(1,nl)
      firsol=-f(2,1)
c ccc  if (iw2.eq.1) write(iwr6,888) (fimir1(i),i=1,nl)
c ccc  if (iw2.eq.1) write(iwr6,888) (fidir1(i),i=1,nl)
c ccc  if (iw2.eq.1) write(iwr6,888) (finir1(i),i=1,nl)
      return
      end
c **********************************************************************
c cray .if [ $concord .eq. 1 ]
c cray nzyx  hornel==
c cray .endif
      subroutine hornel( uu, tt, ind)
c-- calcul des transmissions a partir de coefficients de pade
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/gz1/ih1,ih2,ih3,ihc1,ihc2,ihc3,ico2,io3
c
      common /padel/ ga(16,6),gb(16,6),npad(16)
      dimension uu(8),tt(8)
      do 9 ia=1,8
      tt(ia)=1.0
9     continue
      do 10 ia=1,8
      if ((ia.eq.1).and.(ih1 .ne.1)) go to 10
      if ((ia.eq.2).and.(ico2.ne.1)) go to 10
      if ((ia.eq.3).and.(io3 .ne.1)) go to 10
      if ((ia.eq.4).and.(ihc3.ne.1)) go to 10
      if ((ia.eq.5).and.(ih2 .ne.1)) go to 10
      if ((ia.eq.6).and.(ih3 .ne.1)) go to 10
      if ((ia.eq.7).and.(ihc1.ne.1)) go to 10
      if ((ia.eq.8).and.(ihc2.ne.1)) go to 10
c _rad1
c
      l=2*ia-ind
      np=npad(l)
      npp1=np+1
      xn=0.0
      xd=1.0
      z=uu(ia)
      do 1 j=1,np
      i=npp1-j
      xn=xn*z+ga(l,i)
      xd=xd*z+gb(l,i)
1     continue
      r=xn/xd
      tt(ia)=r
10    continue
c     write(iwr6,888) (uu(l),tt(l),l=1,8)
888   format(1x,8(1x,e10.4,f9.5))
      return
      end
c **********************************************************************
c cray .if [ $concord .eq. 1 ]
c cray nzxy cloud===
c cray .endif
      subroutine cloudy(bsol,cntrb,f,fup,fdn,bint,em0,indcl,mx
     *,nl,ncloud)
      parameter(nir=50)
c _dp implicit double precision (a-h,o-z)
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/ra1/ttra1(nir),hhra1(nir),ozra1(nir),cldra1(nir),
     .           pmbra1(nir),zkmra1(nir)
cnew----------
        common/di1/it1,it2,it3,ic1,ic2,ic3,inu
cnew----------
c
      dimension indcl(nir),cntrb(nir,nir),fup(nir,nir),fdn(nir,nir),
     1 bint(nir),f(2,nir)
      do 110 i=1,nl
      do 110 j=2,nl
      fup(j,i)=fup(1,i)
      fdn(j,i)=fdn(1,i)
110   continue
      do 190 k=1,ncloud
      icloud=indcl(k)
      icp1=icloud+1
      fcldup=bint(icp1)
      fclddn=bint(icloud)
      fm=0.0
      fd=0.0
c ini do 180 i=1,nl
cnew--------
c     do 180 ii=1,4
c     i=1
c     if (ii.eq.2) i=it1
c     if (ii.eq.3) i=it2
c     if (ii.eq.4) i=it3
      do 180 ii=1,2
      i=1
      if (ii.eq.2) i=it3
cnew--------
      im=i-1
      ip=i+1
      if(i.le.icloud) go to 150
      fm=fcldup
      if(i.eq.icp1) go to 141
      do 140 j=icp1,im
      fm=fm+cntrb(i,j)
140   continue
141   fdn(icp1,i)=fdn(1,i)
      fup(icp1,i)=bint(i)+fcldup-fm
      go to 180
150   fd=fclddn
      if(i.eq.icloud) go to 170
      do 160 j=ip,icloud
      fd=fd+cntrb(i,j)
160   continue
170   fdn(icp1,i)=-bint(i)+fclddn-fd
      fup(icp1,i)=fup(1,i)
180   continue
      do 181 i=1,nl
      f(1,i)=fup(icp1,i)
      f(2,i)=fdn(icp1,i)
181   continue
      if(ncloud.gt.1) go to 190
      if (cldra1(icloud).ne.1.0) go to 190
      return
190   continue
      mx1=mx+1
      fm=0.0
      fd=0.0
      do 291 i=1,nl
      ip=i+1
      cloud=1.0
      if(i.gt.mx) go to 280
      ci=cldra1(i)
      fd=fdn(ip,i)*ci
      do 260 j=i,mx
      j1=j+1
      j2=j+2
      cloud=cloud*(1.00-cldra1(j))
      if (j.eq.mx) go to 270
      ccld=cloud*cldra1(j1)
      fd=fd+fdn(j2,i)*ccld
260   continue
270   fd=fd+fdn(1,i)*cloud
      f(2,i)=fd
      go to 290
280   f(2,i)=fdn(1,i)
290   continue
291   continue
      f(1,1)=em0*bsol-(1.0-em0)*f(2,1)
      do 250 i=2,nl
      im=i-1
      cloud=1.0
      if(i.gt.mx1) go to 220
      fm=fup(i,i)*cldra1(im)
      i2=im
      i1=i
      go to 221
220   fm=fup(mx1,i)*cldra1(mx)
      i2=mx
      i1=mx1
221   do 222 j=1,i2
      ij=i1-j
      cloud=cloud*(1.00-cldra1(ij))
      if(ij.eq.1) go to 223
      ccld=cloud*cldra1(ij-1)
      fm=fm+fup(ij,i)*ccld
222   continue
223   fm=fm+fup(1,i)*cloud
      f(1,i)=fm
250   continue
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  infrs===
c cray .endif
      subroutine infras (dt0,jlat,iilw)
c
c +++ bcmir3s.for +++
c
c.......   parametrisation de thompson-warren pour flux au sommet  .....
c                             (jas,1982,39,2667-2680)
c                          de berlyand-berlyand pour flux net en surface
c                             si ta.ge.282.5
c                             (voir fung et al.,1984,rev.geoph.space
c                              phys.,22,177-193)
c                          de maykut      pour flux descendant
c                             si ts.lt.277.5
c.......................................................................
c
      parameter(np=18,nilw=7,nir=50)
c _dp implicit double precision (a-h,o-z)
c
c    **   input   **
c
       common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
       common/ir0/stfir0
       common/ra1/ttra1(nir),hhra1(nir),ozra1(nir),cldra1(nir),
     .            pmbra1(nir),zkmra1(nir)
        common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .             ajosu0(np),ajsu0(np,nilw)
        common/di1/it1,it2,it3,ic1,ic2,ic3,inu
        common/so5/tauso5(nir),omeso5(nir),cgso5(nir),c1iso5(nir)
        common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
        common/ev4/rhev4(np),rhev5(np,12)
        common/ir8/facir8(np)
c
c     **  output  **
c
       common/ir1/fimir1(nir),fidir1(nir),finir1(nir)
c
        dimension bc(3,4),dc(6),ac(4)
        data bc/243.414,-34.7968,10.2790,
     .          2.60065,-1.62064,0.634856,
     .          4.40272e-03,-2.26092e-02,1.12265e-02,
     .         -2.05237e-05,-9.67000e-05,5.62925e-05/
        data dc/-3.1,-0.4146,4.084e-03,-4.44e-03,-80.0,-0.40/
c
        em0  =1.
        cnebu=c1iso5(1)
        ts   =t4ssu0(jlat,iilw)-273.15
        ttn  =ttra1(it2)-273.15
        tscn =ts-ttn
        if (tscn.gt.60.) tscn=60.
c.......................................................................
c
c 1. flux montant au sommet
c
        do 1 i=1,4
        ac(i)=bc(1,i) + bc(2,i)*rhev4(jlat) +
     .                  bc(3,i)*rhev4(jlat)*rhev4(jlat)
    1   continue
        firm1=ac(1)+ac(2)*ts+ac(3)*ts*ts+ac(4)*ts*ts*ts
        firm2=ac(1)+ac(2)*ttn+ac(3)*ttn*ttn+ac(4)*ttn*ttn*ttn
        f1   =dc(1)+dc(2)*tscn+dc(3)*tscn*tscn
        f2   =dc(4)*tscn*(tscn+dc(5))*(rhev4(jlat)+dc(6))
        firm3=f1+f2
c
        firm =firm1 - cnebu*(firm1-firm2+firm3)
c
c 2. flux infrarouge net en surface
c
        ts    = t4ssu0(jlat,iilw)
        tair  = ttra1(1)
        firms = stfir0 * em0 * ts   * ts   * ts   * ts
        firma = stfir0 * em0 * tair * tair * tair * tair
c
        es1 =10.553-2667./tair
        if (tair.ge.273.) es1=9.4051-2354./tair
        es  =10.**es1
c ---   es tension de vapeur saturante en surface , en mb
        er  = es * rhzev1(jlat,iilw)
        facn=1. - facir8(jlat)*cnebu*cnebu
c
c berlyand-berlyand
        firdn1=firma * (0.39 - 0.05*sqrt(er))*facn+4.*stfir0*em0*
     .                  dt0*tair*tair*tair
c
c maykut
        fird  =firma * (0.7855+0.000312*(cnebu*10.)**2.75)
        firdn2=firms - fird
c ccc
c       firdn2=firdn1
c ccc
c
      if (tair.ge.282.5) then
       firdn  = firdn1
      else
       if (tair.ge.277.5) then
        firdn  = ( (tair - 277.5) *firdn1
     .           + (282.5 - tair) *firdn2) / 5.
       else
        firdn  = firdn2
       end if
      end if
c
c ---------------------------------------------------------------------
c
c --- output bcm
c
        finir1(it3)=firm
        finir1(1)  =firdn
        fimir1(it3)=firm
        fidir1(1)  =firms-firdn
c
        return
        end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  c3rad===
c cray .endif
      subroutine c3rad(j,iilw)
c
c +++ bcmir6.for +++
c
      parameter(np=18,nilw=7)
c _dp implicit double precision (a-h,o-z)
c
       common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
       common/cl1/clcl1(np),clhcl1,ztcl1(np),zbcl1(np),
     .                             ptcl1(np),pbcl1(np)
       common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
       common/ir0/stfir0
       common/ir2/fupir2(np),fdnir2(np),fnir2(np),divir2(np)
       common/ir4/c3ir4(np,nilw)
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
c
c ----------------------------------------------------------------------
c
c --- inversion de la relation de budyko du flux i.r. net en surface
c     (cfr. saltzman et ashe, 1976, tellus, p.308)
c
       data c1/.254/,c2/4.95e-3/
c
       es = q0ev1(j)*1.0e3/(0.622*rhzev1(j,iilw))
       pp = c1-c2*es
       fm = stfir0*t4ssu0(j,iilw)**4
       de = fm*pp
       if (de.ge.0) de = amax1(de,1.e-10)
       t1 = fdnir2(j)-fm
       t2 = 1. + t1/de
       c3ir4(j,iilw)=t2/(clcl1(j)**2)
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  aja=====
c cray .endif
      subroutine iniaja
c
c +++ bcmj1.for +++
c
      parameter(np=18,npp=19,nilw=7)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical       lgm18k,taiga,otter,summer,co2var
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/varl02/lgm18k,taiga,otter,summer,co2var
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
      common/p01/szp01(np,nilw),spp01(np,nilw)
      common/p02/szp02(np,nilw)
      dimension secz0(46,4),azc1(np,4)
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .          ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
      dimension ajo0(9),ajwa0(9),ajla0(9),y(np)
      common/ca1/ajca1(np,nilw)
      dimension ajc0(46,4),ajc1(np,4),x(46),z(npp)
c
      dimension ajolgm(18),ajglgm(18),ajnlgm(18),ajflgm(18)
      dimension zjclgm(18),zjglgm(18),zjnlgm(18),zjflgm(18)
c ... Pour simulation du dernier maximum glaciaire
c
       data ajo0  /.78,.72,.62,.56,.475,.43 ,.30 ,.71 ,.935/
c ... fraction oceanique totale
c
       data ajwa0 /.78,.72,.62,.56,.475,.395,.169,.133,.004/
c ... fraction oceanique libre de glaces - moyenne annuelle
c
       data ajla0 /.22,.28,.38,.41,.341,.299,.298,.062,.006/
c ... fraction continentale - moyenne annuelle
c
       data ajc0/
     . 0., 0., 0., 0.,31.,25.,27.,22.,19.,17.,15.,11.,13., 6., 4., 0.,
     .30*0.,
c ... groenland
c
     . 0., 0., 0., 0., 8.,11., 7.,20.,13.,22.,40.,46.,49.,46.,39.,33.,
     .33.,32.,36.,36.,31.,32.,31.,28.,27.,25.,24.,22.,20.,18.,13., 9.,
     . 7., 4., 5., 6., 7., 6., 4., 3., 8.,10.,12.,13.,14.,16.,
c ... amerique du nord
c
     . 0., 0., 0., 0., 0., 9., 7.,17.,21.,41.,64.,84.,87.,86.,83.,77.,
     .33.,64.,74.,74.,74.,70.,63.,62.,55.,56.,53.,47.,46.,43.,43.,42.,
     .39.,37.,31.,22.,18.,15., 9., 7., 3., 5*0.,
c ... eurasie
c
     .16*0.,
     . 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 6., 9.,14.,22.,24.,
     .25.,25.,27.,27.,27.,28.,32.,30.,33.,31.,27.,19.,14.,17./
c ... afrique
c ... source : donnees climap traitees par f.mercier
c
       data secz0/
     .   0.,   0.,   0.,   0.,1397.,1584.,1881.,2146.,2461.,2621.,
     .2484.,2288.,1917.,2206.,1661., 500.,
     .30*0.,
c ... groenland : donnees de Radok et al., 1982, traitees par I.Marsiat
c
c    .   0.,   0.,   0.,   0., 586.,1200.,1522.,1860.,1897.,1936.,
c    .2003.,1875.,1066.,1608.,1235.,   0.,
c    .30*0.,
c ... groenland : donnees Groenland fournies par I.Marsiat (mai 1987)
c
     .   0.,   0.,   0.,   0., 150., 262., 107., 129., 154., 139.,
     . 190., 335., 307., 452., 656., 475., 486., 482., 461., 479.,
     . 439., 412., 564., 721., 842., 925., 947., 805., 656., 540.,
     . 548., 748.,1003.,1487., 952.,1150., 493., 442., 555., 177.,
     . 361., 359., 528., 476., 369., 360.,
c ... amerique du nord
c
     .   0.,   0.,   0.,   0.,   0., 104.,  81., 121., 106., 123.,
     . 187., 275., 312., 320., 363., 288., 486., 387., 359., 482.,
     . 546., 609., 608., 776.,1143.,1140.,1506.,2136.,2192.,1985.,
     .1710.,1250., 658., 551., 508., 587., 514., 461., 396., 336.,
     . 157.,   0.,   0.,   0.,   0.,   0.,
c ... eurasie
c
     .20*0.,
     .   0.,   0.,   0.,   0.,   0.,   0.,   0., 540., 647., 574.,
     . 359., 363., 404., 502., 469., 453., 405., 418., 486., 519.,
     . 571., 582., 565., 586., 369., 659./
c ... afrique
c ... source : donnees climap traitees par f.mercier
c
      data zjnlgm
     ./   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 444.4,1418.2,
     . 2065.6,2505.9,2592.0,2482.2,2231.8,1690.6,1553.1,1185.8,   0.0/
      data ajnlgm
     ./0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0042,0.0560,
     . 0.1417,0.2120,0.2412,0.2153,0.1773,0.1199,0.1199,0.0454,0.0000/
c ... Laurentide -18000 (CLIMAP)
c
      data zjflgm
     ./   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     .    0.0,1244.5,1856.0,2185.8,1836.7,1809.5,1576.1, 969.1,   0.0/
      data ajflgm
     ./0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     . 0.0000,0.0259,0.1023,0.1106,0.1954,0.2333,0.2523,0.0620,0.0000/
c ... Fennoscandienne -18000 (CLIMAP)
c
      data zjglgm
     ./   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     .    0.0,   0.0,   0.0,1133.8,1736.2,2055.8,1877.5,1083.7,  75.0/
      data ajglgm
     ./0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     . 0.0000,0.0000,0.0000,0.0380,0.0903,0.1083,0.1532,0.0921,0.0005/
c ... Groenlandaise -18000 (CLIMAP) 
c
      data zjclgm
     ./ 615.8, 605.3, 587.4, 601.1, 678.1, 951.7,1452.8,1450.5, 986.3,
     .  735.6, 571.0, 470.6, 501.1, 474.6, 179.1,  81.3, 487.5,   0.0/
c ... Continent moyen (sans les calottes) -18000 (CLIMAP)
c
      data ajolgm
     ./0.7606,0.7296,0.7421,0.6968,0.6273,0.5699,0.5458,0.5426,0.4745,
     . 0.4032,0.3593,0.2894,0.1718,0.1319,0.2468,0.4056,0.7986,0.9995/
c ... Ocean -18000 (CLIMAP) - Traitement : H. Gallee, fevrier 1990
c
c
c **********************************************************************
c
       jcalot = 1 + int(45./gs)
c ...  jcalot : latitude minimum pour situer des petits glaciers
c
c      Traitement de l'ocean
c      *********************
       call pol(ajwa0,y,9,np)
c
       if (lgm18k) then
        call pol(ajolgm,ajosu0,18,np)
        do 1 j=1,jp
        y(j) = amin1(ajosu0(j),y(j))
c ...   A modifier pour LGM
 1      continue
       else
        call pol(ajo0,ajosu0,9,np)
       end if
c
       do 10 j=1,jp
c
       ajw = y(j)/ajosu0(j)
       ajw = amax1(ajw,.005) 
c ...  condition fraction minimale de leads
c
       ajasu1(j,2)=ajw*ajosu0(j)
       ajasu1(j,1)=ajosu0(j)-ajasu1(j,2)
       do 10 iilw=1,2
       szp01(j,iilw) = 0.
       szp02(j,iilw) = 0.
   10  continue
c
c **********************************************************************
c
c      Traitement du continent
c      ***********************
       call pol(ajla0,y,9,np)
c
       if (lgm18k) then
c
c ...   18000 BP (cfr. CLIMAP)
c       ----------------------
        do 23 j=1,jp
        ajasu1(j,3)=y(j)
        ajl = y(j)/(1.-ajosu0(j))
        ajl = amin1(ajl,1.0000)
        ajl = amax1(ajl,0.0001)
c ...   condition fraction maximale de continent libre de neige
c                          minimale
c ...   A modifier pour LGM
        ajasu1(j,3)=ajl*(1.-ajosu0(j))
        ajasu1(j,4)=1.-ajosu0(j)-ajasu1(j,3)
 23     continue
c
        call pol(ajglgm,y,18,np)
        do 25 j=1,jp
        ajasu1(j,5) = amax1(0.0,y(j))
 25     continue
c
        call pol(ajnlgm,y,18,np)
        do 26 j=1,jp
        ajasu1(j,6) = amax1(0.0,y(j))
 26     continue
c
        call pol(ajflgm,y,18,np)
        do 27 j=1,jp
        ajasu1(j,7) = amax1(0.0,y(j))
        do 27 ical=5,7
        ajca1(j,ical)=ajasu1(j,ical)
 27     continue
c
c **********************************************************************
c
c ...  Temps Present
c      -------------
       else
c
        do 30 j=1,jp
        ajasu1(j,3)=y(j)
        ajl = y(j)/(1.-ajosu0(j))
        if (ajl.lt..0001) ajl = .0001
c ...   condition fraction minimale de continent libre de neige
c
        ajasu1(j,3)=ajl*(1.-ajosu0(j))
        ajasu1(j,4)=1.-ajosu0(j)-ajasu1(j,3)
        if (ajasu1(j,4).lt.0.0003.and.j.ge.jcalot) then
         daj        =0.0003 - ajasu1(j,4)
         ajasu1(j,3)=ajasu1(j,3) - daj
         ajasu1(j,4)=0.0003
        end if
        if (j.ge.jcalot) then
         ajasu1(j,4)=ajasu1(j,4)-0.0002
         ajasu1(j,6)=0.0001
         ajasu1(j,7)=0.0001
c ...    petits glaciers permettant de calculer l'accumulation
c ...    rem.: 0.0003 est une fraction suffisamment petite 
c              pour ne pas etre superieure a ajasu1(j,3) + ajasu1(j,4)
        else
         ajasu1(j,6) = 0.
         ajasu1(j,7) = 0.
        end if
   30   continue
c
c **********************************************************************
c
        do 41 ical=1,4
        do 42 jj=1,46
        j=47-jj
        x(j)=ajc0(jj,ical)
   42   continue
        call intpol(x,z,46,npp)
        do 43 j=1,jp
        ajc1(j,ical) = .5*(z(j)+z(j+1))
   43   continue
   41   continue
        do 44 j=1,jp
        ajcat = 0.
        do 45 ical=1,4
        ajcat = ajcat + ajc1(j,ical)
   45   continue
        if (ajcat.gt.0.) ajcat = (1.-ajosu0(j)) / ajcat
        do 46 ical=1,4
         ajc1(j,  ical) = ajc1(j,ical) * ajcat
        if (ical.ne.4)
     .  ajca1(j,4+ical) = ajc1(j,ical)
   46   continue
         ajasu1(j,5) = ajca1(j,5)
        if (j.ge.jcalot.and.ajc1(j,2).le.0.) then
c ...   il n'est pas permis aux petits glaciers generant la Laurentide
c       de depasser les limites du continent Nord-Americain.
         ajasu1(j,6) = 0.
         ajasu1(j,4) = ajasu1(j,4) + 0.0001
        end if
        if (j.ge.jcalot.and.ajc1(j,3).le.0.) then
c ...   il n'est pas permis aux petits glaciers generant la Fennoscandienne
c       de depasser les limites du continent Eurasiatique.
         ajasu1(j,7) = 0.
         ajasu1(j,4) = ajasu1(j,4) + 0.0001
        end if
   44   continue
c
       end if
c
c **********************************************************************
c
c --- initialisation des altitudes continentales
c
       if (lgm18k) then
c
c ...  18000 BP (CLIMAP)
c      -----------------
c
        call pol(zjglgm,y,18,np)
        do 55 j=1,jp
        szp01(j,5) = amax1(0.0,y(j))
 55     continue
c
        call pol(zjnlgm,y,18,np)
        do 56 j=1,jp
        szp01(j,6) = amax1(0.0,y(j))
 56     continue
c
        call pol(zjflgm,y,18,np)
        do 57 j=1,jp
        szp01(j,7) = amax1(0.0,y(j))
 57     continue
c
        call pol(zjclgm,y,18,np)
        do 58 j=1,jp
        szp01(j,3) = amax1(0.0,y(j))
        szp01(j,4) = amax1(0.0,y(j))
        do 58 iilw = 3,nilw
        szp02(j,iilw) = szp01(j,iilw)
 58     continue
c
c ...   Temps present
c       -------------
       else
c
        do 70 ical=1,4
        do 71 jj=1,46
        j=47-jj
        x(j)=secz0(jj,ical)
   71   continue
        call intpol(x,z,46,jpp)
        do 72 j=1,jp
        azc1(j,ical)=.5*(z(j)+z(j+1))
   72   continue
   70   continue
c
        do 73 j=1,jp
        szp01(j,5)=azc1(j,1)
        szp02(j,5)=szp01(j,5)
        ajtot = 0.
        ajhtot = 0.
        do 74 ical=2,4
         ajtot =  ajtot + ajc1(j,ical)
        ajhtot = ajhtot + ajc1(j,ical) * azc1(j,ical)
 74     continue
        if (ajtot.ne.0.) zjtot=ajhtot/ajtot
        do 75 iilw=3,4
        szp01(j,iilw)=zjtot
        szp02(j,iilw)=szp01(j,iilw)
 75     continue
c
        szp01(j,6   )= 500.
        szp01(j,7   )=1500.
c ...   altitudes arbitraires pour reperer la ligne d'equilibre
c       elles ne seront plus utilisees si inpu = .true.
c       dans ce cas c'est la lecture de cal.dat qui prevaut
        do 76 iilw=6,7
c       szp01(j,iilw)=azc1(j,iilw-4)
c ...   altitude reelle moyenne zonale du continent sous-jacent.
        szp02(j,iilw)=szp01(j,iilw)
 76     continue
 73     continue
       end if
c
c **********************************************************************
c
      if (.not.season) then
       do 102 j=1,jp
       do 102 iilw=1,4
       ajsu0(j,iilw) = ajasu1(j,iilw)
  102  continue
      end if
c
      do 103 j=1,jp
      do 103 iilw=5,nilw
      ajsu0(j,iilw) = ajasu1(j,iilw)
  103 continue
c
c **********************************************************************
c
      write(iwr1,60)(iilw,iilw=1,7)
   60 format(/'   -- iniaja --',/,' ajca1',7i12)
      write(iwr1,61)(j,(ajca1(j,iilw),iilw=1,7),j=1,jp)
   61 format(i6,7f12.5)
      write(iwr1,62)(iilw,iilw=1,7)
   62 format(' szp01',7i12)
      write(iwr1,63)(j,(szp01(j,iilw),iilw=1,7),j=1,jp)
   63 format(i6,7f12.1)
      write(iwr1,64)(iilw,iilw=1,7)
   64 format(' ajasu1',7i12)
      write(iwr1,61)(j,(ajasu1(j,iilw),iilw=1,7),j=1,jp)
c
      write(iwr6,739)
 739  format(/,' fin --- iniaja ---',//,1x)
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  la1=====
c cray .endif
      subroutine inila1
c
c +++ bcml1.for +++
c
      parameter (np=18)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/sn0/ajsn0(12,np)
c
      dimension ajs0(12,9),x(9),y(np)
c
      data ajs0/36*0.00,
     &.175,.200,.160,.115,.075,.045,.070,.065,.065,.100,.110,.140,
     &.770,.780,.450,.125,.050,.015,0.00,0.00,.010,.020,.105,.555,
     &.940,.940,.870,.695,.270,.090,.030,.020,.020,.190,.795,.895,
     &1.00,1.00,1.00,1.00,.870,.360,.120,.070,.125,.880,1.00,1.00,
     &1.00,1.00,1.00,1.00,.995,.990,.660,.590,.685,1.00,1.00,1.00,
     &1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/
c
      do 1 i=1,12
c
      do 2 j=1,9
      x(j)=ajs0(i,j)
    2 continue
c
      call pol(x,y,9,np)
c
      js = ( np / 9 ) * 3 - 1
      do 3 j=1,js-1
      y(j)=0.
    3 continue
c
      do 4 j=js,np
      if (y(j).lt.0.    ) y(j)=0.
      if (y(j).gt.0.9999) y(j)=0.9999
    4 continue
c
      do 5 j=1,np
      ajsn0(i,j)=y(j)
    5 continue
    1 continue
c
      write(iwr1,60)(i,i=1,12)
   60 format(//'   -- inila1 -- ',//,
     .   6x,' fraction continentale recouverte de neige',/,7x,41('-'),
     . /,6x,'j',4x,12i8)
      write(iwr1,61)(j,(ajsn0(i,j),i=1,12),j=1,np)
   61 format((i7,4x,12f8.4))
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  la2=====
c cray .endif
      subroutine inila2(il,is)
c
c +++ bcml2.for +++
c
      parameter(np=18,npp=19,nilw=7)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical       lgm18k,taiga,otter,summer,co2var
      double precision  hssn2,hspsn2
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/varl02/lgm18k,taiga,otter,summer,co2var
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/pal/apal,fpal,ipal,ipaleo
c
c --- constantes physiques
c
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables surface
c
       common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
       common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .             vopsn3(np,nilw), vosn3(np,nilw),
     .            hspsn3(np,nilw), hssn3(np,nilw),
     .            ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
       common/sn5/rsusn5(np,nilw),rslsn5(np,nilw),rssn5(np,nilw)
       common/sn6/hfasn6(np,nilw),  nsn6(np,nilw),nnsn6(np,nilw)
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .           ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
       common/la1/tcmla1(np),tcdla1(np),csmla1(np),csdla1(np)
       common/la2/tcmla2(np),tcdla2(np)
       common/la3/t4sla3(np,nilw)
       common/ca1/ajca1(np,nilw)
       common/ca2/caca2(np,nilw),daca2(np,nilw)
c
      dimension t4a0(9),t4s0(9),t4j0(9)
      dimension t4a(np),t4s(np),t4j(np),ass0(9),ass(np)
      dimension fc0 (9) ,fd0 (9) ,fc  (np),fd  (np)
      dimension fc00(18),fd00(18),fc18(18),fd18(18)
c
c     data t4a0  /300.1,297.7,293.5,288.3,281.8,276.2,270.1,262.4,255.1/
c ... ancienne initialisation (mv30)
c
      data t4a0/26.3,26.3,23.2,15.9, 8.4,  2.2, -5.5,-12.7,-18.0/
      data t4s0/26.0,24.0,18.0, 8.2,-2.4,-10.6,-22.6,-24.8,-33.9/
      data t4j0/26.1,27.3,27.2,23.5,18.4, 14.4, 12.3,  2.8,  1.5/
c ... source : warren & schneider (1979) moyenne annuelle,janvier,juillet
c
      data ass0/3*0.00,.175,.770,.940,1.00,1.00,1.00/
c ... source : robock             (1980) moyenne de janvier
c
      data fc0 /.0634,.1410,.2540,.3289,.3735,.3339,.3120,.1999,.0448/
c ...      fc0 : fraction de prairie / toundra ...
c
      data fd0 /.0   ,.0042,.0188,.0290,.0914,.2404,.3202,.0   ,.0   /
c ...      fd0 : fraction de foret   / taiga   ...
c
c ... source : robock - mwr p.275 (1980)
c
      data fc00
     ./0.1072,0.2267,0.2602,0.2569,0.3029,0.4122,0.4546,0.4894,0.5349,
     . 0.5561,0.5091,0.4025,0.5263,0.7559,0.5124,0.3420,0.1754,0.0000/
c ...      fc00 : fraction de prairie  / toundra
c
      data fd00 
     ./0.1447,0.0714,0.0349,0.0586,0.0697,0.0418,0.0573,0.0824,0.0981,
     . 0.1678,0.2513,0.3120,0.3246,0.1416,0.0146,0.0000,0.0000,0.0000/
c ...      fd00 : fraction de foret    / taiga
c ... source : CLIMAP climat present - Traitement par cmapeik.f ,fevrier 1990
c
      data fc18
     ./0.3075,0.3326,0.3202,0.3312,0.3853,0.4573,0.4740,0.4923,0.5200,
     . 0.4861,0.4500,0.4259,0.5381,0.5414,0.4944,0.3364,0.0265,0.0000/
c ...      fc18 : fraction de prairie  / toundra
      data fd18
     ./0.0412,0.0441,0.0205,0.0340,0.0487,0.0200,0.0406,0.0458,0.0234,
     . 0.0275,0.0472,0.0516,0.0396,0.0072,0.0000,0.0000,0.0000,0.0000/
c ...      fd18 : fraction de foret    / taiga
c ... source : CLIMAP 18000 BP       - Traitement par cmapeik.f ,fevrier 1990
c
c      taiga = .false.
c ...  cfr bcmctr
c ...  utilisation taiga feedback ou non pour le climat present
c
c      otter = .false.
c ...  cfr bcmctr
c ...  utilisation taiga feedback Otterman
c      ATTENTION : otter est egalement introduit en variable locale 
c                  dans bcmb0, bcml2, bcml3 !
c                  ne pas oublier de faire la modification appropriee !
c
       call pol(t4a0,t4a,9,np)
       call pol(t4s0,t4s,9,np)
       call pol(t4j0,t4j,9,np)
       call pol(ass0,ass,9,np)
       if (lgm18k) then
       call pol(fc18,fc ,18,np)
       call pol(fd18,fd ,18,np)
       else
       call pol(fc0 ,fc ,9,np)
       call pol(fd0 ,fd ,9,np)
       end if
c
      if (.not.inpu) then
c
       do 10 j=1,jp
c
       if (season) then
           t4su0(j     ) = t4s(j) +273.15
c
c sell     hssn2(j,is  ) = ass(j)**2
c sell     hssn3(j,is  ) = hssn2(j,is)
c ...      epaisseur de la neige (m): parametrisation de sellers
c
       else
           t4su0(j     ) = t4a(j) + tfsth4
       end if
       do 11 iilw=il,is
c
       if (season) then
          t4s   (j)      = t4s(j) + tfsth4
          t4ssu0(j,iilw) = t4s(j)
          t4sla3(j,iilw) = t4s(j)
          t4rsu0(j,iilw) = t4s(j)
       else
          t4ssu0(j,iilw) = t4a(j) + tfsth4
          t4sla3(j,iilw) = t4a(j) + tfsth4
          t4rsu0(j,iilw) = t4a(j) + tfsth4
       end if
c
          t4dsu1(j,iilw) = t4a(j) + tfsth4
 11    continue
          tcdla1(j)      = t4a(j) + tfsth4
          tcdla2(j)      = t4j(j) + tfsth4
 10    continue
      else
       do 210 j=1,jp
       t4s(j) = 0.
       ajc    = 0.
       ajsn   = 0.
       do 212 iilw=4,nilw
       t4s(j) = t4s(j) + ajsu0(j,iilw) * t4ssu0(j,iilw)
       if (t4ssu0(j,iilw).ne.0.) ajc    = ajc    + ajsu0(j,iilw)
       ajsn   = ajsn   + ajsu0(j,iilw)
 212   continue
c
       if (ajc.gt.0.) then
        t4s(j) = t4s(j) / ajc
       else
        t4s(j) = t4su0(j)
       end if
c
       ass(j) = ajsn / (1.-ajosu0(j))
 210   continue
      end if
c
      do 310 j=1,jp
c
      if (ass(j).lt.0.) ass(j) = 0.
      if (ass(j).gt.1.) ass(j) = 1.
c
      ajca = 0.
      do 312 iilw=5,nilw
      ajca = ajca + ajsu0(j,iilw)
  312 continue
c
                ajsu0(j,il  ) = (1.-ajosu0(j)) * (1.-ass(j))
                ajsu0(j,is  ) = (1.-ajosu0(j)) *     ass(j)  - ajca
c
      if (ajsu0(j,is  ).lt.0.) then
                ajsu0(j,il  ) = ajsu0(j,il) + ajsu0(j,is)
                ajsu0(j,is  ) = 0.
      end if
c
      do 313 iilw=il,is
      if (ajsu0(j,iilw).gt.0.and.t4ssu0(j,iilw).eq.0.) then
       t4ssu0(j,iilw) = t4s(j)
      end if
       t4sla3(j,iilw) = t4ssu0(j,iilw)
  313 continue
      if (t4ssu0(j,is  ).gt.tfsth4) then
       t4ssu0(j,is)   = tfsth4
       t4sla3(j,is)   = tfsth4
      end if
c
      if (ajsu0(j,il).lt.0.) then
       dajct = - ajsu0(j,il) + .1e-4
       ajct  = 0.
       do 332 iilw = 5,nilw
        ajct = ajct + ajsu0(j,iilw)
 332   continue
       do 333 iilw=5,nilw
       ajsu0(j,iilw) = ajsu0(j,iilw) - dajct * ajsu0(j,iilw) / ajct
 333   continue
       write(iwr6,336) j,ajsu0(j,il)
 336   format(' inila2 j=',i3,3x,'aj land =',e12.4,' < 0')
       ajsu0(j,il) = 0.
       ajsu0(j,is) = 0.1e-4
      end if
c
      hssn2(j,is) = (ajsu0(j,is) / (1.-ajosu0(j)))**2
      hssn3(j,is) = hssn2(j,is)
c ... epaisseur de la neige (m): parametrisation de sellers
c
      daca2(j,is) = 1.d0
 310  continue
c
      if (.not.(taiga.or.inpu)) then
       rsuc = .850
c petz rsuc = .850
       rslc = .700
       rsud = .400
       rsld = .400
       do 410 j=1,jp
       if (fc(j).lt.0.) fc(j) = 0.
       if (fd(j).lt.0.) fd(j) = 0.
       fcd = fc(j) + fd(j)
       if (fcd.le.0.) then
        rsusn5(j,is) = 0.85
        rslsn5(j,is) = 0.70
       else
        rsusn5(j,is) = (fc(j)*rsuc + fd(j)*rsud) / fcd
        rslsn5(j,is) = (fc(j)*rslc + fd(j)*rsld) / fcd
       end if
 410   continue
c
      else
       if (otter) then
        do 420 j=2,jp-1
        if (tcdla2(j-1).lt.268.15) then
         rsus = .850
c petz   rsus = .850
         rsls = .7
        else
         rsus = .4
         rsls = .4
        end if
c
        if (tcdla2(j  ).lt.268.15) then
         rsu = .850
c petz   rsu = .850
         rsl = .7
        else
         rsu = .4
         rsl = .4
        end if
c
        if (tcdla2(j+1).lt.268.15) then
         rsun = .850
c petz   rsun = .850
         rsln = .7
        else
         rsun = .4
         rsln = .4
        end if
c
        rsusn5(j,is) = (rsus + 2.*rsu + rsun) /4.
        rslsn5(j,is) = (rsls + 2.*rsl + rsln) /4.
c
c ...  'ecocline' : delimitation taiga-toundra pour caracteriser albedo
c                   du champ de neige (cfr otterman & al, jam, 1984)
c       albedo max neige sur toundra et calotte : .85 < alb < .7
c       albedo max neige sur taiga              : .4 (cfr Robock, 1980)
 420    continue
        rsusn5(1,is) = rsusn5(2,is)
        rslsn5(1,is) = rslsn5(2,is)
        rsusn5(jp,is) = rsusn5(jp-1,is)
        rslsn5(jp,is) = rslsn5(jp-1,is)
       else
        if (summer) then
         tmax = 290.
         tran =  12.
        else
         tmax = 274.
         tran =  14.
        end if
        call taifdb(tmax,tran,is)
       end if
      end if
c
      do 450 j=1,jp
       rssn5(j,is)=(rsusn5(j,is)-rslsn5(j,is))/5.
        nsn6(j,is)= 0
       nnsn6(j,is)= 0
      hfasn6(j,is)= 0.
 450  continue
c
      write(iwr1,460) (j,fc(j),fd(j),rsusn5(j,is),rslsn5(j,is),j=1,jp)
 460  format(/,14x,'fc',10x,'fd  rs* (cold)  rs* (melt)',
     .    18(/,(i4,4f12.4)))
      return
      end

      subroutine taifdb(tmax,tran,is)
      parameter(np=18,npp=19,nilw=7)
c _dp implicit double precision (a-h,o-z)
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/la2/tcmla2(np),tcdla2(np)
      common/sn5/rsusn5(np,nilw),rslsn5(np,nilw),rssn5(np,nilw)
      rsmx = 0.85
      rsmx0= 0.70
      rsmn = 0.4
      do 430 j=1,jp
      ftun = (tmax - tcdla2(j)) / tran
      ftun = amax1(0.,ftun)
      ftun = amin1(1.,ftun)
      rsusn5(j,is) = ftun * rsmx + (1.-ftun) * rsmn
      rslsn5(j,is) = ftun * rsmx0+ (1.-ftun) * rsmn
c ... Adaptation parametrisation Harvey, 1988, Climatic Change, p.196
c       (parametrisation basee sur la temperature estivale)
 430  continue
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tland===
c cray .endif
      subroutine tland(il,is)
c
c +++ bcml3.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical       lgm18k,taiga,otter,summer,co2var
      logical fulwri
      double precision  dxld,dxsn
      double precision  r1p,r2p,r3p,r1n,r2n,r3n,flux,bil,ajold
c
      dimension ajold(np,nilw),fsnos(np)
      dimension r1p(np),r2p(np),r3p(np)
      dimension r1n(np),r2n(np),r3n(np),flux(np),bil(np)
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/varl02/lgm18k,taiga,otter,summer,co2var
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables flux de chaleur
c
      common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .           pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
      common/hv1/htshv1(6,np),shshv1(np)
c
c --- variables surface
c
      double precision  hssn2,hspsn2
      common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
      common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .           vopsn3(np,nilw), vosn3(np,nilw),
     .           hspsn3(np,nilw), hssn3(np,nilw),
     .           ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/su1/tmsu1(np),t4dsu1(np,nilw),t4msu1(np,nilw),t4hsu1,
     .          ajasu1(np,nilw),ajmsu1(np,nilw),bjmsu1(np,nilw)
      common/su3/hswsu3(np,nilw),cksu3(np,nilw)
      common/la0/cslla0(np),csila0(np),csla0(np),tcla0(np)
      common/la1/tcmla1(np),tcdla1(np),csmla1(np),csdla1(np)
      common/la2/tcmla2(np),tcdla2(np)
      common/la3/t4sla3(np,nilw)
      common/bi0/surbi0(np,nilw),atmbi0(np),topbi0(np)
c
c **********************************************************************
c
c --- initialisation
c
      fulwri = .false.
c
      if (it.eq.0) then
       call inila1
       call inila2(il,is)
c
c **********************************************************************
c
      else
c
       if (season) then
c
c --- energy before computation
c
        do 1 j=1,jp
        r1p(j) = t4sla3(j,il) * ajsu0(j,il) / cslla0(j)
        r2p(j) = t4rsu0(j,il) * ajsu0(j,is) / cslla0(j)
        r3p(j) =  hssn2(j,is) * ajsu0(j,is) * hfsth4 * c33th4
        do 2 n=1,6
        do 2 iilw=il,is
        htshv0(n,j,iilw) = 0.
 2      continue
 1      continue
c
        pr = 1.
c ...   pr = 1. => schema numerique pronostique
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
       else
c
        pr = 0.
c ...   pr = 0. => schema numerique diagnostique
c
       end if
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c --- modele de sol (deardorff, 1978)
c
c --- continent libre de neige
c
       do 202 j=1,jp
       if (ajsu0(j,il).ne.0.) then
        scs = 0.
        sch = 0.
        do 203 n=1,6
        scs = scs + cshv0(n,j,il)
        sch = sch + chhv0(n,j,il)
 203    continue
        t4old = t4sla3(j,il)
        call tsol(t4sla3(j,il),cslla0(j),scs,sch,tcdla1(j),pr,
     .            alpha,beta)
        t4hts = alpha *t4sla3(j,il) + beta *t4old
        if (tsclim) t4hts = t4ssu0(j,il)
        do 204 n=1,6
        htshv0(n,j,il) = chhv0(n,j,il) - cshv0(n,j,il) *t4hts
 204    continue
       else
        t4sla3(j,il) = 0.
       end if
 202   continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c --- continent recouvert de neige et champ de neige
c
       do 212 j= 1,jp
       dtsnow =  0.
c
       if (ajsu0(j,is).ne.0.) then
        scs = 0.
        sch = 0.
        cshv0(6,j,is) = cshv0(4,j,is)*rhsth4
        chhv0(6,j,is) = chhv0(4,j,is)*rhsth4
c ...   pour tenir compte de la partie chaleur de fonte retiree 
c       du manteau neigeux lors de la sublimation
        do 213 n=1,6
        scs = scs + cshv0(n,j,is)
        sch = sch + chhv0(n,j,is)
 213    continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        if (season) then
         t4rsu0(j,is)  = t4sla3(j,is)
c ...    temperature du champ de neige a la fin du pas de temps precedent
c
         cdsh          =  cdsth4 *.6    /  hssn3(j,is)
         schr          =  sch *    cdsh / (scs+cdsh)
         scsr          = cdsh *(1.-cdsh / (scs+cdsh))
         t4old         = t4rsu0(j,il)
         call tsol(t4rsu0(j,il),cslla0(j),scsr,schr,tcdla1(j),pr,
     .             alpha,beta)
         t4rl          = alpha *t4rsu0(j,il) + beta *t4old
         t4sla3(j,is)  = (sch +cdsh *t4rl) / (scs +cdsh)
c
         if (t4sla3(j,is).gt.tfsth4) then
          t4rsu0(j,il)  = t4old
          schr          = cdsh * tfsth4
          scsr          = cdsh
          call tsol(t4rsu0(j,il),cslla0(j),scsr,schr,tcdla1(j),pr,
     .              alpha,beta)
          t4rl          = alpha *t4rsu0(j,il) + beta *t4old
          f1            = (sch +cdsh *t4rl) - (scs +cdsh) *tfsth4
          dtsnow        = t4sla3(j,is) - t4rsu0(j,is)
          if (f1.lt.0.d0) t4sla3(j,is) = tfsth4
          t4hts         = tfsth4
         else
          f1            = 0.
          t4hts         = t4sla3(j,is)
         end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        else
         t4old         = t4sla3(j,is)
         call tsol(t4sla3(j,is),cslla0(j),scs,sch,tcdla1(j),pr,
     .             alpha,beta)
         t4hts         = alpha *t4sla3(j,is) + beta *t4old
        end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        if (tsclim) t4hts = t4ssu0(j,is)
        do 214 n=1,6
        htshv0(n,j,is)= chhv0(n,j,is) - cshv0(n,j,is) *t4hts
 214    continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
       else
        f1           = 0.
        dtsnow       = 0.
         hssn3(j,is) = 0.
        t4sla3(j,is) = 0.
        t4rsu0(j,il) = t4sla3(j,il)
       end if
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
       if (season) then
c
c  -- variation des fractions continent-neige --
c
        ajold(j,il) = ajsu0(j,il)
        ajold(j,is) = ajsu0(j,is)
        ajcont      = ajsu0(j,il) + ajsu0(j,is)
        tcla0(j) = (ajsu0(j,is)*t4sla3(j,is) +ajsu0(j,il)*t4sla3(j,il))
     .           / (ajsu0(j,is)              +ajsu0(j,il)             )
        if (summer) then
c       if (iday.gt.180.and.iday.le.210) 
        if (iday.eq.180) 
     .  tcmla2(j) = tcmla2(j) + tcla0(j)
        else
        tcmla2(j) = tcmla2(j) + tcla0(j)
        end if
c
        call ajland(dtsnow,f1,cslla0(j),ajcont,hfal,hmel,htasn,fsnos,
     .              j,il,is)
c
        if (ajcont.gt.0.) then
         accsn3(j,is) = accsn3(j,is) + hfal
         ablsn3(j,is) = ablsn3(j,is) + hmel
c
         hssn3(j,is) = hssn2(j,is) + .1d0*( 0.01d0 -hssn2(j,is)       )
     .                                   *(10.00d0+t4sla3(j,is)-tfsth4)
         if (hssn3(j,is).lt.0.01d0     ) hssn3(j,is) = 0.01d0
         if (hssn3(j,is).gt.hssn2(j,is)) hssn3(j,is) = hssn2(j,is)
c
          vosn3(j,is) =  vosn3(j,is) + ajsu0(j,is) * hssn2(j,is)
          ajsn3(j,is) =  ajsn3(j,is) + ajsu0(j,is)
c
         if (.not.assn3(j,is).gt.ajsu0(j,is)) then
          assn3(j,is) =  ajsu0(j,is)
         end if
c
         dxld   = ajsu0(j,il)-ajold(j,il)
c
         if (dxld.gt.0.) then
          dxsn   =  0.
         else
          dxsn   = -dxld
          dxld   =  0.
         end if
c
c  -- ajustement temperatures du continent et du champ de neige
c
c   - le champs de neige s'etend -
         if (dxsn.ne.0.) then
          t4rsu0(j,il) = (ajold(j,is)*t4rsu0(j,il)+dxsn*t4sla3(j,il))
     .                 /  ajsu0(j,is)
         end if
c
c   - le champs de neige se retire -
         if (dxld.ne.0.) then
          t4sla3(j,il) = (ajold(j,il)*t4sla3(j,il)+dxld*t4rsu0(j,il))
     .                 /  ajsu0(j,il)
          hswsu3(j,il) = (ajold(j,il)*hswsu3(j,il)+dxld             )
     .                 /  ajsu0(j,il)
c ...     hypothese : le facteur de disponibilite en eau 
c                     sous le manteau neigeux est egal a 1
         end if
c
        end if
c
       end if
c
 212   continue
c
c      if (it.ge.2400) then 
c       if (it.eq.2400) then
c        open(unit=18,status='new',file='bcms.d')
c        rewind 18
c       end if
c       if (mod(it,20).eq.0) then
c        write(18,180)mm
c180     format(i3,48x,'60-65N ! ',51x,'65-70N ! ',
c    .          /,2(5x,'t *',4x,'t r*',5x,'h *',5x,'f *',
c    .              5x,'S0',5x,'IR',5x,'HS',5x,'HL',5x,'FC ! '))
c       end if
c        write(18,181)(t4ssu0(j,is),t4rsu0(j,il),hssn3(j,is),fsnos(j),
c    .                 (htshv0(n,j,is),n=1,5),j = 13,14)
c181     format(2(2f8.2,f8.3,f8.4,5f7.2,' ! '))
c      end if
c      if (it.eq.nts) close(unit=18)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
       if (season) then
c
c --- energy input
c
        do 3 j=1,jp
        flux(j) = 0.
        do 4 iilw=il,is
        do 4 n=1,6
        flux(j) = flux(j) + ajold(j,iilw)*htshv0(n,j,iilw)
 4      continue
        flux(j) = flux(j) * delt
 3      continue
c
c --- energy after  computation
c
        do 7 j=1,jp
        r1n(j) = t4sla3(j,il) * ajsu0(j,il) / cslla0(j)
        r2n(j) = t4rsu0(j,il) * ajsu0(j,is) / cslla0(j)
        r3n(j) =  hssn2(j,is) * ajsu0(j,is) * hfsth4 * c33th4
 7      continue
c
c --- energy balance
c
        do 5 j=1,jp
        bil(j) = (r1n(j) - r1p(j) + r2n(j) - r2p(j) - r3n(j) + r3p(j)
     .   - flux(j)) / delt
c       if (.not.tsclim.and.abs(bil(j)).gt..05) write(iwr6,6)it,j,bil(j)
c ...   attention, pour faire la verification, il ne faut pas oublier 
c       de supprimer la conduction a long terme dans le sol.
 6      format('  -- tland --  * error *  iter no ',i5,3x,'j =',i2,3x,
     .          'bilan =',e12.4)
 5      continue
c
c --- print
c
        if (.not.fulwri) go to 60
        if (iprint.eq.1) go to 61
        if (mts.ne.0)    go to 60
        if (monpr.gt.mm) go to 60
 61     continue
        write(iwr1,62)
 62     format(/'  -- tland --')
        write(iwr1,63)
 63     format(5x,'j',9x,'r1p',9x,'r2p',9x,'r3p',8x,'flux',
     .                 9x,'r1n',9x,'r2n',9x,'r3n',9x,'bil')
        write(iwr1,64)(j,r1p(j),r2p(j),r3p(j),flux(j),
     .                   r1n(j),r2n(j),r3n(j),bil(j),j=1,jp)
 64     format((i6,8e12.4))
 60     continue
c
       end if
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c --- moyennes continentales
c
       do 400 j=1,jp
c
       if (season) then
         csla0(j) = ajsu0(j,is)                / cslla0(j) +
     .              ajsu0(j,il)                / cslla0(j)
        csmla1(j) = csmla1(j) + csla0(j)
        tcmla1(j) = tcmla1(j)
     .            + t4sla3(j,il) *ajsu0(j,il) /cslla0(j)
     .            + t4rsu0(j,il) *ajsu0(j,is) /cslla0(j)
       end if
c
       do 402 iilw=il,is
       do 402 n=1,6
       htshv1(n,j)   =htshv1(n,j)    + ajsu0(j,iilw)*htshv0(n,j,iilw)
       surbi0(j,iilw)=surbi0(j,iilw) + ajsu0(j,iilw)*htshv0(n,j,iilw)
 402   continue
 400   continue
c
       if (.not.tsclim) then
        do 410 j=1,jp
        t4ssu0(j,il) = t4sla3(j,il)
        t4ssu0(j,is) = t4sla3(j,is)
 410    continue
       end if
c
      end if
c
c **********************************************************************
c
c --- print sur fichier bcmls.d
c
      if (.not.mts.eq.0) go to 600
c
      write(4,631)(t4sla3(j,il),j=1,jp)
  631 format((9f8.2,4x,'tlan'))
      write(4,632)(ajsu0(j,il),j=1,jp)
  632 format((9f8.3,4x,'aj l'))
      write(4,637)(hswsu3(j,il),j=1,jp)
  637 format((9f8.3,4x,'w  l'))
      write(4,633)(htshv0(1,j,il),j=1,jp)
  633 format((9f8.2,4x,'h1 l'))
      write(4,634)(htshv0(2,j,il),j=1,jp)
  634 format((9f8.2,4x,'h2 l'))
      write(4,635)(htshv0(3,j,il),j=1,jp)
  635 format((9f8.2,4x,'h3 l'))
      write(4,636)(htshv0(4,j,il),j=1,jp)
  636 format((9f8.2,4x,'h4 l'))
c
      write(4,641)(t4sla3(j,is),j=1,jp)
  641 format((9f8.2,4x,'tsno'))
      write(4,642)(ajsu0(j,is),j=1,jp)
  642 format((9f8.3,4x,'aj s'))
      write(4,648)(hssn2(j,is),j=1,jp)
  648 format((9f8.3,4x,'hs s'))
      write(4,643)(htshv0(1,j,is),j=1,jp)
  643 format((9f8.2,4x,'h1 s'))
      write(4,644)(htshv0(2,j,is),j=1,jp)
  644 format((9f8.2,4x,'h2 s'))
      write(4,645)(htshv0(3,j,is),j=1,jp)
  645 format((9f8.2,4x,'h3 s'))
      write(4,646)(htshv0(4,j,is),j=1,jp)
  646 format((9f8.2,4x,'h4 s'))
      write(4,647)ian,mm
  647 format(' cycle annuel no',i3,5x,'mois',i3)
c
  600 continue
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ajlnd===
c cray .endif
      subroutine ajland(dtsnow,f1,csl,ajcont,hfal,hmel,htasnj,fsnos,
     .                  j,il,is)
c
c +++ bcml4.for
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables flux de chaleur verticaux
c
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables dynamiques
c
      common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .              pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
      dimension htasn(np)
      common/p01/szp01(np,nilw),spp01(np,nilw)
c
c --- variables surface
c
      common/sn0/ajsn0(12,np)
      double precision  hssn2,hspsn2
      common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
      common/sn3/ablsn3(np,nilw),accsn3(np,nilw),
     .           vopsn3(np,nilw), vosn3(np,nilw),
     .           hspsn3(np,nilw), hssn3(np,nilw),
     .           ajpsn3(np,nilw), ajsn3(np,nilw),assn3(np,nilw)
      common/sn6/hfasn6(np,nilw),  nsn6(np,nilw),nnsn6(np,nilw)
      common/ql2/qlql2(np),fsql2(np)
      common/ca2/caca2(np,nilw),daca2(np,nilw)
      common/la3/t4sla3(np,nilw)
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .             ajosu0(np),ajsu0(np,nilw)
      common/su4/t4su4(np),t4asu4(np,nilw)
      dimension fsf(np),fsnos(np),fsnoss(np),hsnmel(np),hsnfal(np)
c
c **********************************************************************
c
c --- fraction zonale du continent considere
c
      jw              = 0
      fsf   (j   )    = 0.
      fsnos (j   )    = 0.
      fsnoss(j   )    = 0.
      hsnmel(j   )    = 0.
      hsnfal(j   )    = 0.
      htasn (j   )    = 0.
c
c **********************************************************************
c
      if (ajcont.gt.0.) then
c
c --- temperatures
c
c      tts =  ( ajsu0(j,is) *(t4asu4(j,is)+t4sla3(j,is))
c    .        + ajsu0(j,il) *(t4asu4(j,il)+t4sla3(j,il)) )
c    .     /  ((ajsu0(j,is)  + ajsu0(j,il))*2.           )
c ... exp. ibm24
c
c      tts =  ( ajsu0(j,is) *              t4sla3(j,is)
c    .        + ajsu0(j,il) *              t4sla3(j,il)  )
c    .     /  ( ajsu0(j,is)  + ajsu0(j,il)               )
c ... exp. ibm25 ...
c
       tts = t4su0(j) - t4su4(j)
     .  + ( ajsu0(j,is) *t4asu4(j,is) + ajsu0(j,il) *t4asu4(j,il))
     .   /( ajsu0(j,is)               + ajsu0(j,il)              )
c
       if (tts.lt.0.75*t4su0(j)) write(iwr6,210) it,j,
     .  il,ajsu0(j,il),t4sla3(j,il),t4asu4(j,il),
     .  is,ajsu0(j,is),t4sla3(j,is),t4asu4(j,is),tts
  210  format(' * error in ajland : it =',i4,'  j =',i2,
     . /,3x,'il =',i2,'  aj = ',e11.3,'  t4s = ',f8.2,'  t4a = ',f8.2,
     . /,3x,'is =',i2,'  aj = ',e11.3,'  t4s = ',f8.2,'  t4a = ',f8.2,
     .                                '  tts = ',f8.2)
c
       jw              = 1
       fsnos (j   )    = ajsu0(j,is) / ajcont
       hspsn2(j,is)    = hssn2(j,is)
        fsql2(j   )    = 0.
c
c ----------------------------------------------------------------------
c
c --- fonte de la neige sur l'ancienne fraction du champ de neige
c
       enmel           = 0.
       enplu           = 0.
c
       if (t4sla3(j,is).gt.tfsth4) then
                          ft =  1.
        if (dtsnow.gt.0.) ft = (t4sla3(j,is) - tfsth4) / dtsnow
        enmel         =  delt * f1 * ft
        hsnmel(  j)    =  enmel  / (c33th4*hfsth4)
        htasn (  j)    =               f1 *(1.-ft) *fsnos(j) *ajcont
c ...   htasn : flux de chaleur renvoye a l'atmosphere et 
c            comptabilise dans verification bilan energetique de surface
        htady2(5,j)    =  htady2(5,j) +f1 *(1.-ft) *fsnos(j) *ajcont
        t4sla3(j,is)   =  tfsth4
       end if
c
       if (  fsnos (j).gt.0.                         ) then
c ...  champ de neige existe deja dans la bande avant l'appel de ajland
c
         dhsnow         =  hssn2(j,is) - hsnmel(j)
        if ((hsnmel(j).gt.hssn2(j,is)).and.(il.ne.is)) then
c ...   fonte totale de la neige sur le champ de neige
c
         hsnplu         = -dhsnow
c
c         enplu         =  hsnplu *hfsth4 *c33th4 *csl
c        t4rsu0(j,il)   =  t4rsu0(j,il) + enplu
c ...    s'il reneige immediatement sur la surface, celle-ci est 
c        'isolee' et ne parvient pas a evacuer la chaleur sensible
c        une meilleure solution est donc le renvoi vers l'atmosphere.
c
          enplu         =  hsnplu *hfsth4 *c33th4 
         htady2(5,j)    =  htady2(5,j) + fsnos(j) *ajcont *enplu /delt
c
          hssn2(j,is)   =  0.d0
         fsnoss(j)      =  0.d0
         t4sla3(j,is)   =  0.d0
          fsql2(j)      =  hssn2(j,is)
        else
         hsnplu         =  0.d0
         if (il.ne.is.and.hssn2(j,is).lt.10.) then
          fsnoss(j)      =  fsnos(j) * (hssn2(j,is)-hsnmel(j))
     .                               / (hssn2(j,is)-hsnmel(j)/2.)
c         fsnoss(j)      =  fsnos(j) * exp(-0.5 *hsnmel(j) /hssn2(j,is))
c ...     Parametrisation de Harvey, J.of Climate, 1988, p.1072
c
c         fsnoss(j)      =  fsnos(j) * (hssn2(j,is)-hsnmel(j)/2.)
c    .                               /  hssn2(j,is)
c ...     Parametrisation personnelle (H.G.)
c
          if (fsnoss(j).lt.0.0001d0) fsnoss(j) = 0.0001d0
         else
          fsnoss(j)      = 1.d0
         end if
         hssn2(j,is) =  fsnos(j) * (hssn2(j,is)-hsnmel(j)) / fsnoss(j)
         fsql2(j)    = hsnmel(j)
        end if
       end if
c
c ----------------------------------------------------------------------
c
c --- fraction de continent sur laquelle tombe la neige
c
       if (is.ne.il) then
                         fsf(j) =  0.
        if (tts.lt.280.) fsf(j) =.05 * ( 280.-tts )
        if (tts.lt.260.) fsf(j) =  1.
c ... parametrisation de harvey
c
c
c --- fraction de neige tombant sur calotte
c
       else
                         fsf(j) =  0.
        if (tts.lt.281.) fsf(j) =      ( 281.-tts ) / 18.
        if (tts.lt.272.) fsf(j) =.5  + ( 272.-tts ) / 38.
        if (tts.lt.253.) fsf(j) =  1.
c ... parametrisation de ledley, 1976 (these)
c
       end if
c
c --- recouvrement
c
       psi          = 1.
c ...  psi          : parametre de recouvrement
c     (psi compris entre 0 et 1)
       ql           = 0.
       ql0          = 0.
c
       if (fsf(j).ne.0.) then
c
c  -- precipitations liquides    (m eau)
        qlh         = htshv0(4,j,is) * delt / hvsth4
c ...   qlh  : diminution epaisseur couche de neige par evaporation
c
        ql0         = qlql2(j)       * delt / 86400.e3
c ...   ql0  : precipitations au large de la calotte
c
c       dh  = szp01(j,is)/1000. - 2.
c       ql0 = ql0* (2. ** (-amax1(0.0,dh)))
c ...  parametrisation de pollard, 1983, j.g.r.
c
c       ql0 = ql0* daca2(j,is)
c ...  moyenne zonale de la parametrisation de pollard, 1983, j.g.r.
c
        if (il.eq.is) then
         ql0 = ql0* caca2(j,is)
        else
         ql0 = ql0* 0.8
        end if
c ...  parametrisation d'oerlemans, 1982, j.of climatology
c
        ql = ql0 + qlh
        if (ql.lt.0.) ql = 0.
c
        hsnfal(  j) = ql / c33th4
          hfal0     = ql0/ c33th4
c ...  masse de neige precipitee  (m neige)
c
c  -- liberation de chaleur latente equivalente dans l'atmosphere
        htasn (  j) =(ql * hfsth4 / delt) *fsf(j) *ajcont + htasn (  j)
        htady2(5,j) =(ql * hfsth4 / delt) *fsf(j) *ajcont + htady2(5,j)
c
        if (ql.le.0.) fsf(j) = 0.
c
        if (fsnoss(j).le.0.) then
                             t4sla3(j,is) = tfsth4
         if (fsnos(j).le.0.) t4rsu0(j,il) = t4sla3(j,il)
                              fsnos(j)    =    fsf(j)
        else
         if (         fsf(j).gt.0.          ) then
          if((1.-psi)*fsf(j).gt.1.-fsnoss(j))
     .     psi = 1. - (1.-fsnoss(j)) / fsf(j)
          if(    psi *fsf(j).gt.   fsnoss(j))
     .     psi =          fsnoss(j)  / fsf(j)
         else
           psi = 1.
         end if
         fsnos(j)    =  fsnoss(j)    + (1.-psi) * fsf(j)
        end if
       else
        fsnos(j)    =  fsnoss(j)
       end if
c
c --- nouvelle epaisseur de neige
c
       if (fsnos(j).ne.0.)
     . hssn2(j,is) = (fsnoss(j) *hssn2(j,is) + hsnfal(j) *fsf(j))
     .             /  fsnos (j)
c
c --- energie fournie a la surface
c
       if (ajsu0(j,is).gt.0.001) then
        htshv0(5,j,il) = 0.
        htshv0(5,j,is) = - htasn(j) / ajsu0(j,is)
c       if (abs(htshv0(5,j,is)).gt.100.) 
c    .   write(iwr6,100)it,j,is,htshv0(5,j,is),htasn(j),ajsu0(j,is)
c100     format(' --- ajland ---',i7,2i3,
c    .          ' hts(5) =',e12.4,' hta* =',e12.4,' aj =',e12.4)
c *** l'ordre de ces 2 instructions est important ***
       else
        htshv0(5,j,is) = 0.
        htshv0(5,j,il) = - htasn(j) / ajsu0(j,il)
c       if (abs(htshv0(5,j,il)).gt.100.) 
c    .   write(iwr6,100)it,j,il,htshv0(5,j,il),htasn(j),ajsu0(j,il)
       end if
c
c --- fraction nouvelle de neige
c
       ajsu0(j,il) = ajcont * ( 1. - fsnos(j))
       ajsu0(j,is) = ajcont *        fsnos(j)
c *** l'ordre de ces 2 instructions est important ***
c
       fsql2(j)    = fsql2(j) * (1.e6/c33th4) * (86400./delt)
c
      end if
c
c **********************************************************************
c
      if (iprint.eq.1) go to 61
      if (mts.ne.0)    go to 60
      if (monpr.gt.mm) go to 60
   61 continue
      if (j .eq.1) write(iwr1,62)il,is,it
   62 format(/,'  -- ajland -- (',i1,',',i1,')',3x,'it = ',i6//,
     . 5x,'j',8x,'ho *',4x,'hsnmel',3x,'fsnow**',4x,'hsnfal',7x,'fsf',
     . 5x,'hta 5',7x,'hs2 *',7x,'hs3 *',4x,'fsnow*',
     . 4x,'ts *',4x,'tb *')
      if (jw.eq.1) write(iwr1,63)j,hspsn2(j,is),hsnmel(j),fsnoss(j),
     . hsnfal(j),fsf(j),htasn(j),hssn2(j,is),hssn3(j,is),fsnos(j),
     . t4sla3(j,is),t4rsu0(j,il)
   63 format((i6,f12.4,4f10.4,f10.2,2f12.4,f10.4,2f8.2))
   60 continue
c
       hfal  = hsnfal(j)
      if (fsnos(j).gt.0.) then
       hfalr = amin1(1.,fsf(j)/fsnos(j))
      else
       hfalr = 0.
      end if
                           hfasn6(j,is) = hfal0* hfalr * 28.512e6/delt
c ... 28.512 e+6 = 3600.*24. *100.     * 330.              / 100.
c ...              24 heures *100.cm/m * rho(neige tassee) / rho(neige f
c
      hmel   = hsnmel(j)
                                      nsn6(j,is)= nsn6(j,is)+delt/86400.
      if (hsnmel(j).eq.0.)            nsn6(j,is)= 0
                                     nnsn6(j,is)=nnsn6(j,is)+delt/86400.
      if (t4ssu0(j,is).lt.tfsth4-10.)nnsn6(j,is)= 0
      htasnj = htasn (j)
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  csland==
c cray .endif
      subroutine csland(j,il,is)
c
c +++ bcml5.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables surface
c
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/la0/cslla0(np),csila0(np),csla0(np),tcla0(np)
       common/la1/tcmla1(np),tcdla1(np),csmla1(np),csdla1(np)
c
c **********************************************************************
c
c --- calcul des capacites calorifiques des secteurs continentaux
c
       csla0(j)= ajsu0(j,il) /cslla0(j) + ajsu0(j,is) /cslla0(j)
      csmla1(j)=csmla1(j) + csla0(j)
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tsol====
c cray .endif
      subroutine tsol (t4n,csol,scs,sch,t4nd,pr,alpha,beta)
c
c +++ bcml6.for +++
c
c _dp implicit double precision (a-h,o-z)
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
      data ct/31104.e+03/
c ...      ct : constante de temps du cycle (s)
c
c --- initialisation du schema de resolution
c
c   equation diagnost.: pr   = 0.    equation pronost.: pr   = 1
c                       beta = 0.                       beta = 0.25
c                       alpha= 1.                       alpha= 0.75
c
      if (pr.gt.0.5) then
       alpha = 0.75d00
       beta  = 0.25d00
c schema purement explicite
c      alpha = 0.
c      beta  = 1.
      else
       alpha = 1.
       beta  = 0.
      end if
c
c --- resolution
c
      pct  = 2.*pi/ct
c     pct  = 0.
c ... pct  = 0. : en vue de la verification du bilan energetique 
c                 du modele de continent
c
      csdt = (csol * scs + pct) * delt
      tind1= delt * csol * sch
      tind2= delt * t4nd * pct
      t4n = (t4n * (pr - beta  * csdt) + tind1 + tind2)
     .    /        (pr + alpha * csdt)
c
      if (pr.lt.0.5) then
       t4n =  sch / scs
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  inics===
c cray .endif
      subroutine inics
c
c +++ bcmlc0.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables surface
c
      common/la0/cslla0(np),csila0(np),csla0(np),tcla0(np)
c
c **********************************************************************
c
c --- calcul des capacites calorifiques des secteurs continentaux
c
      do 10 j=1,jp
c     csila0(j)= .331e-6
c ... valeur proposee par sellers, 1983
c
      csila0(j)= .404e-6
c ... valeur modifiee - equivaut a un coeff.accum.chaleur = 0.5 W.m-2.K-1
c
c     csila0(j)= .202e-6
c ... valeur modifiee - equivaut a un coeff.accum.chaleur = 1.0 W.m-2.K-1
c
c     cslla0(j)= .321e-6
c ... valeur proposee par sellers, 1983
c
      cslla0(j)= .253e-6
c ... valeur modifiee - equivaut a un coeff.accum.chaleur = 0.8 W.m-2.K-1
c
c     cslla0(j)= .202e-6
c ... valeur modifiee - equivaut a un coeff.accum.chaleur = 1.0 W.m-2.K-1
c
   10 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  rflls===
c cray .endif
      subroutine reflls(il,is)
c
c +++ bcmlr.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9,nda=120)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical       fixalb
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      common/th4/tfsth4,hfsth4,hvsth4,cdsth4,c33th4,chfth4,rhsth4
c
c --- variables flux de chaleur verticaux
c
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables vapeur d'eau
c
      common/ql2/qlql2(np),fsql2(np)
c
c --- variables surface
c
      common/lr0/hswlr0(np,nda)
      common/rf0/rgsrf0(np,nilw)
      common/rf1/rgwrf1(np,12),rgarf1(np,nilw)
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .             ajosu0(np),ajsu0(np,nilw)
      common/su3/hswsu3(np,nilw),cksu3(np,nilw)
      double precision  hssn2,hspsn2
      common/sn2/hssn2(np,nilw),hspsn2(np,nilw)
      common/sn5/rsusn5(np,nilw),rslsn5(np,nilw),rssn5(np,nilw)
c
c **********************************************************************
c
c --- initialisation
c
      fixalb = .false.
c
c --- secteur neige     (is)
c
      call reflca(il,is)
c
c **********************************************************************
c
c --- secteur continent (il)
c
      do 20 j=1,jp
c
c  -- moyenne annuelle du contenu en eau du sol
c
       do 10 ida=nda,2,-1
       hswlr0(j,ida) = hswlr0(j,ida-1)
 10    continue
       hswlr0(j,1  ) = hswsu3(j,il)
       nnda = it
       if (nnda.gt.nda) nnda=nda
       hswloc        = 0.d0
       do 11 ida=1,nnda
       hswloc        = hswloc + hswlr0(j,ida)
 11    continue
       hswloc        = hswloc / nnda
c
      if (ajsu0(j,il).ne.0.) then
c
c --- bilan hydrique du continent => w factor
c
       if (hswsu3(j,il).eq.0.  ) hswsu3(j,il) = 1.
       hswsu3(j,il) = ( hswsu3(j,il) +delt *(qlql2(j)/86400.d3) /.15d0 )
     .  / (1. -delt *htshv0(4,j,il) /(hswsu3(j,il) *hvsth4 *.15))
       if (hswsu3(j,il).lt.0.01) hswsu3(j,il) = 0.01
       if (hswsu3(j,il).gt.1.  ) hswsu3(j,il) = 1.
c
c --- heat diffusivity (sa76)
c
       if (hswsu3(j,il).lt..7) then
        cksu3(j,il) = 12.e2+(hswsu3(j,il)-.5)*10.e2
        if (cksu3(j,il).lt.12.e2) cksu3(j,il) = 12.e2
       else
        cksu3(j,il) = 14.e2+(hswsu3(j,il)-.7)*20.e2
       end if
c
c --- albedo
c
       if (fixalb) then
        rgsrf0(j,il) = rgarf1(j,il)
       else
c
c  -- parametrisation de saltzman et ashe, 1976, tellus, p315
c
        if  (hswloc.gt.0.7) then
          rgsrf0(j,il) = .30 - 0.50 *(hswloc -.7)
        else
         if (hswloc.gt.0.5) then
          rgsrf0(j,il) = .35 - 0.25 *(hswloc -.5)
         else
          rgsrf0(j,il) = .35
         end if
        end if
c
       end if
      else
       hswsu3(j,il)=1.
        cksu3(j,il)=0.
       rgsrf0(j,il)=0.
      end if
 20   continue
c
c **********************************************************************
c
      return
      end
c
c bcmo
c
c cray .if [ $concord .eq. 1 ]
c cray nzyx  oc0=====
c cray .endif
      subroutine inioc0
c
c +++ bcmo0.for +++
c
      parameter(np=18,npp=19,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables surface
c
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .                hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
       dimension           zoc (nho),dzoc (nho)
c
          tic0 = 24.d00*3600.d00
c ...           time step (sec)
        iimic0 = delt/tic0
c ...            local total number of time steps
c
       hcric0 = 0.5
c ...           critical height of ice
       amnic0 = 0.005
c ...           minimum ocean extend
c
       data  zoc /5.,10.,15.,20.,30.,50.,75.,100.,150./
       data dzoc /5., 5., 5., 5.,10.,20.,25., 25., 50./
c      data  zoc /5.,10.,15.,20.,30.,50.,75.,100.,150.,250.,350.,500./
c      data dzoc /5., 5., 5., 5.,10.,20.,25., 25., 50.,100.,100.,150./
c ...  attention a l'initialisation des profils en temperature (cfr.bcmo2)
c
      do 1 n=1,nho
         zoc2(n)= zoc(n)
        dzoc2(n)=dzoc(n)
    1 continue
c
      do 5 j=1,jp
           dic2(j)= 30.d0
c ...    profondeur initiale couche melangee oceanique sous la glace
c        remarque importante : cette valeur ne peut en aucun cas etre 
c                 egale ou superieure a 150 m, 
c                 mais doit etre strictement egale a une valeur de zoc
         nhiic4(j)= 5
c ...    nhiic4 : contient numero du niveau du modele oceanique dont 
c             profondeur correspond a celle de dic2
         nhoic4(j)= nhiic4(j)
c ...    nhoic4 : contient numero du niveau du modele oceanique dont
c             profondeur correspond a celle de hoc2 utilisee dans fixdep
    5 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  oc1=====
c cray .endif
      subroutine inioc1
c
c +++ bcmo1.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
      common/oc1/ckwoc1,daoc1(npp),aloc1(npp)
c
c --- variables locales initialisees
c
      dimension d0(10),al0(10)
c
      ckwoc1 = 1.7e5
c
c     data d0    /90., 90.,88.,83.,75.,65.,55.,47.,42.,39./
c ... couche de melange oceanique : data de sellers
c
      data d0    /47.,125.,60.,57.,50.,50.,53.,70.,70.,90./
c ... couche de melange oceanique : data de manabe et stouffer
c
      data al0   /.24,.25,.32,.43,.44 ,.57,.62,.54,.20 ,0./
c ... fraction oceanique totale
c
c --- interpolations
c
       call intpol(al0,aloc1,10,npp)
       call intpol(d0 ,daoc1,10,npp)
c
      write(iwr1,12)
   12 format(//'   -- inioc1 --')
      write(iwr1,13)(daoc1(j),j=1,jpp)
   13 format(/,5x,'profondeur couche melangee (manabe-stouffer) :',
     .      /,(5x,10f8.1))
      write(iwr1,14)(aloc1(j),j=1,jpp)
   14 format(/,5x,'fraction ocean frontiere laterale :',
     .      /,(5x,10f8.3))
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  mixdp===
c cray .endif
c     *****************
      subroutine mixdep
c     *****************
c
c +++ bcmo10.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
      parameter(nol=0)
      parameter(nro=10)
c _dp implicit double precision (a-h,o-z)
c
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical fulwri
c
      double precision  fsol,fnsol,fadv,ustar,energ,xxx
      double precision  enp,enn,flo,den,bio
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables flux de chaleur verticaux
c
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables dynamiques
c
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
c
c --- variables surface
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  acric3,ablic3,amxic3
      common/ic3/acric3(np),ablic3(np),amxic3(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
      common/lat/j
c
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .           hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
      common/oc7/ceoc7
      double precision  ssto,heqo,fadvo
      common/o0/heqo(np),ssto(np),fadvo(np)
      double precision  qoo1,dto1
      common/o1/qoo1(np),dto1(np)
c
      character*6  route
      common/rou1/ route(np,nro)
      common/rou2/iroute(np,nro,2),airout(np,nro),hirout(np,nro),iiirou
c
c **********************************************************************
c
      do 10 iro = nro,2,-1
       route(j,iro)   =  route(j,iro-1)
      airout(j,iro)   = airout(j,iro-1)
      hirout(j,iro)   = hirout(j,iro-1)
      do 10 k=1,2
      iroute(j,iro,k) = iroute(j,iro-1,k)
 10   continue
       route(j,1  )   = 'mixdep'
      airout(j,1  )   = 1.d0-awic2(j)
      hirout(j,1  )   = hiic2(j)
      iroute(j,1,1)   = it
      iroute(j,1,2)   = iiirou
c
c **********************************************************************
      fulwri = .false.
c **********************************************************************
c
      heqo(j)=heqoc2(j)
      ssto(j)=sstoc2(j)
c
c **********************************************************************
c
c     flux solaires en surface
c
      fsol      =   chhv0(1,j,2)
      fsoc2(j)  =   chhv0(1,j,2)
c
c **********************************************************************
c
c     flux non solaires en surface (w.m-2)
c
      fnsol     = ( ablic3(j) + amxic3(j)*awpic2(j)/awic2(j) ) / tic0
      ablic3(j) =   0.
      amxic3(j) =               0.
      do 100 n=2,4
      fnsol     = fnsol + chhv0(n,j,2) - sstoc2(j) * cshv0(n,j,2)
  100 continue
      fnsoc2(j) = fnsol
c
c **********************************************************************
c
c       flux d'advection meridienne calcule a partir des gradients de ss
c       temps delt precedent
c
      fadv      =         chhv0(6,j,2) - sstoc2(j) * cshv0(6,j,2)
      fadvo(j)  = fadv
c     qoo1 (j)  = qooc3
      qoo1 (j)  = qooc3 - fadv *tic0
c
c **********************************************************************
c
c       vitesse de friction
c
      u4mb      =  .5 * (u4dy8(j) + u4dy8(j+1))
      ustar     = ceoc7 * sqrt(abs(u4mb))
      th4dy8(j) = ustar
c
c **********************************************************************
c
c       determination des nouvelles profondeur et temperature de
c       la couche limite oceanique
c
       call cmo (j,fsol,fnsol,fadv,ustar,tic0,energ)
c
c **********************************************************************
c
c       bilan energetique du modele de couche melangee
c
       heqoc2(j) = ( hoc2(j)*sstoc2(j) + energ ) / sstoc2(j)
c
c **********************************************************************
c
c --- correction si sst < t congelation de l'eau de mer (tfth3)
c
          acric3(j) = 0.d0
          dawic4    = 0.d0
          ice1      = 0
          ice2      = 0
c
          lfxd = 0
      if (sstoc2(j).le.tfth3) then
          lfxd = 1
          ice1 = 1
          hoc2(j)   = dic2(j)
          nhoic4(j) = nhiic4(j)
          tws  = tfth3
          twn  = sstoc2(j)
      else
       if (j.gt.1.and.j.lt.jp) then
        if (qoo1(j).lt.0.d0.and.awic2(j+1).le..75d0) then
         lfxd = 1
         ice2 = 1
        end if
       end if 
      end if
c
      if (lfxd.eq.1) then
         xxx=0.d0
       if (nol.ne.0) then
         do 40 icc=1,nol
         xxx=xxx+dzoc2(nho+1-icc)*twoc2(nho+1-icc,j)
c ...    xxx: contenu energetique en-dessous de la couche melangeable 
c
 40      continue
       endif
       if (ice2.eq.1) then
         hoc2(j)   =  zoc2(nho-nol)
         nhoic4(j) =       nho-nol
         sstoc2(j) = (heqoc2(j) *sstoc2(j) - xxx) /hoc2(j)
c ...    homogeneisation temperature sur epaisseur couche melangeable
c
         dawic4 = -qoo1(j) *awic2(j)*(1.d0 -awic2(j+1))
     .          / ( ocath2 *hoc2(j) *(sstoc2(j) -tfth3)
     .            + hfith3          * hcric0           )
         dawic4 = dmax1(dawic4,0.d0)
         twn  =  tfth3 - hcric0 *hfith3 /(hoc2(j) *ocath2)
         tws  = (sstoc2(j) *awic2(j) -twn *dawic4) / (awic2(j) -dawic4)
       end if
c
       acric3(j) = ocath2 *(twn -tfth3) *zoc2(nho-nol)
       sstoc2(j) = tws
       heqoc2(j) = (sstoc2(j)*zoc2(nho-nol)+xxx)/sstoc2(j)
       do 41 i=1,nho-nol
       twoc2(i,j) = sstoc2(j)
 41    continue
      end if
c
c **********************************************************************
c
c       bilan energetique de la routine
c
       enp    = heqo(j)   *ssto(j)      *ocath2
       enn    = heqoc2(j) *sstoc2(j)    *ocath2 *(1.d0 -dawic4) 
     .        +(hoc2(j)   *tfth3   +xxx)*ocath2 *       dawic4 
     .        + acric3(j) *(ice1 + ice2 *dawic4)
c
       flo    = (fsoc2(j) +fnsoc2(j) +fadv)
       den    = (enn - enp) / tic0
       bio    =  flo - den
c
c ifte
       if (.not.fulwri)      go to 50
       if (dabs(bio).gt.0.1) go to 51
       if (iprint.eq.1)      go to 51
       if (mts.ne.0)         go to 50
       if (monpr.gt.mm)      go to 50
c then
   51  continue
       write(iwr1,52)it,j,enn,enp,flo,bio,dawic4
   52  format(/,i5,i3,' abs(',e12.4,' -',e12.4,' ; /tic0 ; -',e12.4,
     .  ') =',e12.4,' >.1 w.m-2    da =',e12.4,'  -- mixdep --')
       write(iwr1,53)
   53  format(8x,' rout.',4x,'it  i0',11x,'ai',11x,'hi')
       write(iwr1,54)
     .(route(j,iro),(iroute(j,iro,k),k=1,2),airout(j,iro),hirout(j,iro),
     .         iro=1,nro)
   54  format((8x,a6,i6,i4,2e13.5))
   50  continue
c
c **********************************************************************
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  cmo=====
c cray .endif
      subroutine cmo (i,fsol,fnsol,fadv,ustar,dtd,energ)
c   ********************************************************
c   * subroutine cmo (i,fsol,fnsol,fadv,ustar,dtd,energ) *
c   *    computes mixed-layer temperature and depth ,      *
c   *      vertical profile of temperature below it        *
c   ********************************************************
c
c +++ bcmo11.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
      parameter(nol=0)
c _dp implicit double precision (a-h,o-z)
c
      integer kn,kn1
c
      logical vdif
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      double precision lambda,monino,lpsl
      double precision  dtd
      double precision  r,d1,d2
      double precision  t,tm,hr,f,delz,aux,algcp,br,dbr,bn,entflu
      double precision  dr,dh,dflu,hrin
      double precision  ust3,d,rup,rlow,delt,tst
      double precision  ratiol,dileng,mdrho,rad
      double precision  g1,g2,cp1,cp3,hsl,hbw,dzz,x,sl
      double precision  ustar,fsol,fnsol,fadv,energ,zzzz
      double precision  hrad ,gsol,gnsol,gadv
      dimension t(nho),x(nho+2,nho+2),sl(nho+2)
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/sdt,tht,it
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision m1th5,m2th5,m3th5,m4th5,m5th5,k1th5,k2th5,k3th5
      double precision  a1th5,a2th5,c2th5,c4th5
      common/th5/m1th5,m2th5,m3th5,m4th5,m5th5,a1th5,a2th5,c2th5,c4th5,
     .k1th5,k2th5,k3th5
c
c --- variables dynamiques
c
      common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
c
c --- variables surface
c
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .                 hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      common/corad/r,d1,d2
c
      vdif=.false .
c ... vertical diffusion of heat is taken into account when vdif=.true.
c
      data algcp/4.905d-07/
      r  =  0.671d00
      d1 =  1.000d00
      d2 = 17.100d00
c
c     initialisation
c
      kmax=nho
      do 16 j=1,kmax
   16 t(j)=twoc2(j,i)
      hr= hoc2(i)
      tm=twoc2(1,i)
      f=fjdy1(i)
c
      kn =kmax-nol
      if (hr .eq.zoc2(kmax-nol) ) go to 1100
 1201 kn =kn-1
 1200 if (hr .lt.zoc2(kn)) go to 1201
 1100 kn1=kn+1
      if (kn1.gt.kmax    ) kn1=kmax
      delz = zoc2(kn1)-hr
c
      ust3 = ustar**3
      if (ust3.le.0.d0) ust3=1.d-20
      gsol = fsol * dtd
      gnsol= fnsol* dtd
      gadv = fadv * dtd
      call bretd (hr,br,dbr)
      bn   = (algcp/rhwth2) * (fnsol + fsol*(1. + br/hr))
      lambda = ustar / f
      if (lambda.eq.0.d00) lambda = 1.d-10
      monino=1.d10
      if(bn.ne.0.d0)monino=ust3/bn
c
      g1=(algcp/rhwth2)*(fsol+fnsol)
      g2=(algcp/rhwth2)*fsol
      lpsl=ratiol(hr,lambda,ust3,g1,g2,br)
      cp1=(2.d0*(1.d0-m5th5)*lpsl+m4th5)/6.d0
      cp3=((m4th5-lpsl)*(m2th5+m3th5)+m5th5*m3th5*lpsl)/3.d0
c
c     select entrainment or retreat
c
      itst = 0
      if((cp3*ust3-cp1*hr*bn).le.0.d0)goto 3300
      itst = itst + 1
c
c     entrainment
c
      if (kn.eq.kmax-nol) go to 3400
      itst = itst + 10
      hsl=dileng(hr,lambda,monino)
      call ephil2(hr,ust3,bn,hsl,hbw,cp1,cp3)
      entflu = hbw * rhwth2 * dtd /9.81d00
      dr = mdrho(tm,t(kn1))
      if (dr.lt.0.d0) dr=0.d0
      dh   = delz
      dflu = hr * dr * dh
 3211 if (entflu.le.dflu) go to 3210
      tm=(tm*hr+t(kn1)*dh)/zoc2(kn1)
      hr=zoc2(kn1)
      kn1=kn1+1
      if (kn1.gt.kmax-nol) then
        entflu=0.d0
      else
        entflu=entflu-dflu
        dr=mdrho(tm,t(kn1))
        if (dr.lt.0.d0) dr=0.d0
        dh=dzoc2(kn1)
        dflu=hr*dr*dh
      endif
      goto 3211
 3210 aux=0.d0
      if (dr.ne.0.d0.and.kn1.le.(kmax-nol)) then
        aux=entflu/(hr*dr)
        tm=(tm*hr+t(kn1)*aux)/(hr+aux)
      endif
      hr=hr+aux
      kn=kn1-1
      if(kn1.gt.kmax)kn1=kmax
      delz=zoc2(kn1)-hr
      goto 3400
c
 3300 continue
c
c     retrait
c
      itst = itst + 100
      hrin = hr
      zzzz=zoc2(1)
      call rphil3(hr,lambda,ust3,g1,g2,br,zzzz)
      if (hr  .eq.hrin) go to 3400
      if (hrin.eq.zoc2(kmax)) go to 3310
c
c     compute tm, t(kn1), so as to conserve potential energy
c
      itst = itst + 1000
      d = hr
      if (zoc2(kn).gt.hr) d = zoc2(kn)
      delt = tm - t(kn1)
c
      aux = hrin * (d + delz) / (d * zoc2(kn1))
      tm  = t(kn1) + aux * delt
c
      aux = hrin * (hrin - d) / (zoc2(kn1) * (zoc2(kn1) - d))
      t(kn1) = t(kn1) + aux * delt
c
c     compute new kn, kn1, delz and set adjusted tm into released layer
c
 3310 continue
      itst = 10000
 3321 if (hr.ge.zoc2(kn)) go to 3320
      t(kn) = tm
      kn    = kn - 1
      go to 3321
 3320 kn1 = kn + 1
      if (kn1.gt.kmax) kn1=kmax
      delz = zoc2(kn1) - hr
c
 3400 continue
      jtst=0
      tst=0.d00
c
c     end of the mixed layer depth computation phase
c
c
c     heat distribution
c
      if (kn.eq.kmax) go to 3410
      jtst=jtst+1000
      rup   = rad(hr)
      hrad  = zoc2(kn1)
      rlow  = rad(hrad)
      t(kn1)= t(kn1) + gsol * (rup - rlow) / (delz  * ocath2)
      tst  = tst    +         rup - rlow
      if (kn1.eq.kmax) go to 3410
      iaux  = kn1 + 1
      do 3406 j = iaux,kmax
      rup = rlow
      hrad= zoc2(j)
      rlow= rad(hrad)
      t(j)= t(j)   + gsol * (rup - rlow) / (dzoc2(j) * ocath2)
      tst   = tst    +         rup - rlow
 3406 continue
 3410 continue
      jtst=jtst+100
      hrad=zoc2(kmax)
      tm = tm + (gnsol + (1.d00-rad(hr)+rad(hrad)) *gsol + gadv)
     .          / (hr*ocath2)
      tst= tst         +  1.d00-rad(hr)+rad(hrad)
c
c     diffusion processes
c
      kord=kmax-kn+1
      if (.not.vdif) goto 3420
      do 3415 ic=1,kord
      sl(ic)=0.d0
      do 3415 jc=1,kord
      x(ic,jc)=0.d0
 3415 continue
      call difu(x,sl,dtd,hr,tm,t,delz,kn,kn1,kmax,kord)
 3420 continue
c
c     convective adjustement of water column if gravitationally unstable
c
      if (kn.ne.kmax) then
        do 3605 jc=1,1000
        iconv=0
c
c       mix mixed layer and layer kn1 if unstable.
c
        if (mdrho(tm,t(kn1)).lt.0.d0) then
          tm=(tm*hr+t(kn1)*delz)/zoc2(kn1)
          t(kn1)=tm
          hr=zoc2(kn1)
          kn1=kn1+1
          iconv=1
 3610     continue
          if (kn1.gt.kmax) goto 3615
          if (mdrho(tm,t(kn1)).lt.0.d0) then
            tm=(tm*hr+t(kn1)*dzoc2(kn1))/zoc2(kn1)
            t(kn1)=tm
            hr=zoc2(kn1)
            kn1=kn1+1
            iconv=1
            goto 3610
          endif
 3615     continue
          kn=kn1-1
          if (kn1.gt.kmax) kn1=kmax
          delz=zoc2(kn1)-hr
        endif
c
c       mix layers below the ml. if unstable.
c
        if (kn1.lt.kmax) then
          kmaxm1=kmax-1
          do 3620 ic=kn1,kmaxm1
          if (mdrho(t(ic),t(ic+1)).lt.0.d0) then
            dzz=dzoc2(ic)
            if (ic.eq.kn1) dzz=delz
            t(ic)=(t(ic)*dzz+t(ic+1)*dzoc2(ic+1))/(dzz+dzoc2(ic+1))
            t(ic+1)=t(ic)
            iconv=1
          endif
 3620     continue
        endif
        if (iconv.eq.0) goto 3625
        if (jc.eq.1000) write(iwr6,100)
 100    format (/,1x,'convective adjustment does not converge after'
     &         ,'1000 iterations')
 3605   continue
 3625   continue
      endif
c
c     update t in the ml.
c
      do 3630 j = 1,kn
 3630 t(j) = tm
      if (hr.gt.zoc2(kmax-nol)) then
        hr=zoc2(kmax-nol)
        kn=kmax-nol
        kn1=kn+1
        delz=zoc2(kn1)-hr
      endif
      do 3706 j = 1,kmax
 3706 twoc2(j,i) = t(j)
      sstoc2(i)  = tm
      hoc2(i)  = hr
c
c     contenu energetique
c
      kn1 = kn + 1
      kn2 = kn + 2
c
      if (kn.ge.kmax.or.hr.ge.zoc2(nho-nol)) then
       energ = 0.d00
      else
       energ =         t(kn1)*delz
       if (kn1.ge.kmax) go to 57
       do 56 kc=kn2,kmax-nol
       energ = energ + t(kc)*dzoc2(kc)
 56    continue
 57    continue
      end if
      if (nol.gt.0) then
       do 58 kc=kmax-nol+1,kmax
       energ = energ + t(kc)*dzoc2(kc)
 58    continue
      end if
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  bredtd==
c cray .endif
      subroutine bretd(h,br,dbr)
c
c          ************************************************
c          *        subroutine bretd (h,br,dbr)           *
c          *   computes brad(h) = br and dbrad/dh = dbr*  *
c          ************************************************
c
c _dp implicit double precision (a-h,o-z)
      double precision  br,dbr,r,d1,d2,ed1,ed2,h,x,y
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common /corad/ r,d1,d2
      x=-h/d1
      if (x.lt.-170.d0) x=-170.d0
      y=-h/d2
      if (y.lt.-170.d0) y=-170.d0
      ed1=dexp(x)
      ed2=dexp(y)
      br=r*(-2.*d1+ed1*(h+2.*d1))+(1-r)*(-2.*d2+ed2*(h+2.*d2))
      dbr=-r*(1.+h/d1)*ed1-(1-r)*(1.+h/d2)*ed2
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  mdrho===
c cray .endif
      double precision  function mdrho(t1,t2)
c
c         *************************
c         * rho(t2,s2)-rho(t1,s1) *
c         *************************
c
c _dp implicit double precision (a-h,o-z)
      double precision  t1,t2
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      mdrho=(t1-t2)* 0.19494d00
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  rad=====
c cray .endif
      double precision  function rad(h)
c
c         ************************************
c         * r(z) : paulson and simpson, 1977 *
c         ************************************
c
c _dp implicit double precision (a-h,o-z)
      double precision  h,r,d1,d2,x,y
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common /corad/r,d1,d2
c
      x=-h/d1
      if (x.lt.-170.d0) x=-170.d0
      y=-h/d2
      if (y.lt.-170.d0) y=-170.d0
      rad=r*dexp(x)+(1.d0-r)*dexp(y)
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ephil2==
c cray .endif
      subroutine ephil2(hr,ust3,bn,hsl,hbw,cp1,cp3)
c
c         *************************************************
c         * subroutine ephil2(hr,ust3,bn,hsl,hbw,cp1,cp3) *
c         *          computes -h*b'w'(-h) = hbw           *
c         *           following gaspar (1985)             *
c         *************************************************
c
      implicit double precision (a-h,o-z)
      double precision m1th5,m2th5,m3th5,m4th5,m5th5,k1th5,k2th5,k3th5
      common/th5/m1th5,m2th5,m3th5,m4th5,m5th5,a1th5,a2th5,c2th5,c4th5
     .,k1th5,k2th5,k3th5
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      sp=(m2th5+m3th5)*ust3-0.5d0*hr*bn
      ap=cp3*ust3-cp1*hr*bn
      aux=0.5d0*ap-cp1*sp
      hbw=(dsqrt(aux*aux+2.d0*c4th5*hsl*hsl*ap*sp)-0.5d0*ap-cp1*sp)/
     *(c4th5*hsl*hsl-cp1)
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  rphil3==
c cray .endif
      subroutine rphil3(hr,lambda,ust3p,g1,g2,br,hmin)
c
c         ******************************************************
c         *  subroutine rphil3(hr,lambda,ust3p,g1,g2,br,hmin)  *
c         * computes the retreat depth following gaspar (1985) *
c         ******************************************************
c
      implicit double precision (a-h,o-z)
      double precision m1th5,m2th5,m3th5,m4th5,m5th5,k1th5,k2th5,k3th5
      double precision lambda
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/th5/m1th5,m2th5,m3th5,m4th5,m5th5,a1th5,a2th5,c2th5,c4th5
     .,k1th5,k2th5,k3th5
c
      if(g1.lt.0.d0)return
      if(hr.eq.hmin)return
      hp=hr
      fp=g1*hp+g2*br+ust3p*(k3th5-k1th5*k2th5/(ratiol(hp,lambda,ust3p,
     *g1,g2,br)+k2th5))
      if(fp.le.0.d0)return
c
10    continue
      ha=hp-10.d0
      if(ha.lt.hmin)ha=hmin
      call bretd(ha,br,dbr)
      fa=g1*ha+g2*br+ust3p*(k3th5-k1th5*k2th5/(ratiol
     *(ha,lambda,ust3p,g1,g2,br)+k2th5))
      if(fa.le.0.d0)goto 20
c
      if(ha.eq.hmin)then
        hr=hmin
        return
      end if
      hp=ha
      fp=fa
      goto 10
c
20    hold=ha
      iter=0
c
30    iter=iter+1
      h=(fp*ha-fa*hp)/(fp-fa)
      call bretd(h,br,dbr)
      f=g1*h+g2*br+ust3p*(k3th5-k1th5*k2th5/(ratiol(h,lambda,ust3p,
     *g1,g2,br)+k2th5))
      if(dabs(h-hold).lt.0.01d0)then
        hr=h
        return
      else
        if(iter.le.10) then
          hold=h
          if(f.le.0.)then
            ha=h
            fa=f
          else
            hp=h
            fp=f
          end if
        else
          write(iwr6,98)
98        format(3x,'regula does not converge after 10 iter ')
          hr=h
          return
        end if
      end if
      goto 30
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  dileng==
c cray .endif
      double precision  function dileng(hr,lambda,monino)
c
c         *********************************************************
c         *           function dileng(hr,lambda,monino)           *
c         * computes h/l = hsl, where l is the dissipation length *
c         *              following gaspar (1985)                  *
c         *********************************************************
c
      implicit double precision (a-h,o-z)
      double precision m1th5,m2th5,m3th5,m4th5,m5th5,k1th5,k2th5,k3th5
      double precision lambda,monino
      common/th5/m1th5,m2th5,m3th5,m4th5,m5th5,a1th5,a2th5,c2th5,c4th5,
     .k1th5,k2th5,k3th5
c
      aux2=1.d0
      aux=hr/(0.4d0*lambda)
      if(aux.gt.1.d0)aux2=aux
      aux=hr/monino
      if(dabs(aux).gt.25.d0)aux=(monino/dabs(monino))*25.d0
      dileng=a1th5+a2th5*aux2*dexp(aux)
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ratiol==
c cray .endif
      double precision  function ratiol(hr,lambda,ust3,g1,g2,br)
c
c         *****************************************************************
c         *            function ratiol(hr,lambda,ust3,g1,g2,br)           *
c         * computes lp/l = lpsl, where lp is the pressure redistribution *
c         *  length and l the dissipation length, following gaspar (1985) *
c         *****************************************************************
c
      implicit double precision (a-h,o-z)
      double precision m1th5,m2th5,m3th5,m4th5,m5th5,k1th5,k2th5,k3th5
      double precision lambda
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/th5/m1th5,m2th5,m3th5,m4th5,m5th5,a1th5,a2th5,c2th5,c4th5,
     .k1th5,k2th5,k3th5
c
      aux2=hr/(0.4d0*lambda)
      if(aux2.le.1.d0) then
        ratiol=1.d0
      else
        aux=(g1*hr+g2*br)/ust3
        if(dabs(aux).gt.25.d0)aux=(aux/dabs(aux))*25.d0
        aux=dexp(aux)
        ratiol=(a1th5+a2th5*aux2*aux)/(a1th5+a2th5*aux)
      end if
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  difu====
c cray .endif
c
c         ****************************************************
c         *   difu(x,sl,dtd,hr,tm,t,delz,kn,kn1,kmax,kord)   *
c         * computes variations of temperature due to verti- *
c         *           cal diffusion processes                *
c         ****************************************************
c
      subroutine difu(x,sl,dtd,hr,tm,t,delz,kn,kn1,kmax,kord)
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      double precision  dtd
      double precision  hr,tm,t,delz,x,sl
      double precision  b,cdif,hfic,zdo,tdo,cn1,cn2,cn3
      double precision  ex1,ex2,ex3,ex4,ex5,dzz
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      dimension t(nho),x(kord,kord),sl(kord)
c
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .        hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
c
      b=0.5d0
c         if b=1,the numerical scheme is implicit.
c         if b=0,the numerical scheme is explicit.
      cdif=2.d-05
c         cdif is the vertical diffusion coefficient of heat.
      hfic=10.d0
c         hfic is the fictive thickness just below the ml.
      zdo=1.d+6
      tdo=276.d0
c         tdo is the temperature of the deep ocean.
c
      nl=kmax
      cn1=dtd*cdif
      cn2=cn1/(hr*hfic)
      cn3=cn1/hfic
c
c     1. case of ocean nearly mixed.
c     ==============================
c
c     we must resolve the system:
c         * x.t=sl
c
      if (kn.eq.(nl-1)) then
c
c     a) computation of auxiliary variables:
c     --------------------------------------
c
        ex1=1.d0-(1.d0-b)*cn2
        ex2=(1.d0-b)*cn2
        ex3=(1.d0-b)*cn3/delz
        ex4=1.d0-((2.d0*(1.d0-b)*cn1)/((delz+zdo)*delz))
     &  -((1.d0-b)*cn3/delz)
        ex5=2.d0*cn1/((delz+zdo)*delz)
c
c     b)computation of matrix components:
c     -----------------------------------
c
        x(1,1)=1.d0+b*cn2
        x(1,2)=-b*cn2
        x(2,1)=-b*cn3/delz
        x(2,2)=1.d0+((2.d0*b*cn1)/((delz+zdo)*delz))
     &  +(b*cn3/delz)
        sl(1)=ex1*tm+ex2*t(nl)
        sl(2)=ex3*tm+ex4*t(nl)+ex5*tdo
c
c     c) resolution of the system of equations:
c     -----------------------------------------
c
        call simq (x,sl,2,ks)
        tm=sl(1)
        t(nl)=sl(2)
      endif
c
c     2) general case.
c     ================
c
c     we have to resolve the system:
c          * x.t=sl
c
      if (kn.lt.(nl-1)) then
c
c     a) computation of matrix order:
c     -------------------------------
c
        nord=nl-kn+1
c
c     b) computation of matrix components:
c     ------------------------------------
c
        x(1,1)=1.d0+b*cn2
        x(1,2)=-b*cn2
        do 20 jc=3,nord
        x(1,jc)=0.d0
 20     continue
        sl(1)=(1.d0-(1.d0-b)*cn2)*tm+(1.d0-b)*cn2*t(kn1)
c
        ex1=2.d0*cn1/((dzoc2(kn1+1)+delz)*delz)
        ex2=b*cn3/delz
        ex3=(1.d0-b)*cn3/delz
        x(2,1)=-ex2
        x(2,2)=1.d0+b*ex1+ex2
        x(2,3)=-b*ex1
        if (nord.gt.3) then
          do 30 jc=4,nord
          x(2,jc)=0.d0
 30       continue
        endif
        sl(2)=ex3*tm+(1.d0-(1.d0-b)*ex1-ex3)
     &  *t(kn1)+(1.d0-b)*ex1*t(kn1+1)
c
        if (kn.lt.(nl-2)) then
          nordm1=nord-1
          do 40 ic=3,nordm1
          icm2=ic-2
          do 50 jc=1,icm2
          x(ic,jc)=0.d0
 50       continue
          dzz=dzoc2(kn+ic-2)
          if (ic.eq.3) dzz=delz
          ex1=2.d0*cn1/((dzoc2(kn+ic-1)+dzz)*dzoc2(kn+ic-1))
          ex2=2.d0*cn1*((1.d0/(dzoc2(kn+ic-1)+dzoc2(kn+ic)))+
     &    (1.d0/(dzoc2(kn+ic-1)+dzz)))/dzoc2(kn+ic-1)
          ex3=2.d0*cn1/((dzoc2(kn+ic-1)+dzoc2(kn+ic))*dzoc2(kn+ic-1))
          x(ic,icm2+1)=-b*ex1
          x(ic,icm2+2)=1.d0+b*ex2
          x(ic,icm2+3)=-b*ex3
          if (kn.lt.(nl-3).and.ic.ne.(nordm1)) then
            icp2=ic+2
            do 60 jc=icp2,nord
            x(ic,jc)=0.d0
 60         continue
          endif
          sl(ic)=(1.d0-b)*ex1*t(kn+ic-2)+(1.d0-(1.d0-b)*ex2)*
     &    t(kn+ic-1)+(1.d0-b)*ex3*t(kn+ic)
 40       continue
        endif
c
        nordm2=nord-2
        do 70 jc=1,nordm2
        x(nord,jc)=0.d0
 70     continue
        dzz=dzoc2(nl-1)
        if (kn.eq.(nl-2)) dzz=delz
        ex1=2.d0*cn1/((dzoc2(nl)+dzz)*dzoc2(nl))
        ex2=2.d0*cn1*((1.d0/(dzoc2(nl)+zdo))+
     &  (1.d0/(dzoc2(nl)+dzz)))/dzoc2(nl)
        ex3=2.d0*cn1/((dzoc2(nl)+zdo)*dzoc2(nl))
        x(nord,nord-1)=-b*ex1
        x(nord,nord)=1.d0+b*ex2
        sl(nord)=(1.d0-b)*ex1*t(nl-1)+(1.d0-(1.d0-b)*ex2)
     &  *t(nl)+ex3*tdo
c
c     c) resolution of the system of equations:
c     -----------------------------------------
c
        call simq (x,sl,nord,ks)
        tm=sl(1)
        do 90 ic=kn1,nl
        t(ic)=sl(ic-kn1+2)
 90     continue
      endif
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  simq====
c cray .endif
      subroutine simq (a,b,no,ks)
c
c         **************************************************
c         *          subroutine simq (a,b,no,ks)           *
c         * computes the solution of a set of simultaneous *
c         *            linear equations ax=b               *
c         **************************************************
c
      implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      dimension a(1),b(1)
c
c     forward solution
c
      tol=0.0d0
      ks=0
      jj=-no
      do 300 j=1,no
      jy=j+1
      jj=jj+no+1
      biga=0.d0
      it=jj-j
      do 301 i=j,no
c
c     search for maximum coefficient in column
c
      ij=it+i
      if(dabs(biga)-dabs(a(ij))) 302,301,301
  302 biga=a(ij)
      imax=i
  301 continue
c
c     test for pivot less than tolerance (singular matrix)
c
      if (dabs(biga)-tol) 303,303,304
  303 ks=1
      write(iwr6,100)
  100 format (//,1x,'singular matrix in subroutine simq')
      return
c
c     interchange rows if necessary
c
  304  i1=j+no*(j-2)
      it=imax-j
      do 305 k=j,no
      i1=i1+no
      i2=i1+it
      save=a(i1)
      a(i1)=a(i2)
      a(i2)=save
c
c     divide equation by leading coefficient
c
  305 a(i1)=a(i1)/biga
      save=b(imax)
      b(imax)=b(j)
      b(j)=save/biga
c
c     eliminate next variable
c
      if (j-no) 306,307,306
  306 iqs=no*(j-1)
      do 300 ix=jy,no
      ixj=iqs+ix
      it=j-ix
      do 310 jx=jy,no
      ixjx=no*(jx-1)+ix
      jjx=ixjx+it
  310 a(ixjx)=a(ixjx)-(a(ixj)*a(jjx))
  300 b(ix)=b(ix)-(b(j)*a(ixj))
c
c     back solution
c
  307 ny=no-1
      it=no*no
      do 80 j=1,ny
      ia=it-j
      ib=no-j
      ic=no
      do 80 k=1,j
      b(ib)=b(ib)-a(ia)*b(ic)
      ia=ia-no
   80 ic=ic-1
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  fixdp===
c cray .endif
c     ***********************
      subroutine fixdep(ices)
c     ***********************
c
c +++ bcmo16.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
      parameter(nol=0)
      parameter(nro=10)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical fulwri
c
      double precision  fsol,fnsol,fadv
      double precision  enp,enn,flo,den,bio,eice1,eice2
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables flux de chaleur verticaux
c
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables dynamiques
c
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .             th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
c
c --- variables surface
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  acric3,ablic3,amxic3
      common/ic3/acric3(np),ablic3(np),amxic3(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
      common/lat/j
c
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
      double precision  ssto,heqo,fadvo
      common/o0/heqo(np),ssto(np),fadvo(np)
      double precision  qoo1,dto1
      common/o1/qoo1(np),dto1(np)
c
      character*6  route
      common/rou1/ route(np,nro)
      common/rou2/iroute(np,nro,2),airout(np,nro),hirout(np,nro),iiirou
c
c **********************************************************************
c
      do 10 iro = nro,2,-1
       route(j,iro)   =  route(j,iro-1)
      airout(j,iro)   = airout(j,iro-1)
      hirout(j,iro)   = hirout(j,iro-1)
      do 10 k=1,2
      iroute(j,iro,k) = iroute(j,iro-1,k)
 10   continue
       route(j,1  )   = 'fixdep'
      airout(j,1  )   = 1.d0-awic2(j)
      hirout(j,1  )   = hiic2(j)
      iroute(j,1,1)   = it
      iroute(j,1,2)   = iiirou
c
c **********************************************************************
c
      fulwri = .false.
c
      heqo(j)=heqoc2(j)
      ssto(j)=sstoc2(j)
c
c **********************************************************************
c
c     flux solaires en surface
c
      fsol      =   chhv0(1,j,2)
      fsoc2(j)  =   chhv0(1,j,2)
c
c **********************************************************************
c
c     flux non solaires en surface (w.m-2)
c
      fnsol     = ( ablic3(j) + amxic3(j)*awpic2(j)/awic2(j) ) / tic0
      ablic3(j) =   0.
      amxic3(j) =               0.
      do 100 n=2,4
      fnsol     = fnsol + chhv0(n,j,2) - sstoc2(j) * cshv0(n,j,2)
  100 continue
      fnsoc2(j) = fnsol
c
c **********************************************************************
c
c       flux d'advection meridienne calcule a partir des gradients de ss
c       temps delt precedent
c
      fadv      =         chhv0(6,j,2) - sstoc2(j) * cshv0(6,j,2)
      fadvo(j)  = fadv
c     qoo1(j)   = qooc3
      qoo1(j)   = qooc3 - fadv *tic0
c
c **********************************************************************
c
c       determination des nouvelles profondeur et temperature de
c       la couche limite oceanique
c
      sstoc2(j) = ssto(j) + ((fsol +fnsol +fadv)*tic0/ocath2)/hoc2(j)
      do 40 n=1,nhoic4(j)
      twoc2(n,j)=sstoc2(j)
 40   continue
c
c **********************************************************************
c
c       bilan energetique du modele de couche melangee
c
       heqoc2(j) = heqoc2(j) *       ssto(j) /sstoc2(j) 
     .            +hoc2  (j) *(1.d0 -ssto(j) /sstoc2(j))
c
c **********************************************************************
c
c --- correction si sst < t congelation de l'eau de mer (tfth3)
c
      acric3(j)=0.d0
      dawic4   =0.d0
      eice1    =0.d0
      eice2    =0.d0
c
      if (sstoc2(j).le.tfth3) then
        acric3(j) =(sstoc2(j) - tfth3) *ocath2 *hoc2(j)
        eice1     = acric3(j)
        heqoc2(j) = heqoc2(j) *       sstoc2(j) /tfth3
     .             +hoc2  (j) *(1.d0 -sstoc2(j) /tfth3)
        sstoc2(j) = tfth3
        do 41 n=1,nhoic4(j)
        twoc2(n,j)=sstoc2(j)
 41     continue
        hoc2(j)   = dic2(j)
        nhoic4(j) = nhiic4(j)
      else
       if (ices.gt.0.and.nhoic4(j).gt.nhiic4(j)) then
        dawic4 = -qoo1(j) *awic2(j)*(1.d0 -awic2(j+1))
     .         / ( ocath2 *hoc2(j) *(sstoc2(j) -tfth3)
     .           + hfith3          * hcric0           )
        dawic4 = dmax1(dawic4,0.d0)
        twn  =  tfth3 - hcric0 *hfith3 /(hoc2(j) *ocath2)
        tws  = (sstoc2(j) *awic2(j) -twn *dawic4) / (awic2(j) -dawic4)
c
        acric3(j) = ocath2 *hoc2(j) *(twn - tfth3)
        eice2     = acric3(j)
        heqoc2(j) = heqoc2(j) *       sstoc2(j) /tws
     .             +hoc2  (j) *(1.d0 -sstoc2(j) /tws)
        sstoc2(j) = tws
        do 42 n=1,nho-nol
        twoc2(n,j)=sstoc2(j)
 42     continue
       end if
      end if
c
c **********************************************************************
c
c       bilan energetique de la routine
c
       enp    = heqo(j)   * ssto(j)          *ocath2         *awic2(j)
       enn    =(heqoc2(j) * sstoc2(j)        *ocath2 +eice1) *awic2(j)
     .        +(hoc2(j)   *(tfth3 -sstoc2(j))*ocath2 +eice2) *dawic4
c
       flo    = (fsoc2(j) +fnsoc2(j) +fadv)                  *awic2(j)
       den    = (enn - enp) / tic0
       bio    =  flo - den
c
c ifte
       if (.not.fulwri)      go to 50
       if (dabs(bio).gt.0.1) go to 51
       if (iprint.eq.1)      go to 51
       if (mts.ne.0)         go to 50
       if (monpr.gt.mm)      go to 50
c then
   51  continue
       write(iwr1,52)it,j,enn,enp,flo,bio,dawic4
   52  format(/,i5,i3,' abs(',e12.4,' -',e12.4,' ; /tic0 ; -',e12.4,
     .  ') =',e12.4,' >.1 w.m-2    da =',e12.4,'  -- fixdep --')
       write(iwr1,53)
   53  format(8x,' rout.',4x,'it  i0',11x,'ai',11x,'hi')
       write(iwr1,54)
     .(route(j,iro),(iroute(j,iro,k),k=1,2),airout(j,iro),hirout(j,iro),
     .         iro=1,nro)
   54  format((8x,a6,i6,i4,2e13.5))
   50  continue
c
c **********************************************************************
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  oc2=====
c cray .endif
      subroutine inioc2
c
c +++ bcmo2.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
      parameter(nol=0)
c
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c --- constantes physiques
c
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables surface
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
c
      dimension t4a0(9),t4s0(9),t4a(np),t4s(np)
      dimension twa0(9),tws0(9),twa(np),tws(np)
      dimension twa1(9),tws1(9)
      dimension ais0(9),ais(np)
      dimension tw0(9,18),h0(18),x(18),y(np)
c
c ----------------------------------------------------------------------
c
      kmax = nho
c
      data t4a0/26.3,26.3,23.2,15.9, 8.4,  2.2, -5.5,-12.7,-18.0/
      data t4s0/26.0,24.0,18.0, 8.2,-2.4,-10.6,-22.6,-24.8,-33.9/
c ... source : warren & schneider (1979) moyenne annuelle et de janvier
c
      data twa0/27.4,27.9,26.1,22.3,15.8,  8.7,  5.7,  3.4,- 0.7/
      data tws0/27.4,27.1,24.7,19.9,12.8,  6.5,  3.9,  2.5,  0.4/
c ... source : oort (1983) s.s.t. moyenne annuelle et de janvier
c
      data ais0/5*0.0,.085,.445,.825,.995/
c ... source : robock(198 ) fraction ocean couverte de glace
c              (moyenne janvier)
c
c _u28data h0 /75.,81.,94.2,116.7,125.,140.8,12*150./
      data h0 /                              18*150./
      data tw0/
     .300.9,300.9,300.9,300.9,300.9,300.9,300.9,296.6,291.9,
     .300.8,300.8,300.8,300.8,300.8,300.8,300.8,298.8,293.6,
     .300.2,300.2,300.2,300.2,300.2,300.2,300.2,299.7,295.0,
     .299.0,299.0,299.0,299.0,299.0,299.0,299.0,299.0,296.1,
     .297.4,297.4,297.4,297.4,297.4,297.4,297.4,297.4,295.9,
     .295.0,295.0,295.0,295.0,295.0,295.0,295.0,295.0,294.8,
     .291.9,291.9,291.9,291.9,291.9,291.9,291.9,291.9,291.9,
     .288.1,288.1,288.1,288.1,288.1,288.1,288.1,288.1,288.1,
     .284.2,284.2,284.2,284.2,284.2,284.2,284.2,284.2,284.2,
     .281.1,281.1,281.1,281.1,281.1,281.1,281.1,281.1,281.1,
     .278.1,278.1,278.1,278.1,278.1,278.1,278.1,278.1,278.1,
     .275.0,275.0,275.0,275.0,275.0,275.0,275.0,275.0,275.0,
     .271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,
     .271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,
     .271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,
     .271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,
     .271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,
     .271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2,271.2/
c ... source : output modele bom (th.fichefet, oct.1988)
c
      do 1 j=1,9
      t4a0(j) = t4a0(j) + 273.15
      t4s0(j) = t4s0(j) + 273.15
      twa0(j) = twa0(j) + 273.15
      tws0(j) = tws0(j) + 273.15
    1 continue
c
      call pol(t4a0,t4a,9,np)
      call pol(t4s0,t4s,9,np)
c
      do 2 j=1,8
      twa1(j) = .5 * (twa0(j) + twa0(j+1))
      tws1(j) = .5 * (tws0(j) + tws0(j+1))
    2 continue
c
      twa1(9) = .5 * (twa0(9) + tfth3)
      tws1(9) = .5 * (tws0(9) + tfth3)
c
      call pol(twa1,twa,9,np)
      call pol(tws1,tws,9,np)
c
      call pol(ais0,ais,9,np)
c
c ----------------------------------------------------------------------
c
      if (season) then
       do 5 i=1,nho
       do 6 j=1,18
c _u28 x(j) = tw0(i,j)
       x(j) = tw0(1,j)
 6     continue
       call pol(x,y,18,np)
       do 7 j=1,jp
       twoc2(i,j) = y(j)
 7     continue
 5     continue
       call pol(h0,y,18,np)
c _ino
c
c  -- cycle saisonnier
       do 10 j=1,jp
       t4ssu0(j,1) = t4s(j)
       tsiic2(j)   = t4ssu0(j,1)
c
c
                            awic2(j) = 1. - ais(j)
       if (ais(j).lt.1.e-4) awic2(j) = 1.
c
       if (awic2(j).lt.amnic0) then
           awic2(j)   = amnic0
       end if
c
c       awic2(j)   = 1.
c ...  absence de glace marine sur tout l'ocean a l'instant initial
c
       awpic2(j)   = awic2(j)
        ajsu0(j,2) = awic2(j) * ajosu0(j)
        ajsu0(j,1) = ajosu0(j) - ajsu0(j,2)
c
       if (ajsu0(j,2).ge.ajosu0(j)) then
         hoc2(j)   = y(j)
         hr        = y(j)
c
        kn =kmax-nol
        if (hr .eq.zoc2(kmax-nol) ) go to 1100
 1201   kn =kn-1
 1200   if (hr .lt.zoc2(kn)) go to 1201
 1100   kn1=kn+1
        if (kn1.gt.kmax    ) kn1=kmax
        delz = zoc2(kn1)-hr
c
         heqoc2(j) = hr * twoc2(1,j) + delz * twoc2(kn1,j)
        if (kn1.lt.kmax) then
         do 8 i=kn1+1,kmax
         heqoc2(j) = heqoc2(j) + dzoc2(i) * twoc2(i,j)
 8       continue
        end if
         heqoc2(j) = heqoc2(j) / twoc2(1,j)
c
        t4ssu0(j,2) = twoc2(1,j)
c _ino  t4ssu0(j,2) = tws(j)
c
         hiic2(j)   = 0.
        hipic2(j)   = 0.
c
       else
        t4ssu0(j,2) = tfth3
         hiic2(j)   = 3.
        hipic2(j)   = 3.
          hoc2(j)   = dic2(j)
        heqoc2(j)   = zoc2(nho)
        do 11 i=1,nho
        twoc2(i,j)  = tfth3
 11     continue
       end if
         twic2(j)   = t4ssu0(j,2)
        twiic2(j)   = tfth3
        sstoc2(j)   = t4ssu0(j,2)
c
 10    continue
c
      else
c  -- moyenne annuelle
c
       do 202 j=1,jp
       t4ssu0(j,1)=t4a(j)
       t4ssu0(j,2)=twa(j)
 202   continue
c
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tocea===
c cray .endif
c     *****************
      subroutine tocean
c     *****************
c
c +++ bcmo3.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
      parameter(nro=10)
c
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      logical fulwri
c
      double precision  r1p,r2p,r3p,r1n,r2n,r3n,r1d,r2d,r3d,bioj
      double precision  r1nhm,r2nhm,r3nhm,r1phm,r2phm,r3phm,r1,r2,r3
      double precision  flvice,flvsea,flvihm,flvshm,flvtot,biohm
      double precision  remen,remenb
      double precision  htsi,htsw,awl,ail,aiimic,toc,hoc
      dimension htsi(6,np),htsw(6,np),awl(np),toc(np),hoc(np)
      dimension flvice(np),flvsea(np)
      dimension r1p(np),r2p(np),r3p(np)
      dimension r1n(np),r2n(np),r3n(np)
      dimension r1d(np),r2d(np),r3d(np),bioj(np)
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables dynamiques
c
       common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
c
c --- variables flux de chaleur verticaux
c
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
      common/hv1/htshv1(6,np),shshv1(np)
c
c --- variables surface
c
      common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .          ajosu0(np),ajsu0(np,nilw)
      common/su2/hocsu2(np),hicsu2(np),awsu2(np),haisu2(np)
c
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  acric3,ablic3,amxic3
      common/ic3/acric3(np),ablic3(np),amxic3(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
c
      common/oc1/ckwoc1,daoc1(npp),aloc1(npp)
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
      dimension fov(np)
      common/oc6/otoc6(npp),otmoc6(npp),flooc6(npp)
      dimension doc(npp)
c
c --- variables bilan
c
       common/bi0/surbi0(np,nilw),atmbi0(np),topbi0(np)
c
c ----------------------------------------------------------------------
c
      common/iab/iabsur,iabfon
      logical iabsur,iabfon
c
      common/lat/j
      common/rem/remen,remenb
c
c
      character*6  route
      common/rou1/ route(np,nro)
      common/rou2/iroute(np,nro,2),airout(np,nro),hirout(np,nro),iiirou
c
c **********************************************************************
c
c --- initialisation
c
      fulwri = .false.
c
      if (it.eq.0) then
       call inioc0
       call inioc1
       if (.not.inpu) call inioc2
c
c **********************************************************************
c
      else
c
       if (season) then
c
c --- print
c
        if (iprint.eq.1)     go to 601
        if (.not.(mts.eq.0)) go to 600
        if (monpr.gt.mm)     go to 600
 601    continue
        write(iwr1,602)it,ian,mm,iday
 602    format(//'   -- modele glace marine - couche superf. ocean. --',
     .  3x,'it =',i6,3x,i3,'e annee',3x,i2,'e mois',3x,i4,'e jour')
        write(iwr1,603)
 603    format(/3x,'j',5x,'sstoc2',7x,'tw',9(7x),'mld')
        write(iwr1,604)(j,sstoc2(j),(twoc2(i,j),i=1,nho),hoc2(j),j=1,jp)
 604    format((i4,2x,f7.2,9(f7.2),f7.2))
c ...   ne pas oublier de modifier les formats 603 et 604 % nho
c
 600    continue
c
       end if
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
       if (season) then
c --- cycle saisonnier
c
c --- initialisation avant boucle temporelle locale
c
        do 100 j=1,jp
        awl   (j)     =0.
        do 110 n   =1,6
        do 110 iilw=1,2
        htshv0(n,j,iilw)=0.
 110    continue
 100    continue
c
c ----------------------------------------------------------------------
c
        aiimic = dble(iimic0)
c
        do 200 iiic=1,iimic0
        iiirou = iiic
c
c ----------------------------------------------------------------------
c
c --- flux de chaleur meridien oceanique
c
        if (ocen) then
         do 61 j=1,jp
         hoc(j) =           hoc2(j)
         if (nhoic4(j).gt.nhiic4(j))
     .   hoc(j) =                                        dic2(j)
c    .   hoc(j) = awic2(j) *hoc2(j) + (1.d0 - awic2(j)) *dic2(j)
 61      continue
         do 62 j=1,jp
         if (j.eq.1 ) doc(1)   = (3.*hoc(1  )-hoc(2   )) / 2.
         if (j.lt.jp) doc(j+1) = (   hoc(j  )+hoc(j +1)) / 2.
         if (j.eq.jp) doc(jpp) = (3.*hoc(jp )-hoc(jp-1)) / 2.
 62      continue
         call ocean(twic2,doc,iiic)
        end if
c
c ----------------------------------------------------------------------
c
c --- research of ice height and extend for all latitude
c
        do 250 j=1,jp
c
        do 10 iro = nro,2,-1
         route(j,iro)   =  route(j,iro-1)
        airout(j,iro)   = airout(j,iro-1)
        hirout(j,iro)   = hirout(j,iro-1)
        do 10 k=1,2
        iroute(j,iro,k) = iroute(j,iro-1,k)
 10     continue
         route(j,1  )   = 'tocean'
        airout(j,1  )   = 1.d0-awic2(j)
        hirout(j,1  )   = hiic2(j)
        iroute(j,1,1)   = it
        iroute(j,1,2)   = iiirou
c
        twwpl = twic2(j)
c
c  -- calculation of energy of the system before iteration
        r1p(j) = twic2(j) *heqoc2(j) * awic2(j)  *ajosu0(j)
        r2p(j) = ( twiic2(j)*(dic2(j)-c88th3*hiic2(j))
     .           + tfth3             *c88th3*hiic2(j)
     .           + tfth3    *(hoc2(j)  -dic2(j))
     .           + twic2(j) *(heqoc2(j)-hoc2(j))       )
     .                             *(1.-awic2(j)) *ajosu0(j)
        r3p(j) = hiic2(j)          *(1.-awic2(j)) *ajosu0(j)
        fov(j) = 0.
c
c --- calculation of energy balance of ocean
        qooc3   =  0.
        do 260 n=1,6
        qooc3   =  qooc3 + chhv0(n,j,2) - cshv0(n,j,2) * twwpl
 260    continue
        qooc3   =  qooc3 * tic0
c
c --- epaisseur de la glace marine
        if (hiic2(j).eq.0.) goto 254
        call latmix
        fov(j)=fovoc3(j)
        call thickn
        if ((.not.(iabsur)).and.(.not.(iabfon))) goto 255
        call adjvem
        if (hiic2(j).eq.0.) go to 256
 255    continue
c
        if (qooc3.gt.0.) call qzerop
        if (qooc3.le.0.) call qzeron
        if (hiic2(j).eq.0.) goto 256
        goto 258
 254    continue
        call iceocn
        go to 258
 256    continue
        nhoic4(j) = nhiic4(j)
 258    continue
c
c  -- calculation of energy of the system after local iteration
c
        flvice(j) = 0.
        flvsea(j) = 0.
c
        if (tsclim) then
         t41hts = t4ssu0(j,1)
         t42hts = t4ssu0(j,2)
        else
         t41hts = tsiic2(j)
         t42hts = twwpl
        end if
        do 270 n=1,6
        htsi(n,j) = chhv0(n,j,1) - cshv0(n,j,1) * t41hts
        htshv0(n,j,1) = htshv0(n,j,1) + (1.-awpic2(j)) *htsi(n,j)
        flvice(j) = flvice(j) + htsi(n,j)
c
        htsw(n,j) = chhv0(n,j,2) - cshv0(n,j,2) * t42hts
        htshv0(n,j,2) = htshv0(n,j,2) +     awpic2(j)  *htsw(n,j)
        flvsea(j) = flvsea(j) + htsw(n,j)
c
 270    continue
        flvice(j) = flvice(j) * (1.-awpic2(j))*ajosu0(j)
        flvsea(j) = flvsea(j) *     awpic2(j) *ajosu0(j)
c
        r1n(j) =  twic2(j) * heqoc2(j) *    awic2(j)  *ajosu0(j)
        r2n(j) =(twiic2(j) *(dic2(j)-c88th3*hiic2(j))
     .          + tfth3             *c88th3*hiic2(j)
     .          + tfth3    *(hoc2(j)  -dic2(j))
     .          + twic2(j) *(heqoc2(j)-hoc2(j))      )
     .                                 *(1.-awic2(j)) *ajosu0(j)
        r3n(j) =             hiic2(j)  *(1.-awic2(j)) *ajosu0(j)
        r1d(j) = (r1n(j)-r1p(j))*(ocath2/tic0)
        r2d(j) = (r2n(j)-r2p(j))*(ocath2/tic0)
        r3d(j) =-(r3n(j)-r3p(j))*(hfith3/tic0)
        bioj(j)= (r1n(j)+r2n(j)-(r1p(j)+r2p(j)))*(ocath2/tic0)
     .          -(r3n(j)-r3p(j))*(hfith3/tic0) -flvice(j) -flvsea(j)
 250    continue
c
        call dav(r1p,r1phm,cp)
        call dav(r2p,r2phm,cp)
        call dav(r3p,r3phm,cp)
c
        call dav(flvice,flvihm,cp)
        call dav(flvsea,flvshm,cp)
c
        call dav(r1n,r1nhm,cp)
        call dav(r2n,r2nhm,cp)
        call dav(r3n,r3nhm,cp)
c
        r1 = (r1nhm-r1phm)*(ocath2/tic0)
        r2 = (r2nhm-r2phm)*(ocath2/tic0)
        r3 =-(r3nhm-r3phm)*(hfith3/tic0)
c
        flvtot = flvshm + flvihm
c
        biohm  = flvtot - r1 - r2 - r3
c
        if (.not.fulwri)                      go to 650
        if (.not.tsclim.and.abs(biohm).gt..1) go to 651
        if (iprint.eq.1)                      go to 651
        if (.not.(mts.eq.0))                  go to 650
        if (monpr.gt.mm)                      go to 650
 651    continue
        write(iwr1,652)it,iiic,biohm,r1,r2,r3,flvihm,flvshm,r1phm,r1nhm
 652    format(/' iteration',i4,1x,i1,3x,'bilan du modele',13x,':',
     .   e12.4,3x,'r1 =',e12.4,3x,'r2 =',e12.4,3x,'r3 =',e12.4,/,
     .        63x,'ii =',e12.4,3x,'is =',e12.4,/,
     .        63x,'r1p=',e15.7,19x,'r1n=',e15.7)
        write(iwr1,654)
 654    format(/8x,'j',9x,'r1d',9x,'r2d',9x,'r3d',
     .                 6x,'flvsea',6x,'flvice',8x,'bioj')
        write(iwr1,655)(j,r1d(j),r2d(j),r3d(j),
     .                flvsea(j),flvice(j),bioj(j),j=1,jp)
 655    format(i9,6e12.4)
 650    continue
c
        do 251 j=1,jp
c
        awl(j) = awl(j) + awpic2(j)
c ...   aw local
c
        awpic2(j)=awic2(j)
        hipic2(j)=hiic2(j)
 251    continue
c
 200    continue
c
c ----------------------------------------------------------------------
c
c --- assignations apres les boucles temporelles locales
c
        do 300 j=1,jp
c
        t4ssu0(j,1) = tsiic2(j)
        t4ssu0(j,2) =  twic2(j)
         ajsu0(j,1) = ajosu0(j) * (1.-awic2(j))
         ajsu0(j,2) = ajosu0(j) *     awic2(j)
c
        ail = aiimic - awl(j)
c
        if (ail.ne.0.) then
         do 312 n=1,6
         htshv0(n,j,1) = htshv0(n,j,1) / ail
 312     continue
        else
         do 313 n=1,6
         htshv0(n,j,1) = 0.
 313     continue
        end if
c
        do 314 n=1,6
        htshv0(n,j,2)=htshv0(n,j,2) / awl(j)
 314    continue
c
c -- 1)
        do 302 n=1,6
        htshv1(n,j) =htshv1(n,j) +ajosu0(j)*ail   *htshv0(n,j,1) /aiimic
        surbi0(j,1) =surbi0(j,1) +ajosu0(j)*ail   *htshv0(n,j,1) /aiimic
c -- 2)
        htshv1(n,j) =htshv1(n,j) +ajosu0(j)*awl(j)*htshv0(n,j,2) /aiimic
        surbi0(j,2) =surbi0(j,2) +ajosu0(j)*awl(j)*htshv0(n,j,2) /aiimic
 302    continue
c
        hocsu2(j) = hocsu2(j) +  hoc2(j)
        hicsu2(j) = hicsu2(j) + hiic2(j)
         awsu2(j) =  awsu2(j) + awic2(j)
        haisu2(j) = haisu2(j) + awic2(j) *hiic2(j)
c
 300    continue
c
c ----------------------------------------------------------------------
c
c --- assignation des variables generales
c
        if (.not.tsclim) then
         do 322 j=1,jp
c
         if (hiic2(j).ne.0.) then
          t4ssu0(j,1)=tsiic2(j)
         else
          t4ssu0(j,1)= 0.
         end if
          t4ssu0(j,2)= twic2(j)
 322     continue
        end if
c
c ----------------------------------------------------------------------
c
c --- print
c
        if (iprint.eq.1)     go to 631
        if (.not.(mts.eq.0)) go to 630
        if (monpr.gt.mm)     go to 630
 631    continue
c
        write(iwr1,632)
 632    format(/' caracteristiques glace et ocean ',
     .          'en cycle saisonnier :')
c
        write(iwr1,636)
 636    format(/15x,'tsi',9x,'twi',10x,'hi',10x,'aw',
     .           9x,'foc',10x,'tw')
        write(iwr1,637)(j,t4ssu0(j,1),twiic2(j),hiic2(j),awic2(j),
     .                 fov(j),twic2(j),j=1,jp)
 637    format((i6,2f12.2,3f12.3,f12.2))
 630    continue
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
       else
c  --- moyenne annuelle
c
c
c ----------------------------------------------------------------------
c
c  --- flux de chaleur meridien oceanique
c
        if (ocen) then
         do 66 j=1,jp
         toc (j) = t4ssu0(j,2)
 66      continue
         call ocean(toc,daoc1,1)
        end if
c
c ----------------------------------------------------------------------
c
c --- assignation des variables generales
c
        do 400 j=1,jp
        do 401 iilw=1,2
c
        if (ajsu0(j,iilw).ne.0.) then
         scs = 0.
         sch = 0.
         do 402 n=1,6
         scs = scs + cshv0(n,j,iilw)
         sch = sch + chhv0(n,j,iilw)
 402     continue
         if (.not.tsclim) t4ssu0(j,iilw) = sch / scs
c
         do 403 n=1,6
         htshv0(n,j,iilw) =
     .    chhv0(n,j,iilw) - cshv0(n,j,iilw) *t4ssu0(j,iilw)
 403     continue
c
         do 404 n=1,6
         htshv1(n,j)    = htshv1(n,j)   +ajsu0(j,iilw)*htshv0(n,j,iilw)
         surbi0(j,iilw) = surbi0(j,iilw)+ajsu0(j,iilw)*htshv0(n,j,iilw)
 404     continue
        else
         t4ssu0(j,iilw) = 0.
        end if
c
 401    continue
 400    continue
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
       end if
c
      end if
c
c --- print sur fichier bcmio.d
c
      if (.not.mts.eq.0) go to 610
c
       write(3,611)(tsiic2(j),j=1,jp)
 611   format((9f8.2,4x,'tice'))
       write(3,612)(ajsu0(j,1),j=1,jp)
 612   format((9f8.3,4x,'aj i'))
       write(3,613)(htshv0(1,j,1),j=1,jp)
 613   format((9f8.2,4x,'h1 i'))
       write(3,614)(htshv0(2,j,1),j=1,jp)
 614   format((9f8.2,4x,'h2 i'))
       write(3,615)(htshv0(3,j,1),j=1,jp)
 615   format((9f8.2,4x,'h3 i'))
       write(3,616)(htshv0(4,j,1),j=1,jp)
 616   format((9f8.2,4x,'h4 i'))
       write(3,617)(hiic2(j),j=1,jp)
 617   format((9f8.3,4x,'hice'))
c
       write(3,621)(twic2(j),j=1,jp)
 621   format((9f8.2,2x,'twat'))
       write(3,622)(ajsu0(j,2),j=1,jp)
 622   format((9f8.3,2x,'aj w'))
       write(3,618)(hoc2(j),j=1,jp)
 618   format((9f8.2,2x,'hocn'))
       write(3,623)(htshv0(1,j,2),j=1,jp)
 623   format((9f8.2,2x,'h1 w'))
       write(3,624)(htshv0(2,j,2),j=1,jp)
 624   format((9f8.2,2x,'h2 w'))
       write(3,625)(htshv0(3,j,2),j=1,jp)
 625   format((9f8.2,2x,'h3 w'))
       write(3,626)(htshv0(4,j,2),j=1,jp)
 626   format((9f8.2,2x,'h4 w'))
       write(3,627)(ajosu0(j),j=1,jp)
 627   format((9f8.3,4x,'aj o'))
       write(3,628)(flooc6(j),j=1,jp)
 628   format((9e12.4,4x,'flo'))
       write(3,629)ian,mm
 629   format(' cycle annuel no',i3,5x,'mois',i3)
c
 610  continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  thicn===
c cray .endif
c     *****************
      subroutine thickn
c     *****************
c
c +++ bcmo4.for +++
c
c     calculation of surface temperature of ice and determination
c     of the variation of ice thickness due to vertcal heat fluxes
c     at the top and at the bottom of the ice layer
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
      parameter(nro=10)
      logical fulwri
c
c _dp implicit double precision (a-h,o-z)
      double precision  sch,scs,chc,csc,fa,dispen,dspenb
      double precision  r2p,r3p,r2n,r3n,r2d,r3d,enin,fov,rem,denerg
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables flux de chaleur verticaux
c
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables couplage
c
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
c
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
c
      logical iabsur,iabfon
      common/iab/iabsur,iabfon
c
      common/lat/j
      double precision  delhi,delhb
      common/deh/delhi,delhb
      double precision  remen,remenb
      common/rem/remen,remenb
c
      character*6  route
      common/rou1/ route(np,nro)
      common/rou2/iroute(np,nro,2),airout(np,nro),hirout(np,nro),iiirou
c
c **********************************************************************
c
      fulwri = .false.
c
      do 10 iro = nro,2,-1
       route(j,iro)   =  route(j,iro-1)
      airout(j,iro)   = airout(j,iro-1)
      hirout(j,iro)   = hirout(j,iro-1)
      do 10 k=1,2
      iroute(j,iro,k) = iroute(j,iro-1,k)
 10   continue
       route(j,1  )   = 'thickn'
      airout(j,1  )   = 1.d0-awic2(j)
      hirout(j,1  )   = hiic2(j)
      iroute(j,1,1)   = it
      iroute(j,1,2)   = iiirou
c
c ----------------------------------------------------------------------
c
c     energy before computation
c
      r2p = ( twiic2(j) *(dic2(j) -c88th3*hiic2(j))
     .      + tfth3               *c88th3*hiic2(j) )
     .                * (1.-awic2(j))
      r3p =  hiic2(j) * (1.-awic2(j))
      fov =  fovoc3(j)* (1.-awic2(j))
c
c     calculation of the fluxes at the ice surface
c
      sch       = 0.
      scs       = 0.
      do 110 n=1,4
      sch       = sch + chhv0(n,j,1)
      scs       = scs + cshv0(n,j,1)
  110 continue
c
      chc       = (cdith3 / hiic2(j)) *  tfth3
      csc       =  cdith3 / hiic2(j)
      sch       = sch + chc
      scs       = scs + csc
c
c     determination of the new surface temperature
      tsiic2(j) = sch / scs
      fcdoc3    = (cdith3 /hiic2(j)) * (tfth3 -tsiic2(j))
c
      iabsur=.false.
      if (tsiic2(j) .gt.tfith3) go to 101
      remen     = 0.
      go to 104
  101 tsiic2(j) = tfith3
c
      fa        = 0.
      do 111 n=1,4
      fa        = fa  + chhv0(n,j,1) - cshv0(n,j,1) *tsiic2(j)
  111 continue
      fcdoc3    = (cdith3 /hiic2(j)) * (tfth3 -tsiic2(j))
c
c     determination of available energy for ablation at ice surface
      dispen    = (fa+fcdoc3)*tic0
c
c     calculation of the variation of ice thickness due to heat fluxes
c     at the surface
      delhi     = dispen / hfith3
      delhi     =-dmin1(hiic2(j),delhi)
c
c     twiic2 adjustement
c
      twiic2(j) =(twiic2(j) * (dic2(j)-c88th3* hiic2(j))
     .           -tfth3               *c88th3* delhi    )
     .          /             (dic2(j)-c88th3*(hiic2(j)+delhi))
c
      hiic2(j)  = hiic2(j)+delhi
      iabsur    =.true.
c
c     calculation of the remaining energy
      remen     = dispen+delhi*hfith3
      if (delhi.eq.(-hipic2(j))) go to 105
c
c     calculation of ice thickness variation due to vertical heat
c     fluxes at the bottom
  104 continue
      iabfon    =.false.
      dspenb    =(fovoc3(j)-fcdoc3)*tic0
      delhb     = dspenb / hfith3
      if (delhb.lt.0.) go to 103
      iabfon    =.true.
      delhb     = dmin1(hiic2(j),delhb)
      delhb     =-delhb
c
c     twiic2 adjustement
c
      twiic2(j) =(twiic2(j) * (dic2(j)-c88th3* hiic2(j))
     .           -tfth3               *c88th3* delhb    )
     .           /            (dic2(j)-c88th3*(hiic2(j)+delhb))
c
      hiic2(j)  = hiic2(j) + delhb
      remenb    = dspenb+delhb*hfith3
      go to 102
  103 delhb     =-delhb
c
c     twiic2 adjustement
c
      twiic2(j) =(twiic2(j) * (dic2(j)-c88th3* hiic2(j))
     .           -tfth3               *c88th3* delhb    )
     .          /             (dic2(j)-c88th3*(hiic2(j)+delhb))
c
      hiic2(j)  = hiic2(j) + delhb
      remenb    = 0.
  102 continue
      go to 106
  105 continue
      remen = remen + (fovoc3(j)-fcdoc3) *tic0
  106 continue
      fovoc3(j) = 0.
c
c     energy after  computation
c
      r2n =(twiic2(j) * (dic2(j) -c88th3*hiic2(j))
     .    + tfth3                *c88th3*hiic2(j) )
     .                * (1.-awic2(j))
      r3n =  hiic2(j) * (1.-awic2(j))
      r2d = ocath2 *(r2n -r2p) /tic0
      r3d = hfith3 *(r3n -r3p) /tic0
      rem =(remen +remenb) *(1.-awic2(j)) /tic0
c
      enin   = ((sch-chc) -(scs-csc)*tsiic2(j)) *(1.-awic2(j))
      denerg = r2d -r3d -enin -fov +rem
c ... denerg : doit etre nul ; erreur admise : .001w.m-2
c
      if (fulwri.and.dabs(denerg).gt..001) then
       write(iwr1,62)it,j,r2d,r3d,enin,fov,rem,denerg
 62    format(/,i5,i3,' abs(',15x,e12.4,' -',e12.4,' -',e12.4,' -',
     .  e12.4,' +',e12.4,') =',e12.4,' >.001 w.m-2   -- thickn --')
       write(iwr1,63)
 63    format(8x,' rout.',4x,'it  i0',11x,'ai',11x,'hi')
       write(iwr1,64)
     .(route(j,iro),(iroute(j,iro,k),k=1,2),airout(j,iro),hirout(j,iro),
     .         iro=1,nro)
 64    format((8x,a6,i6,i4,2e13.5))
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  adjvm===
c cray .endif
c     *****************
      subroutine adjvem
c     *****************
c
c +++ bcmo5.for +++
c
c     temperature adjustement for vertical melting
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
      parameter(nro=10)
      logical fulwri
c _dp implicit double precision (a-h,o-z)
      double precision  cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      double precision  chhv0,cshv0,htshv0
      double precision  tic0,hcric0,amnic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      double precision  acric3,ablic3,amxic3
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      double precision  delhi,delhb
      double precision  remen,remenb
      double precision  deltwi,hoop,amx
      double precision  r1p,r2p,r3p,r1n,r3n,r1d,r2d,r3d,ren,flo,denerg
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/time2/delt,tht,it
c
c --- constantes physiques
c
       common/th2/cpwth2,rhwth2,ocath2
       common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables flux de chaleur vertical
c
       common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables surface
c
       common/ic0/tic0,hcric0,amnic0,iimic0
       common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .        hiic2(np),hipic2(np),awic2(np),awpic2(np)
       common/ic3/acric3(np),ablic3(np),amxic3(np)
       common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .        hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
c
      common/lat/j
      common/deh/delhi,delhb
      common/rem/remen,remenb
c
c
      character*6  route
      common/rou1/ route(np,nro)
      common/rou2/iroute(np,nro,2),airout(np,nro),hirout(np,nro),iiirou
c
c **********************************************************************
c
      fulwri = .false.
c
      do 10 iro = nro,2,-1
       route(j,iro)   =  route(j,iro-1)
      airout(j,iro)   = airout(j,iro-1)
      hirout(j,iro)   = hirout(j,iro-1)
      do 10 k=1,2
      iroute(j,iro,k) = iroute(j,iro-1,k)
 10   continue
       route(j,1  )   = 'adjvem'
      airout(j,1  )   = 1.d0-awic2(j)
      hirout(j,1  )   = hiic2(j)
      iroute(j,1,1)   = it
      iroute(j,1,2)   = iiirou
c
c ----------------------------------------------------------------------
c
c     energy before computation
      r1p =    twic2(j) * hoc2(j)                *      awic2(j) *ocath2
      r2p = ( twiic2(j) *(dic2(j)-c88th3*hiic2(j))
     .      + tfth3     *(hoc2(j)  -dic2(j))    )*(1.d0-awic2(j))*ocath2
      r3p =   hiic2(j)                           *(1.d0-awic2(j))*hfith3
      hoop= 0.d0
      amx = 0.d0
c
      if (hiic2(j).le.0.d0) then
       deltwi    = (remen +remenb) /(dic2(j)*ocath2)
       twiic2(j) =               twiic2(j) + deltwi
       twiic2(j) =       awic2(j) *  twic2(j)
     .           + (1.d0-awic2(j))*(twiic2(j) *dic2(j) 
     .                        +tfth3 *(hoc2(j)-dic2(j))) /hoc2(j)
c
       ablic3(j) = - (twic2(j)-twiic2(j)) *hoc2(j) *ocath2
       awic2(j)  =  1.d0
c
       do 22 n=1,6
       hoop = hoop + chhv0(n,j,2) - twic2(j) *cshv0(n,j,2)
 22    continue
       ablic3(j) = ablic3(j) - (1.d0-awpic2(j)) *hoop *tic0
       amx       = amxic3(j)                          /tic0
       ices      = 0
       call fixdep(ices)
       twic2(j)  = sstoc2(j)
      end if
c
c     energy after  computation
      r1n =  twic2(j) * hoc2(j)                  *    awic2(j) *ocath2
      r2n =(twiic2(j) *(dic2(j)-c88th3*hiic2(j)) 
     .    + tfth3     *(hoc2(j)  -dic2(j))      )*(1.d0-awic2(j))*ocath2
      r3n =  hiic2(j)                            *(1.d0-awic2(j))*hfith3
      r1d =(r1n -r1p)/tic0
      r2d =(r2n -r2p)/tic0
      r3d =(r3n -r3p)/tic0
      ren =(remen +remenb)*(1.-awpic2(j))/tic0
      flo =(amx   +hoop  )*    awpic2(j)
      denerg = r1d +r2d -r3d -ren -flo
c
      if (fulwri.and.dabs(denerg).gt..1) then
       write(iwr1,62)it,j,r1d,r2d,r3d,ren,flo,denerg
 62    format(/,i5,i3,' abs(',e12.4,' +',e12.4,' -',e12.4,' -',e12.4,
     .     ' -',e12.4,') =',e12.4,' >.1 w.m-2   -- adjvem --')
       write(iwr1,63)
 63    format(8x,' rout.',4x,'it  i0',11x,'ai',11x,'hi')
       write(iwr1,64)
     .(route(j,iro),(iroute(j,iro,k),k=1,2),airout(j,iro),hirout(j,iro),
     .         iro=1,nro)
 64    format((8x,a6,i6,i4,2e13.5))
      end if
c
      remen  = 0.
      remenb = 0.
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  qzerp===
c cray .endif
c     *****************
      subroutine qzerop
c     *****************
c
c +++ bcmo6.for +++
c
c     computation of the variation of ocean temperature and sea-ice
c     extend for positive qo
c
      parameter(np=18,npp=19,nhd=2,nho=9)
      parameter(nro=10)
      logical fulwri
c _dp implicit double precision (a-h,o-z)
      double precision  delta,astar,heat1,heat2,remana
      double precision  r1n,r2n,r3n,r1p,r2p,r3p,r1d,r2d,r3d,denerg
      double precision  abl,amx,awf,awn,deltaf
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/time2/delt,tht,it
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables couplage
c
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  acric3,ablic3,amxic3
      common/ic3/acric3(np),ablic3(np),amxic3(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
c
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
c
      common/lat/j
c
c
      character*6  route
      common/rou1/ route(np,nro)
      common/rou2/iroute(np,nro,2),airout(np,nro),hirout(np,nro),iiirou
c
c **********************************************************************
c
      fulwri = .false.
c
      do 10 iro = nro,2,-1
       route(j,iro)   =  route(j,iro-1)
      airout(j,iro)   = airout(j,iro-1)
      hirout(j,iro)   = hirout(j,iro-1)
      do 10 k=1,2
      iroute(j,iro,k) = iroute(j,iro-1,k)
 10   continue
       route(j,1  )   = 'qzerop'
      airout(j,1  )   = 1.d0-awic2(j)
      hirout(j,1  )   = hiic2(j)
      iroute(j,1,1)   = it
      iroute(j,1,2)   = iiirou
c
c ----------------------------------------------------------------------
c
c     energy before computation
      r1p =  twic2(j) *  heqoc2(j)              *     awic2(j) *ocath2
      r2p =(twiic2(j) * (dic2(j) -c88th3*hiic2(j))
     .     +tfth3     *           c88th3*hiic2(j)
     .     +tfth3     * (hoc2(j)  -dic2(j))
     .     +twic2(j)  * (heqoc2(j)-hoc2(j))       )
     .                                          * (1.-awic2(j))*ocath2
      r3p =  hiic2(j)                           * (1.-awic2(j))*hfith3
      amx = amxic3(j) * awic2(j) /tic0
c
c     calculation of the new ice extend due to lateral melting
      delta  =((1.-awic2(j))*awic2(j)*qooc3)/(hfith3*hiic2(j))
      deltaf =     delta
c
      if (nhoic4(j).gt.nhiic4(j)) then
       if (sstoc2(j-1).gt.sstoc2(j)) then
         awf = -0.5d0 + (sstoc2(j-1) -tfth3) / (sstoc2(j-1) -sstoc2(j))
        if  (awf.lt.awic2(j)) then
                          deltaf = 0.d0
        else
         awn = awic2(j) + delta
         if (awf.lt.awn)  deltaf = awf - awic2(j)
        end if
       end if
      end if
c
      awic2(j)=awic2(j) + deltaf
      astar   =awic2(j)
      if (awic2(j).gt.1.d0) awic2(j)=1.d0
c
c ... correction ablic3 : dillution on new open ocean extend
c
c     calculation of flux heat used to warm lateral cold water
      heat1 =(twiic2(j) *(dic2(j) -c88th3*hiic2(j))
     .       +tfth3               *c88th3*hiic2(j)
     .       +tfth3     *(hoc2(j) -dic2(j))        )
     .                                           *(awic2(j)-awpic2(j))
     .       + twic2(j) * hoc2(j)                *          awpic2(j)
      heat2 =  twic2(j) * hoc2(j)                * awic2(j)
c
c     if a>1 then a=1 and the remaining flux energy goes to ocean
       remana    =          (delta-deltaf) *hiic2(j) *hfith3
      if (astar.gt.1.d0) then
       remana    = remana + (astar-1.d0  ) *hiic2(j) *hfith3
      end if
c
c     correction of the flux available for open ocean
      ablic3(j) = qooc3 * (awpic2(j)*awpic2(j)/awic2(j) -1.)
     .          + ( (heat1 -heat2) *ocath2 + remana ) / awic2(j)
c
c     computation of ocean temperature (due to vertical heat fluxes)
      ices = 0
      call fixdep(ices)
      twic2 (j)=sstoc2(j)
c
      if (astar.ge.1.) then
       twiic2(j)=twic2(j)
       hiic2(j)=0.
      end if
c
c     energy after  computation
      r1n =  twic2(j) *  heqoc2(j)              *     awic2(j) *ocath2
      r2n =(twiic2(j) * (dic2(j) -c88th3*hiic2(j))
     .     + tfth3               *c88th3*hiic2(j)
     .     + tfth3    * (hoc2(j)  -dic2(j))
     .     + twic2(j) * (heqoc2(j)-hoc2(j))       )
     .                                          * (1.-awic2(j))*ocath2
      r3n =  hiic2(j)                           * (1.-awic2(j))*hfith3
      r1d = (r1n - r1p)/tic0
      r2d = (r2n - r2p)/tic0
      r3d = (r3n - r3p)/tic0
      abl = ablic3(j) * awic2(j) /tic0
      hoo = qooc3     *awpic2(j) /tic0  + amx
      denerg = r1d +r2d -r3d -hoo -abl
c
      if (fulwri.and.dabs(denerg).gt..1) then
       write(iwr1,62)it,j,r1d,r2d,r3d,hoo,abl,denerg
 62    format(/,i5,i3,' abs(',e12.4,' +',e12.4,' -',e12.4,' -',e12.4,
     .  ' -',e12.4,') =',e12.4,' >.1 w.m-2','   -- qzerop --')
       write(iwr1,63)
 63    format(8x,' rout.',4x,'it  i0',11x,'ai',11x,'hi')
       write(iwr1,64)
     .(route(j,iro),(iroute(j,iro,k),k=1,2),airout(j,iro),hirout(j,iro),
     .         iro=1,nro)
 64    format((8x,a6,i6,i4,2e13.5))
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  qzern===
c cray .endif
c     *****************
      subroutine qzeron
c     *****************
c
c +++ bcmo7.for +++
c
c     computation of the variation of ocean temperature and sea-ice
c     extend for negative qo
c
      parameter(np=18,npp=19,nhd=2,nho=9)
      parameter(nro=10)
      logical fulwri
c _dp implicit double precision (a-h,o-z)
      double precision  delta,hdelta,hicn,astar,hiacr,twista,hoc
      double precision  r1d,r2d,r3d,r1p,r2p,r3p,r1n,r2n,r3n,fsea,denerg
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/time2/delt,tht,it
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables couplage
c
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  acric3,ablic3,amxic3
      common/ic3/acric3(np),ablic3(np),amxic3(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
c
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
c
      common/lat/j
c
c
      character*6  route
      common/rou1/ route(np,nro)
      common/rou2/iroute(np,nro,2),airout(np,nro),hirout(np,nro),iiirou
c
c **********************************************************************
c
      fulwri = .false.
c
      do 10 iro = nro,2,-1
       route(j,iro)   =  route(j,iro-1)
      airout(j,iro)   = airout(j,iro-1)
      hirout(j,iro)   = hirout(j,iro-1)
      do 10 k=1,2
      iroute(j,iro,k) = iroute(j,iro-1,k)
 10   continue
       route(j,1  )   = 'qzeron'
      airout(j,1  )   = 1.d0-awic2(j)
      hirout(j,1  )   = hiic2(j)
      iroute(j,1,1)   = it
      iroute(j,1,2)   = iiirou
c
c ----------------------------------------------------------------------
c
c     energy before computation
      r1p =  twic2(j) *  hoc2(j)                *       awic2(j) *ocath2
      hoc =              hoc2(j)
c ... hoc :  pour prevoir le passage de hoc2 =150m a hoc2 =30m ds fixdep
c
      r2p =(twiic2(j) * (dic2(j) -c88th3*hiic2(j))
     .     + tfth3               *c88th3*hiic2(j)
     .     + tfth3    * (hoc2(j)  -dic2(j))       )
     .                                          * (1.d0-awic2(j))*ocath2
      r3p =  hiic2(j)                           * (1.d0-awic2(j))*hfith3
      fsea=(qooc3 + amxic3(j))                  *       awic2(j) /tic0
c
c     computation of ocean temperature (due to vertical heat fluxes)
      ices = 1
      call fixdep(ices)
      twic2(j)=sstoc2(j)
c
c     if tw < 271.2 then tw=271.2 and calculation of the new ice extend
c     due to lateral accretion
c
      if (acric3(j).lt.0.d0) then
       if (nhoic4(j).gt.nhiic4(j)) then
         delta     =-dawic4
        hdelta     = hcric0
c ...  cas de la progression du front de glace vers l'equateur sans que 
c      la sst moyenne de la bande ait atteint le point de congelation. 
c
       else
         delta     = awic2(j) * acric3(j) / (hfith3 *hiic2(j))
        hdelta     = hiic2(j)
       end if
         acric3(j) = 0.d0
         awic2(j)  = awic2(j) + delta
         astar     = awic2(j)
         if (awic2(j).lt.amnic0) awic2(j)=amnic0
         hicn =((1.d0-awpic2(j))*hiic2(j) +(awpic2(j) -awic2(j))*hdelta)
     .        /( 1.d0- awic2(j)                                        )
c
c     calculation of the variation of ocean temperature under ice
c     due to lateral accretion
       twiic2(j) =
     .  (twiic2(j) *(dic2(j) -c88th3*hiic2(j)) *(1.d0    -awpic2(j))
     .  + tfth3    *(dic2(j) -c88th3*hdelta  ) *(awpic2(j)-awic2(j)) )
     .           / ( (1.d0-awic2(j)) *(dic2(j) - c88th3*hicn) )
        hiic2(j) = hicn
c
c     if a < amnic0 then a=amnic0 and the deficit of energy cools the
c     water under the ice
       if (astar.gt.amnic0) go to 2
       twiic2(j)=twiic2(j) -( hfith3*hdelta +(twiic2(j)-tfth3) 
     &    *(dic2(j) - c88th3*hdelta)   *ocath2 ) *(awic2(j)-astar)
     & / ( (dic2(j) - c88th3*hiic2(j)) *ocath2   *(1.d0 -awic2(j)) )
 2     continue
c
c     if twi < 271.2 then twi=271.2 and calculation of the new
c     ice thickness
       if (twiic2(j).gt.tfth3) go to 3
       twista=twiic2(j)
       twiic2(j)=tfth3
        hiacr   =
     .  (twiic2(j)-twista) *(dic2(j)-c88th3*hiic2(j)) * (ocath2/hfith3)
c
c     twiic2 adjustement
c
       twiic2(j) =(twiic2(j) *(dic2(j)-c88th3*hiic2(j))
     .            - tfth3             *c88th3*hiacr    )
     .           /            (dic2(j)-c88th3*(hiic2(j)+hiacr))
        hiic2(j) =  hiic2(j)+hiacr
 3     continue
c
      end if
c
c     energy after  computation
      r1n =( twic2(j) *  hoc2(j)               
     .     + tfth3    * (hoc     -hoc2(j)) )   *       awic2(j)  *ocath2
      r2n =(twiic2(j) * (dic2(j) -c88th3*hiic2(j))
     .     + tfth3               *c88th3*hiic2(j)
     .     + tfth3    * (hoc     -dic2(j))       )
     .                                         * (1.d0-awic2(j)) *ocath2
      r3n =  hiic2(j)                          * (1.d0-awic2(j)) *hfith3
      r1d = (r1n - r1p) /tic0
      r2d = (r2n - r2p) /tic0
      r3d = (r3n - r3p) /tic0
c
      denerg = r1d +r2d -r3d - fsea
c ... denerg : bilan energetique : doit etre nul 
c              (erreur permise < .1w.m-2)
c
      if (fulwri.and.dabs(denerg).gt..1) then
       write(iwr1,62)it,j,r1d,r2d,r3d,fsea,denerg,dawic4
 62    format(/,i5,i3,' abs(',e12.4,' +',e12.4,' -',e12.4,' -',e12.4,
     .     ') =',e12.4,' >.1 w.m-2   da =',e12.4,'   -- qzeron --',i2)
       write(iwr1,63)
 63    format(8x,' rout.',4x,'it  i0',11x,'ai',11x,'hi')
       write(iwr1,64)
     .(route(j,iro),(iroute(j,iro,k),k=1,2),airout(j,iro),hirout(j,iro),
     .         iro=1,nro)
 64    format((8x,a6,i6,i4,2e13.5))
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  latmix==
c cray .endif
c     *****************
      subroutine latmix
c     *****************
c
c +++ bcmo8.for
c
c     computation of the variation of water temperatures due to lateral
c     mixing
c     computation of the oceanic heat flux at ice bottom
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      double precision  cdif,twn,rap
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c --- constantes physiques
c
      double precision  cpwth2,rhwth2,ocath2
      common/th2/cpwth2,rhwth2,ocath2
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables couplage
c
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  acric3,ablic3,amxic3
      common/ic3/acric3(np),ablic3(np),amxic3(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
c
      double precision  ssto,heqo,fadvo
      common/o0/heqo(np),ssto(np),fadvo(np)
      double precision  qoo1,dto1
      common/o1/qoo1(np),dto1(np)
      common/lat/j
c
c ----------------------------------------------------------------------
c
      if   (qoo1(j).lt.0.d0.and.nhoic4(j).gt.nhiic4(j)) then
       dto1(j) = twic2(j) - twiic2(j)
       rap     = 0.d0
      else
       rap     = 1.d0
       if  (                    nhoic4(j).le.nhiic4(j)) dto1(j) = 0.d0
      end if
c
      cdif= .25d0 *rap
c ... cdif: facteur de melange
c
      twin =twiic2(j)+cdif*      awic2(j) *( twic2(j)-twiic2(j)-dto1(j))
      if (twin.lt.tfth3) then
       dto1(j) = - (tfth3    -twiic2(j)) /(cdif *awic2(j)) 
     .           + (twic2(j) -twiic2(j))
       twin    =    tfth3
      end if
      twn  = twic2(j)+cdif*(1.d0-awic2(j))*(twiic2(j)- twic2(j)+dto1(j))
     .                         *  (dic2(j)-c88th3*hiic2(j)) /dic2(j)
      amxic3(j)=(twn-twic2(j)) *   dic2(j) *ocath2
      twiic2(j)= twin
c
      if (twiic2(j).gt.tfth3.and.awic2(j).lt.1.) then
       fovoc3(j)=(dic2(j)-c88th3*hiic2(j)) *(twiic2(j)-tfth3)*ocath2
     .          / tic0
       twiic2(j)=tfth3
      else
       fovoc3(j)=0.
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  iceocn==
c cray .endif
c     *****************
      subroutine iceocn
c     *****************
c
c +++ bcmo9.for +++
c
c     this subroutine treats the case of an ice-free ocean
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      double precision  remen,remenb
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c --- constantes physiques
c
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
c --- variables flux de chaleur verticaux
c
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables surface
c
      double precision  tic0,hcric0,amnic0
      common/ic0/tic0,hcric0,amnic0,iimic0
      double precision  tsiic2,twiic2,twic2,dic2
      double precision  hiic2,hipic2,awic2,awpic2
      common/ic2/tsiic2(np),twiic2(np),twic2(np),dic2(np),
     .       hiic2(np),hipic2(np),awic2(np),awpic2(np)
      double precision  acric3,ablic3,amxic3
      common/ic3/acric3(np),ablic3(np),amxic3(np)
      double precision  dawic4
      common/ic4/dawic4,jic4,nhiic4(np),nhoic4(np)
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .       hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
      double precision  qooc3,fovoc3,fcdoc3
      common/oc3/qooc3,fovoc3(np),fcdoc3
c
      common/lat/j
      common/rem/remen,remenb
c
c ----------------------------------------------------------------------
c
c --- calculation of the variation of ocean temperature
c
      call mixdep
      twic2(j) = sstoc2(j)
c
c     if tw < tfth3 then tw=tfth3 and creation of ice
c
      if (acric3(j).lt.0.d0  ) then
       if (sstoc2(j).gt.tfth3) then
        awic2 (j)= 1.d0 - dawic4
        hiic2 (j)= - acric3(j) / hfith3
c ... cas de la progression du front de glace vers l'equateur sans que 
c     la sst moyenne de la bande ait atteint le point de congelation. 
c
       else
        awic2(j) = 1.d0 + ( acric3(j) / (hfith3 *hcric0) )
        hiic2(j) = hcric0
       end if
        acric3(j)= 0.d0
        tsiic2(j)= tfth3
        twiic2(j)= tfth3
      else
       hiic2(j) = 0.
       awic2(j) = 1.
       twiic2(j)= twic2(j)
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ocean===
c cray .endif
      subroutine ocean(toc,doc,loc)
c
c +++ bcmoc.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      double precision  chhv0,cshv0,htshv0
      double precision  cpwth2,rhwth2,ocath2
      double precision  twoc2,sstoc2,fsoc2,fnsoc2,hoc2,heqoc2
      double precision  tic0,hcric0,amnic0
c
      double precision  dscm,pi2r,uabs,c1,choc
      double precision  toc,f2m,ah,ai,xw,akw,dk,rss
      double precision  fl,flv,voc,flx,flxv,flxd
c
      dimension toc(np),tocc(np)
      dimension f2m(np),ah(npp),ai(npp),xw(npp)
      dimension akw(npp),dk(npp),rss(np)
      dimension fl (npp),flv (npp),voc(npp)
      dimension flx(npp),flxv(npp),aflv(npp),bflv(npp)
      dimension flxd(np),          aflx(npp),bflx(npp)
      dimension ax(np),bx(np),cx(np),dx(np)
c
      dimension doc(npp)
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
       common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
       common/gr2/argr2(np),atgr2
c
c --- constantes physiques
c
       common/th2/cpwth2,rhwth2,ocath2
c
c --- variables dynamiques
c
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy1/fjdy1(np),fjjdy1(npp),betdy1(npp)
       common/dy8/u4dy8(npp),u4mdy8(npp),
     .            th4dy8(np),tsody8(npp),tsidy8(npp),tsdy8(npp)
c
c --- variables flux de chaleur
c
       common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
c
c --- variables surface
c
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
c
       common/ic0/tic0,hcric0,amnic0,iimic0
       common/oc1/ckwoc1,daoc1(npp),aloc1(npp)
       common/oc2/twoc2(nho,np),sstoc2(np),fsoc2(np),fnsoc2(np),
     .        hoc2(np),heqoc2(np),zoc2(nho),dzoc2(nho)
       common/oc6/otoc6(npp),otmoc6(npp),flooc6(npp)
c
c ----------------------------------------------------------------------
c
c --- initialisation
c
          dscm =  ca * cp
          pi2r =  ca * pi * 2.
c
         aekmn = .2
c ...    aekmn : fraction oceanique libre minimum sur laquelle circule l
c                courant d'ekman
c
c --- conditions limites
c
      do 10 j=1,jpp,jp
          ah(j) = 0.
          ai(j) = 0.
c
          xw(j) = 0.
          fl(j) = 0.
          flv(j)= 0.
         aflv(j)= 0.
         bflv(j)= 0.
          voc(j)= 0.
          akw(j) = 0.
          flx(j) = 0.
         aflx(j) = 0.
         bflx(j) = 0.
          flxv(j)= 0.
          dk(j) = 0
   10 continue
c
c ----------------------------------------------------------------------
c
c --- flux de chaleur oceanique
c
      do 20 j=1,jp
       f2m(j) = ajsu0(j,1) / ajosu0(j)
   20 continue
c
      do 30 j=2,jp
       ah(j) = 0.5 * (ajsu0(j,2)+ajsu0(j-1,2))
       ai(j) = 0.5 * (f2m  (j-1)+f2m  (j)    )
c
       xw(j) = ah(j)  * pi2r  * c10gr1(j)
c
c --- flux diffusif oceanique
c
       akw(j) = ckwoc1 * aloc1(j) * (1.-ai(j))
       fl (j) = -akw(j)*doc(j)*(toc(j)-toc(j-1))/dscm
c ...  unites :   (m2.k.s-1)
c
       flx (j) = fl(j)  * xw(j) * ocath2
c ...  conversion (m2.k.s-1) * (m) * (j.m-3.k-1) -> (w)
c
c --- courant d'ekman
c
       uabs = abs (u4dy8(j))
c ifte
      if (.not.(               ekm.and.uabs.gt.0.) ) go to 300
c then
                      c1 = 3.728e-5 *  .5 /sqrt(fjdy1(j)*uabs)
      if (uabs.gt.6.) c1 = 1.522e-5 *  .5 /sqrt(fjdy1(j)     )
c ... modification parametrisation =>  .5 /sqrt(fjdy1...
c                       a la place de 1.0
       voc(j) = -c1 * u4dy8(j)
c  ifte
      if (.not.voc(j).gt.0.) go to 320
c  then
       flv(j) = toc(j-1) *voc(j)
      aflv(j) =           voc(j)
      bflv(j) = 0.
      go to 321
c  else
  320 continue
       flv(j) = toc(j  ) *voc(j)
      aflv(j) = 0.
      bflv(j) =           voc(j)
  321 continue
c  endif
c
       flv(j) = flv(j) *doc(j) *exp(-doc(j)/39.)
      aflv(j) =aflv(j) *doc(j) *exp(-doc(j)/39.)
      bflv(j) =bflv(j) *doc(j) *exp(-doc(j)/39.)
       if (ajsu0(j-1,2)/ajosu0(j-1).lt.aekmn) then
         flv(j)=0.
        aflv(j)=0.
        bflv(j)=0.
       end if
       if (ajsu0(j  ,2)/ajosu0(j  ).lt.aekmn) then
         flv(j)=0.
        aflv(j)=0.
        bflv(j)=0.
       end if
      go to 301
c else
  300 continue
       flv(j) = 0.
      aflv(j) = 0.
      bflv(j) = 0.
  301 continue
       flxv(j)   = flv(j) * xw(j) * ocath2
c ...  conversion en (w)
c
c --- flux oceanique total
c
       flooc6(j) = flx(j) + flxv(j)
c
      if (season) then
       otmoc6(j) = otmoc6(j) + flooc6(j)*tic0/delt
      else
       otmoc6(j) = otmoc6(j) + flooc6(j)
      end if
c
   30 continue
c
      do 40 j=1,jp
      flxd(j) = (-flooc6(j)+flooc6(j+1)) / argr2(j)
c ... conversion en (w.m-2)
c ----------------------------------------------------------------------
c
c --- flux de chaleur oceanique transforme par unite surface oceanique
c
          rss(j)  =   dscm * ca * (s10gr1(j+1)-s10gr1(j))
          dk(j)   =   akw(j) * doc(j) * ah(j) * c10gr1(j)
   40 continue
c
      do 62 j=1,jp
      if ( .not. ocen .or. ajsu0(j,2).le.0.)      go to 60
       cshv0(6,j,2)= (ocath2*(dk(j+1)+dk(j))/rss(j))/ajsu0(j,2)
      if (j.eq.1.or.j.eq.jp) go to 61
       chhv0(6,j,2)= (ocath2*(dk(j+1)*toc(j+1)+dk(j)*toc(j-1))/rss(j))
     .             / ajsu0(j,2)
 61   if (j.eq.1) chhv0(6,j,2)=(ocath2*dk(j+1)*toc(j+1)/rss(j))
     .             / ajsu0(j,2)
      if (j.eq.jp)chhv0(6,j,2)=(ocath2*dk(j  )*toc(j-1)/rss(j))
     .             / ajsu0(j,2)
 60    continue
 62   continue
c
c --- contribution du courant d'ekman
c     resolution par un schema numerique semi-implicite
c
      alpha = 0.5
      beta  = 1.0 - alpha
      do 630 j=1,jp
      fflx    = tic0 /(hoc2(j) *ajsu0(j,2) *rss(j) *2.*pi/cp)
      aflx(j) = fflx * xw(j)
      bflx(j) = fflx * xw(j+1)
 630  continue
      do 631 j=1,jp
      cx0     = aflx(j) *aflv(j)
      bx0     = aflx(j) *bflv(j) - bflx(j) *aflv(j+1)
      ax0     =                    bflx(j) *bflv(j+1)
      cx  (j) =   -cx0 *alpha
      bx  (j) =                 -bx0 *alpha                      +1.
      ax  (j) =                              ax0 *alpha
      if  (j.le.1 ) then
        dx  (j) = (              bx0*toc(j) -ax0*toc(j+1)) *beta +toc(j)
      else
       if (j.lt.jp) then
        dx  (j) = (cx0*toc(j-1) +bx0*toc(j) -ax0*toc(j+1)) *beta +toc(j)
       else
        dx  (j) = (cx0*toc(j-1) +bx0*toc(j)              ) *beta +toc(j)
       end if
      end if
      tocc(j) = toc(j)
 631  continue
      logequ  = 1
      call pq(ax,bx,cx,dx,tocc,logequ)
      do 632 j=2,jp
      flxv(j) = (aflv(j)*(alpha*tocc(j-1)+beta*toc(j-1))
     .          +bflv(j)*(alpha*tocc(j)  +beta*toc(j  ))) *xw(j) *ocath2
 632  continue
c
c     write(6,636)
c636  format(/,' BCMOC - RESOLUTION EKMAN VIA SCHEMA SEMI-IMPLICITE',
c    .       /,' **************************************************',
c    .       /,'    j',9x,'toc',8x,'tocc')
c     write(6,637)(j,toc(j),tocc(j),j=1,jp)
c637  format((i5,2f12.3))
c
      do 63 j=1,jp
         choc = 0.
      if (ajsu0(j,2)/ajosu0(j).ge.aekmn)
     .   choc = (flxv(j) -flxv(j+1)) /(ajsu0(j,2) *rss(j) *2.*pi/cp)
         chhv0(6,j,2) = chhv0(6,j,2) + choc
 63   continue
c ----------------------------------------------------------------------
c
c --- print
c
      if (loc.ne.1)    go to 50
      if (iprint.eq.1) go to 51
      if (mts.ne.0.or.monpr.gt.mm) go to 50
   51 continue
      write(iwr1,52)ckwoc1
   52 format (//'   -- ocean --    ckw = ',e12.5,//
     .5x,'j',5x,'t water',8x,'area',3x,'m.l.depth',5x,'k diff.',
     .2x,'fd.(m2.k.s-1)',3x,'fd.(w)',4x,'v0 ekman',3x,'ekman (w)',
     .2x,'fl.tot.(w)',4x,'div.flux')
      write(iwr1,53)(j,toc(j),argr2(j),hoc2(j),akw(j),fl(j),flx(j),
     .              voc(j),flxv(j),flooc6(j),flxd(j),j=1,jp)
   53 format ((i6,f12.2,e12.4,f12.2,7e12.4))
   50 continue
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  rflio===
c cray .endif
      subroutine reflio
c
c +++ bcmor.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables rayonnement
c
       common/so1/sotso1(np),cmuso1(np),bhso1(np,6)
       common/so11/zso11(np),czso11(np)
       common/cl1/clcl1(np),clhcl1,ztcl1(np),zbcl1(np),
     .                             ptcl1(np),pbcl1(np)
c
c --- variables surface
c
       common/rf0/rgsrf0(np,nilw)
       common/rf1/rgwrf1(np,12),rgarf1(np,nilw)
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/su3/hswsu3(np,nilw),cksu3(np,nilw)
c
c **********************************************************************
c
      do 20 j=1,jp
c
c --- glace marine
c
      if (ajsu0(j,1).ne.0.) then
c
c --- w factor
c
           hswsu3(j,1) = .2 + .16*(t4ssu0(j,1)-268.05)
       if (hswsu3(j,1).gt.1. )   hswsu3(j,1)=1.
       if (hswsu3(j,1).lt.0.2)   hswsu3(j,1)=0.2
c
c --- heat diffusivity
c
       cksu3(j,1) = 16.e2
c
c --- albedo
c
       if (season) then
            rgsrf0(j,1)   =0.75-0.035*(t4ssu0(j,1)-263.05)
        if (rgsrf0(j,1).gt.0.75)   rgsrf0(j,1)=0.75
        if (rgsrf0(j,1).lt.0.40)   rgsrf0(j,1)=0.40
c
        if (czso11(j).gt.0.) then
         delta = 0.01766 / czso11(j)  -0.0221
         if (delta.gt..5) delta = .5
         if (delta.lt.0.) delta = 0.
         rgsrf0(j,1) = rgsrf0(j,1) * (1.+delta)
         if (rgsrf0(j,1).gt..9999) rgsrf0(j,1) = .9999
c ...  correction cosinus de la distance zenithale (robock, 1980, mwr, p.274)
c
        end if
       end if
      else
       hswsu3(j,1)=0.
        cksu3(j,1)=0.
       rgsrf0(j,1)=0.
      end if
c
c **********************************************************************
c
c --- ocean
c
c
c --- w factor
c
      hswsu3(j,2) = 1.
c
c --- heat diffusivity
c
       cksu3(j,2) = 1000.e2
c
c --- albedo
c
      rgsrf0(j,2) = .07                   *     clcl1(j)
     .  +(.05 / (1.1*cmuso1(j)**1.4 +.15))* (1.-clcl1(j))
   20 continue
c
c **********************************************************************
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  pol=====
c cray .endif
      subroutine pol(a,b,na,nb)
c
c +++ bcmpl.for +++
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      dimension a(na),b(nb)
      xt  = 1.
      xa0 = xt / (2.*real(na))
      dxa = xt /     real(na)
       xa = xa0
      xb0 = xt / (2.*real(nb))
      dxb = xt /     real(nb)
       xb = xb0
c
      if (nb.ge.na) then
       ja = 0
       jb = 0
c
c
c do until
 110   continue
       ja = ja + 1
       xa = xa0 + dxa*real(ja-1)
       xa1= xa0 + dxa*real(ja)
       da = a(ja+1)-a(ja)
c
c  do until
 1100  continue
       jb = jb + 1
       xb = xb0 + dxb*real(jb-1)
       dx = xb-xa
       b(jb) = a(ja) + dx * da / dxa
       if (.not.(xb.ge.xa1) ) go to 1100
c
       if (.not.(ja.ge.na-1)) go to 110
c
c  do while
 1110  continue
       if (jb.ge.nb) go to 1111
       jb = jb + 1
       xb = xb0 + dxb*real(jb-1)
       dx = xb-xa
       b(jb) = a(ja) + dx * da / dxa
       go to 1110
 1111  continue
c
c
      else
       write(iwr6,6)na,nb
 6     format(//,' mauvaise resolution horizontale : ',
     .          ' na =',i3,3x,'nb =',i3)
      end if
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ra1=====
c cray .endif
      subroutine inira1
c
c +++ bcmr1s.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- constantes physiques
c
       common/th1/cpath1,rath1,akth1
c
c --- variables nuages
c
       common/cl1/clcl1(np),clhcl1,ztcl1(np),zbcl1(np),
     .                             ptcl1(np),pbcl1(np)
       dimension pt0(9),pb0(9),zc0(9),zb0(9)
c
       common/cl2/h2ocl2(np),taucl2(np)
       dimension h2ocl0(9),tauc0(9)
c
       common/cl3/clacl3(np),clvcl3(np,12)
       dimension cl0(9),clouv0(18,12),facbe0(18),facber(np)
       logical london
c
c --- variables vapeur d'eau
c
       common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
       common/ev3/c1ev3(np,12)
c...................
        common/ev4/rhev4(np),rhev5(np,12)
        common/ir8/facir8(np)
c...................
       dimension c10(9),c10v(9,12),c10oa(9)
c      dimension rh00(9)
       dimension rh0(18,nilw)
        dimension clark(9)
c..................
        dimension rhtrov(9,2),rhtro(9)
c..................
       logical c10log
c
c --- variables gaz
c
       common/gz0/uggz0(np),ucgz0(np)
       dimension ug0(9),uc0(9)
       common/gz1/ih1,ih2,ih3,ihc1,ihc2,ihc3,ico2,io3
       common/gz2/ppmgz2
c
c --- variables aerosols
c
        common/as0/tas0(np,nilw,3),oas0(np,nilw,3),gas0(np,nilw,3)
        common/as1/tas1(np,nilw,3,12),oas1(np,nilw,3,12),
     .             gas1(np,nilw,3,12)
        common/as2/tas2(np,nilw,3),oas2(np,nilw,3),gas2(np,nilw,3)
c
       dimension tas100(9,3,12),oas100(9,3,12),gas100(9,3,12)
       dimension tas200(9,3,12),oas200(9,3,12),gas200(9,3,12)
        dimension tas110(9,3)   ,oas110(9,3)   ,gas110(9,3)
        dimension tas210(9,3)   ,oas210(9,3)   ,gas210(9,3)
        dimension oas101(9,3,6),oas102(9,3,6),oas201(9,3,6),
     . oas202(9,3,6),gas101(9,3,6),gas102(9,3,6),gas201(9,3,6),
     . gas202(9,3,6)
c
c --- variables reflectivite
c
       common/rf0/rgsrf0(np,nilw)
       common/rf1/rgwrf1(np,12),rgarf1(np,nilw)
       dimension rsh0(9,3),rsl0(9),rsi0(9),rsh(np,4)
c
       common/fac/facsc,facco2,faccl
c
c --- variables locales
c
       dimension x9(9),x18(18),y(np)
        dimension y1(9),y2(9),y3(9),y4(9),y5(9),y6(9)
        dimension y12(18),y22(18),y32(18),y42(18),y52(18),y62(18)
        dimension yy1(18,3,12),yy2(18,3,12),yy3(18,3,12),
     .            yy4(18,3,12),yy5(18,3,12),yy6(18,3,12)
        dimension zz1(18,3),zz2(18,3),zz3(18,3),
     .            zz4(18,3),zz5(18,3),zz6(18,3)
        dimension rhy1(18),rhy2(9,12),rhy3(9),rhy4(18)
       dimension clouva(108),clouvb(108)
       equivalence (clouv0(1,1),clouva(1))
       equivalence (clouv0(1,7),clouvb(1))
c
c ----------------------------------------------------------------------
c
c
      data  cl0         /.512,.438,.408,.472,.572,.635,.635,.612,.548/
      data london       /.true./
c ... nebulosites : donnees de london
c
c     data pt0/ .48e5,.47e5,.49e5,.49e5,.54e5,.58e5,.57e5,.58e5,.60e5/
c     data pb0/ .62e5,.60e5,.60e5,.62e5,.68e5,.73e5,.72e5,.73e5,.73e5/
c ... reference   : listing    oa78
c
      data pt0/ .58e5,.57e5,.59e5,.59e5,.64e5,.68e5,.67e5,.68e5,.70e5/
      data pb0/ .72e5,.70e5,.70e5,.72e5,.78e5,.83e5,.82e5,.83e5,.83e5/
c ... reference   : oa78 + 100mb (nuages plus bas)
c
       data zc0/6.2,6.3,6.1,5.7,5.0,4.4,4.3,4.2,4.0/
       data zb0/4.1,4.5,4.5,4.0,3.2,2.7,2.6,2.5,2.5/
c ... reference   : article    oa78
c
      data h2ocl0/9*0.08/
      data tauc0/7.6,7.1,6.6,6.5,7.0,7.5,8.0,7.0,6.0/
c     data tauc0/9.0,8.0,7.0,6.5,7.0,7.5,8.0,7.0,6.0/
c     data tauc0/10.,8.0,7.0,6.5,7.0,7.5,8.0,7.0,6.0/
c     data tauc0/10.,8.0,7.0,6.5,7.0,7.5,8.0,7.0,6.5/
c     data tauc0/10.,8.0,7.0,6.5,6.8,7.1,7.2,7.0,6.6/
c _p0 data tauc0/10.,8.0,7.0,6.5,6.8,7.1,7.4,7.5,7.5/
c ... reference   : gallee, de plus en plus ose
c
c     data tauc0/7. ,7. ,7. ,7. ,7. ,7. ,7. ,7. ,7. /
c ... reference   : article    oa78
c
c     data tauc0/7.5,7.0,6.0,5.5,6.0,6.5,7.0,6.5,6.0/
c     data tauc0/7.6,7.1,6.6,6.5,6.8,7.1,7.2,7.0,6.6/
c ... reference   : article    chou & al,1981
c
      data clouva
     ./0.623,0.599,0.543,0.479,0.449,0.454,0.491,0.549,0.599,
     . 0.642,0.669,0.682,0.687,0.662,0.631,0.600,0.552,0.528,
c   0.565     1
     . 0.618,0.599,0.540,0.472,0.444,0.451,0.496,0.562,0.614,
     . 0.645,0.661,0.670,0.668,0.640,0.605,0.577,0.541,0.513,
c   0.561     2
     . 0.612,0.606,0.560,0.488,0.438,0.439,0.486,0.550,0.606,
     . 0.645,0.664,0.674,0.670,0.636,0.599,0.558,0.518,0.577,
c   0.561     3
     . 0.615,0.632,0.596,0.520,0.459,0.444,0.473,0.538,0.607,
     . 0.656,0.687,0.704,0.698,0.661,0.630,0.614,0.584,0.588,
c   0.580     4
     . 0.581,0.622,0.624,0.558,0.481,0.454,0.473,0.524,0.590,
     . 0.651,0.692,0.725,0.738,0.719,0.721,0.746,0.753,0.781,
c   0.603     5
     . 0.562,0.631,0.664,0.619,0.547,0.491,0.471,0.494,0.549,
     . 0.625,0.690,0.729,0.737,0.723,0.727,0.761,0.792,0.899/
c   0.617     6
       data clouvb /
     . 0.547,0.617,0.668,0.658,0.596,0.523,0.476,0.473,0.512,
     . 0.579,0.662,0.726,0.735,0.718,0.727,0.769,0.816,0.926,
c   0.619     7
     . 0.555,0.616,0.669,0.653,0.584,0.510,0.456,0.447,0.479,
     . 0.547,0.638,0.711,0.736,0.736,0.753,0.798,0.841,0.954,
c   0.609     8
     . 0.578,0.612,0.638,0.622,0.558,0.475,0.431,0.431,0.468,
     . 0.539,0.619,0.703,0.759,0.777,0.791,0.818,0.836,0.908,
c   0.597     9
     . 0.584,0.605,0.612,0.568,0.496,0.439,0.430,0.465,0.507,
     . 0.562,0.635,0.716,0.768,0.777,0.776,0.785,0.785,0.795,
c   0.586    10
     . 0.596,0.603,0.581,0.522,0.472,0.447,0.457,0.496,0.550,
     . 0.619,0.684,0.731,0.748,0.725,0.694,0.668,0.634,0.609,
c   0.578    11
     . 0.613,0.603,0.557,0.495,0.467,0.461,0.481,0.536,0.598,
     . 0.649,0.691,0.717,0.720,0.691,0.644,0.591,0.536,0.584/
c   0.575    12
c ... nebulosites : donnees de berlyand et strokina
c
       data facbe0/
     ..880,.850,.820,.850,.880,.885,.890,.895,.900,
     ..895,.890,.890,.890,.900,.910,.895,.880,.820/
c ... facteur correctif des nebulosites de berlyand -> london (moyennes
c
       data ug0/0.25,0.26,0.275,0.30,0.33,0.35,0.36,0.365,0.37/
       data uc0/0.23,0.24,0.255,0.28,0.31,0.33,0.34,0.345,0.35/
c ----------------------------------------------------------------------
c
c ... donnees pour les aerosols
c ... 9 latitudes, ocean (1) ou continent (2), pbl-troposphere-stratosph
c
c ... moyennes mensuelles (0)
c
      data tas100  /
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 ,
     * 9*0.20 , 9*0.025 , 9*0.007 /
c     ./0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007/
c optical depth    aerosol    continental sector
c
      data tas200 /
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 ,
     *  9*0.05 , 9*0.025 , 9*0.007 /
c     ./0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,
c     . 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007/
c stratosphere
c optical depth    aerosol   oceanic sector
c
      data oas101 /
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 /
c     ./0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975/
       data oas102 /
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 ,
     *   9*0.872212 , 9*0.872212 , 9*0.997975 /
c     ./0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975/
c single scattering    aerosol   continental sector
c
      data oas201 /
     * 9*0.982545 , 9*0.872212 , 9*0.997975 ,
     * 9*0.982545 , 9*0.872212 , 9*0.997975 ,
     * 9*0.982545 , 9*0.872212 , 9*0.997975 ,
     * 9*0.982545 , 9*0.872212 , 9*0.997975 ,
     * 9*0.982545 , 9*0.872212 , 9*0.997975 ,
     * 9*0.982545 , 9*0.872212 , 9*0.997975 /
c     ./0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975/
       data oas202 /
     *  9*0.982545 , 9*0.872212 , 9*0.997975 ,
     *  9*0.982545 , 9*0.872212 , 9*0.997975 ,
     *  9*0.982545 , 9*0.872212 , 9*0.997975 ,
     *  9*0.982545 , 9*0.872212 , 9*0.997975 ,
     *  9*0.982545 , 9*0.872212 , 9*0.997975 ,
     *  9*0.982545 , 9*0.872212 , 9*0.997975 /
c     ./0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975,
c     . 0.982545,0.982545,0.982545,0.982545,0.982545,
c     . 0.982545,0.982545,0.982545,0.982545,
c     . 0.872212,0.872212,0.872212,0.872212,0.872212,
c     . 0.872212,0.872212,0.872212,0.872212,
c     . 0.997975,0.997975,0.997975,0.997975,0.997975,
c     . 0.997975,0.997975,0.997975,0.997975/
c single scattering     aerosol    oceanic sector
c
      data gas101 /
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 /
c     ./0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246/
       data gas102 /
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 ,
     * 9*0.647596 , 9*0.647596 , 9*0.624246 /
c     ./0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246/
cc assymetry factor     aerosol    continental sector
c
      data gas201 /
     * 9*0.739002 , 9*0.647596 , 9*0.624246 ,
     * 9*0.739002 , 9*0.647596 , 9*0.624246 ,
     * 9*0.739002 , 9*0.647596 , 9*0.624246 ,
     * 9*0.739002 , 9*0.647596 , 9*0.624246 ,
     * 9*0.739002 , 9*0.647596 , 9*0.624246 ,
     * 9*0.739002 , 9*0.647596 , 9*0.624246 /
c     ./0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246/
       data gas202 /
     *  9*0.739002 , 9*0.647596 , 9*0.624246 ,
     *  9*0.739002 , 9*0.647596 , 9*0.624246 ,
     *  9*0.739002 , 9*0.647596 , 9*0.624246 ,
     *  9*0.739002 , 9*0.647596 , 9*0.624246 ,
     *  9*0.739002 , 9*0.647596 , 9*0.624246 ,
     *  9*0.739002 , 9*0.647596 , 9*0.624246 /
c     ./0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246,
c     . 0.739002,0.739002,0.739002,0.739002,0.739002,
c     . 0.739002,0.739002,0.739002,0.739002,
c     . 0.647596,0.647596,0.647596,0.647596,0.647596,
c     . 0.647596,0.647596,0.647596,0.647596,
c     . 0.624246,0.624246,0.624246,0.624246,0.624246,
c     . 0.624246,0.624246,0.624246,0.624246/
c assymetry factor     aerosol     oceanic sector
c
c ... moyennes annuelles (1)
c
      data tas110
     ./0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,
c pbl
     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c troposphere
     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007/
c stratosphere
c optical depth    aerosol    continental sector
c
      data tas210
     ./0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
c pbl
     . 0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,
c troposphere
     . 0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007/
c stratosphere
c optical depth    aerosol   oceanic sector
c
      data oas110
     ./0.872212,0.872212,0.872212,0.872212,0.872212,
     . 0.872212,0.872212,0.872212,0.872212,
c pbl
     . 0.872212,0.872212,0.872212,0.872212,0.872212,
     . 0.872212,0.872212,0.872212,0.872212,
c troposphere
     . 0.997975,0.997975,0.997975,0.997975,0.997975,
     . 0.997975,0.997975,0.997975,0.997975/
c stratosphere
c single scattering    aerosol   continental sector
c
      data oas210
     ./0.982545,0.982545,0.982545,0.982545,0.982545,
     . 0.982545,0.982545,0.982545,0.982545,
c pbl
     . 0.872212,0.872212,0.872212,0.872212,0.872212,
     . 0.872212,0.872212,0.872212,0.872212,
c troposphere
     . 0.997975,0.997975,0.997975,0.997975,0.997975,
     . 0.997975,0.997975,0.997975,0.997975/
c stratosphere
c single scattering     aerosol    oceanic sector
c
      data gas110
     ./0.647596,0.647596,0.647596,0.647596,0.647596,
     . 0.647596,0.647596,0.647596,0.647596,
c pbl
     . 0.647596,0.647596,0.647596,0.647596,0.647596,
     . 0.647596,0.647596,0.647596,0.647596,
c troposphere
     . 0.624246,0.624246,0.624246,0.624246,0.624246,
     . 0.624246,0.624246,0.624246,0.624246/
c stratosphere
c assymetry factor     aerosol    continental sector
c
      data gas210
     ./0.739002,0.739002,0.739002,0.739002,0.739002,
     . 0.739002,0.739002,0.739002,0.739002,
c pbl
     . 0.647596,0.647596,0.647596,0.647596,0.647596,
     . 0.647596,0.647596,0.647596,0.647596,
c troposphere
     . 0.624246,0.624246,0.624246,0.624246,0.624246,
     . 0.624246,0.624246,0.624246,0.624246/
c stratosphere
c assymetry factor     aerosol     oceanic sector
c
c ---------------------------------------------------------------------
c
      data  c10oa       /2.92,2.91,3.13,3.00,2.78,2.80,2.42,2.04,1.63/
c ...       c10  : profil vertical de vapeur d'eau (cfr data oa78)
c
      data  c10         /3.11,3.53,3.64,3.68,3.16,2.75,2.40,1.88,1.61/
      data  c10v        /3.43,4.12,4.20,4.00,3.66,3.56,2.65,0.94,0.50,
     .                   3.43,4.12,4.20,4.05,3.83,3.47,2.65,0.94,0.50,
     .                   3.34,4.03,3.96,3.85,3.57,3.35,2.94,1.32,0.50,
     .                   3.14,3.75,3.84,3.66,3.34,3.01,2.75,1.43,0.51,
     .                   2.98,3.39,3.56,3.54,2.91,2.75,2.53,1.94,1.76,
     .                   2.96,3.23,3.41,3.54,2.52,2.30,2.20,2.03,1.94,
     .                   2.94,3.18,3.27,3.33,2.64,2.11,2.14,1.91,1.79,
     .                   2.93,3.06,3.18,3.48,2.92,2.40,2.36,2.22,2.10,
     .                   3.00,3.17,3.31,3.65,3.27,2.87,2.69,2.36,2.17,
     .                   2.97,3.31,3.69,3.84,3.49,3.02,2.92,2.33,2.32,
     .                   3.08,3.70,4.06,4.10,3.58,3.08,2.75,1.82,0.94,
     .                   3.26,3.90,4.10,4.10,3.73,3.51,2.94,1.38,0.51/
c ...       c10v : profil vertical humidite specifique calcule d'apres
c ...              oort,83 en utilisant les data a 1000mb et 700mb .
c
c      data  rh00        /.831,.756,.708,.721,.752,.811,.877,.917,.912/
c ... donnees moyennes zonales oa78
c
       data  rh0 /.79,.80,.79,.77,.77,.77,.78,.78,.79,
     .            .81,.82,.83,.83,.82,.80,.81,.81,.80,
c ... glace marine
     .           .79,.80,.79,.77,.77,.77,.78,.78,.79,
     .           .81,.82,.83,.83,.82,.80,.81,.81,.80,
c ... ocean
c    .           .77,.71,.60,.53,.52,.54,.57,.58,.62,
c    .           .78,.75,.77,.78,.79,.80,.80,.80,.80,
c ... continent
c    .           .77,.71,.60,.53,.52,.54,.57,.58,.62,
c    .           .78,.75,.77,.78,.79,.80,.80,.80,.80/
c ... champ de neige continental
c ... these harvey
c
     .           .79,.80,.73,.66,.63,.59,.63,.66,.69,
     .           .72,.75,.77,.78,.79,.80,.80,.80,.80,
c ... continent
     .           .79,.80,.73,.66,.63,.59,.63,.66,.69,
     .           .72,.75,.77,.78,.79,.80,.80,.80,.80,
c ... champ de neige continental
     .           .79,.80,.73,.66,.63,.59,.63,.66,.69,
     .           .72,.75,.77,.78,.79,.80,.80,.80,.80,
c ... calotte groenland
     .           .79,.80,.73,.66,.63,.59,.63,.66,.69,
     .           .72,.75,.77,.78,.79,.80,.80,.80,.80,
c ... calotte laurentide
     .           .79,.80,.73,.66,.63,.59,.63,.66,.69,
     .           .72,.75,.77,.78,.79,.80,.80,.80,.80/
c ... calotte finno-scandinave
c ... correction these harvey :
c         entre  0 et 10 n (id.ocean)
c      et entre 10 et 50 n (calibrage sur oa78 moyenne ocean - continent
c
c....................
      data rhtro /0.54,0.45,0.44,0.48,0.54,
     .            0.56,0.58,0.54,0.51/
      data rhtrov/0.46,0.38,0.42,0.50,0.58,
     .            0.59,0.60,0.55,0.51,
     .            0.58,0.52,0.46,0.46,0.51,
     .            0.56,0.56,0.53,0.51/
      data clark /0.53,0.58,0.62,0.665,0.71,0.745,0.78,0.82,0.85/
c....................
c
      data rsh0/.06,.06,.07,.09,.11,.15,.17,.18,.18,
     .          .06,.06,.07,.08,.09,.12,.13,.14,.14,
     .          .06,.06,.06,.06,.07,.08,.09,.10,.10/
      data rsl0/.08,.13,.18,.16,.15,.16,.16,.16,.16/
      data rsi0/4*0.4,.43,.49,.56,.64,.72/
c
c ----------------------------------------------------------------------
c
       ih1    = 1
       ih2    = 1
       ih3    = 1
       ihc1   = 1
c _rad1
c
       ihc2   = 1
       ihc3   = 1
       ico2   = 1
       io3    = 1
c
c      facsc  = 0.
c      facco2 = 1.
c ...  cfr. bcmctr
       faccl  = 1.
       facaer = 1.
c
c     c3    =0.6
c ... gradient thermique vertical
c
c     co2stp= 1.94e-3
      if (.not.inpu) ppmgz2= 330.*facco2
c
c --- interpolations
c
       call pol(pt0,ptcl1,9,np)
       call pol(pb0,pbcl1,9,np)
       call pol(zc0,ztcl1,9,np)
       call pol(zb0,zbcl1,9,np)
c
       call pol(facbe0,facber,18,np)
       do 1 mm=1,12
       do 2 j =1,18
       x18(j)=clouv0(j,mm)
    2  continue
       call pol(x18,y,18,np)
       do 3 j =1,jp
       clvcl3(j,mm)=y(j)
       if (london) clvcl3(j,mm) = clvcl3(j,mm)*facber(j)
    3  continue
    1  continue
c
       if (london) then
        call pol(cl0,clacl3,9,np)
       else
        do 122 j=1,jp
        clacl3(j)=0.
        do 123 mm=1,12
        clacl3(j)=clacl3(j)+clvcl3(j,mm)
 123    continue
        clacl3(j)=clacl3(j)/12.
 122    continue
       end if
c
       call pol(h2ocl0,h2ocl2,9,np)
       call pol(tauc0, taucl2,9,np)
c
       call pol(ug0,uggz0,9,np)
       call pol(uc0,ucgz0,9,np)
c
c      write(iwr1,128)
  128  format(/,3x,'  aerosols -- moyennes annuelles --  data  --',/,
     .     4x,'lat',4x,'tau',8x,'omega',4x,'    g',//)
c     write (iwr1,127)(j,tas110(j,1),oas110(j,1),gas110(j,1),j=1,9)
c     write (iwr1,126)(j,tas110(j,2),oas110(j,2),gas110(j,2),j=1,9)
c     write (iwr1,125)(j,tas110(j,3),oas110(j,3),gas110(j,3),j=1,9)
c     write (iwr1,127)(j,tas210(j,1),oas210(j,1),gas210(j,1),j=1,9)
c     write (iwr1,126)(j,tas210(j,2),oas210(j,2),gas210(j,2),j=1,9)
c     write (iwr1,125)(j,tas210(j,3),oas210(j,3),gas210(j,3),j=1,9)
  127 format(2x,'it = 1',/,(2x,i3,3x,f9.6,3x,f9.6,3x,f9.6))
  126 format(2x,'it = 2',/,(2x,i3,3x,f9.6,3x,f9.6,3x,f9.6))
  125 format(2x,'it = 3',/,(2x,i3,3x,f9.6,3x,f9.6,3x,f9.6))
c
c
c--  aerosols  et humidite relative tropospherique moyenne ---------
c
c --- moyennes annuelles -------------------------------------------
c
        do 30 it=1,3
        do 31 j =1,9
        y1(j)  = tas110(j,it)*facaer
        y2(j)  = tas210(j,it)*facaer
        y3(j)  = oas110(j,it)
        y4(j)  = oas210(j,it)
        y5(j)  = gas110(j,it)
        y6(j)  = gas210(j,it)
   31   continue
        call pol(y1,y12,9,np)
        call pol(y2,y22,9,np)
        call pol(y3,y32,9,np)
        call pol(y4,y42,9,np)
        call pol(y5,y52,9,np)
        call pol(y6,y62,9,np)
        do 32 j=1,np
        zz1(j,it) = y12(j)
        zz2(j,it) = y22(j)
        zz3(j,it) = y32(j)
        zz4(j,it) = y42(j)
        zz5(j,it) = y52(j)
        zz6(j,it) = y62(j)
   32   continue
   30   continue
        do 33 it=1,3
        do 33 j =1,np
        do 33 iilw=1,nilw
        if (iilw.gt.2) goto 34
        tas2(j,iilw,it) = zz2(j,it)
        oas2(j,iilw,it) = zz4(j,it)
        gas2(j,iilw,it) = zz6(j,it)
        goto 33
   34   tas2(j,iilw,it) = zz1(j,it)
        oas2(j,iilw,it) = zz3(j,it)
        gas2(j,iilw,it) = zz5(j,it)
   33   continue
c..........................................
        call pol (rhtro,rhy1,9,np)
        call pol (clark,facir8,9,np)
c..........................................
c
c --- moyennes mensuelles -------------------------------------------
c
        do 370 j =1,9
        do 370 it=1,3
        do 371 mm=1,6
        oas100(j,it,mm)=oas101(j,it,mm)
        oas200(j,it,mm)=oas201(j,it,mm)
        gas100(j,it,mm)=gas101(j,it,mm)
        gas200(j,it,mm)=gas201(j,it,mm)
  371   continue
        do 372 mm=7,12
        oas100(j,it,mm)=oas102(j,it,mm-6)
        oas200(j,it,mm)=oas202(j,it,mm-6)
        gas100(j,it,mm)=gas102(j,it,mm-6)
        gas200(j,it,mm)=gas202(j,it,mm-6)
  372   continue
  370   continue
c
        do 35 mm=1,12
        do 35 it=1,3
        do 36 j =1,9
        y1(j)  = tas100(j,it,mm)*facaer
        y2(j)  = tas200(j,it,mm)*facaer
        y3(j)  = oas100(j,it,mm)
        y4(j)  = oas200(j,it,mm)
        y5(j)  = gas100(j,it,mm)
        y6(j)  = gas200(j,it,mm)
   36   continue
        call pol(y1,y12,9,np)
        call pol(y2,y22,9,np)
        call pol(y3,y32,9,np)
        call pol(y4,y42,9,np)
        call pol(y5,y52,9,np)
        call pol(y6,y62,9,np)
        do 37 j=1,np
        yy1(j,it,mm) = y12(j)
        yy2(j,it,mm) = y22(j)
        yy3(j,it,mm) = y32(j)
        yy4(j,it,mm) = y42(j)
        yy5(j,it,mm) = y52(j)
        yy6(j,it,mm) = y62(j)
   37   continue
   35   continue
        do 38 mm=1,12
        do 38 it=1,3
        do 38 j =1,np
        do 38 iilw=1,nilw
        if (iilw.gt.2) goto 39
        tas1(j,iilw,it,mm) = yy2(j,it,mm)
        oas1(j,iilw,it,mm) = yy4(j,it,mm)
        gas1(j,iilw,it,mm) = yy6(j,it,mm)
        goto 38
   39  continue
       if (iilw.gt.4) goto 40
       tas1(j,iilw,it,mm) = yy1(j,it,mm)
       oas1(j,iilw,it,mm) = yy3(j,it,mm)
       gas1(j,iilw,it,mm) = yy5(j,it,mm)
       goto 38
   40  continue
       tas1(j,iilw,it,mm) = 0.
        oas1(j,iilw,it,mm) = yy3(j,it,mm)
        gas1(j,iilw,it,mm) = yy5(j,it,mm)
   38   continue
c................................................
        do 41 j=1,9
        fac=(rhtrov(j,2)-rhtrov(j,1))/6
        do 41 mm=1,12
        rhy2(j,mm)=rhtrov(j,2)+fac*(mm-7)
   41   continue
        do 42 mm=1,12
        do 43 j=1,9
        rhy3(j)=rhy2(j,mm)
   43   continue
        call pol(rhy3,rhy4,9,np)
        do 44 j=1,np
        rhev5(j,mm)=rhy4(j)
   44   continue
   42   continue
c.............................................
c
      if (.not.season) then
       do 52   it=1,3
       do 52    j=1,np
       do 53 iilw=1,4
       tas0(j,iilw,it) = tas2(j,iilw,it)
       oas0(j,iilw,it) = oas2(j,iilw,it)
       gas0(j,iilw,it) = gas2(j,iilw,it)
  53   continue
       do 54 iilw=5,nilw
       tas0(j,iilw,it) = 0.
       oas0(j,iilw,it) = oas2(j,iilw,it)
       gas0(j,iilw,it) = gas2(j,iilw,it)
  54   continue
  52   continue
c..........................................
       do 55 j=1,np
       rhev4(j)=rhy1(j)
  55   continue
c...........................................
c
      end if
c
c --------------------------------------------------------------------
c
      c10log = .true.
c
      if (c10log) then
       call pol(c10  ,c1ev1,9,np)
      else
       call pol(c10oa,c1ev1,9,np)
      end if
c
       do 170 mm=1,12
       do 171 j=1,9
       x9(j) = c10v(j,mm)
  171  continue
       call pol(x9,y,9,np)
       do 172 j=1,np
       c1ev3(j,mm) = y(j)
  172  continue
  170  continue
c
       do 173 iilw=1,nilw
       do 174 j=1,18
       x18(j)=rh0(j,iilw)
  174  continue
c      call pol(rh00,y,9,np)
       call pol(x18,y,18,np)
       do 175 j=1,jp
       rhev1(j,iilw)=y(j)
  175  continue
  173  continue
c
c
      do 130 i=1,3
      do 131 j=1,9
      x9(j)=rsh0(j,i)
  131 continue
      call pol(x9,y,9,np)
      do 132 j=1,jp
      rsh(j,i)=y(j)
      if (i.eq.2) rgarf1(j,2)=y(j)
  132 continue
  130 continue
      do 135 j=1,jp
      rsh(j,4)=rsh(j,2)
      do 136 mm=1,12
                   ii  = 1  + (mm-1) /3
                   ii1 = ii + 1
      if (ii.eq.4) ii1 = 1
                   im  = mm - (ii-1) *3 - 1
      rgwrf1(j,mm) = rsh(j,ii) + (im/3.) *(rsh(j,ii1)-rsh(j,ii))
  136 continue
  135 continue
c
      call pol(rsl0,y,9,np)
      do 140 j=1,jp
      rgarf1(j,3)=y(j)
  140 continue
c
      call pol(rsi0,y,9,np)
      do 150 j=1,jp
      rgarf1(j,1)=y(j)
      rgarf1(j,4)=y(j)
       do 151 iilw=4,nilw
       rgarf1(j,iilw)=y(j)
  151  continue
  150 continue
c
c **********************************************************************
c
       do 163 j=1,jp
       do 163 iilw=1,nilw
       rgsrf0(j,iilw)=rgarf1(j,iilw)
 163   continue
c
      if (.not.season) then
       do 162 j=1,jp
       clcl1(j)=clacl3(j)
 162   continue
       call av(clcl1,clhcl1,cp)
      end if
c
c **********************************************************************
c
      write(iwr1,10)london
   10 format(//'   -- inirad --   london = ',l1,//,
     . 6x,'j',11x,'cl',11x,'pt',11x,'pb',11x,
     .            'zt',11x,'zb',8x,'h2ocl',9x,'tauc')
      write(iwr1,11)(j,clacl3(j),ptcl1(j),pbcl1(j),ztcl1(j),zbcl1(j),
     .              h2ocl2(j),taucl2(j),j=1,jp)
   11 format((i7,7(1x,e12.5)))
c
      write(iwr1,12)(mm,mm=1,12)
   12 format(/7x,'clvcl3 (nebulosites, donnees berlyand & strokina)',
     .      /12i8)
      write(iwr1,13)((clvcl3(j,mm),mm=1,12),j=1,jp)
   13 format(12f8.2)
c
c      write(iwr1,28)
c   28 format(/,3x,'  aerosols   --  moyennes annuelles  --',/,
c     .    4x,'lat',4x,'tau',8x,'omega',4x,'    g',//)
c      do 29 ii=1,nilw
c      write(iwr1,27) ii
c   27 format(/,5x,'secteur = ',i3,/)
c      write(iwr1,22)(j,tas0(j,ii,1),oas0(j,ii,1),gas0(j,ii,1),j=1,np)
c      write(iwr1,23)(j,tas0(j,ii,2),oas0(j,ii,2),gas0(j,ii,2),j=1,np)
c      write(iwr1,24)(j,tas0(j,ii,3),oas0(j,ii,3),gas0(j,ii,3),j=1,np)
c   22 format(2x,'it = 1',/,(2x,i3,3x,f9.6,3x,f9.6,3x,f9.6))
c   23 format(2x,'it = 2',/,(2x,i3,3x,f9.6,3x,f9.6,3x,f9.6))
c   24 format(2x,'it = 3',/,(2x,i3,3x,f9.6,3x,f9.6,3x,f9.6))
c   29 continue
c
      write(iwr1,14)(mm,mm=1,12)
   14 format(/7x,'rgwrf1 (reflectivites ocean, donnees oa78)',
     .      /12i8)
      write(iwr1,15)((rgwrf1(j,mm),mm=1,12),j=1,jp)
   15 format(12f8.2)
c
      write(iwr1,16)(i,i=1,4)
   16 format(/6x,'j',5x,'rsi',5x,'rsw',5x,'rsl',5x,'rss',
     .                  '  rsc g.',3x,'laur.','  fn-sc.',
     .        3x,'rsw',i2,3i8)
      write(iwr1,17)(j,(rgarf1(j,iilw),iilw=1,nilw),
     .              (rsh(j,i),i=1,4),j=1,jp)
   17 format((i7,11f8.4))
c
      write(iwr1,20)(iilw,iilw=1,nilw)
   20 format(/6x,'j',10x,'ug',10x,'uc',2x,'rh',7(i8,4x),6x,'c1')
      write(iwr1,21)(j,uggz0(j),ucgz0(j),(rhev1(j,iilw),iilw=1,nilw),
     .              c1ev1(j),j=1,jp)
   21 format((i7,10(e12.4)))
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  ras=====
c cray .endif
      subroutine iniras(dt0,j,iilw,ntrop,nl5)
c
c +++ bcmr2.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9,nir=50)
c _dp implicit double precision (a-h,o-z)
      logical fulwri
c
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c --- constantes physiques
      common/time2/ delt,tht,itime
c
      common/th1/cpath1,rath1,akth1
c
c --- variables dynamiques
c
       common/p00/p4p00,p3p00,p2p00,p1p00
       common/p01/szp01(np,nilw),spp01(np,nilw)
       common/dy0/comdy0,sigdy0,gdy0,f0dy0,q2dy0,crfdy0,betdy0,cedy0
       common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .            pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
c
c --- variables surface
c
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajsu0(np,nilw)
       common/su4/t4su4(np),t4asu4(np,nilw)
c
c --- discretisation atmosphere en 3 couches
c
        common/di1/ it1,it2,it3,ic1,ic2,ic3,inu
        common/di2/ j13(3),j23(3)
c
c --- variables radiatives
c
        common/as0/tas0(np,nilw,3),oas0(np,nilw,3),gas0(np,nilw,3)
       common/so4/pso4(nir),wo3so4(nir)
       common/so5/tauso5(nir),omeso5(nir),cgso5(nir),c1iso5(nir)
       common/so6/taeso6(nir),oaeso6(nir),gaeso6(nir)
c
       common/cl1/clcl1(np),clhcl1,ztcl1(np),zbcl1(np),
     .                             ptcl1(np),pbcl1(np)
c ...  ptcl1 : cloud top pressure (pa)
       common/cl2/h2ocl2(np),taucl2(np)
c
       common/gz0/uggz0(np),ucgz0(np)
       common/ra1/ttra1(nir),hhra1(nir),ozra1(nir),cldra1(nir),
     .            pmbra1(nir),zkmra1(nir)
c ... zkmra1 : altitude (km)
c
c --- variables h2o
c
       common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
c
c --- variables locales
c
      dimension ddp(50)
c
c ----------------------------------------------------------------------
c
      fulwri =.false.
c
      reso   = 100.e+02
c --  reso est un parametre pour choix de la discretisation verticale
c --  reso entre 5000 et 10000 pa
      pbnebu = pbcl1(j)
      ptnebu = ptcl1(j)
c
c --- initialisation des pressions (pa)
c
      p0=1.
c ... pour infra , la pression p0=0. au sommet du modele est prohibee
c
c --- pression tropopause (pa)
c
       ptro=1.e2*(200.-110.*tanh(.05*(t4su0(j)-285.)))
c      ptro=.2e5
c
c --- initialisation des temperatures en surface
c
       tsur=t4ssu0(j,iilw)
c ...  tsur:temperature de la surface (k)
c ----------------------------------------------------------------------
c
c --- calcul de la discretisation verticale pour calculs radiatifs
c
c --- hypothese  :  pression de la base du nuage inferieure a 850 mb
c
c-1  calcul de la pression au sommet de la pbl ,pour calcul profil de t
c
      toppbl  =  spp01(j,iilw) - 15000.*(spp01(j,iilw)/100000.)*
     1           (1. - ((100000.-spp01(j,iilw))/spp01(j,iilw)))
      deppbl  =  spp01(j,iilw) - toppbl
      deppbl  =  amax1(3000.,deppbl)
      toppbl  =  spp01(j,iilw) - deppbl
c --  toppbl=850 pour spp01=1000mb  ----> toppbl=570 pour spp01=600mb
c
c-2  ajustement (si necessaire) de la pression de la base et du sommet
c           de la couche nuageuse effective
      if ((toppbl-pbnebu).lt.1000.)  then
       pbnebi  =  pbnebu
       pbnebu  =  pbnebi - 1000. - (pbnebi - toppbl)
       ptnebu  =  ptnebu - 1000. - (pbnebi - toppbl)
      end if
c
c-3     boundary layer   spp01(j,iilw) ---> toppbl  ( 2 ou 3 couches)
c
      if ((spp01(j,iilw)-toppbl).ge.10000.) then
c-3.1   3 couches entre surface et sommet pbl
       ipbl   = 4
       nc1    = 3
       ddp1   = spp01(j,iilw) - toppbl
       ddp1   = ddp1/nc1
       ddp(1) = ddp1
       ddp(2) = ddp1
       ddp(3) = ddp1
c
c-3.2   :  toppbl -  nuage  : nombre variable de couches
c
       nc2    = int((toppbl-pbnebu)/reso) + 1
       do 72 i=nc1+1,nc1+nc2
       ddp(i) = (toppbl-pbnebu) / nc2
 72    continue
      else
c-3.1   2 couches entre surface et sommet pbl
       ipbl   = 3
       nc1    = 2
       ddp1   = spp01(j,iilw) - toppbl
       ddp1   = ddp1/nc1
       ddp(1) = ddp1
       ddp(2) = ddp1
c
c-3.2    : toppbl -  nuage  : nombre variable de couches
c
       nc2    = int((toppbl-pbnebu)/reso) + 1
       do 73 i=nc1+1,nc1+nc2
       ddp(i) = (toppbl-pbnebu) / nc2
 73    continue
c
      end if
cendif
c
c-4    :   couche nuageuse   (1 couche)
c
       nc3      = 1
       inu      = nc1+nc2+1
       ddp(inu) = pbnebu - ptnebu
c
c-5    :   sommet nuage  -  tropopause  : nombre variable de couches
c
      if (ptnebu.gt.ptro) then
c
       nc4    = int((ptnebu-ptro)/reso) + 1
       do 2 i=inu+1,inu+nc4
       ddp(i) = (ptnebu-ptro) / nc4
 2     continue
      else
c
       ptro  =  ptnebu
       nc4   =  0
c
      end if
c
c-6    :  stratosphere   : nombre variable de couches
c
       nc5    = int((ptro-p0)/reso) + 1
       istra  = inu+nc4
       do 3 i=istra+1,istra+nc5
       ddp(i) = (ptro-p0) / nc5
 3     continue
c
c --- calcul du nombre total de couches ... nc
c         et du nombre total de niveaux ... nl
       nc = nc1 + nc2 + nc3 + nc4 + nc5
       nl = nc  + 1
c
c --- determination des niveaux principaux
c
c ----------     sommet pbl (toppbl)      --------------------
       nl1 = nc1 + 1
c ----------       base du nuage          --------------------
       nl2 = nl1 + nc2
c ----------       sommet du nuage        --------------------
       nl3 = nl2 + nc3
c ----------       tropopause             --------------------
       nl4 = nl3 + nc4
       ntrop = nl4
c ----------       sommet du modele       --------------------
        nl5 = nl4 + nc5
c ----------------------------------------------------------------------
c
c ---  parametres pour atmosphere 3 couches
c
        it1=inu
        it2=it1+1
        it3=nl5
        ic1=it1-1
        ic2=it2-1
        ic3=it3-1
c
        j13(1)=1
        j23(1)=it1
        j13(2)=it1
        j23(2)=it2
        j13(3)=it2
        j23(3)=it3
c -----------------------------------------------------------------
c
       if (mod(j-1,20).eq.0.and.fulwri) then
        write(iwr1,5)j,iilw,spp01(j,iilw),pbnebu,ptnebu,ptro,p0,ipbl,
     .                      toppbl
 5      format(/,'   -- inira2 --',
     .         /' j     = ',i8  ,4x,' iilw  = ',i3,5x,' p4    = ',f8.0,
     .         /' pbcl1 = ',f8.0,4x,' ptcl1 = ',f8.0,
     .         /' ptro  = ',f8.0,4x,' p0    = ',f8.0,
     .         /' ipbl  = ',i8  ,4x,' toppbl= ',f8.1)
        write(iwr1,6)nc,nl,nl1,nl2,nl3,nl4,nl5
 6      format(/,' nc  = ',i3,3x,'nl  = ',i3,
     .  /' nl1 = ',i3,3x,'nl2 = ',i3,3x,'nl3 = ',i3,3x,'nl4 = ',i3,
     .   3x,'nl5 = ',i3)
       end if
c ----------------------------------------------------------------------
c
c --- calcul des pressions aux niveaux ... pmb (pa)
c
       pmbra1(1)=spp01(j,iilw)
       do 10 i=2,nl5
       pmbra1(i)=pmbra1(i-1)-ddp(i-1)
 10    continue
c ----------------------------------------------------------------------
c
c --- calcul des temperatures et des altitudes pour les niveaux ... tt,z
c  -- troposphere
       do 20 i=1,nl4
       call ttzz(pmbra1(i),t2dy2(j),z2,ttra1(i),zz)
       zkmra1(i)=zz    /1000.
c ...   conversion altitudes (m)  -> (km) pour profil ozone
 20    continue
c
       dtsol  =  tsur  -  ttra1(1)
       gamma  =  dtsol / (zkmra1(ipbl)-zkmra1(1))
       t4asu4(j,iilw) = ttra1(1)
c ...  t4asu4  : temperature de l air ramenee de 500 mb a la surface
c
       do 21 i=1,ipbl-1
       ttra1(i) = ttra1(i) + gamma * ( zkmra1(ipbl) - zkmra1(i))
c ...  ttra1(1) : temperature de l'air en surface (= tsur)
c
 21    continue
       dt0    =  tsur - ttra1(1)
c ...  dt0      : difference de tsur et temp. air en surface (= 0 ici)
c
       zztro= zz
c  -- stratosphere
       do 22 i=nl4+1,nl5
       ttra1(i) =ttra1(nl4)
       zz       = zztro
     . - (rath1*ttra1(nl4)/gdy0)*alog(pmbra1(i)/pmbra1(nl4))
       zkmra1(i)=zz/1000.
c ...   conversion altitudes (m)  -> (km) pour profil ozone
 22    continue
c
c --- conversion pressions (pa) -> (mb) pour infra *
c
       do 23 i=1,nl5
       pmbra1(i)=pmbra1(i)/100.
       pso4(i)=pmbra1(i)/pmbra1(1)
 23    continue
c ----------------------------------------------------------------------
c
c --- calcul des nebulosites dans les couches ... cld
c
       do 30 i=1,nc
c ini  omeso5(i)=.9981
       omeso5(i)= 0.9989 - 0.004*exp(-0.15*taucl2(j))
c--  parametrisation  fouquart-bonnel --------------------------------
c ---------------------------------------------------------------------
       cgso5(i)=.85
c ---------
       cldra1(i)=0.
       tauso5(i)=0.
c ---------
       if (i.eq.inu) cldra1(i)=clcl1(j)
       if (i.eq.inu) tauso5(i)=taucl2(j)
   30  continue
c
       r23 = 0.
       do 31 i=1,nc
       l=nl5-i
       r23=1.-(1.-cldra1(l))*(1.-r23)
       c1iso5(l)=r23
c ...  c1iso5   =total cloudiness above the level l  --- random overlap
c
   31  continue
c
c ----------------------------------------------------------------------
c
c --- calcul des humidites specifiques aux niveaux ... hh (g/g)
c
                                   tzon = t4su0(j)
      if (iilw.le.4) then 
       if (tzon.lt.t4asu4(j,iilw)) tzon = t4asu4(j,iilw)
      else
                                   tzon = t4asu4(j,iilw)
      end if
c ... utilisation de tsur -moyenne zonale- pour calculer qs
c     correction :   tair -moyenne zonale- pour les cas d'inversion thermique
c
      if (tzon.ge.273.) esl = 9.4051-2354./tzon
      if (tzon.lt.273.) esl = 10.553-2667./tzon
          es    = 100.* 10.**esl
c ...     es    : tension de vapeur saturante (pa)
          q0ev1(j) = .622*rhzev1(j,iilw)*es/(pmbra1(1)*100.)
          q0ev1(j) =  q0ev1(j)/(1.+q0ev1(j))
          hhra1(1) =  q0ev1(j)
       i = 2
 401  continue
      if (      i .gt.     nl5) go to 400
      if (ttra1(i).le.ttra1(1)) go to 400
       hhra1(i) = hhra1(1)
       i = i + 1
       go to 401
 400   continue
       if (i.le.nl5) then
        inl = i
        inl1= i-1
        do 40 i=inl,nl5
        hhra1(i)=hhra1(inl1)*(pmbra1(i)/pmbra1(inl1))**c1ev1(j)
        if (hhra1(i).lt.3.e-6) hhra1(i)=3.e-6
 40     continue
       end if
c ----------------------------------------------------------------------
c
c --- calcul du profil d'ozone ... oz (cm stp) au-dessus du niveau coura
c     d'apres lacis-hansen 1974
c
       aoz    = uggz0(j)
       boz    = 20.
       coz    =  5.
       zoz1   =(aoz+(aoz*exp(-boz/coz)))
       do 50 i=1,nc
       zoz2   =(zkmra1(i)-boz)/coz
       ozra1(i)= zoz1/(1.+exp(zoz2))
   50  continue
       ozra1(nl5)=  0.
c
       do 51 i=1,nc
       wo3so4(i)=ozra1(i)-ozra1(i+1)
   51  continue
c
c ----  donnees pour aerosols   ---------------------------------------
c
c ------ 1 : pbl                ---------------------------------------
        do 52 i=1,nc1
        oaeso6(i) = oas0  (j,iilw,1)
        gaeso6(i) = gas0  (j,iilw,1)
        tae       = tas0  (j,iilw,1) / nc1
        taeso6(i) = tae*0.730719
        if (iilw.le.2) taeso6(i)=tae*0.912819
   52   continue
c ------ 2 : troposphere ( jusque tropopause)      --------------------
        do 54 i=nl1,nl4-1
        oaeso6(i) = oas0  (j,iilw,2)
        gaeso6(i) = gas0  (j,iilw,2)
        tae       = tas0  (j,iilw,2) / (nl4-nl1)
        taeso6(i) = tae*0.730719
   54   continue
c ------ 3 : stratosphere ( au dessus tropopause )  -------------------
        do 55 i=nl4,nc
        oaeso6(i) = oas0  (j,iilw,3)
        gaeso6(i) = gas0  (j,iilw,3)
        tae       = tas0  (j,iilw,3) / (nc-nl4+1)
        taeso6(i) = tae*0.682188
   55   continue
c ---------------------------------------------------------------------
c
       if (mod(j-1,20).eq.0.and.fulwri) then
        write(iwr1,7)(zkmra1(i),pmbra1(i),ttra1(i),hhra1(i),ozra1(i),
     .             i=1,nl5)
        write(iwr1,8)(cldra1(i),tauso5(i),c1iso5(i),taeso6(i),oaeso6(i),
     .             i=1,nc)
    7   format(/,5x,'zkm',5x,'pmb',7x,'  t  ',7x,' h2o ',7x,'  o3  ',
     .         /,(1x,f7.2,2x,f8.2,4x,f7.3,4x,e9.3,4x,e9.3))
    8   format(/,3x,'nebu',7x,'taun',5x,'overlap',2x,'tauaer',3x,
     .         'omeso6',/,(1x,f6.3,4x,f7.3,4x,f7.4,2x,f8.5,2x,f8.5))
       end if
c ----------------------------------------------------------------------
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  clmrd===
c cray .endif
      subroutine clmrad
c
c +++ bcmr4s.for +++
c
c     routine interpolation temporelle pour les nuages et
c       les profils verticaux de vapeur d'eau et aerosols
c       a partir des moyennes mensuelles
c       pour modele utilise en cycle saisonnier
c
      parameter(np=18,npp=19,nilw=7,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables radiatives
c
       common/cl1/clcl1(np),clhcl1,ztcl1(np),zbcl1(np),
     .                             ptcl1(np),pbcl1(np)
       common/cl3/clacl3(np),clvcl3(np,12)
       common/fac/facsc,facco2,faccl
       common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
       common/ev3/c1ev3(np,12)
        common/ev4/rhev4(np),rhev5(np,12)
       logical c1seas
c
c --- variables radiatives (suite : aerosols)
c
      common/as0/tas0(np,nilw,3),oas0(np,nilw,3),gas0(np,nilw,3)
      common/as1/tas1(np,nilw,3,12),oas1(np,nilw,3,12),
     .           gas1(np,nilw,3,12)
c
c **********************************************************************
c
      c1seas = .true.
c
c --- parametres temporels de l'interpolation
c
      ita=int(12.*tsm)
      itsm=int(tsm)
      itm=mod(it,ita)
      itj=mod(itm,itsm)
      i  =1+itm/itsm
      ip1=mod(i+1,12)
      if(ip1.eq.0) ip1=12
c
c --- interpolations
c
      do 1 j   =1,jp
c
c --- interpolation des aerosols
      do 2 iilw=1,nilw
      do 2 ihd =1,3
      tas0(j,iilw,ihd)=tas1(j,iilw,ihd,i)        +
     .                itj*(tas1(j,iilw,ihd,ip1)-tas1(j,iilw,ihd,i)) /tsm
      oas0(j,iilw,ihd)=oas1(j,iilw,ihd,i)        +
     .                itj*(oas1(j,iilw,ihd,ip1)-oas1(j,iilw,ihd,i)) /tsm
      gas0(j,iilw,ihd)=gas1(j,iilw,ihd,i)        +
     .                itj*(gas1(j,iilw,ihd,ip1)-gas1(j,iilw,ihd,i)) /tsm
    2 continue
c
c --- interpolation des observations des nuages
      clcl1(j)=clvcl3(j,i)+itj*(clvcl3(j,ip1)-clvcl3(j,i))/tsm
      clcl1(j)=clcl1(j)*faccl
c
c --- interpolation des observations des profils de vapeur d'eau
      if (c1seas)
     .c1ev1(j)=c1ev3 (j,i)+itj*(c1ev3 (j,ip1)-c1ev3 (j,i))/tsm
c
c --- interpolation des humidites relatives tropospheriques moyennes
c
      rhev4(j)=rhev5 (j,i)+itj*(rhev5(j,ip1)-rhev5(j,i))/tsm
c
    1 continue
c
      call av(clcl1,clhcl1,cp)
c
c **********************************************************************
c
c --- print
c
      if (iprint.eq.1) go to 601
      if (mts.ne.0) go to 600
      if (monpr.gt.mm) go to 600
  601 continue
      write(iwr1,602)(clcl1(j),j=1,jp)
  602 format(//'   -- clmrad --',//,' cloudiness :',2(/,9f8.3))
      write(iwr1,603)(c1ev1(j),j=1,jp)
  603 format(                     /,' h2o profile:',2(/,9f8.3))
      write(iwr1,604)(tas0(j,1,1),j=1,jp)
  604 format(//'   -- aerosol --',//,' tau : iilw=1,it=1 :',2(/,9f8.3))
      write(iwr1,605)(oas0(j,1,1),j=1,jp)
  605 format(                     /,'  omega             :',2(/,9f8.3))
      write(iwr1,606)(gas0(j,1,1),j=1,jp)
  606 format(                     /,'  g                 :',2(/,9f8.3))
  600 continue
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  evarad==
c cray .endif
      subroutine evarad(j,iilw)
c
c +++ bcmr5.for +++
c
      parameter(np=18,npp=19,nilw=7)
c _dp implicit double precision (a-h,o-z)
      double precision  tfth3,tfith3,hfith3,cdith3,c88th3
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
      common/th3/tfth3,tfith3,hfith3,cdith3,c88th3
c
      common/p00/p4p00,p3p00,p2p00,p1p00
      common/p01/szp01(np,nilw),spp01(np,nilw)
c
       common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
       common/su0/t4su0(np),t4ssu0(np,nilw),t4rsu0(np,nilw),
     .           ajosu0(np),ajvsu0(np,nilw)
       common/su3/hswsu3(np,nilw),cksu3(np,nilw)
c
c      if (t4ssu0(j,iilw).ge.273.) then
c       esl = 9.4051-2354./t4ssu0(j,iilw)
c       es  = 10.**esl
c ...   es  : tension de vapeur saturante (mb)
c
c       e10 = .85 * hswsu3(j,iilw) * (t4ssu0(j,iilw) - 273.) +3.35
c ...   e10 : tension de vapeur reelle    (mb)
c             (saltzman et ashe, 1976, tellus, p310)
c
c                                 rhzev1(j,iilw) = e10/es
c       if (rhzev1(j,iilw).gt.1.) rhzev1(j,iilw) = 1.
c
c      else
c       esl = 10.553-2667./t4ssu0(j,iilw)
c       es  = 10.**esl
c ...   es  : tension de vapeur saturante (mb)
c
c       rhzev1(j,iilw) = rhev1(j,iilw) *(spp01(j,iilw) /p4p00 -.02) /.98
c       e10 = rhzev1(j,iilw) *es
c       if (t4ssu0(j,iilw).lt.tfth3) hswsu3(j,iilw) = 0.2
c      end if
c
c     rhzev1(j,iilw) = rhev1(j,iilw) *(spp01(j,iilw) /p4p00 +1.0) /2.0
      rhzev1(j,iilw) = rhev1(j,iilw) *(spp01(j,iilw) /p4p00 -.02) /.98
c     rhzev1(j,iilw) = rhev1(j,iilw)
c
       return
       end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  iniso===
c cray .endif
      subroutine iniso (facco2)
c
c +++ bcms0.for +++
c
c **********************************************************************
c
c ***     solar radiation data     ****
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common /eti200/ uinf(3),aap(3,6),bbp(3,6)
c
c mVax open(unit=15,status='old',file='bcms01.dat')
       open(unit=15,status='old',file='bcms01')
       rewind 15
c
      if (facco2.ne.2.) then
       read (15,1000) uinf(1),uinf(2),uinf(3)
c      write(iwr6,1000) uinf(1),uinf(2),uinf(3)
       read (15,1001) ((aap(i,j),j=1,6),i=1,3)
c      write(iwr6,1001) ((aap(i,j),j=1,6),i=1,3)
       read (15,1001) ((bbp(i,j),j=1,6),i=1,3)
c      write(iwr6,1001) ((bbp(i,j),j=1,6),i=1,3)
 1000  format (3(f6.3))
 1001  format (e15.9)
c
       close(unit=15)
c
      else
c
c mVax open(unit=15,status='old',file='bcms02.dat')
       open(unit=15,status='old',file='bcms02')
       rewind 15
c
       read (15,1000) uinf(1),uinf(2),uinf(3)
c      write(iwr6,1000) uinf(1),uinf(2),uinf(3)
       read (15,1001) ((aap(i,j),j=1,6),i=1,3)
c      write(iwr6,1001) ((aap(i,j),j=1,6),i=1,3)
       read (15,1001) ((bbp(i,j),j=1,6),i=1,3)
c      write(iwr6,1001) ((bbp(i,j),j=1,6),i=1,3)
c
       close(unit=15)
c
      end if
c
      write(iwr6,2000)
2000  format (///,'    end of the solar radiative data lecture',//)
c
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  insol===
c cray .endif
      subroutine insol
c
c +++ bcms2.for +++
c
      parameter(np=18,npp=19,nhd=2,nho=9)
c _dp implicit double precision (a-h,o-z)
      double precision annee
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/pal/apal,fpal,ipal,ipaleo
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
c
c --- variables radiatives
c
       common/fac/facsc,facco2,faccl
       common/so1/sotso1(np),cmuso1(np),bhso1(np,6)
       common/so11/zso11(np),czso11(np)
c
       dimension phi(np),h(np),f(np)
       dimension fa(np),cmua(np),fa1(np),bha(np,6)
c
c ----------------------------------------------------------------------
c
c --- solar input
c
      sc = 1368. + 1368. * 0.01 * facsc
c
      annee = apal
c pour exp. co2
c     annee = real(ipal)
c pour exp. co2
c ... annee = epoque consideree   en milliers d'annees
c ...         annee est negatif pour passe (-1.-2,etc...)
c ...         annee=0.   pour present
c
      if (it.eq.0) call celest(dble(annee))
c
c ----------------------------------------------------------------------
c
      pi2   =  pi  / 2.
      delphi=  pi2 / jp
c
      if (.not.season.and.it.eq.0) then
       iday1=1
       iday2=360
c
c  -- initialisation des variables auxiliaires
       do 102 j=1,jp
       fa(j)  =0.
       cmua(j)=0.
       fa1(j) =0.
       do 103 nh=1,6
       bha(j,nh)=0.
 103   continue
 102   continue
      else
       iday1=iday
       iday2=iday
      end if
c
      if (season.or.it.eq.0) then
c
c --- seasonal cycle
c
       do 202 idayi=iday1,iday2
       do 203 j=1,jp
       phi(j) = delphi * (j-.5)
c
       call solext(sc,phi(j),h(j),del,idayi,j)
c
       zso11(j)    = abs(phi(j) - del)
      czso11(j)    = cos(zso11(j))
       zso11(j)    = 90. - zso11(j) * 180./pi
c...   zso11       : noon solar elevation angle
c
       bhso1(j,1 ) =      sc * cos(phi(j)) * cos(del)
     .  * ( h(j) - sin(h(j)) *cos(h(j)) ) / pi
       do 250 nh=2,6
       bhso1(j,nh) = 2. * sc * cos(phi(j)) * cos(del)
     .  * ( sin(nh*h(j)) *cos(h(j)) - nh *sin(h(j)) *cos(nh*h(j)))
     .  / ( pi * nh * ( nh*nh - 1.))
  250  continue
c
         f(j)      = abs(h(j))/pi
c
         phi(j)    = phi(j) * 180./pi
c
c  -- incrementation des variables auxiliaires
c
       if (.not.season.and.it.eq.0) then
        fa(j)  =fa(j) +f(j)
        fa1(j) =fa1(j)+1.
        cmua(j)=cmua(j)+cmuso1(j)*f(j)
        do 212 nh=1,6
        bha(j,nh)=bha(j,nh)+(bhso1(j,nh)**2)*f(j)
 212    continue
       end if
c
 203   continue
 202   continue
c
       if (.not.season.and.it.eq.0) then
        do 272 j=1,jp
        cmuso1(j)=cmua(j)/fa(j)
        do 273 nh=1,6
        bhso1(j,nh)=sqrt(bha(j,nh)/fa1(j))
 273    continue
 272    continue
c
        write(iwr1,280)(nh,nh=1,6)
 280    format(//,'   -- test moyenne annuelle insol --',/,6x,27('#'),
     .     /,5x,'cosmu',5x,'bh',i3,5i10)
        write(iwr1,281)(cmuso1(j),(bhso1(j,nh),nh=1,6),j=1,jp)
 281    format(18(f10.4,6f10.1,/))
       end if
      end if
c ----------------------------------------------------------------------
c
      if (.not. season) then
c
c --- annual mean
c
       do 302 j=1,jp
       sotso1(j)=.5*sc*cmuso1(j)
 302   continue
      end if
c
c --- print
c
      if (iprint.eq.1) go to 96
      if (mts.ne.0) go to 97
      if (monpr.gt.mm) go to 97
   96 continue
      write(iwr1,98) it,ian,mm,iday,(nh,nh=1,6)
   98 format(//'  -- insol  --',
     . 3x,'it =',i6,3x,i3,'e annee',3x,i2,'e mois',3x,i4,'e jour',
     . //,5x,'j',5x,'phi',2x,'sotso1',2x,'cmuso1',2x,'bh',i4,5i8)
      write(iwr1,99) (j,phi(j),sotso1(j),cmuso1(j),(bhso1(j,nh),nh=1,6),
     .             j=1,jp)
   99 format((i6,f8.3,f8.2,f8.4,6f8.2))
   97 continue
      return
      end
c **********************************************************************
c cray .if [ $concord .eq. 1 ]
c cray nzyx  celst===
c cray .endif
      subroutine celest(annee)
      implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
       common /so10/ perh,ecc,so
c
      dimension ae(19),be(19),ce(19),aob(104),
     . bob(104),cob(104),aop(177),bop(177),cop(177)
c     dimension ap(117),bp(117),cp(117)
c
c   constant
c
      pi=3.14159265358979d0
      pir=pi/180.0d0
      pirr=pir/3600.0d0
c
c   1.earth orbital elements : eccentricity           ecc
c **************************   precessional parameter pre
c                              obliquity              xob
c
c         read amplitude a  mean rate b  phase c
c              they are immediately converted in radians
c
c mVax open(unit=16,status='old',file='bcms2.dat')
       open(unit=16,status='old',file='bcms2')
       rewind  16
c
c   eccentricity
c
      nef=19
      do 1 i=1,nef
      read(16,5000) ae(i),y,z
c     write(iwr1,5000) ae(i),y,z
 5000 format (13x,f11.8,f20.7,f20.6)
      be(i)=y*pirr
      ce(i)=z*pir
    1 continue
c
c   precessional parameter
c
      nem=117
      do 2 i=1,nem
      read(16,5001) 
c     write(iwr1,5111)i 
 5001 format(52x)
 5111 format(2x,i4)
c ... cfr. discussion avec A.B. le 3-fev-1989 
c          (le calcul ci-dessous est imprecis)
c     read(16,5001) ap(i),y,z
c5001 format(10x,f10.7,3x,f13.6,3x,f13.6)
c     bp(i)=y*pirr
c     cp(i)=z*pir
    2 continue
c
c   obliquity
c
      xod=23.320556d0
      nob=104
      do 3 i=1,nob
      read(16,5002) aob(i),y,z
c     write(iwr1,5002) aob(i),y,z
 5002 format(7x,f13.7,2x,f10.6,2x,f10.4)
      bob(i)=y*pirr
      cob(i)=z*pir
    3 continue
c
c   general precession in longitude
c
      xop=3.392506d0
      prm=50.439273d0
      nop=177
      do 31 i=1,177
      read(16,5002) aop(i),y,z
c     write(iwr1,5002) aop(i),y,z
      bop(i)=y*pirr
      cop(i)=z*pir
   31 continue
c
      close(unit=16)
c
      nef=19
      nob=104
      nop=177
      write(iwr1,6900) nef,nob,nop
 6900 format(//,1x,'long term daily insolation',/,1x,
     *'number of terms in',/,6x,'eccentricity',i5,2x,
     *'obliquity',i5,2x,'general precession',i5,//)
c
      t=annee
c
c
c   2.numerical value for ecc pre xob
c ***********************************
c       t is negative for the past   t is in 1000 years
c
      t=t*1000.0d0
      xes=0.0d0
      xec=0.0d0
      do 4 i=1,nef
      arg=be(i)*t+ce(i)
      xes=xes+ae(i)*dsin(arg)
      xec=xec+ae(i)*dcos(arg)
    4 continue
      ecc=dsqrt(xes*xes+xec*xec)
      tra=dabs(xec)
      if(tra.le.1.0d-08) go to 10
      rp=datan(xes/xec)
      if(xec) 11,10,12
   11 rp=rp+pi
      go to 13
   12 if(xes) 14,13,13
   14 rp=rp+2.0d0*pi
      go to 13
   10 if(xes) 15,16,17
   15 rp=1.5d0*pi
      go to 13
   16 rp=0.0d0
      go to 13
   17 rp=pi/2.0d0
   13 perh=rp/pir
c
      prg=prm*t
      do 5 i=1,nop
      arg=bop(i)*t+cop(i)
      prg=prg+aop(i)*dsin(arg)
    5 continue
      prg=prg/3600.0d0+xop
      perh=perh+prg
   54 if(perh) 51,55,53
   51 perh=perh+360.0d0
      go to 54
   53 if(perh.lt.360.0d0) go to 55
      perh=perh-360.0d0
      go to 53
   55 continue
c
      pre=ecc*dsin(perh*pir)
c
      xob=xod
      do 6 i=1,nob
      arg=bob(i)*t+cob(i)
      xob=xob+aob(i)/3600.0d0*dcos(arg)
    6 continue

c
      so=dsin(xob*pir)
      xeq=(datan(4.0d0*pre/pi/so))/pir
c
      ipage=dabs(t/1000.)
c
      write(iwr1,6025) pre,ipage
 6025 format(10x,'pre   = ',f8.5,5x,'page = ',i4)
      ian=t/1000.

c ---------------------------------
c if you want to fix astronomical elements to some
c constant value, then uncomment the following lines
c
c    so = MY_CONSTANT (sine of obliquity)
c    perh = MY_LONG_OF_PERIHELION (in degrees)
c    ecc  = MY_ECCENTRICITY
c ----------------------------------

      write(iwr1,6023) ian,ecc,perh,xob,xeq
      write(iwr6,6023) ian,ecc,perh,xob,xeq
 6023 format(1x,'date = ',i6,3x,'eccen = ',f9.6,3x,'long per = ',f7.2,3x
     *,'obliq = ',f7.3,3x,'cal eq = ',f5.2,/)
c
      return
      end
c **********************************************************************
c cray .if [ $concord .eq. 1 ]
c cray nzyx  solext==
c cray .endif
      subroutine solext(sc,rphi,pihor,sdelta,nd,jlat)
      parameter(np=18)
      implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
       common/so1/sotso1(np),cmuso1(np),bhso1(np,6)
       real sotso1,cmuso1,bhso1,sc,rphi,pihor,sdelta
       common /so10/ perh,ecc,so
c
c   3.daily insolation
c ********************
c
c   option solar date - calendar date
c       daily insolation in w/m2          sotso1(lat)
c       daily averaged zenithal angle     cmuso1(lat)
c
c     input parameters : latitude phi - time t
c *********************
c
      pi=3.14159265358979d0
      pir=pi/180.0d0
c ini step=360.0d0/365.25d0
      step=360.0d0/360.00d0
      test=0.0001d0
c
      ss=sc
      tau=86400.d0
c
      sf=tau*ss/pi
      xl=perh+180.0d0
c
      phi=rphi/pir
c
c
c   4.2 calendar date  ma-ja
c -------------------
c      nd  number of this day in a year of 365 days
c      xlam = true long. sun for mean long. = 0
c      dlamm = mean long. sun for ma-ja
c
      xllp=xl*pir
      xee=ecc*ecc
      xse=dsqrt(1.0d0-xee)
      xlam=(ecc/2.0d0+ecc*xee/8.0d0)*(1.0d0+xse)*dsin(xllp)-xee/4.0d0
     1*(0.5d0+xse)*dsin(2.0d0*xllp)+ecc*xee/8.0d0*(1.0d0/3.0d0+xse)
     2*dsin(3.0d0*xllp)
      xlam=2.0d0*xlam/pir
      dlamm=xlam+(nd-80)*step
      anm=dlamm-xl
      ranm=anm*pir
      xee=xee*ecc
      ranv=ranm+(2.0d0*ecc-xee/4.0d0)*dsin(ranm)+5.0d0/4.0d0*ecc*ecc
     1*dsin(2.0d0*ranm)+13.0d0/12.0d0*xee*dsin(3.0d0*ranm)
      anv=ranv/pir
      tls=anv+xl
c
      call dayins(ecc,xl,so,tls,phi,pir,pi,test,sf,ww,dayl,rdelta,cozm)
c
c ----------------------------------------------------------------------
      sotso1(jlat)=ww/86400.d0
      cmuso1(jlat)=cozm
      pihor       =dayl*pi/24.
      sdelta      =rdelta
c ----------------------------------------------------------------------
c ini write(iwr1,1000) t,nd,phi,sotso1(jlat),cmuso1(jlat),dayl
1000  format (3x,'t =',f8.1,3x,'jour =',i4,3x,'lat =',f7.2,3x,
     *        'rs sommet journalier =',f8.2,3x,'cozm =',f7.4,
     *        3x,'duree du jour =',f6.2)
c
      return
      end
c **********************************************************************
c cray .if [ $concord .eq. 1 ]
c cray nzyx  dayins==
c cray .endif
      subroutine dayins(ecc,xl,so,dlam,phi,pir,pi,test,sf,ww,dayl,
     .                  rdelta,cozm)
      implicit double precision(a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c   output : ww=j/m2 day     dayl=length of day (hours)
c
      rphi=phi*pir
      ranv=(dlam-xl)*pir
      rau=(1.0d0-ecc*ecc)/(1.0d0+ecc*dcos(ranv))
      s=sf/rau/rau
      rlam=dlam*pir
      sd=so*dsin(rlam)
      cd=dsqrt(1.0d0-sd*sd)
      rdelta=datan(sd/cd)
      delta=rdelta/pir
      sp=sd*dsin(rphi)
      cp=cd*dcos(rphi)
      aphi=dabs(phi)
      adelta=dabs(delta)
c
c   singularity for aphi=90 and delta=0
c   particular cases for phi=0  or  delta=0
c
      tt=dabs(aphi-90.0d0)
      if ((tt.le.test).and.(adelta.le.test)) go to 2
      if(adelta.le.test) go to 6
      if(aphi.le.test) go to 7
c
c   label 2 : polar continual night or w=0  dayl=0
c   label 4 : polar continual day
c   label 3 : daily sunrise and sunset
c   label 6 : equinoxes
c   label 7 : equator
c
      at=90.0d0-adelta
      spd=phi*delta
      if (aphi.le.at) go to 3
      if (spd) 2,3,4
    2 dayl=0.00d0
      ww=0.00d0
      cozm=0.d0
      go to 5
    4 dayl=24.00d0
      ww=s*sp*pi
      cozm=sp
      go to 5
    3 tp=-sp/cp
      stp=dsqrt(1.0d0-tp*tp)
      rdayl=dacos(tp)
      dayl=24.0d0*rdayl/pi
      ww=s*(rdayl*sp+cp*stp)
      cozm=sp+(cp*stp/rdayl)
      go to 5
    6 dayl=12.0d0
      ww=s*dcos(rphi)
      cozm=2.d0*dcos(rphi)/pi
      go to 5
    7 dayl=12.0d0
      ww=s*dcos(rdelta)
      cozm=2.d0*dcos(rdelta)/pi
    5 continue
c ---------------------------------------------------------------------
      wwt=ww/86400.d0
      if (wwt.lt.0.1d0)      ww=0.d0
      if (ww.eq.0.d0) goto 8
      if (cozm.lt.1.d-03)  cozm=1.0d-03
8     continue
c ----------------------------------------------------------------------
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  solri===
c cray .endif
         subroutine solari (ysol,jlat,iilw,ntrop,nbniv)
c *126*  subroutine solari (ysol,jlat,iilw,nbniv)
c
c +++ bcms3.for +++
c
      parameter(np=18,nilw=7,nir=50)
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c **** input  ******
c
c --- caracteristiques surface
c
      common/rf0/rgsrf0(np,nilw)
c
c --- vapeur d'eau
c
      common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
c
c --- variables radiatives
c
      common/gz1/ih1,ih2,ih3,ihc1,ihc2,ihc3,ico2,io3
      common/ra1/ttra1(nir),hhra1(nir),ozra1(nir),cldra1(nir),
     .           pmbra1(nir),zkmra1(nir)
      common/so1/sotso1(np),cmuso1(np),bhso1(np,6)
      common/so2/gtso2(np),
     .           atso2(np),ab0so2(np),abcso2(np),albso2(np)
      common/so4/pso4(nir),wo3so4(nir)
      common/so5/tauso5(nir),omeso5(nir),cgso5(nir),c1iso5(nir)
      common/so6/taeso6(nir),oaeso6(nir),gaeso6(nir)
c
c --  local solar
c
      common /eti200/ uinf(3),aap(3,6),bbp(3,6)
c
c ****  output  ****
c
      common /eti510/ fup(nir),fdown(nir),fnets(nir),dfnets(nir)
c
      dimension ud(3,50),aki(2),rj(6,50),
     *rk(6,50),tr(2,50),rl(8),ruef(8)
      dimension udd(50),um(50),tcou(50)
c
      dimension rj1(2,50),rk1(2,50),rmue(50)
      dimension refup(50),trup(50),refz(2,50)
      dimension refae(50),trae(50)
      dimension rayldw(50),raylup(50),to1dw(50),to1up(50),
     *                wdw(50),wup(50),cgdw(50),cgup(50)
c
c  initialisation parametres
c
c ----------------------------------------------------------------------
ccc   rgsrf0(jlat,iilw)= 0.d0
c ----------------------------------------------------------------------
      iatst = 1
      idiffu= 1
cess  idiffu= 9
      iwri  = 2
c
      np1=nbniv+1
c ----------------------------------------------------------------------
      ysol=sqrt(1224.*cmuso1(jlat)*cmuso1(jlat)+1.)/35.
c ----------------------------------------------------------------------
      if (iwri.eq.1) write(iwr6,3005) nbniv,np1,ysol,rgsrf0(jlat,iilw)
3005  format (5x,'****  m=',i3,'   ****',3x,'****  np1=',i3,'   ****',
     * /,5x,'**** ysol=',f8.5,'   ****',3x,'**** albs=',f8.5,'****',/)
c
c  calcul pour contenu d'ozone dans les couches
      udd(np1)=0.0
      do 4001 i=1,nbniv
      l=np1-i
      lp=l+1
      udd(l)=udd(lp)+wo3so4(l)/ysol
4001  continue
      y=0.6024
      um(1)=udd(1)
      do 4002 i=2,np1
      im=i-1
      um(i)=um(im)+wo3so4(im)/y
4002  continue
c
c     ***   initialisation, ground level   ***
c
      refz(1,1)=rgsrf0(jlat,iilw)
      refz(2,1)=0.
      y=1.66
 400  format( 5x,'refz1',9x,'refz2',11x,'tr1',12x,'tr2',11x,'refup',
     *10x,'trup',10x,'refae',10x,'trae',/)
      if (iwri.eq.1) write(iwr6,400)
      do 420 i=2,np1
      j=i-1
      rneb=cldra1(j)
      tr1=0.
      re1=tr1
      re2=0.
      te2=0.
      k=1
c  equivalent zenith angle
 401   format (8(1x,e13.7))
      xmue=(1.-c1iso5(i))/ysol+c1iso5(i)*1.66
      rmue(i)=1./xmue
c ini r=(pso4(j)-pso4(i))
c ---new-----------
      r=(pso4(j)-pso4(i))*(pmbra1(1)/1013.)
c ---new-----------
c  reflectivity of layer j due to rayleigh scattering
      yydw=rmue(i)
      yyup=0.6024
      xeridw=0.0294+0.313*yydw-0.6316*yydw*yydw+0.608*yydw*
     *       yydw*yydw-0.2194*yydw*yydw*yydw*yydw
      xeriup=0.0294+0.313*yyup-0.6316*yyup*yyup+0.608*yyup*
     *       yyup*yyup-0.2194*yyup*yyup*yyup*yyup
c
c -------------- introduction aerosols  --------------------------------
      rayldw(j)=xeridw*r
      raylup(j)=xeriup*r
      if (idiffu.eq.9) rayldw(j)=0.00001
      if (idiffu.eq.9) raylup(j)=0.00001
c
      tautdw=taeso6(j)+rayldw(j)
      gmdw=gaeso6(j)*taeso6(j)/tautdw
      omegdw=(oaeso6(j)*taeso6(j)+rayldw(j)*0.99999)/tautdw
c
      tautup=taeso6(j)+raylup(j)
      gmup=gaeso6(j)*taeso6(j)/tautup
      omegup=(oaeso6(j)*taeso6(j)+raylup(j)*0.99999)/tautup
c
c ----------------------------------------------------------------------
c
c  reflectivity of layer j due to aerosol scattering
c
c     write(iwr6,9661) yydw,yyup,gaeso6(j),taeso6(j),oaeso6(j),rayldw(j),
c    *               raylup(j)
9661  format (////,'    ******  ecriture temporaire   ********',/,
     * 5x,'yydw      =',e12.5,5x,'yyup      =',e12.5,/,
     * 5x,'gaeso6(j) =',e12.5,5x,'taeso6(j) =',e12.5,5x,'oaeso6(j) =',
     * e12.5,/,
     * 5x,'rayldw(j) =',e12.5,5x,'raylup(j) =',e12.5)
c
       call deled2 (gmdw,gmup,omegdw,omegup,tautdw,tautup,
     *             refae(j),trae(j),rmue(i),refup(j),trup(j))
c
      if(rneb.eq.0) goto 430
c
      k=0
c
      to1dw(j)=tauso5(j)+tautdw
      to1up(j)=tauso5(j)+tautup
      wdw(j)=(omeso5(j)*tauso5(j)+oaeso6(j)*taeso6(j)
     *        +rayldw(j))/to1dw(j)
      wup(j)=(omeso5(j)*tauso5(j)+oaeso6(j)*taeso6(j)
     *        +raylup(j))/to1up(j)
      cgdw(j)=(cgso5(j)*omeso5(j)*tauso5(j)+gaeso6(j)*oaeso6(j)
     *          *taeso6(j)) / (wdw(j)*to1dw(j))
      cgup(j)=(cgso5(j)*omeso5(j)*tauso5(j)+gaeso6(j)*oaeso6(j)
     *          *taeso6(j)) / (wup(j)*to1up(j))
c
      call deled2 (cgdw(j),cgup(j),wdw(j),wup(j),to1dw(j),to1up(j),
     *             re1,tr1,rmue(i),re2,te2)
c
c **********************************************************************
430   continue
c **********************************************************************
c
c refz(1,j): reflexion of layer j
c            for ref. of underlaying layer =0
c refz(2,j): r#0
c
      refz(1,i)=(1.-rneb)*(refae(j)+refz(1,j)*trae(j)*
     *trup(j)/(1.-refup(j)*refz(1,j)))
     *   +  rneb * (re1+(refz(1,j)*tr1*te2/(1.-refz(1,j)*re2)))
      tr(1,j)=(trae(j)*(1.-rneb)/(1.-refup(j)*refz(1,j)))
     *   +  rneb * (tr1/(1.-refz(1,j)*re2))
      refz(2,i)=(1.-rneb)*(refae(j)+k*refz(2,j)*
     *trae(j)*trup(j)/(1.-refup(j)*refz(2,j)))
     *   +  rneb * re1
      tr(2,j)=((1.-rneb)*trae(j)/(1.-k*refup(j)*refz(2,j)))
     *   +  rneb * tr1
c
c
      if (iwri.eq.1) write(iwr6,401) refz(1,i),refz(2,i),tr(1,j),
     *                tr(2,j),refup(j),trup(j),refae(j),trae(j)
c
420      continue
c
      do 440 j=1,2
      if (iwri.eq.1) write(iwr6,450) j
 450  format (//,' j =',i6,/)
      if (iwri.eq.1) write(iwr6,451)
 451  format (1x,'     pression      rj1(down)      rk1(up)')
c  rj : downward
c  rk : upward
      rj1(j,np1)=1.
      rk1(j,np1)=refz(j,np1)
      do 440 i=1,nbniv
      l=np1-i
      lp=l+1
      rj1(j,l)=rj1(j,lp)*tr(j,l)
      re1=rj1(j,l)
      rk1(j,l)=re1*refz(j,l)
      if (iwri.eq.1) write(iwr6,401)pso4(l), rj1(j,l),rk1(j,l)
440   continue
c
      do 9 i=1,np1
      rj(1,i)=rj1(1,i)
      rj(2,i)=rj1(2,i)
      rk(1,i)=rk1(1,i)
      rk(2,i)=rk1(2,i)
      if (iwri.eq.1) write(iwr6,1000) rj(1,i),rj(2,i),rk(1,i),rk(2,i)
1000  format (4(5x,e13.6))
    9 continue
      do 8 i=1,np1
      if (iwri.eq.1) write(iwr6,3001) hhra1(i),pso4(i)
3001    format (5x,e15.6,3x,e15.6)
    8 continue
      if (iwri.eq.1) write(iwr6,3003)(hhra1(i),pso4(i),ttra1(i),
     .                                i=1,nbniv)
3003  format (//,10x,'hhra1',10x,'pso4',10x,'ttra1',/,(2x,e12.4,
     *2x,e12.4,2x,e12.4))
c   for downward radiation : ud
c   for upward   radiation : um
c   interactions between o3 absorption and
c     scattering  are neglected
c   amount of absorber in each layer
c   for h2o : ud(1,i)
c   for o2,co2 : ud(2,i)
      ud(3,np1)=0.00
      ud(1,np1)=0.00
      ud(2,np1)=0.00
      do 20 i=1,nbniv
      l=np1-i
      ud(3,l)=udd(l)
   20 continue
      h2otot=0.
      co2tot=0.
c
      to=273.
      g=9.81
      ala=c1ev1(jlat)+1.9
      fach2o=10.*1013.*hhra1(1)/(g*ala)
c
      do 40 i=1,nbniv
      ip=i+1
      tcou(i)=0.5*(ttra1(i)+ttra1(ip))
      fctw=(to/tcou(i))**0.45
      fco2=(to/tcou(i))**1.375
      if (iatst.ne.1) ud(1,i)=.5*543.62*(hhra1(i)+hhra1(ip))*
     *(pso4(i)**1.9-pso4(ip)**1.9)*fctw
      if (iatst.eq.1) ud(1,i)=fach2o*fctw*(pso4(i)**ala-pso4(ip)**ala)
      h2otot=h2otot+ud(1,i)
      ud(2,i)=150.678*fco2*(pso4(i)**1.75-pso4(ip)**1.75)
      co2tot=co2tot+ud(2,i)
 40   continue
  306 format( 7x,'uh2o',10x,'uco2',11x,'udo3',11x,'um03',10x,
     *           'h2otot',10x,'co2tot')
      if (iwri.eq.1) write(iwr6,306)
      do 300 i=1,np1
      if (iwri.eq.1) write(iwr6,305)
     * ud(1,i),ud(2,i),ud(3,i),um(i),h2otot,co2tot
300   continue
  305 format(8(2x,e13.7))
      if (iwri.eq.1) write(iwr6,1003)
 1003 format (/,3x,'**** nebulosite totale au dessus du niveau i ****')
      if (iwri.eq.1) write(iwr6,1004)
 1004 format (/,1x,'** epaisseurs opt.    q h2o')
      do 310 i=1,nbniv
      if (iwri.eq.1) write(iwr6,305) tauso5(i),hhra1(i)
310   continue
      if (iwri.eq.1) write(iwr6,1005)
 1005 format (/)
      if (iwri.eq.1) write(iwr6,305) (c1iso5(i),i=1,np1)
      if (iwri.eq.1) write(iwr6,1005)
c  grey absorption coefficients for simulation absorption
c  by (1)h2o ,(2)co2,o2
      aki(1)=.457
      aki(2)=.00636
      y=1.66
      n=2
c ******calcul des flux avec absorption simulee
      do 90 iab=1,2
      if (iwri.eq.1) write(iwr6,1008) iab
 1008 format (/,1x,'iab =',i6,/)
      rki=aki(iab)
      if ((iab.eq.1).and.(ih1.ne.1)) rki=0.0
      if ((iab.eq.2).and.(ico2.ne.1)) rki=0.0
c -----------
c ccc  rki = 0.d0
c -----------
      refz(1,1)=rgsrf0(jlat,iilw)
      refz(2,1)=0.
      do 100 i=2,np1
      j=i-1
      rneb=cldra1(j)
      aa=ud(iab,j)
      k=1
      s= exp(-rki*aa*y)
      g= exp(-rki*aa/rmue(i))
      re1=0.
      tr1=0.
      re2=0.
      te2=0.
c *****************************************************************
      if(rneb.eq.0) goto 110
c *****************************************************************
c
      k=0
c
      adw=to1dw(j)*wdw(j)
      aup=to1up(j)*wup(j)
      cdw=cgdw(j)*wdw(j)*to1dw(j)
      cup=cgup(j)*wup(j)*to1up(j)
      to1dw(j)=rki*aa+to1dw(j)
      to1up(j)=rki*aa+to1up(j)
      wdw(j)=adw/to1dw(j)
      wup(j)=aup/to1up(j)
      cgdw(j)=cdw/(wdw(j)*to1dw(j))
      cgup(j)=cup/(wup(j)*to1up(j))
c
      call deled2 (cgdw(j),cgup(j),wdw(j),wup(j),to1dw(j),to1up(j),
     *             re1,tr1,rmue(i),re2,te2)
c
c **********************************************************************
110   continue
c **********************************************************************
c
c refz(1,j): reflexion of layer j
c            for ref. of underlaying layer =0
c refz(2,j): r#0
c
      refz(1,i)=(1.-rneb)*g*s*(refae(j)+refz(1,j)*trae(j)*
     *trup(j)/(1.-refup(j)*refz(1,j)))
     *   +  rneb * (re1+(refz(1,j)*tr1*te2/(1.-refz(1,j)*re2)))
      tr(1,j)=(trae(j)*g*(1.-rneb)/(1.-refup(j)*refz(1,j)))
     *   +  rneb * (tr1/(1.-refz(1,j)*re2))
      refz(2,i)=(1.-rneb)*g*s*(refae(j)+k*refz(2,j)*
     *trae(j)*trup(j)/(1.-refup(j)*refz(2,j)))
     *   +  rneb * re1
      tr(2,j)=((1.-rneb)*g*trae(j)/(1.-k*refup(j)*refz(2,j)))
     *   +  rneb * tr1
c
  100 continue
c
      if (iwri.eq.1) write(iwr6,1009) n,iab
 1009 format (3x,'n =',i6,3x,'iab =',i6,/)
      if (iwri.eq.1) write(iwr6,307)
307   format (5x,'refz1',10x,'refz2',12x,'tr1',13x,'tr2',/)
      do 350 i=2,np1
      j=i-1
      if (iwri.eq.1) write(iwr6,305) refz(1,i),refz(2,i),tr(1,j),tr(2,j)
350     continue
      if (iwri.eq.1) write(iwr6,1010)
 1010 format (/)
      do 120 k=1,2
      if (iwri.eq.1) write(iwr6,1011) k
 1011 format (/,' k =',i6,/)
      if (iwri.eq.1) write(iwr6,1012)
 1012 format (1x,'    pression       rj(down)       rk(up)')
      n=n+1
      rj(n,np1)=1.
      rk(n,np1)=refz(k,np1)
      do 120 i=1,nbniv
      l=np1-i
      lp=l+1
      rj(n,l)=rj(n,lp)*tr(k,l)
      re1=rj(n,l)
      rk(n,l)=re1*refz(k,l)
      if (iwri.eq.1) write(iwr6,305)pso4(l), rj(n,l),rk(n,l)
120   continue
c
   90 continue
c
c now the fluxes up(rk) and down(rj) are calculated without
c     gaz absorption n=1 and 2
c    with h2o absorption : n=3,4
c    with co2 absorption : n=5,6
c    n=1,3,5 : underlaying reflectivity =0
c      2,4,6 : #0
c  calculation of the equivalent amount of absorber
c       ueff= -ln(f(k)/f(0))/k
c calcul des quantites efficaces d absorbant
      ee=1.e-10
      do 130 i=1,np1
      do 130 j=1,5,2
      jp=j+1
      rj(j,i)=rj(j,i)-rj(jp,i)+ee
  130 rk(j,i)=rk(j,i)-rk(jp,i)+ee
      do 140 i=1,np1
      do 140 j=2,6,2
      rj(j,i)=ee+rj(j,i)
  140 rk(j,i)=ee+rk(j,i)
      if (iwri.eq.1) write(iwr6,1013)
 1013 format (/,1x,'*** 1: sans reflexion de la couche sous jacente')
      if (iwri.eq.1) write(iwr6,1014)
 1014 format (1x,'*** 2: avec reflexion de la couche sous jacente',/)
 161  format (30x,'u.eff down',25x,'*',30x,'u.eff up')
 162  format (8x,'h2o:1',11x,'h2o:2',11x,'co2:1',11x,'co2:2',
     *4x,'*',7x,'h2o:1',11x,'h2o:2',11x,'co2:1',11x,'co2:2')
 164  format ( 1h ,65('-'),'*',62('-'))
 163  format(1h ,4e16.4,1x,'*',4e16.4)
 165  format(65x,'*')
      if (iwri.eq.1) write(iwr6,161)
      if (iwri.eq.1) write(iwr6,164)
      if (iwri.eq.1) write(iwr6,165)
      if (iwri.eq.1) write(iwr6,162)
      if (iwri.eq.1) write(iwr6,165)
      do 150 i=1,np1
      k=1
      do 160 j=1,2
      rki=aki(j)
      do 160 n=1,2
c  effective amount
      w=alog(rj(n,i)/rj(n+2*j,i))/rki
c  calculation of the transmission
      call tttt(j,w,r1)
      if ((j.eq.1).and.(ih1.ne.1)) r1=1.0
      if ((j.eq.2).and.(ico2.ne.1)) r1=1.0
c ------------
c ccc  r1 = 1.d0
c ------------
      rl(k)=r1
      ruef(k)=w
      w=alog(rk(n,i)/rk(n+2*j,i))/rki
      call tttt(j,w,r1)
      if ((j.eq.1).and.(ih1.ne.1)) r1=1.0
      if ((j.eq.2).and.(ico2.ne.1)) r1=1.0
c ------------
c ccc  r1 = 1.d0
c ------------
      rl(4+k)=r1
      ruef(k+4)=w
      k=k+1
160   continue
166   continue
      if (iwri.eq.1) write(iwr6,163)(ruef(j),j=1,8)
c  upward and downward fluxes with h2o and co2 absorption
      fdown(i)=rj(1,i)*rl(1)*rl(3)+rj(2,i)*rl(2)*rl(4)
  150 fup(i)=rk(1,i)*rl(5)*rl(7)+rk(2,i)*rl(6)*rl(8)
      if (iwri.eq.1) write(iwr6,2014)
 2014 format (//,5x,'fdown(h2o,co2)',6x,'fup(h2o,co2)')
      do 2015 i=1,np1
      if (iwri.eq.1) write(iwr6,2013)  fdown(i),fup(i)
 2013 format (5x,e13.7,5x,e13.7)
 2015 continue
c abs o3
      iab=3
      do 170 i=1,np1
      w=ud(3,i)
      call tttt(iab,w,r1)
      if (io3.ne.1) r1=1.0
c ------------
c ccc  r1 = 1.d0
c ------------
      fdown(i)=r1*fdown(i)
      w=um(i)
      call tttt(iab,w,r1)
      if (io3.ne.1) r1=1.0
c ------------
c ccc  r1 = 1.d0
c ------------
  170 fup(i)=r1*fup(i)
c
c  net fluxes
      do 180 i=1,np1
  180 fnets(i)=fdown(i)-fup(i)
c
c ----------------------------------------------------------------------
c
c --  output => bcm
c
      gtso2(jlat)  = fnets(1)
      atso2(jlat)  = fnets(np1)-fnets(1)
cess  atso2(jlat)  = fnets(ntrop)-fnets(1)
      albso2(jlat) = fup(np1)/fdown(np1)
c ----------------------------------------------------------------------
      fupdwn=100.*(fup(np1)/fdown(np1))
      if (iwri.eq.1) write(iwr6,185) fupdwn
  185 format(' albedo: ',e13.7)
      if (iwri.eq.1) write(iwr6,192)
  192 format(3x,'i',6x,'pso4(i+1)',6x,'pso4(i)',10x,'f up',
     *10x,'f down',10x,'f net',10x,'dfnets')
      if (iwri.ne.1) goto 1020
      write(iwr6,195) np1,fup(np1),fdown(np1),fnets(np1)
1020  continue
c
      dfnets(np1)=0.
      do 190 i=1,nbniv
      l=np1-i
      lp=l+1
      r2=pso4(lp)
      r3=pso4(l)
      dfnets(l)=fnets(lp)-fnets(l)
      if (iwri.ne.1) goto 1021
      write(iwr6,195) l,r2,r3,fup(l),fdown(l),
     *                     fnets(l),dfnets(l)
1021  continue
190   continue
  195 format(1h ,i6,6(2x,e13.7))
      return
      end
c **********************************************************************
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tttt2===
c cray .endif
c computation of gaz transmission by pade approximants
c   sauf pour h2o    exponential sum  (lacis-hansen)
      subroutine tttt2(j,w,r1)
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common /eti200/ uinf(3),aap(3,6),bbp(3,6)
c
      dimension akn(8),pkn(8)
      data akn /4.0e-05,2.0e-03,3.5e-02,3.77e-01,1.95e0,9.40e0,
     *          44.6e0,190.e0/
      data pkn /0.6470e0,0.0698e0,0.1443e0,0.0584e0,0.0335e0,
     *          0.0225e0,0.0158e0,0.0087e0/
c
c
      if (j.eq.1) goto 500
c
      r1=aap(j,6)
      r2=bbp(j,6)
      do 20 i=1,5
      k=6-i
  20  r1=r1*w+aap(j,k)
      do 30 i=1,5
      k=6-i
  30  r2=r2*w+bbp(j,k)
      if((abs(r2).lt.1.e-15).or.(abs(r2).gt.1.e+15)) goto 31
      r1=(r1/r2)*(1-uinf(j))+uinf(j)
c
      goto 32
500   continue
      r1=0.d0
      do 50 i=1,8
      x1=-akn(i)*w
      if (x1.lt.-150.) x1=-150.
      r1=r1+pkn(i)*exp(x1)
50    continue
c
      goto 32
  31  write(iwr6,1016) r2,r1,j
 1016 format (/,2x,'r2 =',e15.8,2x,'r1 =',e15.8,2x,'j =',i6,/)
      stop 1017
 32   return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  sola3===
c cray .endif
      subroutine sola3g (ysol,jlat,iilw,m)
c
c +++ bcms4s.for +++
c
      parameter(np=18,nilw=7,nir=50)
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
c --- caracteristiques surface
c
      common/rf0/rgsrf0(np,nilw)
c
c --- vapeur d'eau
c
      common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
c
c --- ozone
c
      common/gz0/uggz0(np),ucgz0(np)
c
c --- discretistion 3 couches
c
      common/di1/it1,it2,it3,ic1,ic2,ic3,inu
      common/di2/j13(3),j23(3)
c
c --- variables radiatives
c
      common/ra1/ttra1(nir),hhra1(nir),ozra1(nir),cldra1(nir),
     .           pmbra1(nir),zkmra1(nir)
      common/so1/sotso1(np),cmuso1(np),bhso1(np,6)
      common/so2/gtso2(np),
     .           atso2(np),ab0so2(np),abcso2(np),albso2(np)
      common/so4/pso4(nir),wo3so4(nir)
      common/so5/tauso5(nir),omeso5(nir),cgso5(nir),c1iso5(nir)
      common/so6/taeso6(nir),oaeso6(nir),gaeso6(nir)
      common/eti200/ uinf(3),aap(3,6),bbp(3,6)
c
      dimension h2od(3),wat(10),co(10),o3(2),th2o(10),tco2(10),to3(2)
      dimension tcou(nir)
c ----------------------------------------------------------------------
c ccc  rgsrf0(jlat,iilw)= 0.d0
c ----------------------------------------------------------------------
c
      g=9.81d0
      ysol=sqrt(1224.*cmuso1(jlat)*cmuso1(jlat) + 1.) / 35.
c choix parametre pour simplification par diffusivite
c ini-------------
      dif =1.50
      difp=1./dif
      difnu=1.50
c ini-------------
      dif =1.66
      difp=1./dif
      difnu=1.66
c  calculs preliminaires pour partie ciel nuageux
c cnebu = nebulosite totale vue de la surface (recouvrement aleatoire)
      taun=tauso5(inu)
      cnebu=c1iso5(1)
      omegn=omeso5(inu)
      cgn=cgso5(inu)
c  calculs preliminaires pour ciel clair et nuageux
c  aerosols
      aertot=0.0000001
      omegae=0.0000001
      gae   =0.
      aerto3=0.0000001
      omega3=0.0000001
      gae3  =0.
      aerto1=0.0000001
      omega1=0.0000001
      gae1  =0.
      do 203 i=1,ic3
      aertot=aertot+taeso6(i)
      omegae=omegae+oaeso6(i)*taeso6(i)
      gae   =gae+gaeso6(i)*oaeso6(i)*taeso6(i)
203   continue
      omegae=omegae/aertot
      gae   =gae/(omegae*aertot)
      do 204 i=1,ic1
      aerto3=aerto3+taeso6(i)
      omega3=omega3+oaeso6(i)*taeso6(i)
      gae3  =gae3+gaeso6(i)*oaeso6(i)*taeso6(i)
204   continue
      omega3=omega3/aerto3
      gae3  =gae3/(omega3*aerto3)
      do 205 i=ic2+1,ic3
      aerto1=aerto1+taeso6(i)
      omega1=omega1+oaeso6(i)*taeso6(i)
      gae1  =gae1+gaeso6(i)*oaeso6(i)*taeso6(i)
205   continue
      omega1=omega1/aerto1
      gae1  =gae1/(omega1*aerto1)
c rayleigh
      am  =1./ysol
      yd  =ysol
      yu  =difp
      xd  =0.0294+0.313*yd-0.6316*yd*yd+0.608*yd*yd*yd-
     .     0.2194*yd*yd*yd*yd
      xu  =0.0294+0.313*yu-0.6316*yu*yu+0.608*yu*yu*yu-
     .     0.2194*yu*yu*yu*yu
      traydw=xd*pmbra1(1   )/1013.
c -------
c ccc  traydw=0.000001d0
c -------
      tautdw=aertot+traydw
      omdw  =(traydw + omegae*aertot)/tautdw
      trayup=xu*pmbra1(1   )/1013.
c -------
c ccc  trayup=0.000001d0
c -------
      tautup=aertot+trayup
      omup  =(trayup + omegae*aertot)/tautup
      beaedw=(2.-3.*yd*gae)/4.
      beaedw=(aertot*omegae*beaedw + traydw*0.5) /
     .        (aertot*omegae + traydw)
      if (beaedw.lt.0.) beaedw=0.01
      beaeup=(2.-3.*yu*gae)/4.
      beaeup=(aertot*omegae*beaeup + trayup*0.5) /
     .        (aertot*omegae + trayup)
      if (beaeup.lt.0.) beaeup=0.01
      refdw = beaedw*omdw*(tautdw/yd) /
     .        (1.+(1.-omdw+beaedw*omdw)*(tautdw/yd))
      trdw  = 1. /
     .        (1.+(1.-omdw+beaedw*omdw)*(tautdw/yd))
      refup = beaeup*omup*(tautup/yu) /
     .        (1.+(1.-omup+beaeup*omup)*(tautup/yu))
      trup  = 1. /
     .        (1.+(1.-omup+beaeup*omup)*(tautup/yu))
c -------  ajoute par rapport a bcmokn2 fortran t1 -------------
      trayd1 = xd*(pmbra1(it2)-pmbra1(it3))/1013.
      trayu1 = xu*(pmbra1(it2)-pmbra1(it3))/1013.
      trayu3 = xu*(pmbra1(1  )-pmbra1(it1))/1013.
      trayd3 = trayu3
      tautd1=aerto1+trayd1
      omd1  =(trayd1 + omega1*aerto1)/tautd1
      tautu1=aerto1+trayu1
      omu1  =(trayu1 + omega1*aerto1)/tautu1
      tautd3=aerto3+trayd3
      omd3  =(trayd3 + omega3*aerto3)/tautd3
      tautu3=aerto3+trayu3
      omu3  =(trayu3 + omega3*aerto3)/tautu3
      beaed1=(2.-3.*yd*gae1)/4.
      beaed1=(aerto1*omega1*beaed1 + trayd1*0.5) /
     .        (aerto1*omega1 + trayd1)
      if (beaed1.lt.0.) beaed1=0.01
      beaeu1=(2.-3.*yu*gae1)/4.
      beaeu1=(aerto1*omega1*beaeu1 + trayu1*0.5) /
     .        (aerto1*omega1 + trayu1)
      if (beaeu1.lt.0.) beaeu1=0.01
      beaed3=(2.-3.*yu*gae3)/4.
      beaed3=(aerto3*omega3*beaed3 + trayd3*0.5) /
     .        (aerto3*omega3 + trayd3)
      if (beaed3.lt.0.) beaed3=0.01
      beaeu3=(2.-3.*yu*gae3)/4.
      beaeu3=(aerto3*omega3*beaeu3 + trayu3*0.5) /
     .        (aerto3*omega3 + trayu3)
      if (beaeu3.lt.0.) beaeu3=0.01
      refd1 = beaed1*omd1*(tautd1/yd) /
     .        (1.+(1.-omd1+beaed1*omd1)*(tautd1/yd))
      trd1  = 1. /
     .        (1.+(1.-omd1+beaed1*omd1)*(tautd1/yd))
      refu1 = beaeu1*omu1*(tautu1/yu) /
     .        (1.+(1.-omu1+beaeu1*omu1)*(tautu1/yu))
      tru1  = 1. /
     .        (1.+(1.-omu1+beaeu1*omu1)*(tautu1/yu))
      refd3 = beaed3*omd3*(tautd3/yu) /
     .        (1.+(1.-omd3+beaed3*omd3)*(tautd3/yu))
      trd3  = 1. /
     .        (1.+(1.-omd3+beaed3*omd3)*(tautd3/yu))
      refu3 = beaeu3*omu3*(tautu3/yu) /
     .        (1.+(1.-omu3+beaeu3*omu3)*(tautu3/yu))
      tru3  = 1. /
     .        (1.+(1.-omu3+beaeu3*omu3)*(tautu3/yu))
c -------  fin de ajoute par rapport a bcmokn2 fortran t1 -------------
c h2o
      to=273.
      po=1013.
      ala=c1ev1(jlat)+1.9
      fach2o=10.*po*hhra1(1)/(g*ala)
      do 300 iw=1,3
      wvap=0.
      j1=j13(iw)
      j2=j23(iw)
      j2m1=j2-1
      do 301 i=j1,j2m1
      tcou(i)=0.5*(ttra1(i)+ttra1(i+1))
      wvap=wvap+fach2o*(to/tcou(i))**0.45*
     .           (pso4(i)**ala-pso4(i+1)**ala)
c --------
c ccc  wvap = 0.d0
c --------
301   continue
      h2od(iw)=wvap
300   continue
      h2o=h2od(3)+h2od(2)+h2od(1)
      wat(1) =h2o*am
      wat(2) =wat(1)+h2o*dif
      wat(3) =h2od(3)*am
      wat(4) =wat(3)+h2od(3)*dif
c ---------------------------------------------------------------------
      wat(5) =wat(3)+h2od(2)*difnu
cess  wat(5) =wat(3)+h2od(2)*2.0d0
c ---------------------------------------------------------------------
      wat(6) =wat(5)+h2od(1)*dif
      wat(7) =wat(6)+h2od(1)*dif
c ---------------------------------------------------------------------
      wat(8) =wat(7)+h2od(2)*difnu
cess  wat(8) =wat(7)+h2od(2)*2.0d0
c ---------------------------------------------------------------------
      do 10 iw=1,8
      r1=aap(1,6)
      r2=bbp(1,6)
      do 20 i=1,5
      k=6-i
20    r1=r1*wat(iw)+aap(1,k)
      do 30 i=1,5
      k=6-i
30    r2=r2*wat(iw)+bbp(1,k)
      th2o(iw)=(r1/r2)*(1.-uinf(1))+uinf(1)
c -----
c ccc  th2o(iw)= 1.d0
c -----
10    continue
      th2otd=th2o(1)
      th2otu=th2o(2)/th2o(1)
      th2o1d=th2o(3)
      th2o1u=th2o(4)/th2o(3)
      th2o2d=th2o(5)/th2o(3)
      th2o3d=th2o(6)/th2o(5)
      th2o3u=th2o(7)/th2o(6)
      th2o2u=th2o(8)/th2o(7)
c co2
      facop =pso4(1  )**1.75
      facot =((ttra1(1  )+ttra1(it3))*0.5)**1.375
      facot =to**1.375/facot
      facot3=((ttra1(it2)+ttra1(it3))*0.5)**1.375
      facot3=to**1.375/facot3
      facot2=((ttra1(it1)+ttra1(it2))*0.5)**1.375
      facot2=to**1.375/facot2
      facot1=((ttra1(1  )+ttra1(it1))*0.5)**1.375
      facot1=to**1.375/facot1
      facop3=pso4(it2)**1.75-pso4(it3)**1.75
      facop2=pso4(it1)**1.75-pso4(it2)**1.75
      facop1=pso4(1  )**1.75-pso4(it1)**1.75
      co(1) =150.678*facot*facop*am
      co(2) =co(1) + 150.678*facot*facop*dif
      co(3) =150.678*facot3*facop3*am
      co(4) =co(3) + 150.678*facot3*facop3*dif
      co(5) =co(3) + 150.678*facot2*facop2*difnu
      co(6) =co(5) + 150.678*facot1*facop1*dif
      co(7) =co(6) + 150.678*facot1*facop1*dif
      co(8) =co(7) + 150.678*facot2*facop2*difnu
      do 40 ic=1,8
c -------
c cccc co(ic) = 0.0d0
c -------
      r1=aap(2,6)
      r2=bbp(2,6)
      do 50 i=1,5
      k=6-i
50    r1=r1*co(ic)+aap(2,k)
      do 60 i=1,5
      k=6-i
60    r2=r2*co(ic)+bbp(2,k)
      tco2(ic)=(r1/r2)*(1.-uinf(2))+uinf(2)
c -----
c ccc  tco2(ic)= 1.d0
c -----
40    continue
      tco2td=tco2(1)
      tco2tu=tco2(2)/tco2(1)
      tco21d=tco2(3)
      tco21u=tco2(4)/tco2(3)
c     tco22d=tco2(5)/tco2(3)
      tco23d=tco2(6)/tco2(5)
      tco23u=tco2(7)/tco2(6)
c     tco22u=tco2(8)/tco2(7)
c o3
c _CT o3(1) =uggz0(jlat)*am
      o3(1) =ozra1(   1)*am
c _CT o3(2) =o3(1)+uggz0(jlat)*dif
      o3(2) =o3(1)+ozra1(   1)*dif
      do 70 io=1,2
c ------
c ccc  o3(io) = 0.d0
c ------
      r1=aap(3,6)
      r2=bbp(3,6)
      do 80 i=1,5
      k=6-i
80    r1=r1*o3(io)+aap(3,k)
      do 90 i=1,5
      k=6-i
90    r2=r2*o3(io)+bbp(3,k)
      to3(io)=(r1/r2)*(1.-uinf(3))+uinf(3)
c -----
c ccc  to3(io) = 1.d0
c -----
70    continue
      to3td=to3(1)
      to3tu=to3(2)/to3(1)
c ciel clair
      td    =th2otd*tco2td*to3td
      tu    =th2otu*tco2tu*to3tu
c ciel nuageux
c rayleigh et aerosol entierement sous le nuage
c     tudg=to3td*th2o1d*tco21d
c     tddg=th2o3d*tco23d
c     tuug=to3tu*th2o1u*tco21u
c     tdug=th2o3u*tco23u
c ----------------------------------------------------------------------
      tud=to3td*th2o1d*tco21d*trd1
      tdd=th2o3d*tco23d*trd3
      tuu=to3tu*th2o1u*tco21u*tru1
      tdu=th2o3u*tco23u*tru3
c ----------------------------------------------------------------------
c ----old---------------------------------------------------------------
c     gp  =cgn/(1.+cgn)
c     taup=(1.-cgn*cgn)*taun
c     g1 =0.75*(1.-gp)
c     g3d=0.25*(2.-3.*gp*ysol)
c     g3u=0.25*(2.-3.*gp*difp)
c     deno  =1.+g1*taup
c     alphad=(g1*taup+(g3d-g1*ysol)*(1.-exp(-taup/ysol)))/deno
c     alphau=(g1*taup+(g3u-g1*difp)*(1.-exp(-taup/difp)))/deno
c     ts    =(1.-alphad)*th2o2d
c     tsp    =(1.-alphau)*th2o2u
c ----old---------------------------------------------------------------
c ---- new 1  ---  avec deled2  ----------------------------------------
      call deled2 (cgn,cgn,omegn,omegn,taun,taun,re1,tr1,ysol,re2,te2)
      alphad=re1
      alphau=re2
      ts =tr1*th2o2d
      tsp=te2*th2o2u
c ---- new 1  ----------------------------------------------------------
c ---- new 2  ---  avec stephens 1978  ---------------------------------
c     g3d=0.25*(2.-3.*cgn*ysol)
c     g3u=0.25*(2.-3.*cgn*difp)
c     xstep1= 1.-omegn+2.*g3d*omegn
c     xstep2= 1.-omegn+2.*g3u*omegn
c     ustepd=xstep1/(1.-omegn)
c     ustepu=xstep2/(1.-omegn)
c     teffd =sqrt((1.-omegn)*xstep1)*(taun/ysol)
c     teffu =sqrt((1.-omegn)*xstep2)*(taun/difp)
c     rd    =((sqrt(ustepd)+1.)**2)*exp(teffd)-((sqrt(ustepd)-1.)**2)
c    1        *exp(-teffd)
c     ru    =((sqrt(ustepu)+1.)**2)*exp(teffu)-((sqrt(ustepu)-1.)**2)
c    1        *exp(-teffu)
c     alphad=(((ustepd-1.)*exp(teffd))/rd) - exp(-teffd)
c     alphau=(((ustepu-1.)*exp(teffu))/ru) - exp(-teffu)
c     ts =    th2o2d * 4.* sqrt(ustepd) / rd
c     tsp=    th2o2u * 4.* sqrt(ustepu) / ru
c ---- new 2  ----------------------------------------------------------
c calcul des transmissions et reflexions
c ini glob  = (1. - cnebu)*( td * trdw /
c ini.       (1.-rgsrf0(jlat,iilw)*tu*tu*refup))   +
c ini.         cnebu * (tud*ts*tdd) /
c ini.       (1.-rgsrf0(jlat,iilw)*alphau*tdd*tdu)
c ini abss  = (1. - rgsrf0(jlat,iilw))*glob
c ini albp =(1.-cnebu)*( td*refdw*tu    +
c ini.                            (rgsrf0(jlat,iilw)*td*trdw*tu*trup /
c ini.       (1.-rgsrf0(jlat,iilw)*tu*tu*refup)))   +
c ini.      cnebu * (tud*alphad*tuu   +
c ini.                  (rgsrf0(jlat,iilw)*tud*ts*tdd*tdu*tsp*tuu    /
c ini.       (1.-rgsrf0(jlat,iilw)*tdd*tdu*alphau)))
c ---- nouvelle combinaison des couches  -------------------------------
c     tudg=to3td*th2o1d*tco21d
c     tddg=th2o3d*tco23d
c     tuug=to3tu*th2o1u*tco21u
c     tdug=th2o3u*tco23u
      r4  =  rgsrf0(jlat,iilw)
      r3  =  refd3                 + r4 *tdd *tdu / (1. -refu3  *r4)
      r2  =  alphad * th2o2d       + r3 *ts  *tsp / (1. -alphau *r3)
      r1  =  to3td * to3tu * refd1 + r2 *tud *tuu / (1. -refu1  *r2)
c ini albp =(1.-cnebu)*( td*refdw*tu    +
c ini.                            (rgsrf0(jlat,iilw)*td*trdw*tu*trup /
c ini.       (1.-rgsrf0(jlat,iilw)*tu*tu*refup)))   +
c ini.      cnebu * r1
c ----------------------------------------------------------------------
      albp =(1.-cnebu)*( td *  refdw            +
     .                            (rgsrf0(jlat,iilw)*td*trdw*tu*trup /
     .       (1.-rgsrf0(jlat,iilw)*tu*tu*refup)))   +
     .      cnebu * r1
c ----------------------------------------------------------------------
      t1 = 1.
      t2 = tud  / ( 1. - refu1  * r2 )
      t3 = ts   / ( 1. - alphau * r3 )
      t4 = tdd  / ( 1. - refu3  * rgsrf0(jlat,iilw) )
      glob  = (1. - cnebu)*( td * trdw /
     .       (1.-rgsrf0(jlat,iilw)*tu*tu*refup))   +
     .         cnebu * (t1 * t2 * t3 * t4 )
      abss  = (1. - rgsrf0(jlat,iilw))*glob
c ----------------------------------------------------------------------
c
c ---  output =>  bcm
c
      gtso2 (jlat) = abss
      albso2(jlat) = albp
      atso2 (jlat) = 1. - gtso2(jlat) - albso2(jlat)
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  deled2==
c cray .endif
      subroutine deled2(gd,gu,wc1d,wc1u,to1d,to1u,re1,tr1,rmu0,re2,tr2)
c
c +++ bcms5.for +++
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      ffd=gd*gd
      gpd=gd/(1.+gd)
      topd=(1.-wc1d*ffd)*to1d
      wcpd=(1.-ffd)*wc1d/(1.-wc1d*ffd)
c
      ffu=gu*gu
      gpu=gu/(1.+gu)
      topu=(1.-wc1u*ffu)*to1u
      wcpu=(1.-ffu)*wc1u/(1.-wc1u*ffu)
c
      gam1d=0.25*(7.-wcpd*(4.+3.*gpd))
      gam2d=-0.25*(1.-wcpd*(4.-3.*gpd))
      akd  =sqrt(gam1d*gam1d-gam2d*gam2d)
      fac1 =akd*topd
      if (fac1.gt.20.) fac1 =20.
      dexppd=exp(fac1)
      dexpmd=exp(-fac1)
      x2d   =(akd+gam1d)*dexppd+(akd-gam1d)*dexpmd
c
      gam1u=0.25*(7.-wcpu*(4.+3.*gpu))
      gam2u=-0.25*(1.-wcpu*(4.-3.*gpu))
      aku  =sqrt(gam1u*gam1u-gam2u*gam2u)
      fac1 =aku*topu
      if (fac1.gt.20.) fac1 =20.
      dexppu=exp(fac1)
      dexpmu=exp(-fac1)
      x2u   =(aku+gam1u)*dexppu+(aku-gam1u)*dexpmu
c
c-- 'direct' incident radiation
c
      rmu   =rmu0
      akr   =akd*rmu
      fac1  =topd/rmu
      if (fac1.gt.20.) fac1  =20.
      expp  =exp(fac1)
      expm  =exp(-fac1)
      gam3  =0.25*(2.-3.*gpd*rmu)
      gam4  =1.-gam3
      alph1 =gam1d*gam4+gam2d*gam3
      alph2 =gam1d*gam3+gam2d*gam4
      x1    =(1.-akr*akr)
      x3    =(1.-akr)*(alph2+akd*gam3)*dexppd
      x4    =(1.+akr)*(alph2-akd*gam3)*dexpmd
      x5    =2.*akd*(gam3-alph2*rmu)*expm
      x6    =(1.+akr)*(alph1+akd*gam4)*dexppd
      x7    =(1.-akr)*(alph1-akd*gam4)*dexpmd
      x8    =2.*akd*(gam4+alph1*rmu)*expp
      re1   =(wcpd/(x1*x2d))*(x3-x4-x5)
      tr1   =expm*(1.-(wcpd/(x1*x2d))*(x6-x7-x8))
c
      re1   =re1+1.e-10
      tr1   =tr1+1.e-10
c
c-- diffuse incident  radiation (rmu0=0.6024)
c
      rmu   =0.6024
      akr   =aku*rmu
      fac1  =topu/rmu
      if (fac1.gt.20.) fac1  =20.
      expp  =exp(fac1)
      expm  =exp(-fac1)
      gam3  =0.25*(2.-3.*gpu*rmu)
      gam4  =1.-gam3
      alph1 =gam1u*gam4+gam2u*gam3
      alph2 =gam1u*gam3+gam2u*gam4
      x1    =(1.-akr*akr)
      x3    =(1.-akr)*(alph2+aku*gam3)*dexppu
      x4    =(1.+akr)*(alph2-aku*gam3)*dexpmu
      x5    =2.*aku*(gam3-alph2*rmu)*expm
      x6    =(1.+akr)*(alph1+aku*gam4)*dexppu
      x7    =(1.-akr)*(alph1-aku*gam4)*dexpmu
      x8    =2.*aku*(gam4+alph1*rmu)*expp
      re2   =(wcpu/(x1*x2u))*(x3-x4-x5)
      tr2   =expm*(1.-(wcpu/(x1*x2u))*(x6-x7-x8))
c
      re2   =re2+1.e-10
      tr2   =tr2+1.e-10
c
      return
      end
c ***************    tttt    *****************
c cray .if [ $concord .eq. 1 ]
c cray nzyx  tttt====
c cray .endif
c computation of gaz transmission by pade approximants
      subroutine tttt(j,w,r1)
c
c _dp implicit double precision (a-h,o-z)
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common /eti200/ uinf(3),aap(3,6),bbp(3,6)
c
      r1=aap(j,6)
      r2=bbp(j,6)
      do 20 i=1,5
      k=6-i
  20  r1=r1*w+aap(j,k)
      do 30 i=1,5
      k=6-i
  30  r2=r2*w+bbp(j,k)
      if(abs(r2).lt.1.e-30.or.abs(r2).gt.1.e+30) goto 31
      r1=(r1/r2)*(1-uinf(j))+uinf(j)
      goto 32
  31  write(iwr6,1016) r2,r1,j
 1016 format (/,2x,'r2 =',e15.8,2x,'r1 =',e15.8,2x,'j =',i6,/)
      stop 1017
 32   return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  hsens===
c cray .endif
      subroutine hsens(t4n,hsb,hsc,hse,hsf,j,iilw)
c
c +++ bcmt3.for +++
c
      parameter(np=18,npp=19,nilw=7,nhd=2)
c _dp implicit double precision (a-h,o-z)
      logical fulwri
      logical       season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
c
      common/varl00/season,ocen,ekm,tsclim,taclim,simpir,inpu,turbu
      common/varl01/iprint
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
c
      common/time1/tsm,mo,momt,monpr,mm,mts,ian,iday,nts
      common/time2/delt,tht,it
      common/geom1/ca,cp,pi,dg,dels,ds2,gs,gs2,jp,jpp
      common/gr1/deggr1(np),dgcgr1(npp),s10gr1(npp),c10gr1(npp)
c
c --- variables dynamiques
c
      common/dy2/t2mdy2,t2dy2(np),t2ody2(np),dt2dy2(npp),th2dy2(np),
     .           pstdy2(np),htdy2(np),shady2(np),htady2(6,np)
      common/dy5/akdy5(npp,nhd),akody5(npp,nhd)
      common/z00/z2z00(np),z3z00(np),t3z00(np)
      common/p00/p4p00,p3p00,p2p00,p1p00
      common/p01/szp01(np,nilw),spp01(np,nilw)
c
c --- variables thermodynamiques
c
      common/th1/cpath1,rath1,akth1
      double precision  chhv0,cshv0,htshv0
      common/hv0/chhv0(6,np,nilw),cshv0(6,np,nilw),htshv0(6,np,nilw)
      common/so1/sotso1(np),cmuso1(np),bhso1(np,6)
      common/ir4/c3ir4(np,nilw)
      common/cl1/clcl1(np),clhcl1,ztcl1(np),zbcl1(np),
     .                            ptcl1(np),pbcl1(np)
      common/ev1/c1ev1(np),q0ev1(np),rhev1(np,nilw),rhzev1(np,nilw)
c
c --- variables surface
c
      common/su3/hswsu3(np,nilw),cksu3(np,nilw)
c
      fulwri =.true.
c
c --- calcul flux de chaleur sensible (param. de saltzman        68)
c
c     hs0=-(4.02*(t4n-t2dy2(j))-95.02)
c
c --- calcul flux de chaleur sensible (param. de saltzman & ashe 76)
c                            (w.m-2) ( --> (w.cm-2) en fin de routine )
c
      cgr=    4.8e-3
c ... cgr= gam-gamc
      rcp= 1004. *spp01(j,iilw) /p4p00
c ... rcp= rho*cp
      eps= 1200.
      ptpbl  = spp01(j,iilw) - (p4p00 - p3p00)
      call ttzz(ptpbl,t2dy2(j),z2,ttpbl,ztpbl)
      akm = ( 1. + eps * (t4n-ttpbl) / (ztpbl-szp01(j,iilw)) )
      if (t4n.le.ttpbl) akm=1.
      a  = - rcp * akm * cgr
      b  =   rcp * akm / ztpbl
      hs1= a + b * ( t4n - ttpbl )
c
c  -- dt2dy et amu
      dt2dy   = dt2dy2(j) / dels
c
      if (sotso1(j).eq.0.) then
       amu = 0.
      else
       amu = chhv0(1,j,iilw)/sotso1(j)
      end if
c
      a3 = 4.62
      a4 = 316.
      c1 =  .254
      c2 = 4.95e-5
      c7 = rhzev1(j,iilw) * 610.78
c
      if (t4n.lt.273.16) then
       c6 = c7 * 8.239e-2
      else
       c6 = c7 * 7.277e-2
      end if
      c8 = a4 *  c2 * c6
      c9 = a3 * (c1 - c2*c7)
c
      c   =   rcp * eps / ( ztpbl**2 )
      avu = b * ( 1. + hse * hswsu3(j,iilw) )
      gam = avu + ( c9 - c8 ) * ( 1. - c3ir4(j,iilw) *(clcl1(j)**2) )
c
      ck2 = cksu3(j,iilw)  * 2.28e-3
c...  ck2 = cksu3(j,iilw)  * sqrt(omg2/2.)
c
      bet = avu / sqrt( (gam+ck2)**2 + ck2**2 )
c     stz2=((akody5(j,1)+akody5(j+1,1)+akody5(j,2)+akody5(j+1,2))/8.e-6)
c    .            *(dt2dy**2)
      stz2= 1.e12 *(dt2dy**2)
      hs2 = c * ( (1. - bet)**2 ) * stz2
c
      som = 0.
      do 30 n=1,6
      ck1 = cksu3(j,iilw)*sqrt(n*3.64e-5)
      som = som+(bhso1(j,n)**2) / ((gam+ck1)**2 + ck1**2)
   30 continue
      som = som * (amu**2)/2.
      hs3 = som * c
c
      hsc = a - b *(ttpbl -t2dy2(j)) + hs2 + hs3
      hsb = b
c
c  -- parametrisation de saltzman, 1968
c     hsc = - 95.02
c     hsb =    4.02
c
c  -- tuning d'oa78
c     conv = 1.e4 * 4.1855 / 86400.
c     hsc  = -160.1  * conv
c     hsb  =    8.31 * conv
c
c
      if (.not.fulwri) go to 40
      if (iprint.eq.1) go to 41
      if (mts.ne.0)    go to 40
      if (monpr.gt.mm) go to 40
 41   continue
      if (j.eq.1.and.iilw.eq.2) 
     . write(iwr1,60)it
 60    format(/,' -- bcmt3 --  it =',i5,/,
     . '  j  i  z surf  t surf   z pbl   t pbl',5x,'k',5x,'w',
     . 4x,'c3',8x,'a',8x,'b',6x,'hs1',6x,'hs2',6x,'hs3',6x,'hst')
      if (           iilw.eq.5) then
       hst = hs1 + hs2 + hs3
       write(iwr1,61)j,iilw,szp01(j,iilw),t4n,ztpbl,ttpbl,akm,
     . hswsu3(j,iilw),c3ir4(j,iilw),a,b,hs1,hs2,hs3,hst
 61    format(2i3,2(f8.1,f8.2),3f6.2,6f9.2)
      end if
 40   continue
      return
      end
c cray .if [ $concord .eq. 1 ]
c cray nzyx  hlat====
c cray .endif
      subroutine hlat(u4j,u4j1,rh,t4,c3,c4,j,iilw)
c
c +++ bcmt4.for +++
c
      parameter(np=18,nilw=7)
c _dp implicit double precision (a-h,o-z)
c
c     this subroutine calculates the flux of latent heat
c     at the surface
c
      real li,lw,ll
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/th1/cpath1,rath1,akth1
      common/p01/szp01(np,nilw),spp01(np,nilw)
      dimension vg(nilw)
      data vg/2.5,2.5,1.5,2.,2.,2.,2./
c
      li=2.834e+06
c ...    li=latent heat of sublimation (j/kg)
      lw=2.5e+06
c ...    lw=latent heat of vaporisation (j/kg)
      ratio=0.622
c ...    ratio=ratio of molecular weight of water vapor to
c ...    that of dry air
c
      rh75=.75
c ... valeur humidite relative proposee par saltzman,1980
c
      ai=21.87456
      bi=7.66
      aw=17.2694
      bw=35.86
c
      if (t4.le.273.16) then
c     calculation of saturation vapor pressure:
c     -----> over ice:
       es=610.78*exp(ai*(t4-273.16)/(t4-bi))
c
       ll=li
      else
c     calculation of saturation vapor pressure:
c     -----> over ocean:
       es=610.78*exp(aw*(t4-273.16)/(t4-bw))
c
       ll=lw
      end if
c
c     calculation of specific humidity at surface:
      qs=(ratio*es)/spp01(j,iilw)
c
c     calculation of specific humidity at 10 m:
      q10=qs*rh
c
c     calculation of air density:
      rhoa=spp01(j,iilw)/(t4*rath1)
c
      aa = 2.e-3
      bb = 2.e-3
c
      c3 = (ratio**2) * ( (ll**2) / (rath1*cpath1) ) * es
     .   / ( (t4**2) * spp01(j,iilw) )
      c4 = - ll * rhoa * (aa + bb *vg(iilw)) * ratio * es * (1.-rh75)
     .   / spp01(j,iilw)
c
c  -- valeurs dans saltzman,68
c     c3 =   1.27
c     c4 = -38.81
c
      return
      end
c
c     ------------------------------------
      function palco2 (apal)
c     ------------------------------------
c     written by M. Crucifix : 
c     reads co2 from EPICA data
c     ------------------------------------
      real  palco2
      real apal
      real  co2ice(800),f1
      integer i1
c     estimates past CO2 
      open (6001,file="co2ice_800.dat",status="old")
      read (6001,*) co2ice
      i1 = int(apal+800)
      f1 = apal+800-i1
      palco2 = co2ice(i1)+f1*(co2ice(i1+1)-co2ice(i1))
      write(*,*) 'palco2=',palco2

      end

c     ------------------------------------
      function co2inter (ppmgz2, volum, perh, so, ecc)
c     ------------------------------------
c     written by M. Crucifix : 
c     input parameters: 
c      ppmgz2 : CO2 concentration at previous timestep 
c      volum  : total ice volume, en 1e15m3 
c      perh   : longitude of perihelion,  expressed in degrees
c               so need to multiply by pi / 180. to get radians
c      so     : sine of obliquity
c      ecc    : eccentricity
c     parameter : 
c      tau   :  reponse time in timesteps
c     output: 
c      co2inter : CO2 concentration for the present time step
c     ------------------------------------
      real   tau
      real   co2inter, ppmgz2, volum, perh, so, ecc
      real   co2star

      parameter (tau = 1.)

      co2star = 280. - 90. / 51. * volum 
      co2inter = ppmgz2 + (1./tau) * (co2star - ppmgz2)
      return 
      end


