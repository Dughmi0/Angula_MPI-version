c ============================================================
c  ANGULA - MPI VERSION
c  Parallelization over central molecules (imf loop).
c  Each MPI rank processes a subset of central molecules
c  within each configuration. Results are gathered to rank 0
c  at the end of each configuration via MPI_Gather.
c  Best when you have few but very large configurations.
c
c  Compile:  mpif90 -O2 -o angula_mpiB angula_mpi_optionB.for
c  Run:      mpirun -np 8 ./angula_mpiB control.con
c ============================================================

      include 'calc.for'
      include 'energy.for'
      include 'geteuler.for'
      include 'help.for'
      include 'molselec.for'
      include 'ordering.for' 
      include 'plot.for'
      include 'pesc2.for'
      include 'pev2.for' 
      implicit none
c ------ MPI header ------
      include 'mpif.h'
c ------------------------

      
      real Xaxis(3),Yaxis(3),Zaxis(3),centre(3)
      real box(3),rmin,rmax,dr,
     & xc1,yc1,zc1,xc2,yc2,zc2,d,dumr,
     & dvinf,dvsup,rv,dx,dy,dz,
     & modcc1,modcc2,
     & dxc,dyc,dzc,dc,dcsup,dcinf,drc,dangmin,alfamin,alfamax,
     $ xcc,ycc,zcc,pvx1,pvy1,pvz1,pvx2,pvy2,pvz2,alfaxyccl,
     & angulo,corr,psc,nmolvec,pi,
     & cosclcl,modclcl,xclv,yclv,zclv,alfaclclmin,dum,
     $ cosalfaxy,pv,dadd(10),rg,angbiv(100),
     $ rgini1,rgini2,rgi,rgj,stang,totalg,ggbiv,
     & simupesc,rpsc,rpscold,dist(100)
      real*8 eulera(1:5),dcmin
      character*50 dumline
      character*40 nomin,nomout,limitar,nomcontrol,nomcfg(10000)
     & ,nompesc(100),nompev(100),nomeuler(5),namebiv(5),dumtext,
     & anglefile,dumchar,nombre,nomd(100),nommixture,nommolselec,
     & nomHB
      character*20 nom,ext
      integer np,nt,nm,i,j,k,c,a(10),ai,aj,ip,jp,iprogram,iii,iCM,jCM,
     & central,vecino,ndv,idv,ipin,im,jm,ia,ja,imax,nadmin,ii,jj,
     $ iangmin1,iangmin2,pp,ia2,nconf,npesc,npev,nvcadd,
     $ iadd(100),
     $ molpesc(50,4),atpesc(50,4),molpev(50,2,4),atpev(50,2,4),
     $ moladd(50,2,2),atadd(50,2,4),
     $ oppesc(50),oppev(50),opadd(50),ic,m1,m2,m3,m4,m11,m12,m13,m14,
     $ m21,m22,m23,m24,
     & nstapesc(50),nstapev(50),ias,tpesc(50),tpev(50),
     & ibiv,ib,iangb1,iangb2,
     & nangb1,nangb2,nibiv(5),njbiv(5),opbiv1(5),opbiv2(5),filetype,
     & npchk,igaus,ngaus,opcos(100),iangle,ntt,itt,opeuler(1:5),
     & ae10,ae11,ae12,ae13,ae20,ae21,ae22,ae23,
     & euler,confini,conffin,inom,
     & teuler(1:5),nvec,opd,moli,molf,calculate,icon,
     & numd,mold1(100),mold2(100),atd1(100),atd2(100),fixmol,reverse,
     & nmi,nmj,imf,bigerror,checkPBC
      integer natpm(0:1000),nmtype,nmpt(0:100),
     & molmin(100),molmax(0:100),imt,imm
	 
      integer iataHB(10),jataHB(10),kataHB(10),
     &imolHB(10),jmolHB(10),kmolHB(10),condHB,ncondHB,opHB,
     &iatdHB(10),jatdHB(10),iHB,iHBw
      real angHB(10),dHB(10),angHBc,dHBc
      real anginf(2),angsup(2)
      integer angop(2),tg,tgc,tgv,in
      integer numscal,opQ,iq,opico,iico,opcost,icost,iscal,ieLJ,ieQQ
      real scalar(100)
      character*40 action,dialog,nameplot,nameplot1,nameplot2,arg(0:5)
      integer plotop,ntotang,opnom,idum,jdum
      integer*4 iargc,narg,interactive
      integer mixture,molmixi,molmixj,nmolmix
      integer togc,togv,dumi,sep,plotxyz
      integer gir,anglegiri,anglegirj,chkerr,
     & nmolimolj,ffo
      integer opdx,opdy,opdz,moldx1,moldx2,atdx1,atdx2,
     & moldy1,moldy2,atdy1,atdy2,moldz1,moldz2,atdz1,atdz2,
     & idx,idy,idz
      integer nmolshow,imolshow(100),iat,jat,jmin,imin,nstat
      integer opimin,opjmin 
      real resc,Q,cost,ico,ELJ,Ee
      integer imcost,jmcost,iimin,ijmin
      real angtest
      logical moldist
      logical exists
      character*10 nomv(200)

c ****** ENERGY CALCULATION   
      character*40 nomenergy
      real charge(10,1000),E
      integer energyop,iE
      integer il,jl,nerren,erren,nerreu,erreu
	  
c **************hughe arrays
      integer dimat,dimmol,dimEmix,dimEat
      parameter (dimat=300,dimmol=30000,dimEmix=4,dimEat=300)
      real r(3,dimat,dimmol)
      character*3 nomat(dimat,dimmol)
      integer molmix(dimmol),natt(dimmol),order(dimmol,200)
      integer molimolj(2,dimmol)
      real LJ(2,dimEmix,dimEmix,dimEat,dimEat)

c ------ MPI variables ------
      integer myrank, nprocs, ierr
      integer nadmin_local, nadmin_global
c
c  Output buffer for gather:
c  Each (im,jm) pair produces one line: 2 integers + up to 102 reals.
c  We pack them as reals: im_r, jm_r, dc, dist(:), scalar(:), angbiv(:)
c  maxline = 2 + 1 + 100 + 100 + 100 = 303 — use 400 for safety.
      integer, parameter :: MAXLINE=400
      integer, parameter :: MAXPAIRS=500000  ! max pairs per rank per conf
      real,   allocatable :: outbuf_local(:,:)   ! local pair lines
      integer npairs_local                    ! # lines written by this rank
      real,   allocatable :: outbuf_recv(:)   ! gathered packed buffer on root
      integer, allocatable :: npairs_all(:)   ! lines per rank
      integer, allocatable :: recvcounts(:), displs(:)
      integer iline, ipair, irank, nline
      integer sendcount, totalrecv, idx0
c ---------------------------

c ============================================================
c  MPI INITIALISATION
c ============================================================
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

      allocate(outbuf_local(MAXLINE,MAXPAIRS))
      allocate(npairs_all(0:nprocs-1))
      allocate(recvcounts(0:nprocs-1))
      allocate(displs(0:nprocs-1))

c ---- Only rank 0 opens the shared log files ----
      if (myrank .eq. 0) then
       open(unit=1,file='Erren.log')
       write(1,*)'Log file for problems with Energy Calculations'
       write(1,*)''
       write(1,*)''
       close(unit=1)
       open(unit=1,file='Erreu.log')
       write(1,*)'Log file for problems with Euler'
       close(unit=1)
      end if

      checkPBC=0
      reverse=0
      opHB=0
      opd=0
      opcos=0
      energyop=0
      opQ=0
      opico=0
      opcost=0
      resc=1.0
      nmolshow=1
      imolshow(1)=1
      chkerr=0
      gir=0
      anglegiri=8
      anglegirj=9
      plotxyz=1
      bigerror=0
      if(chkerr.eq.1) write(6,*) ''

      narg=iargc()
      do i=0,3
      call getarg(i,arg(i+1))
      end do

      if(narg.eq.1.or.narg.eq.2) then
      nomcontrol=arg(2) 
      if(nomcontrol.eq.'c') nomcontrol='contang.con'
      if(nomcontrol.eq.'t') nomcontrol='test.con'
      inquire(file=nomcontrol,exist=exists)
      if(.not.exists) then
      if (myrank .eq. 0) write(6,*)'The file does not exist!'
      end if
      if (myrank .eq. 0)
     &  write(6,*)'You will work with control file: ',nomcontrol
      end if
      if(narg.eq.1) then
      interactive=0
      end if

      PI=3.14159265358979323846264

c ---- Only rank 0 prints the banner ----
      if (myrank .eq. 0) then
      WRITE(6,*)'***********************************************************************************************'
      WRITE(6,*)'*      AAAA   N   N   GGGG   U  U  L       AAAAA 		MMM        MMM    PPPPPPPPP    IIIIII  *'
      WRITE(6,*)'*      A  A   NN  N  G       U  U  L       A   A 		M  M      M  M    PP     PP      II    *'
      WRITE(6,*)'*      AAAA   N N N  G  GG   U  U  L       AAAAA 		M   M   M    M    PPPPPPPPP      II    *'
      WRITE(6,*)'*      A  A   N  NN  G   G   U  U  L       A   A 		M    M M     M    PP             II    *'
      WRITE(6,*)'*      A  A   N   N   GGGG   UUUU  LLLL    A   A 		M     M      M    PP           IIIIII  *'
      WRITE(6,*)'*                                               											   *'
      WRITE(6,*)'*                                                                                             *'
      WRITE(6,*)'*       ANGULA MPI Version                                                                    *'
      WRITE(6,*)'***********************************************************************************************'
      WRITE(6,*)' '
      WRITE(6,*)''
      WRITE(6,*)'      Dr. Luis Carlos Pardo '
      WRITE(6,*)'      Grup de Caracteritzacio de Materials (UPC) '
      WRITE(6,*)'      for bugs or comments send mail to '
      WRITE(6,*)'      luis.carlos.pardo_upc.edu '
      WRITE(6,*)''		
      WRITE(6,*)''
	  WRITE(6,*)''
	  WRITE(6,*)''
	  WRITE(6,*)'		MPI extension by: '
      WRITE(6,*)'       Rashed M. R. Aldughmi '
	  WRITE(6,*)''
      WRITE(6,*)'****************WARNING**********************'
      WRITE(6,*)'Max number of atoms per molecule:  100'
      WRITE(6,*)'Max number of molecules:           20000'
      WRITE(6,*)'****************WARNING**********************'
      WRITE(6,*)'Max number of species Energy calc: 4'
      WRITE(6,*)'Max number of atoms Energy calc:   100'
      WRITE(6,*)'(To change this contact me, or change the code)'
      WRITE(6,*)'*********************************************'
      WRITE(6,*)''
      WRITE(6,*)''
      WRITE(6,*)''
      WRITE(6,"(a25,i4)")'MPI ranks (Option B):  ',nprocs
      WRITE(6,*)''
      end if

111   continue
      if(narg.eq.0) then
      if (myrank .eq. 0) then
      write(6,*)'Angula control file to use: (h=help,dir)'
      read(5,*) nomcontrol
      end if
c ---- broadcast the chosen control file to all ranks ----
      call MPI_Bcast(nomcontrol,40,MPI_CHARACTER,0,
     &               MPI_COMM_WORLD,ierr)

      if(nomcontrol.eq.'dir'.or.nomcontrol.eq.'ls') then
      if (myrank .eq. 0) call system('dir *.con')
      goto 111
      end if
      if(nomcontrol.eq.'cc') nomcontrol='c2cl6.con'
      if(nomcontrol.eq.'c') nomcontrol='contang.con'
      if(nomcontrol.eq.'a') nomcontrol='agua.con'
      if(nomcontrol.eq.'t') nomcontrol='test.con'
      if (nomcontrol.eq.'h'.or.nomcontrol.eq.'help') then
      if (myrank .eq. 0) then
      call help(1)
      call help(2)
      call help(3)
      end if
      goto 111
      end if
      inquire(file=nomcontrol,exist=exists)
      if(.not.exists) then
      if (myrank .eq. 0)
     &  write(6,"(a20,a20)")nomcontrol,'does not exist' 
      goto 111
      end if
      end if
666   continue

      ntotang=0
      open(unit=1,file=nomcontrol)
      iangle=1

      read(1,*) nomin,interactive,sep,mixture,fixmol,chkerr,nstat
	  
      if(fixmol.lt.0) then
      fixmol=abs(fixmol)
      reverse=1
      end if
      if(chkerr.eq.1.and.myrank.eq.0) write(6,"(a10,6i4)") 
     & nomin,interactive,sep,mixture,fixmol,chkerr,nstat
      read(1,*) nomv(1:nstat)
      do i=1,nstat
      nomv(i)=trim(adjustr(nomv(i)))
      end do
      if(chkerr.eq.1.and.myrank.eq.0) write(6,*) nomv(1:nstat)
      if(chkerr.eq.1.and.myrank.eq.0)
     &  write(6,*)'Error check activated'
      
      read(1,*) molmixi,molmixj
      if (myrank .eq. 0)
     &  write(6,"(a33,i5,a24,i5)") 
     &  'CENTRAL molecule will be: ',molmixi,
     &  ' neighbor will be: ',molmixj
	  
      if(mixture.eq.2) then
       read(1,*)nommixture
       if (myrank .eq. 0)
     &   write(6,*) 'FIXED AXES OPTION ACTIVATED: ', nommixture
       inquire(file=nommixture,exist=exists)
        if(.not.exists) then
        if (myrank .eq. 0) then
        write(6,*)nommixture,'Fixed axes does not exist. Puit!'
        write(6,*)
     & 'You must add it after the line with central and neigh molecules'
        end if
        call MPI_Finalize(ierr)
        stop
        end if
       open(unit=22,file=nommixture)
       read(22,*) dumchar
       read(22,*) Centre(1),Centre(2),Centre(3)
       read(22,*) Xaxis(1),Xaxis(2),Xaxis(3)
       read(22,*) Yaxis(1),Yaxis(2),Yaxis(3)
       read(22,*) Zaxis(1),Zaxis(2),Zaxis(3)
       close(unit=22)
      end if
	  
      if(fixmol.ge.1) then
      read(1,*)nommolselec
      if (myrank .eq. 0) then
      write(6,*) ''
      write(6,*) 'Fix molecules option activated'
      end if
      inquire(file=nommolselec,exist=exists)
      if (myrank .eq. 0)
     &  write(6,*) 'Fix molecules file: ',nommolselec
      if(.not.exists) then
       if (myrank .eq. 0) then
       write(6,*)'The following fix molecules file does not exist:'
       write(6,*)nommolselec
       end if
       call MPI_Finalize(ierr)
       stop
      end if
      end if

      read(1,*) anglefile
      if (myrank .eq. 0)
     &  write(6,*)'Angles will be written in: ',anglefile

      read(1,*) dcinf,dcsup,icm,jcm
      if (myrank .eq. 0)
     &  write(6,"(a21,2f12.6,2i5)")'dcinf,dcsup,icm,jcm: ',
     &  dcinf,dcsup,icm,jcm

      if (myrank .eq. 0) then
      write(6,*)''
      write(6,*)'--------------  ATOMS SELECTION -------------'
      write(6,*)'Atoms to calculate distances: ',icm,jcm
      if(dcinf.ge.0) write(6,*)'Distancia CC minima,maxima',
     & dcinf,dcsup
      end if
      
      if(dcinf.lt.0) then
       moli=int(abs(dcinf))
       molf=int(abs(dcsup))
       opd=1	   
       if (myrank .eq. 0) then
       write(6,*)'ORDERING SELECTION'
       write(6,*)'Initial and final molecule: ',moli,molf
       end if
      end if
	  
      if(icm.eq.0.and.jcm.ne.0) then
       if (myrank .eq. 0) write(6,*)'SURFACE ATOM SELECTION'
       opd=2
       icm=1
      end if	
	  
      if(jcm.eq.0.and.icm.ne.0) then
       if (myrank .eq. 0) write(6,*)'ATOM SURFACE SELECTION'
       opd=3
       jcm=1
      end if
	  
      if(icm.eq.0.and.jcm.eq.0) then
       if (myrank .eq. 0) then
       write(6,*)'SURFACE SURFACE SELECTION'
       write(6,*)'SINGLE NEIGHBOUR OPTION'
       end if
       opd=4
       icm=1
       jcm=1
      end if

      read(1,*) numd
      if (myrank .eq. 0) write(6,*)'Distances to calculate',numd
      do i=1,numd
      read(1,*) nomd(i)
      if (myrank .eq. 0) write(6,*) nomd(i)
      read(1,*) mold1(i),mold2(i),atd1(i),atd2(i)
      if (myrank .eq. 0) write(6,*) mold1(i),mold2(i),atd1(i),atd2(i)
      end do
      if (myrank .eq. 0) write(6,*)''

      read(1,*) numscal
      iscal=1
      do i=1,numscal
      read(1,*)dumtext
      opdx=0
      opdy=0
      opdz=0
      if(dumtext(1:2).eq.'dx')then
      idx=iscal
      opdx=1
      if (myrank .eq. 0)
     &  write(6,*)'Will calculate distance projection in x direction'  
      read(1,*) moldx1,moldx2,atdx1,atdx2	  
      end if
      if(dumtext(1:2).eq.'dy')then
      idy=iscal
      opdy=1
      if (myrank .eq. 0)
     &  write(6,*)'Will calculate distance projection in y direction' 
      read(1,*) moldy1,moldy2,atdy1,atdy2	  
      end if
      if(dumtext(1:2).eq.'dz')then
      idz=iscal
      opdz=1
      if (myrank .eq. 0)
     &  write(6,*)'Will calculate distance projection in z direction'  
      read(1,*) moldz1,moldz2,atdz1,atdz2	  
      end if
      opimin=0
      opjmin=0
      if(dumtext(1:2).eq.'im'.or.dumtext(1:3).eq.'iat')then
      opimin=1
      iimin=iscal
      if (myrank .eq. 0) then
      write(6,*)'The atom for central molecule'
      write(6,*)'used to calculate distance'
      write(6,*)'will be displayed'
      end if
      end if
      if(dumtext(1:2).eq.'jm'.or.dumtext(1:3).eq.'jat')then
      opjmin=1
      ijmin=iscal
      if (myrank .eq. 0) then
      write(6,*)'The atom for neighbour molecule'
      write(6,*)'used to calculate distance'
      write(6,*)'will be displayed'
      end if
      end if
      if(dumtext(1:2).eq.'qf'.or.dumtext(1:2).eq.'Qf')then
      opQ=1
      iQ=iscal
      if (myrank .eq. 0) write(6,*)'Qfactor will be calculated'
      end if
      if(dumtext(1:4).eq.'cost')then
      if (myrank .eq. 0) write(6,*)'Costheta will be calculated'
      opcost=1
      icost=iscal
      read(1,*) imcost,jmcost
      if (myrank .eq. 0) then
      open(unit=77,file='costheta.dat')
      write(77,"(3a5,a15)") 'imol','im','jm','psc'
      close(unit=77)
      end if
      end if
      if(dumtext(1:3).eq.'ico') opico=1
      if(dumtext(1:3).eq.'ico') iico=i
      if(dumtext(1:2).eq.'HB'.or.dumtext(1:2).eq.'hb') then
       opHB=1
       iHBw=iscal
       read(1,*) nomHB
       inquire(file=nomHB,exist=exists)
        if(.not.exists) then
        if (myrank .eq. 0) then
        write(6,*)'The following HB file does not exist!'
        write(6,*)'*******',nomHB
        end if
        call MPI_Finalize(ierr)
        stop
       end if
       if (myrank .eq. 0) write(6,*)'Opening HB file:',nomHB
       open (unit=22,file=nomHB)
       read(22,*) dumtext
       read(22,*) NcondHB
       do ihb=1,ncondHB
       read(22,*) dumtext
       read(22,*) iatdHB(ihb),jatdHB(ihb),dHB(ihb)
       read(22,*) imolHB(ihb),jmolHB(ihb),kmolHB(ihb)
       read(22,*) iataHB(ihb),jataHB(ihb),kataHB(ihb),angHB(ihb)
       end do
       close(unit=22)
       if (myrank .eq. 0) write(6,*)'Closing HB file:',nomHB
      end if
      if(dumtext(1:6).eq.'energy') then
      read(1,*) nomenergy
      if (myrank .eq. 0) write(6,*) 'Energy will be calculated'
      iE=iscal
      iELJ=iscal+1
      iEQQ=iscal+2
      iscal=iscal+2
      energyop=1
      inquire(file=nomenergy,exist=exists)
      if(.not.exists) then
      if (myrank .eq. 0) then
      write(6,*)'The following energy file does not exist!'
      write(6,*)'*******',nomenergy
      end if
      call MPI_Finalize(ierr)
      stop
      end if
      open(unit=22,file=nomenergy)
      read(22,*) dumchar
      if (dumchar(1:7).eq.'closest') energyop=2
      read(22,*) nmolmix,ffo
      if (myrank .eq. 0) then
      if(ffo.eq.1) write(6,*)'LJ parameters will be 6 12'
      if(ffo.eq.2) write(6,*)'LJ parameters will be sigma eps'
      write(6,*)'Units are those of the input parameters ;-)'
      end if
      do ii=1,nmolmix
      read(22,*) natpm(ii)
      end do	
      do ii=1,nmolmix
       do jj=ii,nmolmix      
       read (22,*) dumchar
        do il=1,natpm(ii)
         do jl=1,natpm(jj)
         read(22,*) idum,jdum,LJ(1,ii,jj,il,jl),LJ(2,ii,jj,il,jl)
         LJ(1,jj,ii,jl,il)=LJ(1,ii,jj,il,jl)
         LJ(2,jj,ii,jl,il)=LJ(2,ii,jj,il,jl)
         end do
        end do
       end do
      end do
      read(22,*)dumchar
      do ii=1,nmolmix
       read(22,*)dumchar
      do il=1,natpm(ii)
      read(22,*) idum,charge(ii,il)
      end do
      end do
      close(unit=22)
      end if
      iscal=iscal+1
      end do

      numscal=iscal-1
      if (myrank .eq. 0)
     &  write(6,*)'I will calculate scalars: ',numscal

      read(1,*) npesc
      if (myrank .eq. 0) write(6,*)'Num, Pesc to calculate',npesc
      if(npesc.gt.100) then
       if (myrank .eq. 0) write(6,*) 'Too many PSC!'
       call MPI_Finalize(ierr)
       stop
      end if
      ic=1
      do ii=1,npesc
      read(1,*) nompesc(ii),oppesc(ii)
      if(oppesc(ii).eq.1) opcos(iangle)=1
      iangle=iangle+1
      if (myrank .eq. 0)
     &  write(6,"(a20,2x,a20,1x,i5)")'Statistics on PESC',
     &  nompesc(ii),oppesc(ii)
      if(oppesc(ii).ne.1.and.oppesc(ii).ne.2) then
       if (myrank .eq. 0) write(6,*) 'ERROR deg/cos'
       call MPI_Finalize(ierr)
       stop
      end if
      read(1,*) molpesc(ic,1),molpesc(ic,2),atpesc(ic,1),atpesc(ic,2)
      read(1,*) molpesc(ic,3),molpesc(ic,4),atpesc(ic,3),atpesc(ic,4)
      ic=ic+1
      end do
      ntotang=ntotang+ic-1

      read(1,*) npev
      if (myrank .eq. 0) write(6,*)'Num, Pev to calculate',npev
      if(npev.gt.100) then
       if (myrank .eq. 0) write(6,*) 'Too many PV!'
       call MPI_Finalize(ierr)
       stop
      end if
      ic=1
      do ii=1,npev
      read(1,*) nompev(ii),oppev(ii)
      opcos(iangle)=oppev(ii)
      iangle=iangle+1
      if (myrank .eq. 0)
     &  write(6,"(a20,2x,a20,1x,i2)")'Statistics on PEV',
     &  nompev(ii),oppev(ii)
      if(oppev(ii).gt.4) then
       if (myrank .eq. 0) write(6,*) 'ERROR deg/cos'
       call MPI_Finalize(ierr)
       stop
      end if
      read(1,*) molpev(ic,1,1),molpev(ic,1,2),
     &      molpev(ic,1,3),molpev(ic,1,4),
     &  atpev(ic,1,1),atpev(ic,1,2),atpev(ic,1,3),atpev(ic,1,4)
      read(1,*) molpev(ic,2,1),molpev(ic,2,2),
     &      molpev(ic,2,3),molpev(ic,2,4),
     &  atpev(ic,2,1),atpev(ic,2,2),atpev(ic,2,3),atpev(ic,2,4)
      ic=ic+1
      end do
      ntotang=ntotang+ic-1

      read(1,*) euler
      if(euler.gt.0) then
      if (myrank .eq. 0) then
      write(6,*)''
      write(6,*)'I will calculate euler angles'
      end if
      read(1,*) ae10,ae11,ae12,ae13,ae20,ae21,ae22,ae23
      if(mixture.eq.2) then
      if (myrank .eq. 0) then
      write(6,*)'AXIS OPTION!!!!!'
      end if
      ae10=1
      ae11=2
      ae12=3
      ae13=4
      end if
      read(1,*) nomeuler(1),opeuler(1)
      opcos(iangle)=opeuler(1)
      if(opeuler(1).eq.2)  opcos(iangle)=3
      iangle=iangle+1
      read(1,*) nomeuler(2),opeuler(2)
      opcos(iangle)=opeuler(2)
      if(opeuler(2).eq.2)  opcos(iangle)=3
      iangle=iangle+1
      read(1,*) nomeuler(3),opeuler(3)
      opcos(iangle)=opeuler(3)
      if(opeuler(3).eq.2)  opcos(iangle)=3
      iangle=iangle+1
      read(1,*) nomeuler(4),opeuler(4)
      opcos(iangle)=opeuler(4)
      if(opeuler(4).eq.2)  opcos(iangle)=3
      iangle=iangle+1
      read(1,*) nomeuler(5),opeuler(5)
      opcos(iangle)=opeuler(5)
      if(opeuler(5).eq.2)  opcos(iangle)=3
      ntotang=ntotang+5
      end if
      
      if (myrank .eq. 0)
     &  write(6,*)'TOTAL number of angles: ',ntotang
      if (myrank .eq. 0) write(6,*)''

      write(6,*)
      write(6,*)'-----------------------------------------------------'
      read(1,*) dumtext
      read(1,*)confini,conffin,filetype
      if (myrank .eq. 0) then
      write(6,*)'Starting conf',confini
      write(6,*)'ending conf',conffin
      end if
      Nconf=conffin-confini+1
      if (myrank .eq. 0) write(6,*)'Nconf= ',Nconf
      if(filetype.eq.1.and.myrank.eq.0)
     &  write(6,*)'File type: NEW ANGULA FORMAT'
      if(filetype.ge.2) then
       if (myrank .eq. 0) then
       write(6,*)''
       write(6,*)'Only ANGULA cfg files are accepted'
       end if
       call MPI_Finalize(ierr)
       stop
      end if

      opnom=0
      do i=1,Nconf
      read(1,*) nomcfg(i)
      if(nomcfg(1).eq.'automatic') opnom=1
      if(opnom.eq.1) goto 555
      end do

555   continue

      if(opnom.eq.1) then
      read(1,*) nom,ext
      nom=adjustr(nom)
      ext=adjustl(ext)
      do i=1,Nconf
      inom=(i-1)+confini
      nombre='      '
      if(i.le.9)
     &  write(nombre,"(a20,i1,a1,a4)") nom,inom,'.',ext
      if(inom.ge.10.and.i.le.99)
     &  write(nombre,"(a20,i2,a1,a4)") nom,inom,'.',ext
      if(inom.ge.100.and.i.le.999)
     &  write(nombre,"(a20,i3,a1,a4)") nom,inom,'.',ext
      if(inom.ge.1000.and.i.le.9999)
     &  write(nombre,"(a20,i4,a1,a4)") nom,inom,'.',ext
      nombre=adjustl(nombre)
      nomcfg(i)=nombre
      end do
      end if

      close(unit=1)
	  
      if (myrank .eq. 0) then
      write(6,*)''
      write(6,*)'Labelling molecules...'
      end if
	  
      open (unit=1,file=nomcfg(1))
      read(1,*) nm,nmtype
      read(1,*) box(1:3)
      read(1,*) nmpt(1:nmtype)
      read(1,*) natpm(1:nmtype)
      close(unit=1)
      iii=1	
      if (myrank .eq. 0)
     &  write(6,*)'Number of different moleculer types:',nmtype	   
      do imt=1,nmtype
       if (myrank .eq. 0) then
       write(6,*)''
       write(6,"(a20,i6)")'*** Molecule type:  ',imt
       write(6,"(a20,i6)")'Number of atoms: ',natpm(imt)
       write(6,"(a20,i6)")'First molecule: ',iii
       end if
       do im=1,nmpt(imt)
        molmix(iii)=imt
        natt(iii)=natpm(imt)
        iii=iii+1		
       end do      
       if (myrank .eq. 0) write(6,"(a20,i6)")'Last molecule: ',iii-1
      end do
      if (myrank .eq. 0) then
      write(6,*)'End of labeling molecules'	       
      write(6,*)''	       
      end if

      nmolvec=0
      if(nstat.ne.numd+numscal+ntotang) then
      if (myrank .eq. 0) then
      write(6,*) '********************'
      write(6,*) '*ERROR, nstat wrong*'
      write(6,*) '********************'
      write(6,*) ' Nstat = ',nstat
      write(6,*) ' Nd+Nang = ',numd+numscal+ntotang
      end if
      call MPI_Finalize(ierr)
      stop
      end if

c ---- Only rank 0 opens/writes the single output file ----
      if (myrank .eq. 0) then
      open(unit=9,file=anglefile)
      write(9,"(a17,a25,a19,i5)") 
     &'Angle file from: ',nomcfg(1),' # configurations: ',nconf
      write(9,*) numd+numscal+1,ntotang
      write(9,"(200a12)") 'central','neighbor','d(cond)',nomv(1:nstat)
      end if

c ============================================================
c  MAIN LOOP — ALL RANKS process ALL configurations.
c  Within each config, the imf (central molecule) loop is
c  split round-robin across ranks.
c  After each config, rank 0 gathers and writes results.
c ============================================================
      do iprogram=1,nconf

      if (myrank .eq. 0) then
      write(9,"(a20,i10)")'*************Conf=  ',iprogram
      open(unit=33,file='Erreu.log',access='append')
      write(33,*)''
      write(33,*)'Configuration: ',iprogram
      close(unit=33)
      end if

      if(fixmol.ge.1) 
     & call molselec
     & (chkerr,iprogram,reverse,molimolj,nmolimolj,nommolselec)
      
      inquire(file=nomcfg(iprogram),exist=exists)
      if (myrank .eq. 0) write(6,*)'Opening:    ',nomcfg(iprogram)
      if(.not.exists.and.myrank.eq.0) then
      write(6,*)'File does not exist:',nomcfg(iprogram)
      end if

      if(filetype.eq.1) then
      open (unit=1,file=nomcfg(iprogram))
      read(1,*) nm,nmtype
      read(1,*) box(1:3)
      read(1,*) nmpt(1:nmtype)
      read(1,*) natpm(1:nmtype)
      im=1
      do imt=1,nmtype
       do imm=1,nmpt(imt)
        read(1,"(a50)") dumline
        do ia=1,natpm(imt)
         read(1,*,err=16) nomat(ia,im),r(1,ia,im),r(2,ia,im),r(3,ia,im)
         do c=1,3
         if(r(c,ia,im).gt.box(c)) then
          r(c,ia,im)=r(c,ia,im)-box(c)
          checkpbc=1
         end if
         end do
         do c=1,3
         if(r(c,ia,im).lt.0) then
          r(c,ia,im)=box(c)+r(c,ia,im)
          checkpbc=1
         end if
         end do
         if(imm.ge.2) then
          if(nomat(ia,im).ne.nomat(ia,im-1)) then
          write(6,*)''
          write(6,*)'CHANGE OF ATOM NAME WITHIN A MOLECULE'
          write(6,*)'This is usually a very bad error!'
          write(6,"(3i9,a9)") imt,imm,ia,nomat(ia,im)
          write(6,*)''
          bigerror=1
          end if		
         end if
        end do
        im=im+1
       end do
      end do
      if(checkpbc.eq.1.and.iprogram.eq.nconf.and.myrank.eq.0) then
      write(6,*)'*** Periodic Boundary Conditions corrected ***'
      end if	  
      do i=1,nmtype
      natpm(0)=max(natpm(0),natpm(i))
      end do
      close (unit=1)
      end if

      if (filetype.ge.2) then
      if (myrank .eq. 0) write(6,*) 'Incorrect configuration file'
      call MPI_Finalize(ierr)
      stop
      end if

      if(mixture.eq.2) then
       do c=1,3
       r(c,1,nm+1)=centre(c)   
       r(c,2,nm+1)=Xaxis(c)   
       r(c,3,nm+1)=Yaxis(c)
       r(c,4,nm+1)=Zaxis(c)
       end do
      end if

c ---- Reset local pair buffer for this configuration ----
      nadmin_local=0
      npairs_local=0

      nmi=nm
      nmj=nm
      if(mixture.eq.2) nmi=1

c ============================================================
c  OPTION B: Split the imf loop across ranks (round-robin)
c ============================================================
      do imf=1+myrank, nmi, nprocs

      if(mixture.ne.2) im=imf
      if(mixture.eq.2) im=nm+1

      if(opQ.eq.1) then
      nvec=max(4,molf)	  
      call Qfactor(r,icm,jcm,nvec,im,nm,box,order,Q,
     & mixture,molmixi,molmixj,molmix)
      scalar(iQ)=Q
      end if
      if(opd.eq.1) then
      nvec=molf
      call ordering(r,icm,jcm,nvec,im,nm,box,order,
     &mixture,molmixi,molmixj,molmix)      
      end if
      if(opcost.eq.1) then
      call costeta(r,icm,jcm,imcost,jmcost,im,nm,box,order,cost,
     & mixture,molmixi,molmixj,molmix)
      scalar(icost)=cost
      end if

      do jm=1,nmj  

      if(im.eq.jm.and.mixture.ne.2) goto 50
      if (fixmol.eq.1) then
       do i=1,nmolimolj
       molimolj(2,i)=0
       if(im.eq.molimolj(1,i)) goto 99
       end do
       goto 51
      else if (fixmol.eq.2) then
       do i=1,nmolimolj
       molimolj(1,i)=0
       if(jm.eq.molimolj(2,i)) goto 99
       end do
       goto 50
      else if(fixmol.eq.3) then
       do i=1,nmolimolj
       if(im.eq.molimolj(1,i).and.jm.eq.molimolj(2,i)) then
       xc1=r(1,icm,im)
       yc1=r(2,icm,im)
       zc1=r(3,icm,im)
       xc2=r(1,jcm,jm)
       yc2=r(2,jcm,jm)
       zc2=r(3,jcm,jm)
       xcc=(xc2-xc1)
       ycc=(yc2-yc1)
       zcc=(zc2-zc1)
       xcc=xcc-box(1)*nint(xcc/box(1))
       ycc=ycc-box(2)*nint(ycc/box(2))
       zcc=zcc-box(3)*nint(zcc/box(3))
       dc=sqrt(xcc**2+ycc**2+zcc**2)
       calculate=1
       goto 11
       end if
       end do
       goto 50
      end if
99      continue
      
      if(mixture.ne.2) then
      if(molmix(jm).ne.molmixj) goto 50
      if(molmix(im).ne.molmixi) goto 50
      end if
      if(mixture.eq.2) then
      if(molmix(jm).ne.molmixj) goto 50
      end if

      calculate=0
      xc1=r(1,icm,im)
      yc1=r(2,icm,im)
      zc1=r(3,icm,im)
      xc2=r(1,jcm,jm)
      yc2=r(2,jcm,jm)
      zc2=r(3,jcm,jm)
      xcc=(xc2-xc1)
      ycc=(yc2-yc1)
      zcc=(zc2-zc1)
      xcc=xcc-box(1)*nint(xcc/box(1))
      ycc=ycc-box(2)*nint(ycc/box(2))
      zcc=zcc-box(3)*nint(zcc/box(3))
      dc=sqrt(xcc**2+ycc**2+zcc**2)
      if(opd.eq.0) then 
       if(dc.le.dcsup.and.dc.ge.dcinf) calculate=1
      end if   
      if(opd.eq.1) then 
       do icon=moli,molf
        if(order(im,icon).eq.jm) calculate=1
       end do
      end if   
      if(opd.eq.2) then 
       dcmin=1000000.d0	 
       do iat=1,natt(im)
        if(nomat(iat,im).eq.'Crd') goto 333
         xc1=r(1,iat,im)
         yc1=r(2,iat,im)
         zc1=r(3,iat,im)
         xc2=r(1,jcm,jm)
         yc2=r(2,jcm,jm)
         zc2=r(3,jcm,jm)
         xcc=(xc2-xc1)
         ycc=(yc2-yc1)
         zcc=(zc2-zc1)
         xcc=xcc-box(1)*nint(xcc/box(1))
         ycc=ycc-box(2)*nint(ycc/box(2))
         zcc=zcc-box(3)*nint(zcc/box(3))
         dc=sqrt(xcc**2+ycc**2+zcc**2)
         if(dc.lt.dcmin) dcmin=dc
333      continue	   
        end do	
        dc=dcmin
        if(dcmin.le.dcsup.and.dcmin.ge.dcinf) calculate=1
       end if 
      if(opd.eq.3) then 
       dcmin=1000000.d0
       jmin=666	   
       do jat=1,natt(jm)
        if(nomat(jat,jm).eq.'Crd') goto 334
         xc1=r(1,icm,im)
         yc1=r(2,icm,im)
         zc1=r(3,icm,im)
         xc2=r(1,jat,jm)
         yc2=r(2,jat,jm)
         zc2=r(3,jat,jm)
         xcc=(xc2-xc1)
         ycc=(yc2-yc1)
         zcc=(zc2-zc1)
         xcc=xcc-box(1)*nint(xcc/box(1))
         ycc=ycc-box(2)*nint(ycc/box(2))
         zcc=zcc-box(3)*nint(zcc/box(3))
         dc=sqrt(xcc**2+ycc**2+zcc**2)
         if(dc.lt.dcmin) then
         dcmin=dc
         jmin=jat
         end if
334      continue	   
        end do	
        dc=dcmin        
        if(dcmin.le.dcsup.and.dcmin.ge.dcinf) calculate=1
       end if
      if(opd.eq.4) then 
       dcmin=1000000.d0
       jmin=666	   
       do jat=1,natt(jm)
        do iat=1,natt(im)
        if(nomat(jat,jm).eq.'Crd'.or.nomat(iat,im).eq.'Crd') goto 335
         xc1=r(1,iat,im)
         yc1=r(2,iat,im)
         zc1=r(3,iat,im)
         xc2=r(1,jat,jm)
         yc2=r(2,jat,jm)
         zc2=r(3,jat,jm)
         xcc=(xc2-xc1)
         ycc=(yc2-yc1)
         zcc=(zc2-zc1)
         xcc=xcc-box(1)*nint(xcc/box(1))
         ycc=ycc-box(2)*nint(ycc/box(2))
         zcc=zcc-box(3)*nint(zcc/box(3))
         dc=sqrt(xcc**2+ycc**2+zcc**2)
         if(dc.lt.dcmin) then
         dcmin=dc
         jmin=jat
         imin=iat
         end if
335      continue	   
         end do	
        end do
        dc=dcmin     
        if(dcmin.le.dcsup.and.dcmin.ge.dcinf) then
        calculate=1
        scalar(iimin)=dfloat(imin)
        scalar(ijmin)=dfloat(jmin)		
        end if
       end if

      if(mixture.eq.2) calculate=1

11      continue
      if(calculate.eq.1) then

      nadmin_local=nadmin_local+1
      ibiv=1

      if(opdx.eq.1) then
      if(moldx1.eq.1) xc1=r(1,atdx1,im)
      if(moldx1.eq.2) xc1=r(1,atdx1,jm)
      if(moldx2.eq.1) xc2=r(1,atdx2,im)	  
      if(moldx2.eq.2) xc2=r(1,atdx2,jm)
      if(mixture.eq.2)xc1=r(1,atdx1,nm+1)
      xcc=(xc2-xc1)
      scalar(idx)=xcc	  
      end if
      if(opdy.eq.1) then
      if(moldy1.eq.1) yc1=r(2,atdy1,im)
      if(moldy1.eq.2) yc1=r(2,atdy1,jm)
      if(moldy2.eq.1) yc2=r(2,atdy2,im)	  
      if(moldy2.eq.2) yc2=r(2,atdy2,jm)
      if(mixture.eq.2)  yc1=r(2,atdy1,nm+1)
       ycc=(yc2-yc1)
       scalar(idy)=ycc	  
      end if
      if(opdz.eq.1) then
      if(moldz1.eq.1) zc1=r(3,atdz1,im)
      if(moldz1.eq.2) zc1=r(3,atdz1,jm)
      if(moldz2.eq.1) zc2=r(3,atdz2,im)	  
      if(moldz2.eq.2) zc2=r(3,atdz2,jm)
      if(mixture.eq.2) zc1=r(3,atdz1,nm+1)
       zcc=zc2-zc1
       scalar(idz)=zcc	  
      end if

      dist=0
      do ii=1,numd
      if(mold1(ii).eq.1) then
      xc1=r(1,atd1(ii),im)
      yc1=r(2,atd1(ii),im)
      zc1=r(3,atd1(ii),im)
      end if
      if(mold1(ii).eq.2) then
      xc1=r(1,atd1(ii),jm)
      yc1=r(2,atd1(ii),jm)
      zc1=r(3,atd1(ii),jm)
      end if
      if(mold2(ii).eq.1) then
      xc2=r(1,atd2(ii),im)
      yc2=r(2,atd2(ii),im)
      zc2=r(3,atd2(ii),im)
      end if
      if(mold2(ii).eq.2) then
      xc2=r(1,atd2(ii),jm)
      yc2=r(2,atd2(ii),jm)
      zc2=r(3,atd2(ii),jm)
      end if      
      if(mixture.eq.2) then
      if(mold1(ii).eq.1) then
      xc1=r(1,1,nm+1)
      yc1=r(2,1,nm+1)
      zc1=r(3,1,nm+1)
      end if
      if(mold2(ii).eq.1) then
      xc1=r(1,1,nm+1)
      yc1=r(2,1,nm+1)
      zc1=r(3,1,nm+1)
      end if
      end if
      xcc=(xc2-xc1)
      ycc=(yc2-yc1)
      zcc=(zc2-zc1)
      xcc=xcc-box(1)*nint(xcc/box(1))
      ycc=ycc-box(2)*nint(ycc/box(2))
      zcc=zcc-box(3)*nint(zcc/box(3))
      dist(ii)=sqrt(xcc**2+ycc**2+zcc**2)
      end do

      if(energyop.eq.1.or.energyop.eq.2) then
      call Energy(r,box,natt,im,jm,LJ,charge,E,ELJ,Ee,
     &  mixture,molmix,molmixi,molmixj,ffo,chkerr,
     &  iprogram,energyop,imin,jmin,erren)
      if(erren.eq.1) nerren=nerren+1
      scalar(iE)=E
      scalar(ieLJ)=ELJ
      scalar(ieQQ)=Ee
      end if

      if(opHB.eq.1) then
      scalar(iHBw)=0.0
      do iHB=1,ncondHB
      xc1=r(1,iatdHB(iHB),im)
      yc1=r(2,iatdHB(iHB),im)
      zc1=r(3,iatdHB(iHB),im)	  
      xc2=r(1,jatdHB(iHB),jm)
      yc2=r(2,jatdHB(iHB),jm)
      zc2=r(3,jatdHB(iHB),jm)
      xcc=(xc2-xc1)
      ycc=(yc2-yc1)
      zcc=(zc2-zc1)
      xcc=xcc-box(1)*nint(xcc/box(1))
      ycc=ycc-box(2)*nint(ycc/box(2))
      zcc=zcc-box(3)*nint(zcc/box(3))
      dHBc=sqrt(xcc**2+ycc**2+zcc**2)
      if(imolHB(iHB).eq.1) m1=im
      if(imolHB(iHB).eq.2) m1=jm
      if(jmolHB(iHB).eq.1) m2=im
      if(jmolHB(iHB).eq.2) m2=jm
      if(imolHB(iHB).eq.1) m3=im
      if(imolHB(iHB).eq.2) m3=jm
      if(kmolHB(iHB).eq.1) m4=im
      if(kmolHB(iHB).eq.2) m4=jm
      call pesc(r,box,m1,iataHB(iHB),
     &      m2,jataHB(iHB),
     &      m3,iataHB(iHB),
     &      m4,kataHB(iHB),psc)
      if(psc.gt.1.01.or.psc.lt.-1.01) write(6,*)'ERROR EN PESC2'
      angHBc=abs(acos(psc)*180.00/pi)
      if(psc.eq.180.00) psc=179.9999
      if(angHBc.le.angHB(iHB).and.dHBc.le.dHB(iHB)) then
      scalar(iHBw)=0.5*dHBc/dhb(ihB)+0.5*angHBc/angHB(iHB)
      end if
      end do
      end if

      ic=1      
      do ii=1,npesc
      if(molpesc(ic,1).eq.1) m1=im
      if(molpesc(ic,1).eq.2) m1=jm
      if(molpesc(ic,2).eq.1) m2=im
      if(molpesc(ic,2).eq.2) m2=jm
      if(molpesc(ic,3).eq.1) m3=im
      if(molpesc(ic,3).eq.2) m3=jm
      if(molpesc(ic,4).eq.1) m4=im
      if(molpesc(ic,4).eq.2) m4=jm
      call pesc(r,box,m1,atpesc(ic,1),
     &      m2,atpesc(ic,2),
     &      m3,atpesc(ic,3),
     &      m4,atpesc(ic,4),psc)
      if(oppesc(ii).eq.2) then
      if(psc.gt.1.0.or.psc.lt.-1.0) write(6,*)'ERROR EN PESC2'
      psc=acos(psc)*180.00/pi
      if(psc.eq.180.00) psc=179.9999
      end if
10      continue
      tpesc(ii)=tpesc(ii)+1
      ic=ic+1
      angbiv(ibiv)=psc
      ibiv=ibiv+1
      end do

      ic=1
      do ii=1,npev
      if(molpev(ic,1,1).eq.1) m11=im
      if(molpev(ic,1,1).eq.2) m11=jm
      if(molpev(ic,1,2).eq.1) m12=im
      if(molpev(ic,1,2).eq.2) m12=jm
      if(molpev(ic,1,3).eq.1) m13=im
      if(molpev(ic,1,3).eq.2) m13=jm
      if(molpev(ic,1,4).eq.1) m14=im
      if(molpev(ic,1,4).eq.2) m14=jm
      if(molpev(ic,2,1).eq.1) m21=im
      if(molpev(ic,2,1).eq.2) m21=jm
      if(molpev(ic,2,2).eq.1) m22=im
      if(molpev(ic,2,2).eq.2) m22=jm
      if(molpev(ic,2,3).eq.1) m23=im
      if(molpev(ic,2,3).eq.2) m23=jm
      if(molpev(ic,2,4).eq.1) m24=im
      if(molpev(ic,2,4).eq.2) m24=jm
      if(oppev(ii).le.2) then
      call pev(r,box,m11,atpev(ic,1,1),m12,atpev(ic,1,2),
     &      m13,atpev(ic,1,3),m14,atpev(ic,1,4),
     &      m21,atpev(ic,2,1),m22,atpev(ic,2,2),
     &      m23,atpev(ic,2,3),m24,atpev(ic,2,4),pv)
      end if	 
      if(pv.eq.666.) goto 2 
      if(oppev(ii).eq.2) then
      if(pv.gt.1.00000.or.pv.lt.-1.00000) then
      write(6,*)'ERROR EN PV. PV=',pv
      call pev(r,box,m11,atpev(ic,1,1),m12,atpev(ic,1,2),
     &      m13,atpev(ic,1,3),m14,atpev(ic,1,4),
     &      m21,atpev(ic,2,1),m22,atpev(ic,2,2),
     &      m23,atpev(ic,2,3),m24,atpev(ic,2,4),pv)
      end if
      if(pv.gt.1.0.or.pv.lt.-1.0) write(6,*)'ERROR EN PEV4'
      pv=acos(pv)*180.00/pi
      if(pv.eq.180.0) pv=179.9999
      end if
2      continue
      if(oppev(ii).eq.4) then	  
      call dih(r,box,m11,atpev(ic,1,1),m12,atpev(ic,1,2),
     &      m13,atpev(ic,1,3),m14,atpev(ic,1,4),
     &      m21,atpev(ic,2,1),m22,atpev(ic,2,2),
     &      m23,atpev(ic,2,3),m24,atpev(ic,2,4),pv)
      end if
20      continue
       tpev(ii)=tpev(ii)+1
       ic=ic+1
       angbiv(ibiv)=pv
       ibiv=ibiv+1
      end do

      if(euler.gt.0) then
      call geteuler(im,jm,r,ae10,ae11,ae12,ae13,
     &      ae20,ae21,ae22,ae23,eulera,box,chkerr,euler,erreu)
      if(erreu.eq.1) nerreu=nerreu+1
      do ic=1,5 
      if(opeuler(ic).eq.2) then
      eulera(ic)=eulera(ic)*180.d0/pi
      end if
      angbiv(ibiv)=eulera(ic)
      ibiv=ibiv+1
      end do
      end if

      ibiv=ibiv-1

c ============================================================
c  OPTION B: Instead of writing directly to the file, pack
c  results into the local output buffer for later gathering.
c ============================================================
      npairs_local=npairs_local+1
      if (npairs_local .gt. MAXPAIRS) then
       write(6,*)'ERROR Option B: MAXPAIRS exceeded! Increase MAXPAIRS'
       call MPI_Finalize(ierr)
       stop
      end if
      nline=0
      nline=nline+1 ; outbuf_local(nline,npairs_local)=float(im)
      nline=nline+1 ; outbuf_local(nline,npairs_local)=float(jm)
      nline=nline+1 ; outbuf_local(nline,npairs_local)=dc
      do ii=1,numd
       nline=nline+1
       outbuf_local(nline,npairs_local)=dist(ii)
      end do
      do ii=1,numscal
       nline=nline+1
       outbuf_local(nline,npairs_local)=scalar(ii)
      end do
      do ii=1,ibiv
       nline=nline+1
       outbuf_local(nline,npairs_local)=angbiv(ii)
      end do

      end if ! calculate=1
50      continue
      end do  ! jm loop
51      continue
      end do  ! imf loop (Option B parallel split)

c ============================================================
c  GATHER results from all ranks to rank 0 for this config
c ============================================================

c --- Gather number of pairs each rank found ---
      call MPI_Gather(npairs_local, 1, MPI_INTEGER,
     &                npairs_all,  1, MPI_INTEGER,
     &                0, MPI_COMM_WORLD, ierr)

c --- Gather only the used part of the output buffers ---
      sendcount=MAXLINE*npairs_local
      if (myrank .eq. 0) then
       displs(0)=0
       totalrecv=0
       do irank=0,nprocs-1
        recvcounts(irank)=MAXLINE*npairs_all(irank)
        if (irank.gt.0) displs(irank)=displs(irank-1)+recvcounts(irank-1)
        totalrecv=totalrecv+recvcounts(irank)
       end do
       if (totalrecv.gt.0) then
        allocate(outbuf_recv(totalrecv))
       else
        allocate(outbuf_recv(1))
       end if
      end if
      call MPI_Gatherv(outbuf_local, sendcount, MPI_REAL,
     &                 outbuf_recv, recvcounts, displs, MPI_REAL,
     &                 0, MPI_COMM_WORLD, ierr)

c --- Gather nadmin for nmolvec accumulation ---
      call MPI_Reduce(nadmin_local, nadmin_global, 1,
     &                MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

c --- Rank 0 writes all collected pairs for this configuration ---
      if (myrank .eq. 0) then
       do irank=0,nprocs-1
        do ipair=1,npairs_all(irank)
         idx0=displs(irank)+(ipair-1)*MAXLINE
         im=int(outbuf_recv(idx0+1))
         jm=int(outbuf_recv(idx0+2))
         dc=outbuf_recv(idx0+3)
         do ii=1,numd
          dist(ii)=outbuf_recv(idx0+3+ii)
         end do
         do ii=1,numscal
          scalar(ii)=outbuf_recv(idx0+3+numd+ii)
         end do
         do ii=1,ntotang
          angbiv(ii)=outbuf_recv(idx0+3+numd+numscal+ii)
         end do
         write(9,"(2i12,100f12.6)") 
     &       im,jm,dc,dist(1:numd),scalar(1:numscal),angbiv(1:ntotang)
        end do
       end do
       nmolvec=2*float(nadmin_global)/float(nm)+nmolvec
       deallocate(outbuf_recv)
      end if

c --- Synchronise all ranks before next configuration ---
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      end do
c ============================================================
c  END OF OPTION B MAIN LOOP
c ============================================================

      if (myrank .eq. 0) then
      write(9,"(a20,i10)")'*************Conf=  ',-1
      close(unit=9)
      close(unit=77)

      write(6,*)''
      write(6,"(a20,1x,f10.5,a5,1x,f10.5,a2,1x,f10.5)")
     & 'Neighbours between',dcinf,'and',dcsup,'='
     $      ,nmolvec/(float(nconf))/2.0
      if(nerren.gt.1) then
       write(6,*)'***WARNING!!!! ENERGY ERRORS***'
       write(6,*)'Number of warnings: ', nerren
       write(6,*)'Open ErrEn.log to know where the errors are'	   
      end if
      end if

      if(interactive.eq.0) goto 999

60      continue
      if (myrank .ne. 0) goto 999
      write(6,"(a8)",advance='no')'angula> '
      read(5,*) dialog
      if(dialog.eq.'end'.or.dialog.eq.'fin') goto 999
      if(dialog.eq.'h'.or.dialog.eq.'help') then
      write(6,*)'Type helpd or hd for angle definition'
      call help(1)
      goto 60
      end if
      if(dialog.eq.'hd'.or.dialog.eq.'helpd') then
      call help(2)
      goto 60
      end if
      if(dialog.eq.'nf'.or.dialog.eq.'changefile') then
      write(6,*)'new file:'
      read(5,*) nomcontrol
      goto 666
      end if
      if(dialog.eq.'go'.or.dialog.eq.'c'.or.dialog.eq.'cp') goto 666
      if(dialog.eq.'e') then
      write(action,"(a5,a20)")'edit ',nomcontrol
      write(6,*) action
      call system(action)
      goto 60
      end if
      write(6,*)'NARF! I do not understand'
      goto 60

999   continue

      deallocate(outbuf_local)
      deallocate(npairs_all)
      deallocate(recvcounts)
      deallocate(displs)

c ============================================================
c  MPI FINALISATION
c ============================================================
      call MPI_Finalize(ierr)

      if (myrank .eq. 0) then
      write(6,*)''
      write(6,*)'I have finished the calculations!'
      write(6,*)'*********************************'
      end if

      goto 13
      if(bigerror.eq.1) then
      write(6,*)''
      write(6,*)'THERE IS A BIG ERROR IN THE CALCULATION!!'
      write(6,*)'Check your traductor, please!!!!!'
      write(6,*)''
      end if
      write(6,*)''
      write(6,*)''
16    continue
       write(6,*)'PROBLEM reading your conf'
       write(6,*)'Last atom succesfully read:'
       write(6,"(3a10)")'type',' molecule ','atom'
       write(6,"(3i10)") imt,imm,ia
13    continue
      stop
      end



      SUBROUTINE pause(message)
      CHARACTER(LEN=*) :: message
      if(trim(message)=='') then
      write (*,*) "Press ENTER to continue!"
      else
      write (*,*) message
      end if
      read (*,*)
      END SUBROUTINE
