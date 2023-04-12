      Program TestingShift1
c-------------------------------------------------------------
c     Descriptions
c-------------------------------------------------------------
c	Zero column contains the values of thresholds 
	Real W(0:1000,1000),V(0:1000,1:2),Hid1(1000)
	Real Del1(1000),Del0(1000),Wde(1000)
	Real TInp(1000),TOut(2)
	Real RNY(45000),Slu(45000),Snen(45000)
c	Real Din3(100),Sne(4,45000)
	Real Din3(10000),Sne(4,45000)

c
        Integer NG(1:45000,0:1000)
	Integer Ir1(10000),Ir2(10000),Inn(10000),IC(20)
	Integer In2(30000),In1(30000)
	Integer Ir11(45000),Ir(45000),Inr1(45000),Inr2(45000)
c	Integer In3(100),Iraz(100)
	Integer In3(10000),Iraz(10000)
	Integer Ir12(45000),Unen(45000),Lnen(45000)
	Integer Ushift(4),Lshift(4),Mshift(4,100)

c      
        Logical Pd(10000),SignRep
c
	Character*1 Seq1(10000),Seq11(45000),Sig(10000),Form(4,100),
     *Dd3(5000),Seq(45000),SS,Seq12(45000),PrSign,Stemp,Plus
	Character*6 Desc(10000),Dd2(30000)
	Character*5000 PdbTemp

c-------------------------------------------------------------
c     Input argument: number of sequences 
c-------------------------------------------------------------        

c	Character*6 :: numchar        
c	If(command_argument_count().ne.1)then        
c	Write(*,*)'arg1: number of sequences required.'        
c	Stop        
c	End if        
c	Call get_command_argument(1,numchar)        
c	Read(numchar,*)Lpdb                
        
c-------------------------------------------------------------
c     General Parameters of NN
c-------------------------------------------------------------
c     Speed of training
      Ceta=0.75      
c     Error estimation 
      Cerror=1e-3     
c     Maximal number of iterations
      Icounter=1000      
c     Maximal initial value of weights
      CInitWeight=0.2
c     Activity of neuron 
      BiasNeuron=1.0
c     Coefficient of inertia
      XIner=1.0
c     Number of neurons in the hidden layer
      Nhid=3
c     Number of entries           
      NInp=39                                            
c     Number of outputs
      NOut=1
c-------------------------------------------------------------
c     Reading weights 
c------------------------------------------------------------- 
c       command line version

c       Open(51,File='Weights1',Status='Old')
c	Open(52,File='Weights2',Status='Old')   

c       online version  (db server)

        Open(51,File='Weights1',Status='Old')
	Open(52,File='Weights2',Status='Old')   

	Read(51,*)W
	Read(52,*)V
c-------------------------------------------------------------
c       Applications for single sequence
c-------------------------------------------------------------	
	Por2=0.101
        Do 501 ijk = 1,45000
  501 	Snen(ijk)=0.0

c       command line version

c        Open(54,File='STDIN',Status='Old')
c	Open(55,file='TestingResult.dat',status='unknown')
c	Open(58,file='STDOUT',status='unknown')
c        Open(59,file='Tempor.dat',status='unknown')

c       online version (db server)

        Open(54,File='STDIN',Status='Old')
	Open(55,file='TestingResult.dat',status='unknown')
	Open(58,file='STDOUT',status='unknown')
        Open(59,file='Tempor.dat',status='unknown')


c       Computing number of input sequences        
        Lpdb = 0 
	Do
	Read(54,*, End=100)        
        Lpdb = Lpdb + 1 
	End do
 100    Rewind 54        

c     reading of titles of sequences
        Do 71 li=1,Lpdb 
	Read(54,'(A6)')Dd2(li)
	In2(li)=Index(Dd2(li),'>')
   71 Continue
c
	NRPdb=0
	NVPdb=0
c-------------------------------------------------------------------	
c       beginning of the seraching cycle 
c       li - counter for sequences	
c-------------------------------------------------------------------
	Do 72 li=1,Lpdb   
	ii1=In2(li)
	ii2=1
   	Rewind 54
	If(li.ge.1.5)then
	do 74 lii=1,li-1 
   74   Read(54,*)
	End if 
	Read(54,372)PdbTemp
c sequence starts 8 characters after the ">" symbol ^M        
        PdbTemp=PdbTemp(ii1+7:) 
   82   Continue 
	Stemp=PdbTemp(ii2:ii2)
        If(Stemp.eq.' ')go to 73
        Dd3(ii2)=Stemp
c max protein sequence size: 5000. change value to increase maximal limit:         
  372   Format(A5000)
	ii2=ii2+1
	Go to 82
   73 Continue
	In1(li)=ii2-NInp-1
	Ile=ii2-1

c--------------------------------------------------------------------	
c       A NEW PART DEVELOPED IN JANUARY 2009 
c       calculating of ANN outputs  
c       li=1 number of analyzing sequences  
c       Ile  length of the analyzing sequence 
c       Ninp=39 number of residues in a window
c       In1(li)=Ile-NInp number of windows on a sequence
c       PARAMETERS FOR SHIFTING: 
c       Nshift=4 total number of windows (shifted and non-shifted)
c       Ushift shift for upper window (residues 21-39)
c       Lshift shift for lower window (residues 1-19 )
c       Kshift=20 is always included
c--------------------------------------------------------------------
c23456789
c       parameters for shifting windows
        Nshift=4
        Kshift=20
	DATA Ushift/0,1,0,1/
	DATA Lshift/0,0,-1,-1/
	Do 775 ish=1,Nshift
	Do 776 j=1,Kshift-1
  776	Mshift(ish,j)=j+Lshift(ish)
        Mshift(ish,Kshift)=Kshift
        Do 777 j=Kshift+1,Ninp
  777   Mshift(ish,j)=j+Ushift(ish)
  775   Continue	
c-------------
c        Do 780 j=1,Ninp
c        Write(*,'(4i8)')Mshift(1,j),Mshift(2,j),Mshift(3,j),Mshift(4,j)
c  780   Continue		
c-------------
c       input data #ii2 
	Do 75 ii2=1,In1(li)
	Do 778 ish=1,Nshift
	Do 76 j=1,Ninp
	jsh=Mshift(ish,j)
	jn=(ii2-1)+jsh
	Form(ish,j)=DD3(jn)
	Call Codir1(DD3(jn),IC,Kn)
	Do 77 jj=1,20
   77	TInp(20*(j-1)+jj)=1.0*IC(jj)
   76	Continue     
c       output value Sne
   	Do 78 j1=1,NHid
  	WTemp=BiasNeuron*W(0,j1)
	Do 79 j=1,20*NInp
   79	WTemp=WTemp+TInp(j)*W(j,j1)
	Hid1(j1)=Sigmoid(1.0*WTemp)
   78	Continue
	Do 80 k=1,NOut
	VTemp=BiasNeuron*V(0,k)
	Do 81 j1=1,NHid
   81   VTemp=VTemp+Hid1(j1)*V(j1,k)
	TOut(k)=Sigmoid(1.0*VTemp)
   80	Vnen=TOut(k)
	Sne(ish,ii2)=Vnen
  778   Continue	
c       maximum around shifted windows  
        Snen(ii2)=0
	Unen(ii2)=0
	Lnen(ii2)=0
	Do 779 ish=1,Nshift
	If(Sne(ish,ii2).gt.Snen(ii2))then
	Snen(ii2)=Sne(ish,ii2)
	Unen(ii2)=Ushift(ish)
	Lnen(ii2)=Lshift(ish)
	End if
  779   Continue	
   75   Continue
c       last residues   
c        Write(*,'(5i8)')Mshift(1,Ninp),Mshift(2,Ninp),Mshift(3,Ninp),
c     *	Mshift(4,Ninp),Mshift(5,Ninp)
c        Write(*,'(5A8)')Form(1,Ninp),Form(2,Ninp),Form(3,Ninp),
c     * 	Form(4,Ninp),Form(5,Ninp)      
c--------------------------------------------------------------------	
c       Snen(ii2) outcome for each input window ii2
c       THE END OF A NEW PART DEVELOPED IN JANUARY 2009 
c--------------------------------------------------------------------

	DNSR=0.0
	Dob=0.0
	Nsr=0

c        Do 502 ijk=1,100
        Do 502 ijk=1,10000
 
	In3(ijk)=0
	Iraz(ijk)=0
	Din3(ijk)=0.0
  502   Continue 
c

        Do 83 ii3=1,In1(li)
	Sn00=Snen(ii3)
	If(Sn00.ge.Por2)then
	Nsr=Nsr+1
	Dob=Sn00
	Dnsr=Dnsr+Dob
        In3(Nsr)=ii3+int(0.5*(NInp-1))-1
	Din3(Nsr)=Dob
	end if
   83	Continue

c       identifying repeats from non-zero ANN outputs
c       Nposit - minimal number of values exceeding threshold 
c       Distr - minimal distance between those values
c       PrSign=' ' -> no repeats; PrSign='R' -> repeat

        Nposit=1
        Distpr=13.5
        SignRep=.false.
	PrSign=' '
        if(Nsr.ge.Nposit-0.5)then
	Nvpdb=Nvpdb+1
        Do 94 ii4=2,Nsr
	Iraz(ii4)=In3(ii4)-In3(ii4-1)
	Dsign=abs(1.0*(Iraz(ii4)-Iraz(ii4-1)))
	If(Dsign.le.Distpr.and.Iraz(ii4-1).ge.15)then
	SignRep=.true.
	PrSign='R'
	Nrpdb=Nrpdb+1
	Idistan=Iraz(ii4-1)
	Istart=In3(ii4-2)
	Goto 95
	end if
   94 Continue
      end if
   95 Continue


c--------------------------------------------------------------------------   
c       printing the information about the current sequence li in the screen
c--------------------------------------------------------------------------
c       and to the file 'Tempor.dat'
c       LI - counter for the sequence
c       ILE - the length of the sequence 
c       IN1(LI) - the length of the sequence li MINUS the size of the window (39+1)
c       NSR  - number of values exceeding threshold
c       DNSR - sum of values exceeding threshold
c       NRPDB - number of identified repeats
c       NVPDB - 
c--------------------------------------------------------------------------
	Write(*,173)li,Dd2(li),Ile,Nsr,Dnsr,Nrpdb,Nvpdb
	Write(59,173)li,Dd2(li),Ile,Nsr,Dnsr,Nrpdb,Nvpdb         
  173	Format(1X,I7,1X,A6,2I7,F9.3,2I7)
c       filling in files Pdb_Res1 and Pdb_Res2
c       if some repeats are identified (NSR > NPOSIT)
c       IPRIN=18 - maximal number of displaying repeats

c	If(Nsr.ge.Nposit-0.5)then	

        Iprin=Nsr
	If(Nsr.ge.20)Iprin=20
c
c       printing of summary in Pdb_Res1	
c
c	Write(55,190)Dd2(li),Ile,PrSign,
c     *  (In3(i),'-',Din3(i),i=1,Iprin)
c  190   Format(A7,3X,I5,3X,A1,1X,10(I5,A1,F3.2,1X))
c  
c       printing of details in Pdb_Res1	
c
	Write(55,183)li,Dd2(li),Ile,Nsr,Dnsr
	Write(55,283)
	Write(55,284)(In3(i),i=1,Iprin)
	Write(55,285)(Din3(i),i=1,Iprin)  
	If(SignRep)then
	Write(55,*)
	Write(55,186)Idistan,Istart
	End if
	Write(55,189)  
  183   Format(1x,' # Sequence= ',i4,' ( ',A8, ' ) ',5x,
     *' Length= ',I4,5x,' Number of NN responds= ',I4,5x,
     *' Total sum of NN responds= ',F4.1)
  184   Format(1x,' Positions ',18I6)
  185   Format(1x,' Responds  ',18F6.1)
  283   Format(1x,' Positions / NN Responds : ')
  284   Format(20I6)
  285   Format(20F6.1)
  186   Format(1x,' THE REPEAT WITH DISTANCE ',I6,
     *' AND INITIAL POSITION ',I6,' HAS BEEN IDENTIFIED')
c
c	Filling in Pdb_Res2
c
	Write(58,'(a1,2a6)')'>',Dd2(li),PrSign
	iw1=int(0.5*(NInp-1))
	iw2=iw1+in1(li)+1
	Do 89 ki=1,Ile
	scc=0.0
	icc=0
	lcc=0
	if(ki.ge.iw1.and.ki.le.iw2)then
	scc=Snen(ki-int(0.5*(NInp-1))+1)
	icc=Unen(ki-int(0.5*(NInp-1))+1)
	lcc=Lnen(ki-int(0.5*(NInp-1))+1)
	End if
	Write(58,'(A1,F4.2,I4,I4)')dd3(ki),scc,lcc,icc
   89	Continue
c	     

c        End if  

   72	Continue
  199	Continue 
c-------------------------------------------------------------
c      printing the final statement on the screen 
c      after chacking all the sequences
c-------------------------------------------------------------   
       Write(55,189)       
       Write(55,187)Nposit,Nvpdb
       Write(55,188)1,Nrpdb
       Write(55,189)
  187  Format(1x,'Number of sequences with more than ',I6,
     *' non-zero NN replys equals to ',I6)
  188  Format(1x,'Number of sequences with more than ',I6,
     *' repeats equals to            ',I6)
  189 Format(120('-'))
       Stop
       End 
c-------------------------------------------------------------
c     Function of Neuron
c-------------------------------------------------------------            
	Real Function Sigmoid(z)
	Sigmoid=1.0/(1.0+Exp(-1.0*z))
	Sigmoid=0.8*Sigmoid+0.1
	End
c
	Real Function Binar(z)
	If(z.gt.0.)Binar=1.0
	If(z.le.0.)Binar=0.0
	End
c
	Subroutine Codir1(SS,IV,KN)
	Character*1 SS
	Integer IV(20)
	Do 1 l=1,20
    1	IV(l)=0
        If(SS.eq.'A'.OR.SS.eq.'a')KN=1
        If(SS.eq.'C'.OR.SS.eq.'c')KN=2
        If(SS.eq.'D'.OR.SS.eq.'d')KN=3
        If(SS.eq.'E'.OR.SS.eq.'e')KN=4
        If(SS.eq.'F'.OR.SS.eq.'f')KN=5
        If(SS.eq.'G'.OR.SS.eq.'g')KN=6
        If(SS.eq.'H'.OR.SS.eq.'h')KN=7
        If(SS.eq.'I'.OR.SS.eq.'i')KN=8
        If(SS.eq.'K'.OR.SS.eq.'k')KN=9
        If(SS.eq.'L'.OR.SS.eq.'l')KN=10
        If(SS.eq.'M'.OR.SS.eq.'m')KN=11
        If(SS.eq.'N'.OR.SS.eq.'n')KN=12
        If(SS.eq.'P'.OR.SS.eq.'p')KN=13
        If(SS.eq.'Q'.OR.SS.eq.'q')KN=14
        If(SS.eq.'R'.OR.SS.eq.'r')KN=15
        If(SS.eq.'S'.OR.SS.eq.'s')KN=16
        If(SS.eq.'T'.OR.SS.eq.'t')KN=17
        If(SS.eq.'V'.OR.SS.eq.'v')KN=18
        If(SS.eq.'W'.OR.SS.eq.'w')KN=19
        If(SS.eq.'Y'.OR.SS.eq.'y')KN=20
	IV(KN)=1
	End
c
	Subroutine Codir2(SS,IV,KN)
	Character*1 SS
	Integer IV(20)
	Do 1 l=1,20
    1	IV(l)=0
	Select case(KN)
	Case(1)
	SS='A'
	Case(2)
	SS='C'
	Case(3)
	SS='D'
	Case(4)
 	SS='E'
	Case(5)
 	SS='F'
	Case(6)
 	SS='G'
	Case(7)
 	SS='H'
	Case(8)
 	SS='I'
	Case(9)
 	SS='K'
	Case(10)
  	SS='L'
 	Case(11)
   	SS='M'
	Case(12)
  	SS='N'
	Case(13)
  	SS='P'
	Case(14)
  	SS='Q'
	Case(15)
  	SS='R'
	Case(16)
  	SS='S'
	Case(17)
  	SS='T'
	Case(18)
  	SS='V'
	Case(19)
  	SS='W'
	Case(20)
  	SS='Y'
	End Select
	IV(KN)=1
	End
c
