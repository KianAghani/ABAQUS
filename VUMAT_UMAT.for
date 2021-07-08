      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
	  
	  REAL*8 E,xnu
	  
       dimension stress(ndir+nshr), statev(nstatev),
     & ddsdde(ndir+nshr,ndir+nshr),stran(ndir+nshr), dstran(ndir+nshr),TIME(2), predef(1)

      
	    E=PROPS(1)
	    xnu=PROPS(2)

        ntens  = ndir + nshr
        TIME(1)=totalTime
        TIME(2)=stepTime
      
      do 20 J=1,6
        stress(J)=0.0
        dstran(J)=0.0
20    continue
	  
	  do km = 1,nblock

      do i = 1, ndir
		stress(i) = stressOld(km,i)              
		dstran(i) = strainInc(km,i)
      end do 
      stress(4) = stressOld(km,4)
      dstran(4) = 2.*strainInc(km,4)   
      if (nshr .gt. 1) then      
		stress(5) = stressOld(km,6)
		dstran(5) = 2.*strainInc(km,6)
		stress(6) = stressOld(km,5)
		dstran(6) = 2.*strainInc(km,5)
      end if
       
******************************************************
** Sending vumat                                     *
******************************************************
        
        call xmat(stress,dstran,E,xnu,TIME,ndir,nshr)
     
******************************************************
** Getting vumat                                     *
******************************************************
         
	  do i = 1, ndir
		stressnew(km,i) = stress(i)               
      end do 
      stressnew(km,4) = stress(4)    
      if (nshr .gt. 1) then      
		stressnew(km,5) = stress(6) 
		stressnew(km,6) = stress(5) 
	  end if
      end do
       

      return
      end
	  
	  
******************************************************
** The sub inter_face                                *
******************************************************	  
	  
        subroutine xmat(stress,dstran,E,xnu,TIME,ndir,nshr)
	 
        REAL*8 E,xnu,A,B,C,G
		INTEGER ndir,nshr
        DIMENSION STRESS(ndir + nshr),DDSDDE(ndir + nshr,ndir + nshr),DSTRAN(ndir + nshr)

			G=E/2./(1.+xnu)
	   	    A=(E*(xnu-1.))/(2.*(xnu**2.)+xnu-1.)
		    B=(-E*xnu)/(2.*(xnu**2.)+xnu-1.)
		    C=G
			
	      DO K1=1, 6
          DO K2=1, 6
			DDSDDE(K2, K1)=0.
          END DO
          END DO

          DO I=1,ndir
          DO J=1,ndir
			DDSDDE(J,I)=B
		  END DO
			DDSDDE(I,I)=A
		  END DO
		  
		  DO L=ndir+1,ndir + nshr
			DDSDDE(L,L)=C
		  END DO

		  DO I=1,ndir + nshr
		  DO J=1,ndir + nshr
		       STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
		  END DO
		  END DO
           
	 
	 return
	 end