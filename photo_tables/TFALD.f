!---addFJX_H2CO.f  -- with user supplied subroutine that supplies X-section x Q-yield
!---generates fast-JX 18-bin X-sections (mprather modified file,
!---Killen 07/2024)

      implicit none

      integer, parameter :: NB_ = 100       
      integer, parameter :: NS_ = 40000    
      integer, parameter :: NJ_ = 18        
      real*8 SRB(15,NJ_)                
      real*8, dimension(NB_+1) :: WBIN   
      real*8, dimension(NB_) :: FBIN, ABIN, XBIN
      real*8, dimension(NJ_) :: FFBIN,AABIN,XXBIN
      real*8 W(NS_),F(NS_)                
      real*8 W1,W2,WW, XT
      real*8 XCF3CHO, QCF3CHO    
      integer IBINJ(NS_)
      integer IJX(NB_), ITT               
      integer NB,I,J,J1,J2,K,NB77,I778
      integer     TTT(2)
      character*1 ISXP
      character*6 TITLNEW
      character*6 TITLTBL(2)
      Character*20 TITL77

      open (1, file='wavel-bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(i5)') NB
        if (NB .gt. NB_) stop
        read(1,'(5x,f8.3)') (WBIN(I), I=1,NB+1) 	
        read(1,*)  					
        read(1,*)  					
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8) 
        read(1,*)					
        read(1,'(5x,i5)') (IJX(I),I=16,NB)		 
      close (1)

      open (2, file='solar-p05nm-UCI.dat', status='OLD')
        read(2,*)  					  
        read(2,*)					   
        read(2,'(f10.4,e10.3)') (W(J),F(J), J=1,NS_)        
      close (2)

        IBINJ(:) = 0
       do I=1,NB
         W1 = WBIN(I)       
         W2 = WBIN(I+1)	    
        do J=1,NS_
          if (W(J) .gt. W1) goto 11 
        enddo
          J = NS_ + 1
   11     J1 = J
        do J=J1,NS_
          if (W(J) .gt. W2) goto 12
        enddo
          J = NS_ + 1
   12     J2 = J-1
        do J=J1,J2
          IBINJ(J) = I  
        enddo			
       enddo                                                     

!!!!this binning does not interpolate and is OK for large bins
!     it has 7% error in the very short wavel S-R bins of pratmo.

!not extrapolated, uses Q-yld at 300K min
        TTT(1) = 298.
        TTT(2) = 298.
        TITLTBL(1) = 'CF3CHO '
        TITLTBL(2) = 'CF3CHOb '

! We have only two temperatures
        do K = 1,2
        XT = TTT(K)

!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0
         I = 0

! primary high-resolution wavelength loop - generate input for pratmo reference J's

        do J = 1,NS_
         I = IBINJ(J)
         if (I .gt. 0) then
          WW = W(J)

           call X_CF3CHO (WW, XT, XCF3CHO) 
           call Q_CF3CHO (WW, XT, QCF3CHO) 

           FBIN(I) = FBIN(I) + F(J)
           ABIN(I) = ABIN(I) + F(J)*XCF3CHO *QCF3CHO
           

         endif
        enddo

        do I=1,NB
         if (FBIN(I) .gt. 0.d0) then
           NB77 = I
           ABIN(I) = ABIN(I)/FBIN(I)
         endif
        enddo

!---get Flsn
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
         FFBIN(:) = 0.d0
         AABIN(:) = 0.d0
         !XXBIN(:) = 0.d0
        do I=16,NB
         J = IJX(I)
         FFBIN(J) = FFBIN(J) + FBIN(I)
         AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)
         enddo

        do I=1,15
         do J=1,NJ_
           FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
           AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)*SRB(I,J)
         enddo
        enddo
        do J=1,NJ_
         if (FFBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/FFBIN(J)
        enddo

! save UCI fast-JX v68 data bins for 'FJX_spec.dat'

        TITLNEW = TITLTBL(1)
        if (K .le. 2) then

         write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITLTBL(1),' ',TTT(K), AABIN

         endif

        enddo
     
      stop
      end

      subroutine X_CF3CHO (WW, TT, XCF3CHO)

c---JPL_2010 standard tables.
c---interpolates H2CO cross section vs. wavelength and temperature.
c---    H2CO is linear, limited to range 220 to 294K (THIS MUST BE
c---    CHECKED)

      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: XCF3CHO !

      integer IW,Z
      real*8 X1,X2,X1p,X2p,FT,FTX,FW

C----JPL 2010 Table 4D-3. 
      real*4, parameter, dimension(5)  :: WNM =
     &[248,
     & 266, 281,
     & 308, 320]


      real*4, parameter, dimension(5)  :: X298K =
     &[0.311,
     & 1.15, 2.28,
     & 2.97, 2.25]

      XCF3CHO = 0.d0
      
      if (WW .lt. WNM(1) .or. WW .gt. WNM(5)) goto 2
      
      do Z=1,5
      if (WW .ge. WNM(Z)) IW=Z
      enddo

      FW = WW - WNM(IW)
      FW = max(0.d0, min(1.d0, FW))

      X2p = X298K(IW)

     
      X2 = (X2p + FW*(X298K(IW+1)-X298K(IW)))*1.0d-20
C-----FTX = (TT - 223.d0) / (298.d0 - 223.d0)
C-----FT = max (0.d0, min (1.d0, FTX))
      XCF3CHO = X2

    2 continue
      return
      end


      subroutine Q_CF3CHO (WW, TT, QCF3CHO)

c---JPL_2010 standard tables.
c---interpolates H2CO quantum yield vs. wavelength 

c---JPL 2006 recommendations
      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: QCF3CHO

      integer IW
      real*8 FW,Q1,FT

C----------JPL 2010 Table 4D-4 pg297.
      real*4, parameter, dimension(5)  :: WNQ =
     &[248,
     & 266, 281,
     & 308, 320]
 
      real*4, parameter, dimension(5)  :: Q298 =
     &[0.358,
     & 0.407, 0.0564,
     & 0.0003, 0.0000]

      if (WW .lt. WNQ(1) .or. WW .gt. WNQ(5)) goto 2

      IW = (WW - 250.d0)
      IW = max(1, min(3, IW))
      FW = WW - WNQ(IW)
      FW = max(0.d0, min(1.d0, FW))
      Q1 = Q298(IW) + FW*(Q298(IW+1)-Q298(IW))
      QCF3CHO = Q1 

   2  continue   
      return
      end
