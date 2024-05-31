!---addFJX_H2CO.f  -- with user supplied subroutine that supplies X-section x Q-yield
!---generates fast-JX 18-bin X-sections (mprather modified file,
!---Sebastianelli 02/2023)

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
      real*8 XCH2O, QCH2O    
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
        TTT(1) = 223.
        TTT(2) = 298.
        TITLTBL(1) = 'H2COa '
        TITLTBL(2) = 'H2CO1 '

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

           call X_CH2O (WW, XT, XCH2O) 
           call Q_CH2O (WW, XT, QCH2O) 

           FBIN(I) = FBIN(I) + F(J)
           ABIN(I) = ABIN(I) + F(J)*XCH2O *QCH2O
           

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

      subroutine X_CH2O (WW, TT, XCH2O)

c---JPL_2010 standard tables.
c---interpolates H2CO cross section vs. wavelength and temperature.
c---    H2CO is linear, limited to range 220 to 294K (THIS MUST BE
c---    CHECKED)

      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: XCH2O !

      integer IW,Z
      real*8 X1,X2,X1p,X2p,FT,FTX,FW

C----JPL 2010 Table 4D-3. 
      real*4, parameter, dimension(1:88)  :: WNM =
     &[251, 252, 253, 254, 255,
     & 256, 257, 258, 259, 260, 261,
     & 262, 263, 264, 265, 266, 267,
     & 268, 269, 270, 271, 272, 273,
     & 274, 275, 276, 277, 278, 279,
     & 280, 281, 282, 283, 284, 285,
     & 286, 287, 288, 289, 290, 291,
     & 292, 293, 294, 295, 296, 297,
     & 298, 299, 300, 301, 302, 303,
     & 304, 305, 306, 307, 308, 309,
     & 310, 311, 312, 313, 314, 315,
     & 316, 317, 318, 319, 320, 321,
     & 322, 323, 324, 325, 326, 327,
     & 328, 329, 330, 331, 332, 333,
     & 334, 335, 336, 337, 338]



      real*4, parameter, dimension(1:88)  :: X220K =
     &[0.2027, 0.3375, 0.2882, 0.3410, 0.4480,
     & 0.6277, 0.4364, 0.3019, 0.6174, 0.6006, 0.6593,
     & 0.5972, 1.0861, 0.9427, 0.5225, 0.5348, 1.3666,
     & 1.2409, 0.9816, 0.9544, 1.9498, 1.4225, 0.8008,
     & 0.6505, 2.1705, 2.6053, 1.5555, 1.0193, 2.4691,
     & 2.3451, 1.5558, 0.9659, 0.7200, 4.3354, 4.0635,
     & 2.0781, 1.1443, 3.2053, 3.2291, 1.1572, 1.8487,
     & 0.7873, 3.1468, 7.2266, 4.0313, 2.4641, 1.3496,
     & 4.2570, 3.1589, 0.9310, 1.6493, 0.8694, 3.0470,
     & 7.2662, 4.7087, 4.2876, 1.7561, 1.3775, 3.2890,
     & 1.7223, 0.4604, 1.1992, 0.9128, 5.6446, 5.5406,
     & 2.5043, 5.8085, 3.1252, 0.9560, 1.1882, 1.6164,
     & 0.7199, 0.3241, 0.8644, 1.5234, 6.9123, 4.3292,
     & 1.1820, 3.1451, 3.8763, 1.3807, 0.3326, 0.2108,
     & 0.1607, 0.0966, 0.1236, 0.3805, 1.9080]

      real*4, parameter, dimension(1:88)  :: X294K =
     &[0.2040, 0.3370, 0.2890, 0.3420, 0.4500,
     & 0.6290, 0.4430, 0.3070, 0.6180, 0.6040, 0.6600,
     & 0.6020, 1.0800, 0.9470, 0.5300, 0.5380, 1.3600,
     & 1.2400, 0.9900, 0.9600, 1.9400, 1.4300, 0.8100,
     & 0.6570, 2.1500, 2.5900, 1.5700, 1.0300, 2.4500,
     & 2.3400, 1.5600, 0.9720, 0.7200, 4.2700, 4.0500,
     & 2.0900, 1.1500, 3.1700, 3.2200, 1.1700, 1.8400,
     & 0.7960, 3.1100, 7.1500, 4.0600, 2.4800, 1.3600,
     & 4.2200, 3.1700, 0.9630, 1.6300, 0.8520, 3.0200,
     & 7.2300, 4.7400, 4.2900, 1.7800, 1.3800, 3.2600,
     & 1.7400, 0.4610, 1.1900, 0.9020, 5.6500, 5.5600,
     & 2.5400, 5.7900, 3.1500, 0.9750, 1.1900, 1.6000,
     & 0.7210, 0.3270, 0.8610, 1.5400, 6.8700, 4.3700,
     & 1.2200, 3.1200, 3.8600, 1.4100, 0.3460, 0.2140,
     & 0.1590, 0.0966, 0.1260, 0.3830, 1.9200]

      XCH2O = 0.d0
      
      if (WW .lt. WNM(1) .or. WW .gt. WNM(88)) goto 2
      
      do Z=1,87
      if (WW .ge. WNM(Z)) IW=Z
      enddo

      FW = WW - WNM(IW)
      FW = max(0.d0, min(1.d0, FW))

      X1p = X220K(IW)
      X2p = X294K(IW)

      X1 = (X1p + FW*(X220K(IW+1)-X220K(IW)))*1.0d-20
      X2 = (X2p + FW*(X294K(IW+1)-X294K(IW)))*1.0d-20
      FTX = (TT - 223.d0) / (298.d0 - 223.d0)
      FT = max (0.d0, min (1.d0, FTX))
      XCH2O = X1 + FT*(X2-X1)

    2 continue
      return
      end


      subroutine Q_CH2O (WW, TT, QCH2O)

c---JPL_2010 standard tables.
c---interpolates H2CO quantum yield vs. wavelength 

c---JPL 2006 recommendations
      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: QCH2O 

      integer IW
      real*8 FW,Q1,FT

C----------JPL 2010 Table 4D-4 pg297.
      real*4, parameter, dimension(88)  :: WNQ =
     &[251.0, 252.0, 253.0, 254.0, 255.0, 256.0,
     & 257.0, 258.0, 259.0, 260.0, 261.0, 262.0, 263.0,
     & 264.0, 265.0, 266.0, 267.0, 268.0, 269.0, 270.0, 271.0,
     & 272.0, 273.0, 274.0, 275.0, 276.0, 277.0, 278.0, 279.0,
     & 280.0, 281.0, 282.0, 283.0, 284.0, 285.0, 286.0, 287.0,
     & 288.0, 289.0, 290.0, 291.0, 292.0, 293.0, 294.0, 295.0,
     & 296.0, 297.0, 298.0, 299.0, 300.0, 301.0, 302.0, 303.0,
     & 304.0, 305.0, 306.0, 307.0, 308.0, 309.0, 310.0, 311.0,
     & 312.0, 313.0, 314.0, 315.0, 316.0, 317.0, 318.0, 319.0, 
     & 320.0, 321.0, 322.0, 323.0, 324.0, 325.0, 326.0, 327.0,
     & 328.0, 329.0, 330.0, 331.0, 332.0, 333.0, 334.0, 335.0,
     & 336.0, 337.0, 338.0]
 
      real*4, parameter, dimension(88)  :: Q298 =
     &[0.308, 0.307, 0.306, 0.305, 0.304, 0.304,
     & 0.303, 0.303, 0.304, 0.307, 0.312, 0.318, 0.325,
     & 0.333, 0.343, 0.354, 0.365, 0.377, 0.390, 0.404, 0.418,
     & 0.433, 0.448, 0.464, 0.479, 0.495, 0.512, 0.528, 0.544,
     & 0.560, 0.576, 0.591, 0.606, 0.620, 0.633, 0.645, 0.657,
     & 0.669, 0.680, 0.690, 0.700, 0.710, 0.718, 0.726, 0.734,
     & 0.740, 0.746, 0.751, 0.755, 0.758, 0.761, 0.762, 0.762,
     & 0.762, 0.760, 0.758, 0.754, 0.749, 0.744, 0.737, 0.729,
     & 0.720, 0.709, 0.698, 0.685, 0.671, 0.656, 0.639, 0.622,
     & 0.603, 0.583, 0.561, 0.539, 0.515, 0.489, 0.463, 0.435,
     & 0.406, 0.375, 0.343, 0.310, 0.276, 0.240, 0.203, 0.165,
     & 0.126, 0.085, 0.043]

      if (WW .lt. WNQ(1) .or. WW .gt. WNQ(88)) goto 2

      IW = (WW - 250.d0)
      IW = max(1, min(87, IW))
      FW = WW - WNQ(IW)
      FW = max(0.d0, min(1.d0, FW))
      Q1 = Q298(IW) + FW*(Q298(IW+1)-Q298(IW))
      QCH2O = Q1 

   2  continue   
      return
      end
