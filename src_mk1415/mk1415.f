c pdeng: this is a simplified version of nt14S.f, which only reads NDXSRF, ZNATDN, and LABELS to prepare type 14, 15 cards of A.NIP3
      PROGRAM MK1415
      IMPLICIT REAL*8(A-H,O-Z)
      integer :: BLK_SIZE
      allocatable :: BLK(:)      
      character*1 :: freeform

      COMMON/UNITS/ NL,NB,NZ                                            
      COMMON /CCCC/ IGOM,NZONE,MREG,NZCL,NCINTI,NCINTJ,NCINTK,          
     1              NINTI,NINTJ,NINTK,IMB1,IMB2,JMB1,JMB2,KMB1,         
     2              KMB2,NBS,NBCS,NIBCS,NZWBB,NTRIAG,NRASS,             
     3              NTZSZ,NAREA,LREGA,NARC1,NARC2,NARC3,NARC4,          
     4              NARC5,NARC6,NARC7,NARC8,NON,NSN,NNS,NAN,            
     5              NZ0NE,NSZ,NBLKAD                                    
CC
      if(iargc() > 0) then
        call getarg(1,freeform)
      else
        freeform = 'n'
      end if
      
      OPEN (UNIT=2,FILE='LABELS',FORM='UNFORMATTED')
      OPEN (UNIT=3,FILE='NDXSRF',FORM='UNFORMATTED')
      OPEN (UNIT=4,FILE='ZNATDN',FORM='UNFORMATTED')
CC
      NL=2                                                              
      NB=3                                                              
      NZ=4                                                              
CC                                                                      
      READ(NB)                                                          
      READ(NB) NON,NSN,NNS,NAN,NZONE,NSZ                                
      REWIND NB                                                         
CC                                                                      
      READ(NZ)                                                          
      READ(NZ) IDUM,IDUM,NTZSZ,NNS,NBLKAD                               
      REWIND NZ                                                         
CC                                                                      
      READ(NL,END=4)                                                    
      READ(NL) NTZSZL                                                   
      REWIND NL                                                         
      IF(NTZSZ.EQ.NTZSZL) GO TO 5                                       
      WRITE(6,102) NTZSZ,NTZSZL                                         
   4  NL=0                                                              
   5  CONTINUE                                                          
CC                                                                      
      NCMP=NTZSZ                                                        
      NISO2=NON                                                         
      L1=NISO2*NCMP                                                     
      L2=NNS*((NTZSZ-1)/NBLKAD+1)                                       
      IPTS1=1                                                           
      IPTS2=IPTS1+L1                                                    
      IPTS3=IPTS2+NISO2                                                 
      IPTS4=IPTS3+NNS*NSN
      IPTS5=IPTS4+NCMP
      IPTS6=IPTS5+NCMP
      IPTS7=IPTS6+NSZ
      IPTS8=IPTS7+NCMP
      IPTS9=IPTS8+NSZ
      IPTS10=IPTS9+NSZ
      LAST=IPTS10+L2

      LAST2 = IPTS3 + NISO2*2 + NCMP
      BLK_SIZE = MAX(LAST, LAST2)
      allocate(BLK(BLK_SIZE), stat=IOS)
      if(IOS /= 0) then
            write(6,'("[mk1415]...Failed to allocate working array, &
     & required memory size =",I12)')
            call abort
      end if
                                                       
      CALL RDMIXR(BLK(IPTS1),BLK(IPTS2),BLK(IPTS3),BLK(IPTS4),          
     1            BLK(IPTS5),BLK(IPTS6),BLK(IPTS7),BLK(IPTS8),          
     2            BLK(IPTS9),BLK(IPTS10),NON,NCMP,NNS,NSZ)              
      IPTS4=IPTS3+NISO2                                                 
      IPTS5=IPTS4+NISO2                                                 
                                                   
      CALL TYP14(BLK(IPTS1),BLK(IPTS2),BLK(IPTS3),BLK(IPTS4),           
     1           BLK(IPTS5),NISO2,NCMP, freeform)     

 102  FORMAT('0NUMBER OF COMPOSITIONS ON ZNATDN =',I6/                  
     X       ' NUMBER OF COMPOSITIONS ON LABELS =',I6/                  
     X       ' THE COMPOSITION NAMES ON LABELS WILL BE IGNORED'//)      
      END
      
c ====================== called subroutines
      SUBROUTINE RDMIXR(CONC,SONME,NOS,VOLZ,VFPA,VLSA,NSPA,             
     1                  NSSA,NZSZ,ADEN,NISO2,NCMP,MMS,MSZ)              
      REAL*8 SONME,CONC,DUM8
      COMMON /CCCC/ IGOM,NZONE,MREG,NZCL,NCINTI,NCINTJ,NCINTK,          
     1              NINTI,NINTJ,NINTK,IMB1,IMB2,JMB1,JMB2,KMB1,         
     2              KMB2,NBS,NBCS,NIBCS,NZWBB,NTRIAG,NRASS,             
     3              NTZSZ,NAREA,LREGA,NARC1,NARC2,NARC3,NARC4,          
     4              NARC5,NARC6,NARC7,NARC8,NON,NSN,NNS,NAN,            
     5              NZ0NE,NSZ,NBLKAD                                    
      DIMENSION SONME(NISO2),CONC(NISO2,NCMP),NOS(MMS,1),               
     1          VOLZ(NZONE),VFPA(NZONE),VLSA(MSZ),NSPA(NZONE),          
     2          NSSA(MSZ),NZSZ(MSZ),ADEN(MMS,1)                         
      NDX=3                                                             
      NDZ=4                                                             
      NWDS=2*NON+NAN+4*NSN                                              
      READ (NDX)                                                        
      READ (NDX)                                                        
      READ (NDX) SONME,(DUM8,I=1,NON),(DUM4,I=1,NWDS),
     1           ((NOS(I,L),I=1,NNS),L=1,NSN)
      IF (NSZ.GT.0) READ (NDX) VOLZ,VFPA,VLSA,NSPA,NSSA,NZSZ            
      IF (NSZ.EQ.0) READ (NDX) VOLZ,VFPA,NSPA                           
      REWIND NDX                                                        
      DO 100 I=1,NON                                                    
      DO 100 J=1,NCMP                                                   
      CONC(I,J)=0.0                                                     
  100 CONTINUE                                                          
      READ (NDZ)                                                        
      READ (NDZ) DUM4,DUM4,NTZSZ,NNS,NBLKAD                             
      JU=0                                                              
      IOFF=(NTZSZ-1)/NBLKAD+1                                           
      DO 160 M=1,NBLKAD                                                 
      JL=JU+1                                                           
      JU=JU+IOFF                                                        
      IF (JU.LE.NTZSZ) GO TO 105                                        
      JU=NTZSZ                                                          
      IOFF=JU-JL+1                                                      
  105 READ (NDZ) ((ADEN(N,K),N=1,NNS),K=1,IOFF)                         
      DO 150 N=1,NNS                                                    
      K=0                                                               
      DO 140 J=JL,JU                                                    
      K=K+1                                                             
      IF (ADEN(N,K).EQ.0.0) GO TO 130                                   
      IF (J.GT.NZONE) GO TO 110                                         
      L=NSPA(J)                                                         
      I=NOS(N,L)                                                        
      CONC(I,J)=CONC(I,J)+ADEN(N,K)*VFPA(J)                             
      GO TO 120                                                         
  110 CONTINUE                                                          
      ISUB=J-NZONE                                                      
      IZ=NZSZ(ISUB)                                                     
      L=NSSA(ISUB)                                                      
      I=NOS(N,L)                                                        
      CONC(I,IZ)=CONC(I,IZ)+ADEN(N,K)*VLSA(ISUB)/VOLZ(IZ)               
  120 CONTINUE                                                          
  130 CONTINUE                                                          
  140 CONTINUE                                                          
  150 CONTINUE                                                          
  160 CONTINUE                                                          
      REWIND NDZ                                                        
      RETURN                                                            
      END
      
      SUBROUTINE TYP14(CONC,SONME,ADEN,ENAME,CNAME,NISO2,NCMP,freeform)          
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 FMAT1,FMAT2
      character*8 regnam
      character*1 freeform
      COMMON/UNITS/ NL,NB,NZ                                            
      DIMENSION CONC(NISO2,1),SONME(1),ADEN(1),ENAME(1),FMAT1(6),       
     1          FMAT2(6),CNAME(NCMP)                                    
      EQUIVALENCE (A1,I1)                                               
      DATA FMAT1/8H((2H14, , 8H  10X,6H, 8H        ,                    
     X           8H,3(A6   , 8H,1PE12.5, 8H)))     /                    
      DATA FMAT2/8H((2H14, , 8H  6X, 6H, 8H        ,                    
     X           8H,3(1X,A6, 8H,1PE12.5, 8H)))     /     
CC                                                                      
      DATA A1/8H        /                                               
CC
      IF(NL.EQ.0) GO TO 4                                               
      READ(NL)                                                          
      READ(NL)                                                          
      READ(NL) CNAME                                                    
      REWIND NL                                                         
      GO TO 6                                                           
   4  DO 5 I=1,NCMP                                                     
      I1=I                                                   
      CALL IN2LIT(I1,CNAME(I))
   5  CONTINUE                                                          
   6  CONTINUE                                                          
CC
      OPEN(9,file='ANIP1415',status='replace',form='formatted')
      do I=1,NCMP
        N=0                               
        do J=1,NISO2                                                   
          IF(CONC(J,I).EQ.0.) cycle                                     
          N=N+1                                                             
          ENAME(N)=SONME(J)                                                 
          ADEN(N)=CONC(J,I)                                                 
        end do                                                        
        IF(N.LE.0) cycle
        if(freeform == 'y' .or. freeform == 'Y') then
            CALL CHRFLT (FMAT2(3),CNAME(I),1)
            WRITE(9,FMAT2) (ENAME(K),ADEN(K),K=1,N)
        else
            CALL CHRFLT (FMAT1(3),CNAME(I),1)
            WRITE(9,FMAT1) (ENAME(K),ADEN(K),K=1,N)
        end if
      end do                                                          
      if(NL > 0) then
        if(freeform=='y' .or. freeform=='Y') then 
          do I = 1, NCMP
            write(9,102) CNAME(I),CNAME(I)
          end do
        else
          do I=1,NCMP                                                    
            write(9,101) CNAME(I),CNAME(I)
          end do
        end if
        close(9)
      end if                                                         
C
 100  format('14          ',A6,3(1X,A6,PE12.5))
 101  FORMAT(2H15,4X,2A6)                                               
 102  FORMAT(2H15,4X,A6,1X,A6)                                               
C                                                                       
CC                                                                      
      RETURN                                                            
      END
      
CDECK CHRFLT
      SUBROUTINE CHRFLT(CHR,FLT,LEN)
C
C     CHRFLT REMOVES CHARACTER DATA FROM FLOATING POINT VARIABLE
C     ARRAYS WITHOUT TYPE CONVERSION.
C
C     THE FUNCTIONS PROVIDED BY FLTCHR AND CHRFLT ARE VIOLATIONS OF THE
C     FORTRAN-77 STANDARD, AND, IN FACT, THE PURPOSE OF THE ROUTINES IS
C     DELIBERATELY TO CIRCUMVENT THE STANDARD IN ORDER TO PERMIT
C     PROGRAMS TO MIX CHARACTER AND NUMERICAL DATA IN THE BPOINTER
C     CONTAINER AND IN STANDARD-FILE-FORMAT RECORDS.
C
C                        LATEST VERSION 10/06/95 AT CHOSUN U.
C
C     FLT               INPUT ARRAY
C     CHR               OUTPUT ARRAY
C     LEN               NUMBER OF ELEMENTS
C
      DOUBLE PRECISION FLT, F
C
      CHARACTER*8 CHR, C
      COMMON /CHQQFL/ C
      COMMON /FLQQCH/ F
C
      DIMENSION CHR(*), FLT(*)
C
      IF( LEN.LE.0 ) RETURN
      DO 10 I=1,LEN
C
      F=FLT(I)
      CALL CFCOPY
      CHR(I)=C
   10 CONTINUE
      RETURN
      END
      
CDECK IN2LIT
      SUBROUTINE IN2LIT(INT,WORD)
C
C  IN2LIT CONVERTS AN INTEGER VARIABLE INTO A NCHAR-CHARACTER LITERAL.
C  ON IBM MACHINES THE INTEGER MUST BE INTEGER*4.  ON ALL MACHINES
C  THE VARIABLE CONTAINING THE LITERAL MUST BE REAL (A DOUBLE REAL
C  WORD ON SHORT-WORD MACHINES).
C
C                       LATEST VERSION 10/06/95 AT CHOSUN U.
C
C     INT               THE INPUT INTEGER
C     WORD              THE OUTPUT LITERAL
C
      CHARACTER*8 TEMP,BLANK,IBLANK
      CHARACTER*1 IDIG(6),JDIG(10)
C
      DOUBLE PRECISION WORD
      EQUIVALENCE (IDIG(1),TEMP)
      EQUIVALENCE (BLANK,IBLANK)
C
      DATA JDIG/'0','1','2','3','4','5','6','7','8','9'/,I1/1/,I10/10/
      DATA BLANK/'      '/
C
C  DETERMINE THE DIGITS OF THE INTEGER (BASE 10).
C
      NCHAR=8
      TEMP=BLANK
C
      DO 5 I=1,NCHAR
      IDIG(I)=IBLANK
    5 CONTINUE
C
      K=INT
      DO 10 I=1,NCHAR
      J=NCHAR+1-I
      N=MOD(K,I10)+I1
      K=K/I10
      IDIG(J)=JDIG(N)
      IF(K.LE.0)GO TO 20
   10 CONTINUE
   20 CONTINUE
C
C  ENCODE THE DIGITS, THEN LEFT JUSTIFY.
C
      CALL FLTCHR(WORD,TEMP,I1)
      CALL SQUEZE(WORD,I1)
      RETURN
      END
      
CDECK SQUEZE
C********************************************************************
      SUBROUTINE SQUEZE(WORD,NUM)
C********************************************************************
C  SQUEZE LEFT JUSTIFIES THE NUM WORDS IN WORD AND SQUEEZES OUT
C  IMBEDDED BLANKS.  GOWEST LEFT JUSTIFIES BUT RETAINS IMBEDDED BLANKS.
C
C                       LATEST VERSION 10/06/95 AT CHOSUN U.
C
C     WORD              INPUT ARRAY
C     NUM               NO. OF WORDS IN INPUT ARRAY
C
      DOUBLE PRECISION WORD
C
      CHARACTER*8 TEMP8
      CHARACTER*1 TEMP1, BLANK
C
      DIMENSION WORD(NUM), TEMP1(8)
      EQUIVALENCE (TEMP1(1), TEMP8)
      DATA BLANK/' '/
C
      IENTRY=1
      GO TO 10
      ENTRY GOWEST (WORD,NUM)
      IENTRY=2
   10 CONTINUE
      IF( NUM.LE.0 ) RETURN
C
C  LOOP OVER WORDS.  CONVERT EACH EIGHT-CHARACTER WORD (TEMP8) TO EIGHT
C  ONE-CHARACTER WORDS(TEMP1).
C
      DO 30 N=1,NUM
      CALL CHRFLT(TEMP8,WORD(N),1)
C
C  FIND THE FIRST NON-BLANK CHARACTER.  IF THE WORD IS ENTIRELY BLANK
C  DO NOTHING.
C
      DO 16 J=1,8
      IF( TEMP1(J).NE.BLANK ) GO TO 18
   16 CONTINUE
      GO TO 30
C
C  LEFT JUSTIFY.
C
   18 CONTINUE
      K=1
      DO 22 I=J,8
      IF( IENTRY.EQ.1 .AND. TEMP1(I).EQ.BLANK ) GO TO 22
      TEMP1(K)=TEMP1(I)
      K=K+1
   22 CONTINUE
C
C  BLANK OUT REMAINING CHARACTERS.
C
      IF( K.GT.8 ) GO TO 28
      DO 24 I=K,8
      TEMP1(I)=BLANK
   24 CONTINUE
   28 CONTINUE
      CALL FLTCHR(WORD(N),TEMP8,1)
   30 CONTINUE
      RETURN
      END
C
      SUBROUTINE CFCOPY
C
C  CFCOPY TRANSFERS  A CHARACTER WORD FROM /FLQQCH/ TO
C  /CHQQFL/.  THE SUBROUTINE IS CALLED BY CHRFLT.
C
      DOUBLE PRECISION C,F
      COMMON /CHQQFL/ C
      COMMON /FLQQCH/ F
      C=F
      RETURN
C
      END
      
CDECK FLTCHR
      SUBROUTINE FLTCHR(FLT,CHR,LEN)
C
C     FLTCHR STORES CHARACTER DATA IN FLOATING POINT VARIABLE ARRAYS
C     WITHOUT TYPE CONVERSION.
C
C     THE FUNCTIONS PROVIDED BY FLTCHR AND CHRFLT ARE VIOLATIONS OF THE
C     FORTRAN-77 STANDARD, AND, IN FACT, THE PURPOSE OF THE ROUTINES IS
C     DELIBERATELY TO CIRCUMVENT THE STANDARD IN ORDER TO PERMIT
C     PROGRAMS TO MIX CHARACTER AND NUMERICAL DATA IN THE BPOINTER
C     CONTAINER AND IN STANDARD-FILE-FORMAT RECORDS.
C
C                       LATEST VERSION 10/06/95
C
C     FLT               OUTPUT ARRAY
C     CHR               INPUT ARRAY
C     LEN               NUMBER OF ELEMENTS
C
      DOUBLE PRECISION FLT, F
C
      CHARACTER*8 CHR, C
      COMMON /CHQQFL/ C
      COMMON /FLQQCH/ F
C
      DIMENSION CHR(*), FLT(*)
C
      IF( LEN.LE.0 ) RETURN
      DO 10 I=1,LEN
C
      C=CHR(I)
      CALL FCCOPY
      FLT(I)=F
C
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE FCCOPY
C
C  FCCOPY TRANSFERS A FLOATING POINT WORD FROM /CHQQFL/ TO
C  /FLQQCH/.  THE SUBROUTINE IS CALLED BY FLTCHR.
C
      DOUBLE PRECISION C,F
      COMMON /CHQQFL/ C
      COMMON /FLQQCH/ F
      F=C
      RETURN
      END
