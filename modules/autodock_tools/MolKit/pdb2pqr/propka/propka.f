C
C      PROPKA: A PROTEIN PKA PREDICTOR
C      VERSION 1.0
C      BY HUI LI
C
C      04/25/2004, IOWA CITY
C***********************************************************
C
C      MODIFIED FOR USE WITH PDB2PQR BY TODD DOLINSKY
C      07/19/2005
C
C***********************************************************
C
C      THIS PROGRAM PREDICTS PROTEIN PKA VALUES 
C      ACCORDING TO THE EMPIRICAL RULES PROPOSED BY:
C
C      HUI LI, ANDREW D. ROBERTSON AND JAN H. JENSEN
C
C***********************************************************
C
C PropKa 1.00: a program for protein pKa predictions
C Portion Copyright (C) 2005 Jan H. Jensen and Hui Li
C
C This program is free software; you can redistribute it and/or
C modify it under the terms of the GNU General Public License as
C published by the Free Software Foundation; either version 2 of
C the License, or  (at your option) any later version.
C
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this program; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston,
C MA 02111-1307 USA
C
C For questions or comments contact Jan-Jensen@uiowa.edu
C
      PROGRAM  PROPKA
C     ****************
C     INITIALIZE THE PDB FILE
C
       CHARACTER HEAD*6, SPACE7*7, HYDRGN*1, SPAC56*56,
     $           PDB*(700000), LINE*70
       DIMENSION HEAD(100000)
  
       NATOM=0
       OPEN (12, FILE='PROPKATMP1', STATUS='OLD', FORM='FORMATTED',
     $       ACCESS='SEQUENTIAL')
       OPEN (13, FILE='PROPKATMP2', STATUS='NEW', FORM='FORMATTED',
     $       ACCESS='SEQUENTIAL')
       DO J=1,100000
         READ (12, '(A6,A7,A1,A)', END=100)
     $        HEAD(1), SPACE7,HYDRGN,SPAC56
         IF ((HEAD(1).EQ.'ATOM  ' .OR. HEAD(1).EQ.'HETATM') .AND.
     $       HYDRGN.NE.'H') THEN
           NATOM=NATOM+1
           WRITE (PDB(((NATOM-1)*70)+1:(NATOM*70)+1), '(A6,A7,A1,A)')
     $        HEAD(1), SPACE7,HYDRGN,SPAC56
           WRITE(13, '(A6,A7,A1,A)')
     $        HEAD(1), SPACE7,HYDRGN,SPAC56
         END IF
       END DO
 100   CONTINUE
       CLOSE (12)
C
       IF(NATOM.LT.2)THEN
         WRITE(6,*)'PLEASE CHECK YOUR INPUT PDB FILE!'
         STOP
       END IF 
C
C      - READ INFORMATION -
C
       CALL RUNPROPKA(NATOM, PDB, 'PROPKATMP3')
       
       END

       SUBROUTINE RUNPROPKA(NATOM, PDB, OUTNAME)
C      *******   ******
C
       CHARACTER HEAD*6,NAMATM*5,NAMRES*3,SPACE4*4,
     $           SPACE2*2,NAMCAR*3,TYPCAR*6,AORB*1,
     $           NAMSDC*3,NAMBKB*3,NAMCOL*3,
     $           NAMPRT*3,
     $           TYPHIS*6,TYPCYS*6,TYPLYS*6,TYPTYR*6,
     $           TYPARG*6,
     $           SPACE7*7, HYDRGN*1, SPAC60*60,
     $           SPACE1*1, CHAIN*1, LCARCH*1,
     $           LTYRCH*1, LCYSCH*1, LARGCH*1,
     $           LHISCH*1, LLYSCH*1
C
       LOGICAL   CONV
C
C-PDB
       DIMENSION HEAD(100000),NAMATM(100000),
     $           NAMRES(100000),NUMRES(100000),
     $           X(100000),Y(100000),Z(100000),
     $           CHAIN(100000),
C-BACKBONE O=C-N-H
     $           NUMPRT(10000),NAMPRT(10000),
     $           XPRTON(10000),YPRTON(10000),ZPRTON(10000),
     $           XNITRN(10000),YNITRN(10000),ZNITRN(10000),
     $           XCARBN(10000),YCARBN(10000),ZCARBN(10000),
     $           XOXYGN(10000),YOXYGN(10000),ZOXYGN(10000),
C - PKA SITES -
C-CARBOXYL
     $           NAMCAR(1000),TYPCAR(1000),PKACAR(1000),
     $           LCARO1(1000),LCARO2(1000),LCARRS(1000),
     $           NSDCAR(1000),NBKCAR(1000),NCLCAR(1000),
     $           PK1CAR(1000),PK2CAR(1000),NHBCAR(1000),
     $           LCARCH(1000),
C-HIS
     $           XHISP1(1000),YHISP1(1000),ZHISP1(1000),
     $           XHISP2(1000),YHISP2(1000),ZHISP2(1000),
     $           LHISCG(1000),LHISND(1000),LHISCE(1000),
     $           LHISNE(1000),LHISCD(1000),LHISRS(1000),
     $           TYPHIS(1000),PKAHIS(1000),
     $           NSDHIS(1000),NBKHIS(1000),NCLHIS(1000),
     $           PK1HIS(1000),PK2HIS(1000),NHBHIS(1000),
     $           LHISCH(1000),
C-CYS
     $           LCYSSG(1000),LCYSRS(1000),
     $           TYPCYS(1000),PKACYS(1000),
     $           NSDCYS(1000),NBKCYS(1000),NCLCYS(1000),
     $           PK1CYS(1000),PK2CYS(1000),LCYSCH(1000),
C-TYR
     $           LTYROH(1000),LTYRRS(1000),
     $           TYPTYR(1000),PKATYR(1000),
     $           NSDTYR(1000),NBKTYR(1000),NCLTYR(1000),
     $           PK1TYR(1000),PK2TYR(1000),LTYRCH(1000),
C-LYS
     $           LLYSNZ(1000),LLYSRS(1000),
     $           TYPLYS(1000),PKALYS(1000),
     $           NSDLYS(1000),NBKLYS(1000),NCLLYS(1000),
     $           PK1LYS(1000),PK2LYS(1000),LLYSCH(1000),
C-ARG
     $           XARGP1(1000),YARGP1(1000),ZARGP1(1000),
     $           XARGP2(1000),YARGP2(1000),ZARGP2(1000),
     $           XARGP3(1000),YARGP3(1000),ZARGP3(1000),
     $           XARGP4(1000),YARGP4(1000),ZARGP4(1000),
     $           XARGP5(1000),YARGP5(1000),ZARGP5(1000),
     $           LARGN1(1000),LARGN2(1000),LARGN3(1000),
     $           LARGCD(1000),LARGCZ(1000),LARGRS(1000),
     $           TYPARG(1000),PKAARG(1000),
     $           NSDARG(1000),NBKARG(1000),NCLARG(1000),
     $           PK1ARG(1000),PK2ARG(1000),LARGCH(1000),
C
C -NONE PKA GROUPS-
C-GLN
     $           LGLNNE(1000),LGLNOE(1000),LGLNCD(1000),
     $           LGLNRS(1000),
     $           XGLNP1(1000),YGLNP1(1000),ZGLNP1(1000),
     $           XGLNP2(1000),YGLNP2(1000),ZGLNP2(1000),
C-ASN
     $           LASNND(1000),LASNOD(1000),LASNCG(1000),
     $           LASNRS(1000),
     $           XASNP1(1000),YASNP1(1000),ZASNP1(1000),
     $           XASNP2(1000),YASNP2(1000),ZASNP2(1000),
C-SER
     $           LSEROH(1000),LSERRS(1000),
C-THR
     $           LTHROH(1000),LTHRRS(1000),
C-TRP
     $           XTRPP1(1000),YTRPP1(1000),ZTRPP1(1000),
     $           LTRPCD(1000),LTRPNE(1000),LTRPCE(1000),
     $           LTRPRS(1000),
C
C-PKA DETERMINANTS
C                1: CAR (ASP/GLU)
C                2: HIS
C                3: CYS
C                4: TYR
C                5: LYS
C                6: ARG
     $           NMASS(10,1000),NLOCAL(10,1000),
     $           NAMSDC(10,1000,30),NUMSDC(10,1000,30),
     $           VALSDC(10,1000,30),
     $           NAMBKB(10,1000,30),NUMBKB(10,1000,30),
     $           VALBKB(10,1000,30),
     $           NAMCOL(10,1000,30),NUMCOL(10,1000,30),
     $           VALCOL(10,1000,30),
     $           TOLSDC(10,1000),TOLBKB(10,1000),TOLCOL(10,1000),
     $           TOLMAS(10,1000),TOLLOC(10,1000)
C
C SUBROUTINE-RELATED VARIABLES
C
C
      INTEGER LINELEN, I, J
      CHARACTER PDB*(*), OUTNAME*(100), LINE*(70)
C
C
C      **********************
C      STEP 1. READ PDB FILES
C      **********************
C
       NPRTON=1
       NCAR=0
       NLYS=0
       NTYR=0
       NSER=0
       NTHR=0
       NARG=0
       NASN=0
       NGLN=0
       NHIS=0
       NCYS=0
       NTRP=0

       OPEN (11, FILE=OUTNAME, STATUS='NEW', FORM='FORMATTED',
     $       ACCESS='SEQUENTIAL') 

       LINELEN=70
       DO I = 1, NATOM
          LINE=PDB(((I-1)*LINELEN)+1:(I*LINELEN)+1)
          READ(LINE,'(A6,I5,A5,A1,A3,A1,A1,I4,A4,F8.3,F8.3,F8.3)') 
     $              HEAD(I), NUMATM, NAMATM(I), AORB,
     $              NAMRES(I), SPACE1, CHAIN(I), 
     $              NUMRES(I), SPACE4, X(I), Y(I), Z(I)
C
         IF(NAMRES(I).EQ.'ASP' .OR. NAMRES(I).EQ.'GLU') THEN
           IF(NAMATM(I).EQ.'  N  ') THEN
             NCAR=NCAR+1
             LCARRS(NCAR)=NUMRES(I)
             NAMCAR(NCAR)=NAMRES(I)
             LCARCH(NCAR)=CHAIN(I)
           END IF
           IF(NAMATM(I).EQ.'  OD1' .OR. NAMATM(I).EQ.'  OE1')
     $                  LCARO1(NCAR)=I
           IF(NAMATM(I).EQ.'  OD2' .OR. NAMATM(I).EQ.'  OE2')
     $                  LCARO2(NCAR)=I
         END IF
         IF(NAMATM(I).EQ.'  OXT')THEN
           NCAR=NCAR+1
           LCARRS(NCAR)=NUMRES(I)
           LCARCH(NCAR)=CHAIN(I)
           NAMCAR(NCAR)='C- '
           LCARO1(NCAR)=I
           DO K=1,50
             IF(NAMATM(I-K).EQ.'  O  ')THEN
               LCARO2(NCAR)=I-K
               GO TO 300
             ELSE IF(NAMATM(I+K).EQ.'  O  ')THEN
               LCARO2(NCAR)=I+K
               GO TO 300
             END IF
           END DO
         END IF
 300     CONTINUE
C
         IF((NAMRES(I).EQ.'LYS'.AND.NAMATM(I).EQ.'  NZ ') .OR.
     $                  (I.EQ.1.AND.NAMATM(I).EQ.'  N  ') )THEN
           NLYS=NLYS+1
           LLYSNZ(NLYS)=I
           LLYSRS(NLYS)=NUMRES(I)
           LLYSCH(NLYS)=CHAIN(I)
         END IF
         IF(NAMRES(I).EQ.'TYR'.AND.NAMATM(I).EQ.'  OH ') THEN
           NTYR=NTYR+1
           LTYROH(NTYR)=I
           LTYRRS(NTYR)=NUMRES(I)
           LTYRCH(NTYR)=CHAIN(I)
         END IF
         IF(NAMRES(I).EQ.'CYS'.AND.NAMATM(I).EQ.'  SG ') THEN
           NCYS=NCYS+1
           LCYSSG(NCYS)=I
           LCYSRS(NCYS)=NUMRES(I)
           LCYSCH(NCYS)=CHAIN(I)
         END IF
         IF(NAMRES(I).EQ.'SER'.AND.NAMATM(I).EQ.'  OG ') THEN
           NSER=NSER+1
           LSEROH(NSER)=I
           LSERRS(NSER)=NUMRES(I)
         END IF
         IF(NAMRES(I).EQ.'THR'.AND.NAMATM(I).EQ.'  OG1') THEN
           NTHR=NTHR+1
           LTHROH(NTHR)=I
           LTHRRS(NTHR)=NUMRES(I)
         END IF
         IF(NAMRES(I).EQ.'GLN')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NGLN=NGLN+1
             LGLNRS(NGLN)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  NE2')LGLNNE(NGLN)=I
           IF(NAMATM(I).EQ.'  OE1')LGLNOE(NGLN)=I
           IF(NAMATM(I).EQ.'  CD ')LGLNCD(NGLN)=I
         END IF
         IF(NAMRES(I).EQ.'ASN')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NASN=NASN+1
             LASNRS(NASN)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  ND2')LASNND(NASN)=I
           IF(NAMATM(I).EQ.'  OD1')LASNOD(NASN)=I
           IF(NAMATM(I).EQ.'  CG ')LASNCG(NASN)=I
         END IF
         IF(NAMRES(I).EQ.'ARG')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NARG=NARG+1
             LARGRS(NARG)=NUMRES(I)
             LARGCH(NARG)=CHAIN(I)
           END IF
           IF(NAMATM(I).EQ.'  NE ')LARGN1(NARG)=I
           IF(NAMATM(I).EQ.'  NH2')LARGN2(NARG)=I
           IF(NAMATM(I).EQ.'  NH1')LARGN3(NARG)=I
           IF(NAMATM(I).EQ.'  CZ ')LARGCZ(NARG)=I
           IF(NAMATM(I).EQ.'  CD ')LARGCD(NARG)=I
         END IF
         IF(NAMRES(I).EQ.'HIS')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NHIS=NHIS+1
             LHISRS(NHIS)=NUMRES(I)
             LHISCH(NHIS)=CHAIN(I)
           END IF
           IF(NAMATM(I).EQ.'  CG ')LHISCG(NHIS)=I
           IF(NAMATM(I).EQ.'  ND1')LHISND(NHIS)=I
           IF(NAMATM(I).EQ.'  CE1')LHISCE(NHIS)=I
           IF(NAMATM(I).EQ.'  NE2')LHISNE(NHIS)=I
           IF(NAMATM(I).EQ.'  CD2')LHISCD(NHIS)=I
         END IF
         IF(NAMRES(I).EQ.'TRP')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NTRP=NTRP+1
             LTRPRS(NTRP)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  CD1')LTRPCD(NTRP)=I
           IF(NAMATM(I).EQ.'  NE1')LTRPNE(NTRP)=I
           IF(NAMATM(I).EQ.'  CE2')LTRPCE(NTRP)=I
         END IF
       END DO
 200   CONTINUE
       CLOSE (10)
C
C
C      WRITE(6,*)'THERE ARE', NCAR ,'CARBOXYL GROUPS IN THE PDB FILE'
C      DO ICAR=1, NCAR
C        WRITE(6,'(A3,I5)') NAMCAR(ICAR),LCARRS(ICAR)
C      END DO
C      WRITE(6,*)'THERE ARE', NLYS ,'LYS GROUPS IN THE PDB FILE'
C      DO ILYS=1, NLYS
C        WRITE(6,'(A3,I5)') 'LYS',LLYSRS(ILYS)
C      END DO
C      WRITE(6,*)'THERE ARE', NTYR ,'TYR GROUPS IN THE PDB FILE'
C      DO ITYR=1, NTYR
C        WRITE(6,'(A3,I5)') 'TYR',LTYRRS(ITYR)
C      END DO
C      WRITE(6,*)'THERE ARE', NSER ,'SER GROUPS IN THE PDB FILE'
C      DO ISER=1, NSER
C        WRITE(6,'(A3,I5)') 'SER',LSERRS(ISER)
C      END DO
C      WRITE(6,*)'THERE ARE', NTHR ,'THR GROUPS IN THE PDB FILE'
C      DO ITHR=1, NTHR
C        WRITE(6,'(A3,I5)') 'THR',LTHRRS(ITHR)
C      END DO
C
C
C
C      **********************
C      STEP 2. ASSIGN H ATOMS
C      **********************
C
C      -- DETERMINE BACKBONE PROTON POSITIONS --
C
       NUMPRT(1)=NUMRES(1)
       NAMPRT(1)=NAMRES(1)
       DO IATOM = 2, NATOM
         IF(NAMATM(IATOM).EQ.'  C  ')THEN
           XC=X(IATOM)
           YC=Y(IATOM)
           ZC=Z(IATOM)
         END IF
         IF(NAMATM(IATOM).EQ.'  O  ')THEN
           XO=X(IATOM)
           YO=Y(IATOM)
           ZO=Z(IATOM)
         END IF
         IF(NAMATM(IATOM).EQ.'  N  '.AND.
     $      NUMRES(IATOM).GT.NUMRES(IATOM-1).AND.
     $      NAMRES(IATOM).NE.'PRO')THEN
           XN=X(IATOM)
           YN=Y(IATOM)
           ZN=Z(IATOM)
           VECNRM=SQRT((XC-XO)**2+(YC-YO)**2+(ZC-ZO)**2)
           XVEC=(XC-XO)/VECNRM
           YVEC=(YC-YO)/VECNRM
           ZVEC=(ZC-ZO)/VECNRM
C          -- NOTE N-H BOND LENGTH = 1 ANGSTROM 
           XH=XN+XVEC
           YH=YN+YVEC
           ZH=ZN+ZVEC
           NPRTON=NPRTON+1
           XPRTON(NPRTON)=XH
           YPRTON(NPRTON)=YH
           ZPRTON(NPRTON)=ZH
           XNITRN(NPRTON)=XN
           YNITRN(NPRTON)=YN
           ZNITRN(NPRTON)=ZN
           XOXYGN(NPRTON)=XO
           YOXYGN(NPRTON)=YO
           ZOXYGN(NPRTON)=ZO
           XCARBN(NPRTON)=XC
           YCARBN(NPRTON)=YC
           ZCARBN(NPRTON)=ZC
           NUMPRT(NPRTON)=NUMRES(IATOM)
           NAMPRT(NPRTON)=NAMRES(IATOM)
         END IF
       END DO
C      WRITE(6,*)NPRTON,'AMIDE PROTONS ASSIGNED'
C
C      -- DETERMINE GLN PROTON POSITIONS --
C
       NGLNP=0
       DO IGLN = 1, NGLN
         XC=X(LGLNCD(IGLN))
         YC=Y(LGLNCD(IGLN))
         ZC=Z(LGLNCD(IGLN))
         XO=X(LGLNOE(IGLN))
         YO=Y(LGLNOE(IGLN))
         ZO=Z(LGLNOE(IGLN))
         XN=X(LGLNNE(IGLN))
         YN=Y(LGLNNE(IGLN))
         ZN=Z(LGLNNE(IGLN))
         VECNRM=SQRT((XC-XO)**2+(YC-YO)**2+(ZC-ZO)**2)
         XVEC=(XC-XO)/VECNRM
         YVEC=(YC-YO)/VECNRM
         ZVEC=(ZC-ZO)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN+XVEC
         YH=YN+YVEC
         ZH=ZN+ZVEC
         NGLNP=NGLNP+1
         XGLNP1(IGLN)=XH
         YGLNP1(IGLN)=YH
         ZGLNP1(IGLN)=ZH
C
         XON=(XO+XN)/2.0
         YON=(YO+YN)/2.0
         ZON=(ZO+ZN)/2.0
         VECNRM=SQRT((XC-XON)**2+(YC-YON)**2+(ZC-ZON)**2)
         XVEC=(XC-XON)/VECNRM
         YVEC=(YC-YON)/VECNRM
         ZVEC=(ZC-ZON)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN-XVEC
         YH=YN-YVEC
         ZH=ZN-ZVEC
         NGLNP=NGLNP+1
         XGLNP2(IGLN)=XH
         YGLNP2(IGLN)=YH
         ZGLNP2(IGLN)=ZH
       END DO
C      WRITE(6,*)NGLNP,'GLN PROTONS ASSIGNED'
C
C      -- DETERMINE ASN PROTON POSITIONS --
C
       NASNP=0
       DO IASN = 1, NASN
         XC=X(LASNCG(IASN))
         YC=Y(LASNCG(IASN))
         ZC=Z(LASNCG(IASN))
         XO=X(LASNOD(IASN))
         YO=Y(LASNOD(IASN))
         ZO=Z(LASNOD(IASN))
         XN=X(LASNND(IASN))
         YN=Y(LASNND(IASN))
         ZN=Z(LASNND(IASN))
         VECNRM=SQRT((XC-XO)**2+(YC-YO)**2+(ZC-ZO)**2)
         XVEC=(XC-XO)/VECNRM
         YVEC=(YC-YO)/VECNRM
         ZVEC=(ZC-ZO)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN+XVEC
         YH=YN+YVEC
         ZH=ZN+ZVEC
         NASNP=NASNP+1
         XASNP1(IASN)=XH
         YASNP1(IASN)=YH
         ZASNP1(IASN)=ZH
C
         XON=(XO+XN)/2.0
         YON=(YO+YN)/2.0
         ZON=(ZO+ZN)/2.0
         VECNRM=SQRT((XC-XON)**2+(YC-YON)**2+(ZC-ZON)**2)
         XVEC=(XC-XON)/VECNRM
         YVEC=(YC-YON)/VECNRM
         ZVEC=(ZC-ZON)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN-XVEC
         YH=YN-YVEC
         ZH=ZN-ZVEC
         NASNP=NASNP+1
         XASNP2(IASN)=XH
         YASNP2(IASN)=YH
         ZASNP2(IASN)=ZH
       END DO
C      WRITE(6,*)NASNP,'ASN PROTONS ASSIGNED'
C
C
C      -- DETERMINE ARG PROTON POSITIONS --
C
       NARGP=0
       DO IARG = 1, NARG
         XCD=X(LARGCD(IARG))
         YCD=Y(LARGCD(IARG))
         ZCD=Z(LARGCD(IARG))
         XCZ=X(LARGCZ(IARG))
         YCZ=Y(LARGCZ(IARG))
         ZCZ=Z(LARGCZ(IARG))
         XN1=X(LARGN1(IARG))
         YN1=Y(LARGN1(IARG))
         ZN1=Z(LARGN1(IARG))
         XN2=X(LARGN2(IARG))
         YN2=Y(LARGN2(IARG))
         ZN2=Z(LARGN2(IARG))
         XN3=X(LARGN3(IARG))
         YN3=Y(LARGN3(IARG))
         ZN3=Z(LARGN3(IARG))
C
         XO=(XCD+XCZ)/2.0
         YO=(YCD+YCZ)/2.0
         ZO=(ZCD+ZCZ)/2.0
         VECNRM=SQRT((XN1-XO)**2+(YN1-YO)**2+(ZN1-ZO)**2)
         XVEC=(XN1-XO)/VECNRM
         YVEC=(YN1-YO)/VECNRM
         ZVEC=(ZN1-ZO)/VECNRM
         NARGP=NARGP+1
         XARGP1(IARG)=XN1+XVEC
         YARGP1(IARG)=YN1+YVEC
         ZARGP1(IARG)=ZN1+ZVEC
         NARGP=NARGP+1
         XARGP2(IARG)=XN2+XVEC
         YARGP2(IARG)=YN2+YVEC
         ZARGP2(IARG)=ZN2+ZVEC
C
         XO=XN1
         YO=YN1
         ZO=ZN1
         VECNRM=SQRT((XCZ-XO)**2+(YCZ-YO)**2+(ZCZ-ZO)**2)
         XVEC=(XCZ-XO)/VECNRM
         YVEC=(YCZ-YO)/VECNRM
         ZVEC=(ZCZ-ZO)/VECNRM
         NARGP=NARGP+1
         XARGP3(IARG)=XN2+XVEC
         YARGP3(IARG)=YN2+YVEC
         ZARGP3(IARG)=ZN2+ZVEC
         NARGP=NARGP+1
         XARGP4(IARG)=XN3+XVEC
         YARGP4(IARG)=YN3+YVEC
         ZARGP4(IARG)=ZN3+ZVEC
C
         XO=XN1
         YO=YN1
         ZO=ZN1
         VECNRM=SQRT((XCD-XO)**2+(YCD-YO)**2+(ZCD-ZO)**2)
         XVEC=(XCD-XO)/VECNRM
         YVEC=(YCD-YO)/VECNRM
         ZVEC=(ZCD-ZO)/VECNRM
         NARGP=NARGP+1
         XARGP5(IARG)=XN3+XVEC
         YARGP5(IARG)=YN3+YVEC
         ZARGP5(IARG)=ZN3+ZVEC
C
       END DO
C      WRITE(6,*)NARGP,'ARG PROTONS ASSIGNED'
C
C
C      -- DETERMINE TRP PROTON POSITIONS --
C
       NTRPP=0
       DO ITRP = 1, NTRP
         XNE=X(LTRPNE(ITRP))
         YNE=Y(LTRPNE(ITRP))
         ZNE=Z(LTRPNE(ITRP))
         XCE=X(LTRPCE(ITRP))
         YCE=Y(LTRPCE(ITRP))
         ZCE=Z(LTRPCE(ITRP))
         XCD=X(LTRPCD(ITRP))
         YCD=Y(LTRPCD(ITRP))
         ZCD=Z(LTRPCD(ITRP))
C
         XO=(XCD+XCE)/2.0
         YO=(YCD+YCE)/2.0
         ZO=(ZCD+ZCE)/2.0
         VECNRM=SQRT((XNE-XO)**2+(YNE-YO)**2+(ZNE-ZO)**2)
         XVEC=(XNE-XO)/VECNRM
         YVEC=(YNE-YO)/VECNRM
         ZVEC=(ZNE-ZO)/VECNRM
         NTRPP=NTRPP+1
         XTRPP1(ITRP)=XNE+XVEC
         YTRPP1(ITRP)=YNE+YVEC
         ZTRPP1(ITRP)=ZNE+ZVEC
       END DO
C      WRITE(6,*)NTRPP,'TRP PROTONS ASSIGNED'
C
C
C      -- DETERMINE HIS PROTON POSITIONS --
C
       NHISP=0
       DO IHIS = 1, NHIS
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XO=(XCG+XCE)/2.0
         YO=(YCG+YCE)/2.0
         ZO=(ZCG+ZCE)/2.0
         VECNRM=SQRT((XND-XO)**2+(YND-YO)**2+(ZND-ZO)**2)
         XVEC=(XND-XO)/VECNRM
         YVEC=(YND-YO)/VECNRM
         ZVEC=(ZND-ZO)/VECNRM
         NHISP=NHISP+1
         XHISP1(IHIS)=XND+XVEC
         YHISP1(IHIS)=YND+YVEC
         ZHISP1(IHIS)=ZND+ZVEC
         XO=(XCD+XCE)/2.0
         YO=(YCD+YCE)/2.0
         ZO=(ZCD+ZCE)/2.0
         VECNRM=SQRT((XNE-XO)**2+(YNE-YO)**2+(ZNE-ZO)**2)
         XVEC=(XNE-XO)/VECNRM
         YVEC=(YNE-YO)/VECNRM
         ZVEC=(ZNE-ZO)/VECNRM
         NHISP=NHISP+1
         XHISP2(IHIS)=XNE+XVEC
         YHISP2(IHIS)=YNE+YVEC
         ZHISP2(IHIS)=ZNE+ZVEC
       END DO
C      WRITE(6,*)NHISP,'HIS PROTONS ASSIGNED'
C
C
       DO I=1,1000
         TYPCAR(I)='SUFACE'
         TYPHIS(I)='SUFACE'
         TYPCYS(I)='SUFACE'
         TYPTYR(I)='SUFACE'
         TYPLYS(I)='SUFACE'
         TYPARG(I)='SUFACE'
       END DO
C
C
C        -- FIND DISULFIDE BONDS --
C
       DO ICYS=1,NCYS
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
         DO JCYS=1,NCYS
           IF(JCYS.NE.ICYS.AND.
     $        LCYSRS(ICYS).NE.LCYSRS(JCYS)-3 .AND.
     $        LCYSRS(ICYS).NE.LCYSRS(JCYS)+3   )THEN
             XS=X(LCYSSG(JCYS))
             YS=Y(LCYSSG(JCYS))
             ZS=Z(LCYSSG(JCYS))
             DIS=SQRT((XSG-XS)**2+(YSG-YS)**2+(ZSG-ZS)**2)
             IF(DIS.LT.2.50)THEN
               TYPCYS(ICYS)='BONDED'
             END IF
           END IF
         END DO
       END DO
C
C
C
C      INTRINSIC pKa VALUES WHICH ARE pH INDEPENDENT
C
       DO I=1,1000
         NHBCAR(I)=0
         NHBHIS(I)=0
         NSDCAR(I)=0
         NBKCAR(I)=0
         NCLCAR(I)=0
         NSDHIS(I)=0
         NBKHIS(I)=0
         NCLHIS(I)=0
         NSDCYS(I)=0
         NBKCYS(I)=0
         NCLCYS(I)=0
         NSDTYR(I)=0
         NBKTYR(I)=0
         NCLTYR(I)=0
         NSDLYS(I)=0
         NBKLYS(I)=0
         NCLLYS(I)=0
         NSDARG(I)=0
         NBKARG(I)=0
         NCLARG(I)=0
         DO J=1,10
           NMASS(J,I)=0
           NLOCAL(J,I)=0
           TOLBKB(J,I)=0.0
           TOLSDC(J,I)=0.0
           TOLCOL(J,I)=0.0
           TOLMAS(J,I)=0.0
           TOLLOC(J,I)=0.0
           DO K=1,30
             NAMBKB(J,I,K)='000'
             NUMBKB(J,I,K)=0
             VALBKB(J,I,K)=0.0
             NAMSDC(J,I,K)='000'
             NUMSDC(J,I,K)=0
             VALSDC(J,I,K)=0.0
             NAMCOL(J,I,K)='000'
             NUMCOL(J,I,K)=0
             VALCOL(J,I,K)=0.0
           END DO
         END DO
       END DO
C
       DO ICAR=1, NCAR
         IF(NAMCAR(ICAR).EQ.'C- ')PK1CAR(ICAR)=3.20
         IF(NAMCAR(ICAR).EQ.'ASP')PK1CAR(ICAR)=3.80
         IF(NAMCAR(ICAR).EQ.'GLU')PK1CAR(ICAR)=4.50
       END DO
       DO IHIS=1, NHIS
         PK1HIS(IHIS)=6.50
       END DO
       DO ICYS=1, NCYS
         PK1CYS(ICYS)=99.99
         IF(TYPCYS(ICYS).NE.'BONDED')PK1CYS(ICYS)=9.00
       END DO
       DO ITYR=1, NTYR
         PK1TYR(ITYR)=10.00
       END DO
       PK1LYS(1)=8.0
       DO ILYS=2, NLYS
         PK1LYS(ILYS)=10.50
       END DO
       DO IARG=1, NARG
         PK1ARG(IARG)=12.50
       END DO
C
C
C      *******************
C      STEP 3. DESOLVATION
C      *******************
C
C      ASP/GLU
C
       DO ICAR=1, NCAR
         XO1=X(LCARO1(ICAR))
         YO1=Y(LCARO1(ICAR))
         ZO1=Z(LCARO1(ICAR))
         XO2=X(LCARO2(ICAR))
         YO2=Y(LCARO2(ICAR))
         ZO2=Z(LCARO2(ICAR))
         XO=(XO1+XO2)/2.0
         YO=(YO1+YO2)/2.0
         ZO=(ZO1+ZO2)/2.0
C
         FMASS=0.010
         FLOCAL=0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(1,ICAR)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'ATOM  ' .AND.
     $        NUMRES(IATOM).NE.LCARRS(ICAR))THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF((ABS(XO-XDS).LT.DMASS .AND.
     $         ABS(YO-YDS).LT.DMASS .AND.
     $         ABS(ZO-ZDS).LT.DMASS     )     ) THEN
             DIS=SQRT((XO-XDS)**2+(YO-YDS)**2+(ZO-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(1,ICAR)=NMASS(1,ICAR)+1
             IF(DIS.LT.DLOCL)NLOCAL(1,ICAR)=NLOCAL(1,ICAR)+1
           END IF
           END IF
         END DO
         IF(NMASS(1,ICAR).GT.400)TYPCAR(ICAR)='BURIED'
         TOLMAS(1,ICAR)=MAX(0.0,FMASS*(NMASS(1,ICAR)-400))
         TOLLOC(1,ICAR)=FLOCAL*NLOCAL(1,ICAR)
         PK1CAR(ICAR)=PK1CAR(ICAR)+TOLMAS(1,ICAR)+TOLLOC(1,ICAR)
       END DO
C
C      HIS
C
       DO IHIS=1, NHIS
         XH1=XHISP1(IHIS)
         YH1=YHISP1(IHIS)
         ZH1=ZHISP1(IHIS)
         XH2=XHISP2(IHIS)
         YH2=YHISP2(IHIS)
         ZH2=ZHISP2(IHIS)
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XCT=(XCG+XND+XCE+XNE+XCD)/5.0
         YCT=(YCG+YND+YCE+YNE+YCD)/5.0
         ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL1=4.00
         DLOCL2=6.00
         DMASS=15.50
         NMASS(2,IHIS)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'ATOM  ' .AND.
     $        NUMRES(IATOM).NE.LHISRS(IHIS))THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF((ABS(XCT-XDS).LT.DMASS .AND.
     $         ABS(YCT-YDS).LT.DMASS .AND.
     $         ABS(ZCT-ZDS).LT.DMASS     )     ) THEN
             DIS=SQRT((XCT-XDS)**2+(YCT-YDS)**2+(ZCT-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(2,IHIS)=NMASS(2,IHIS)+1
             IF(DIS.LT.DLOCL1)NLOCAL(2,IHIS)=NLOCAL(2,IHIS)+1
             IF(DIS.LT.DLOCL2)NLOCAL(10,IHIS)=NLOCAL(10,IHIS)+1
           END IF
           END IF
         END DO
         IF(NMASS(2,IHIS).GT.400)TYPHIS(IHIS)='BURIED'
         TOLMAS(2,IHIS)=MIN(0.0,FMASS*(NMASS(2,IHIS)-400))
         IF(NMASS(2,IHIS).GT.400)NLOCAL(2,IHIS)=NLOCAL(10,IHIS)
         TOLLOC(2,IHIS)=FLOCAL*NLOCAL(2,IHIS)
         PK1HIS(IHIS)=PK1HIS(IHIS)+TOLMAS(2,IHIS)+TOLLOC(2,IHIS)
       END DO
C
C      CYS
C
       DO ICYS=1, NCYS
       IF(TYPCYS(ICYS).NE.'BONDED')THEN
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
C
         FMASS=0.010
         FLOCAL=0.070
         DLOCL=3.50
         DMASS=15.50
         NMASS(3,ICYS)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'ATOM  ' .AND.
     $        NUMRES(IATOM).NE.LCYSRS(ICYS))THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XSG-XDS).LT.DMASS .AND.
     $        ABS(YSG-YDS).LT.DMASS .AND.
     $        ABS(ZSG-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XSG-XDS)**2+(YSG-YDS)**2+(ZSG-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(3,ICYS)=NMASS(3,ICYS)+1
             IF(DIS.LT.DLOCL)NLOCAL(3,ICYS)=NLOCAL(3,ICYS)+1
           END IF
           END IF
         END DO
         IF(NMASS(3,ICYS).GT.400)TYPCYS(ICYS)='BURIED'
         TOLMAS(3,ICYS)=MAX(0.0,FMASS*(NMASS(3,ICYS)-400))
         TOLLOC(3,ICYS)=FLOCAL*NLOCAL(3,ICYS)
         PK1CYS(ICYS)=PK1CYS(ICYS)+TOLMAS(3,ICYS)+TOLLOC(3,ICYS)
       END IF
       END DO
C
C      TYR
C
       DO ITYR=1, NTYR
         XOH=X(LTYROH(ITYR))
         YOH=Y(LTYROH(ITYR))
         ZOH=Z(LTYROH(ITYR))
C
         FMASS=0.010
         FLOCAL=0.070
         DLOCL=3.50
         DMASS=15.50
         NMASS(4,ITYR)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'ATOM  ' .AND.
     $        NUMRES(IATOM).NE.LTYRRS(ITYR))THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XOH-XDS).LT.DMASS .AND.
     $        ABS(YOH-YDS).LT.DMASS .AND.
     $        ABS(ZOH-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XOH-XDS)**2+(YOH-YDS)**2+(ZOH-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(4,ITYR)=NMASS(4,ITYR)+1
             IF(DIS.LT.DLOCL)NLOCAL(4,ITYR)=NLOCAL(4,ITYR)+1
           END IF
           END IF
         END DO
         IF(NMASS(4,ITYR).GT.400)TYPTYR(ITYR)='BURIED'
         TOLMAS(4,ITYR)=MAX(0.0,FMASS*(NMASS(4,ITYR)-400))
         TOLLOC(4,ITYR)=FLOCAL*NLOCAL(4,ITYR)
         PK1TYR(ITYR)=PK1TYR(ITYR)+TOLMAS(4,ITYR)+TOLLOC(4,ITYR)
       END DO
C
C      LYS
C
       DO ILYS=1, NLYS
         XNZ=X(LLYSNZ(ILYS))
         YNZ=Y(LLYSNZ(ILYS))
         ZNZ=Z(LLYSNZ(ILYS))
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(5,ILYS)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'ATOM  ' .AND.
     $        NUMRES(IATOM).NE.LLYSRS(ILYS))THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XNZ-XDS).LT.DMASS .AND.
     $        ABS(YNZ-YDS).LT.DMASS .AND.
     $        ABS(ZNZ-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XNZ-XDS)**2+(YNZ-YDS)**2+(ZNZ-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(5,ILYS)=NMASS(5,ILYS)+1
             IF(DIS.LT.DLOCL)NLOCAL(5,ILYS)=NLOCAL(5,ILYS)+1
           END IF
           END IF
         END DO
         IF(NMASS(5,ILYS).GT.400)TYPLYS(ILYS)='BURIED'
         TOLMAS(5,ILYS)=MIN(0.0,FMASS*(NMASS(5,ILYS)-400))
         TOLLOC(5,ILYS)=FLOCAL*NLOCAL(5,ILYS)
         PK1LYS(ILYS)=PK1LYS(ILYS)+TOLMAS(5,ILYS)+TOLLOC(5,ILYS)
       END DO
C
C      ARG
C
       DO IARG=1, NARG
         X1=X(LARGN1(IARG))
         Y1=Y(LARGN1(IARG))
         Z1=Z(LARGN1(IARG))
         X2=X(LARGN2(IARG))
         Y2=Y(LARGN2(IARG))
         Z2=Z(LARGN2(IARG))
         X3=X(LARGN3(IARG))
         Y3=Y(LARGN3(IARG))
         Z3=Z(LARGN3(IARG))
         XCZ=X(LARGCZ(IARG))
         YCZ=Y(LARGCZ(IARG))
         ZCZ=Z(LARGCZ(IARG))
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL=5.00
         DMASS=15.50
         NMASS(6,IARG)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'ATOM  ' .AND.
     $        NUMRES(IATOM).NE.LARGRS(IARG))THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XCZ-XDS).LT.DMASS .AND.
     $        ABS(YCZ-YDS).LT.DMASS .AND.
     $        ABS(ZCZ-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XCZ-XDS)**2+(YCZ-YDS)**2+(ZCZ-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(6,IARG)=NMASS(6,IARG)+1
             IF(DIS.LT.DLOCL)NLOCAL(6,IARG)=NLOCAL(6,IARG)+1
           END IF
           END IF
         END DO
         IF(NMASS(6,IARG).GT.400)TYPARG(IARG)='BURIED'
         TOLMAS(6,IARG)=MIN(0.0,FMASS*(NMASS(6,IARG)-400))
         TOLLOC(6,IARG)=FLOCAL*NLOCAL(6,IARG)
         PK1ARG(IARG)=PK1ARG(IARG)+TOLMAS(6,IARG)+TOLLOC(6,IARG)
       END DO
C
C
C

C      **********************
C      STEP 4. ASP/GLU 
C      **********************
C
C
       DO ICAR=1, NCAR
         XO1=X(LCARO1(ICAR))
         YO1=Y(LCARO1(ICAR))
         ZO1=Z(LCARO1(ICAR))
         XO2=X(LCARO2(ICAR))
         YO2=Y(LCARO2(ICAR))
         ZO2=Z(LCARO2(ICAR))
         XO=(XO1+XO2)/2.0
         YO=(YO1+YO2)/2.0
         ZO=(ZO1+ZO2)/2.0
C
C
C        -- 2. SIDECHAIN INTERACTION --
C
C           - FIND SER -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XO1-XSER).LT.DIS2 .AND.
     $         ABS(YO1-YSER).LT.DIS2 .AND.
     $         ABS(ZO1-ZSER).LT.DIS2     ) .OR.
     $        (ABS(XO2-XSER).LT.DIS2 .AND.
     $         ABS(YO2-YSER).LT.DIS2 .AND.
     $         ABS(ZO2-ZSER).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XSER)**2+(YO1-YSER)**2+(ZO1-ZSER)**2)
             DISO2O=SQRT((XO2-XSER)**2+(YO2-YSER)**2+(ZO2-ZSER)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='SER'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FOH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XO1-XTHR).LT.DIS2 .AND.
     $         ABS(YO1-YTHR).LT.DIS2 .AND.
     $         ABS(ZO1-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XO2-XTHR).LT.DIS2 .AND.
     $         ABS(YO2-YTHR).LT.DIS2 .AND.
     $         ABS(ZO2-ZTHR).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XTHR)**2+(YO1-YTHR)**2+(ZO1-ZTHR)**2)
             DISO2O=SQRT((XO2-XTHR)**2+(YO2-YTHR)**2+(ZO2-ZTHR)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='THR'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FOH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND ASN -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN1=XASNP1(IASN)
           YASN1=YASNP1(IASN)
           ZASN1=ZASNP1(IASN)
           XASN2=XASNP2(IASN)
           YASN2=YASNP2(IASN)
           ZASN2=ZASNP2(IASN)
           IF((ABS(XO1-XASN1).LT.DIS2.AND.
     $         ABS(YO1-YASN1).LT.DIS2.AND.
     $         ABS(ZO1-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XASN1).LT.DIS2.AND.
     $         ABS(YO2-YASN1).LT.DIS2.AND.
     $         ABS(ZO2-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XASN2).LT.DIS2.AND.
     $         ABS(YO1-YASN2).LT.DIS2.AND.
     $         ABS(ZO1-ZASN2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XASN2).LT.DIS2.AND.
     $         ABS(YO2-YASN2).LT.DIS2.AND.
     $         ABS(ZO2-ZASN2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XASN1)**2+(YO1-YASN1)**2+(ZO1-ZASN1)**2)
             DISO21=SQRT((XO2-XASN1)**2+(YO2-YASN1)**2+(ZO2-ZASN1)**2)
             DISO12=SQRT((XO1-XASN2)**2+(YO1-YASN2)**2+(ZO1-ZASN2)**2)
             DISO22=SQRT((XO2-XASN2)**2+(YO2-YASN2)**2+(ZO2-ZASN2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='ASN'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN1=XGLNP1(IGLN)
           YGLN1=YGLNP1(IGLN)
           ZGLN1=ZGLNP1(IGLN)
           XGLN2=XGLNP2(IGLN)
           YGLN2=YGLNP2(IGLN)
           ZGLN2=ZGLNP2(IGLN)
           IF((ABS(XO1-XGLN1).LT.DIS2.AND.
     $         ABS(YO1-YGLN1).LT.DIS2.AND.
     $         ABS(ZO1-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XGLN1).LT.DIS2.AND.
     $         ABS(YO2-YGLN1).LT.DIS2.AND.
     $         ABS(ZO2-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XGLN2).LT.DIS2.AND.
     $         ABS(YO1-YGLN2).LT.DIS2.AND.
     $         ABS(ZO1-ZGLN2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XGLN2).LT.DIS2.AND.
     $         ABS(YO2-YGLN2).LT.DIS2.AND.
     $         ABS(ZO2-ZGLN2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XGLN1)**2+(YO1-YGLN1)**2+(ZO1-ZGLN1)**2)
             DISO21=SQRT((XO2-XGLN1)**2+(YO2-YGLN1)**2+(ZO2-ZGLN1)**2)
             DISO12=SQRT((XO1-XGLN2)**2+(YO1-YGLN2)**2+(ZO1-ZGLN2)**2)
             DISO22=SQRT((XO2-XGLN2)**2+(YO2-YGLN2)**2+(ZO2-ZGLN2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='GLN'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C           - FIND TRP -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ITRP=1,NTRP
           XTRP=XTRPP1(ITRP)
           YTRP=YTRPP1(ITRP)
           ZTRP=ZTRPP1(ITRP)
           IF((ABS(XO1-XTRP).LT.DIS2 .AND.
     $         ABS(YO1-YTRP).LT.DIS2 .AND.
     $         ABS(ZO1-ZTRP).LT.DIS2     ) .OR.
     $        (ABS(XO2-XTRP).LT.DIS2 .AND.
     $         ABS(YO2-YTRP).LT.DIS2 .AND.
     $         ABS(ZO2-ZTRP).LT.DIS2     )     ) THEN
             DISO1P=SQRT((XO1-XTRP)**2+(YO1-YTRP)**2+(ZO1-ZTRP)**2)
             DISO2P=SQRT((XO2-XTRP)**2+(YO2-YTRP)**2+(ZO2-ZTRP)**2)
             DIS=MIN(DISO1P,DISO2P)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='TRP'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LTRPRS(ITRP)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND HIS H-BONDING
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IHIS=1,NHIS
           XHIS1=XHISP1(IHIS)
           YHIS1=YHISP1(IHIS)
           ZHIS1=ZHISP1(IHIS)
           XHIS2=XHISP2(IHIS)
           YHIS2=YHISP2(IHIS)
           ZHIS2=ZHISP2(IHIS)
           IF((ABS(XO1-XHIS1).LT.DIS2.AND.
     $         ABS(YO1-YHIS1).LT.DIS2.AND.
     $         ABS(ZO1-ZHIS1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XHIS1).LT.DIS2.AND.
     $         ABS(YO2-YHIS1).LT.DIS2.AND.
     $         ABS(ZO2-ZHIS1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XHIS2).LT.DIS2.AND.
     $         ABS(YO1-YHIS2).LT.DIS2.AND.
     $         ABS(ZO1-ZHIS2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XHIS2).LT.DIS2.AND.
     $         ABS(YO2-YHIS2).LT.DIS2.AND.
     $         ABS(ZO2-ZHIS2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XHIS1)**2+(YO1-YHIS1)**2+(ZO1-ZHIS1)**2)
             DISO21=SQRT((XO2-XHIS1)**2+(YO2-YHIS1)**2+(ZO2-ZHIS1)**2)
             DISO12=SQRT((XO1-XHIS2)**2+(YO1-YHIS2)**2+(ZO1-ZHIS2)**2)
             DISO22=SQRT((XO2-XHIS2)**2+(YO2-YHIS2)**2+(ZO2-ZHIS2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='HIS'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
C
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))=NAMCAR(ICAR)
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=-FNH*MIN(1.0,VALUE)
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
               IF(NMASS(1,ICAR)+NMASS(2,IHIS).GT.900 .OR.
     $     (NMASS(1,ICAR).GT.400.AND.NMASS(2,IHIS).GT.400))THEN
                 VALSDC(1,ICAR,NSDCAR(ICAR))=-1.60
                 VALSDC(2,IHIS,NSDHIS(IHIS))=+1.60
                 NHBCAR(ICAR)=NHBCAR(ICAR) + 1
                 NHBHIS(IHIS)=NHBHIS(IHIS) + 1
                 PK1CAR(ICAR)=PK1CAR(ICAR)-6.00
                 PK1HIS(IHIS)=PK1HIS(IHIS)+6.00
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C           - FIND CYS-SH -
C
         FSH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED')THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XO1-XS).LT.DIS2.AND.
     $         ABS(YO1-YS).LT.DIS2.AND.
     $         ABS(ZO1-ZS).LT.DIS2    ) .OR.
     $        (ABS(XO2-XS).LT.DIS2.AND.
     $         ABS(YO2-YS).LT.DIS2.AND.
     $         ABS(ZO2-ZS).LT.DIS2    )     ) THEN
             DISO1S=SQRT((XO1-XS)**2+(YO1-YS)**2+(ZO1-ZS)**2)
             DISO2S=SQRT((XO2-XS)**2+(YO2-YS)**2+(ZO2-ZS)**2)
             DIS=MIN(DISO1S,DISO2S)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='CYS'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FSH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
C
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))=NAMCAR(ICAR)
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=-FSH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END IF
         END DO
C
C
C           - FIND TYR-OH -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ITYR=1,NTYR
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF((ABS(XO1-XOH).LT.DIS2.AND.
     $         ABS(YO1-YOH).LT.DIS2.AND.
     $         ABS(ZO1-ZOH).LT.DIS2    ) .OR.
     $        (ABS(XO2-XOH).LT.DIS2.AND.
     $         ABS(YO2-YOH).LT.DIS2.AND.
     $         ABS(ZO2-ZOH).LT.DIS2    )     ) THEN
             DISO1O=SQRT((XO1-XOH)**2+(YO1-YOH)**2+(ZO1-ZOH)**2)
             DISO2O=SQRT((XO2-XOH)**2+(YO2-YOH)**2+(ZO2-ZOH)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='TYR'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FOH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
C
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))=NAMCAR(ICAR)
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=-FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND LYS H-BONDING
C
         FNH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ILYS=1,NLYS
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           FNH=-0.80
           DIS2=4.00
           IF(ILYS.EQ.1)FNH=-1.20
           IF(ILYS.EQ.1)DIS2=4.50
           IF((ABS(XO1-XN).LT.DIS2.AND.
     $         ABS(YO1-YN).LT.DIS2.AND.
     $         ABS(ZO1-ZN).LT.DIS2    ) .OR.
     $        (ABS(XO2-XN).LT.DIS2.AND.
     $         ABS(YO2-YN).LT.DIS2.AND.
     $         ABS(ZO2-ZN).LT.DIS2    )     ) THEN
             DISO1N=SQRT((XO1-XN)**2+(YO1-YN)**2+(ZO1-ZN)**2)
             DISO2N=SQRT((XO2-XN)**2+(YO2-YN)**2+(ZO2-ZN)**2)
             DIS=MIN(DISO1N,DISO2N)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='LYS'
               IF(ILYS.EQ.1)NAMSDC(1,ICAR,NSDCAR(ICAR))='N+ '
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND ARG H-BONDING
C
         FNH=-0.80
         DIS1=2.00
         DIS2=4.00
         DO IARG=1,NARG
           XH1=XARGP1(IARG)
           YH1=YARGP1(IARG)
           ZH1=ZARGP1(IARG)
           XH2=XARGP2(IARG)
           YH2=YARGP2(IARG)
           ZH2=ZARGP2(IARG)
           XH3=XARGP3(IARG)
           YH3=YARGP3(IARG)
           ZH3=ZARGP3(IARG)
           XH4=XARGP4(IARG)
           YH4=YARGP4(IARG)
           ZH4=ZARGP4(IARG)
           XH5=XARGP5(IARG)
           YH5=YARGP5(IARG)
           ZH5=ZARGP5(IARG)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZN2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZN3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZN4).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH5).LT.DIS2.AND.
     $         ABS(YO1-YH5).LT.DIS2.AND.
     $         ABS(ZO1-ZH5).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH5).LT.DIS2.AND.
     $         ABS(YO2-YH5).LT.DIS2.AND.
     $         ABS(ZO2-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS15=SQRT((XO1-XH5)**2+(YO1-YH5)**2+(ZO1-ZH5)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS25=SQRT((XO2-XH5)**2+(YO2-YH5)**2+(ZO2-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             DIS=MIN(DIS,DIS25)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='ARG'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
C
C     -- THERE ARE POSSIBLY 2 H-BONDS BETWEEN ARG AND CARBOXYL
               IF((DIS11.LT.2.2 .AND. DIS22.LT.2.2).OR.
     $            (DIS21.LT.2.2 .AND. DIS12.LT.2.2).OR.
     $            (DIS13.LT.2.2 .AND. DIS24.LT.2.2).OR.
     $            (DIS23.LT.2.2 .AND. DIS14.LT.2.2)    )THEN
                 VALSDC(1,ICAR,NSDCAR(ICAR))=-2.40
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C
C        -- 3. FIND BACKBONE HYDROGEN BONDING --
C
         FBKB=-1.20
         DIS1=2.0
         DIS2=3.5
         DO I=2,NPRTON
           XP=XPRTON(I)
           YP=YPRTON(I)
           ZP=ZPRTON(I)
           XN=XNITRN(I)
           YN=YNITRN(I)
           ZN=ZNITRN(I)
           IF((ABS(XO1-XP).LT.DIS2.AND.
     $         ABS(YO1-YP).LT.DIS2.AND.
     $         ABS(ZO1-ZP).LT.DIS2    ) .OR.
     $        (ABS(XO2-XP).LT.DIS2.AND.
     $         ABS(YO2-YP).LT.DIS2.AND.
     $         ABS(ZO2-ZP).LT.DIS2    )     ) THEN
             DISO1P=SQRT((XO1-XP)**2+(YO1-YP)**2+(ZO1-ZP)**2)
             DISO2P=SQRT((XO2-XP)**2+(YO2-YP)**2+(ZO2-ZP)**2)
             DIS=MIN(DISO1P,DISO2P)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPO1=-(XP-XO1)/DISO1P
             YVPO1=-(YP-YO1)/DISO1P
             ZVPO1=-(ZP-ZO1)/DISO1P
             AGPO1=XVNP*XVPO1 + YVNP*YVPO1 + ZVNP*ZVPO1
             XVPO2=-(XP-XO2)/DISO2P
             YVPO2=-(YP-YO2)/DISO2P
             ZVPO2=-(ZP-ZO2)/DISO2P
             AGPO2=XVNP*XVPO2 + YVNP*YVPO2 + ZVNP*ZVPO2
             AGPO=MAX(AGPO1,AGPO2)
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
               NBKCAR(ICAR)=NBKCAR(ICAR)+1
               NAMBKB(1,ICAR,NBKCAR(ICAR))=NAMPRT(I)
               NUMBKB(1,ICAR,NBKCAR(ICAR))=NUMPRT(I)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALBKB(1,ICAR,NBKCAR(ICAR))=FBKB*VALUE*AGPO
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALBKB(1,ICAR,NBKCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C        -- 4. IONIZABLE INTERACTION --
C
C
C           - FIND CYS(-)
C
         FCOUL=+2.40
         DIS1=4.00
         DIS2=7.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED')THEN
           IF(NMASS(1,ICAR)+NMASS(3,ICYS).GT.900 .OR.
     $     (NMASS(1,ICAR).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XO-XS).LT.DIS2.AND.
     $         ABS(YO-YS).LT.DIS2.AND.
     $         ABS(ZO-ZS).LT.DIS2    )     ) THEN
             DIS=SQRT((XO-XS)**2+(YO-YS)**2+(ZO-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))=NAMCAR(ICAR)
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND TYR(-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO ITYR=1,NTYR
C          IF(NMASS(1,ICAR)+NMASS(4,ITYR).GT.900 .OR.
C    $     (NMASS(1,ICAR).GT.400.AND.NMASS(4,ITYR).GT.400))THEN
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF(ABS(XO-XOH).LT.DIS2.AND.
     $        ABS(YO-YOH).LT.DIS2.AND.
     $        ABS(ZO-ZOH).LT.DIS2    ) THEN
             DIS=SQRT((XO-XOH)**2+(YO-YOH)**2+(ZO-ZOH)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))=NAMCAR(ICAR)
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
             END IF
           END IF
C          END IF
         END DO
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(1,ICAR)+NMASS(5,ILYS).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XO-XN).LT.DIS2.AND.
     $        ABS(YO-YN).LT.DIS2.AND.
     $        ABS(ZO-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XO-XN)**2+(YO-YN)**2+(ZO-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCAR(ICAR)=NCLCAR(ICAR)+1
               NAMCOL(1,ICAR,NCLCAR(ICAR))='LYS'
               IF(ILYS.EQ.1)NAMCOL(1,ICAR,NCLCAR(ICAR))='N+ '
               NUMCOL(1,ICAR,NCLCAR(ICAR))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
C
               NCLLYS(ILYS)=NCLLYS(ILYS)+1
               NAMCOL(5,ILYS,NCLLYS(ILYS))=NAMCAR(ICAR)
               NUMCOL(5,ILYS,NCLLYS(ILYS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(5,ILYS,NCLLYS(ILYS))=-FCOUL*MIN(1.0,VALUE)
               PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(1,ICAR)+NMASS(6,IARG).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF (ABS(XO-XCZ).LT.DIS2.AND.
     $         ABS(YO-YCZ).LT.DIS2.AND.
     $         ABS(ZO-ZCZ).LT.DIS2    )THEN
             DIS=SQRT((XO-XCZ)**2+(YO-YCZ)**2+(ZO-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCAR(ICAR)=NCLCAR(ICAR)+1
               NAMCOL(1,ICAR,NCLCAR(ICAR))='ARG'
               NUMCOL(1,ICAR,NCLCAR(ICAR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
C
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))=NAMCAR(ICAR)
               NUMCOL(6,IARG,NCLARG(IARG))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=-FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END DO
C
      END DO
C
C
C      **********************
C      STEP 5. HIS
C      **********************
C
C
       DO IHIS=1, NHIS
         XH1=XHISP1(IHIS)
         YH1=YHISP1(IHIS)
         ZH1=ZHISP1(IHIS)
         XH2=XHISP2(IHIS)
         YH2=YHISP2(IHIS)
         ZH2=ZHISP2(IHIS)
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XCT=(XCG+XND+XCE+XNE+XCD)/5.0
         YCT=(YCG+YND+YCE+YNE+YCD)/5.0
         ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
C
C
C           - FIND CYS-HIS H-BONDING
C
C        IT SEEMS THAT CYS-S...H...N-HIS SYSTEM
C        TENDS TO BE CYS-S(-)...(+)H-N-HIS
C        PK(CYS) < PK(HIS)
C
         FSN=1.60
         DIS1=3.00
         DIS2=4.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED')THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XH1-XS).LT.DIS2.AND.
     $         ABS(YH1-YS).LT.DIS2.AND.
     $         ABS(ZH1-ZS).LT.DIS2    ) .OR.
     $        (ABS(XH2-XS).LT.DIS2.AND.
     $         ABS(YH2-YS).LT.DIS2.AND.
     $         ABS(ZH2-ZS).LT.DIS2    )     ) THEN
             DISH1S=SQRT((XH1-XS)**2+(YH1-YS)**2+(ZH1-ZS)**2)
             DISH2S=SQRT((XH2-XS)**2+(YH2-YS)**2+(ZH2-ZS)**2)
             DIS=MIN(DISH1S,DISH2S)
             IF(DIS.LT.DIS2)THEN
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))='CYS'
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=FSN*MIN(1.0,VALUE)
C
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='HIS'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=-FSN*MIN(1.0,VALUE)
C
C              STRONG H-BONDING IF BURIED
               IF(NMASS(2,IHIS)+NMASS(3,ICYS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
                 VALSDC(2,IHIS,NSDHIS(IHIS))=+3.60
                 VALSDC(3,ICYS,NSDCYS(ICYS))=-3.60
               END IF
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END IF
         END DO
C
C
C
C           - FIND TYR-OH -
C             HIS-NH...(-)O-TYR DECREASES TYR'S PK
C
         FOH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ITYR=1,NTYR
           XO=X(LTYROH(ITYR))
           YO=Y(LTYROH(ITYR))
           ZO=Z(LTYROH(ITYR))
           IF((ABS(XH1-XO).LT.DIS2.AND.
     $         ABS(YH1-YO).LT.DIS2.AND.
     $         ABS(ZH1-ZO).LT.DIS2    ) .OR.
     $        (ABS(XH2-XO).LT.DIS2.AND.
     $         ABS(YH2-YO).LT.DIS2.AND.
     $         ABS(ZH2-ZO).LT.DIS2    )     ) THEN
             DISH1O=SQRT((XH1-XO)**2+(YH1-YO)**2+(ZH1-ZO)**2)
             DISH2O=SQRT((XH2-XO)**2+(YH2-YO)**2+(ZH2-ZO)**2)
             DIS=MIN(DISH1O,DISH2O)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='HIS'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND ASN -
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF((ABS(XH1-XASN).LT.DIS2.AND.
     $         ABS(YH1-YASN).LT.DIS2.AND.
     $         ABS(ZH1-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XASN).LT.DIS2.AND.
     $         ABS(YH2-YASN).LT.DIS2.AND.
     $         ABS(ZH2-ZASN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XASN)**2+(YH1-YASN)**2+(ZH1-ZASN)**2)
             DISH2=SQRT((XH2-XASN)**2+(YH2-YASN)**2+(ZH2-ZASN)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))='ASN'
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=FOH*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF((ABS(XH1-XGLN).LT.DIS2.AND.
     $         ABS(YH1-YGLN).LT.DIS2.AND.
     $         ABS(ZH1-ZGLN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XGLN).LT.DIS2.AND.
     $         ABS(YH2-YGLN).LT.DIS2.AND.
     $         ABS(ZH2-ZGLN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XGLN)**2+(YH1-YGLN)**2+(ZH1-ZGLN)**2)
             DISH2=SQRT((XH2-XGLN)**2+(YH2-YGLN)**2+(ZH2-ZGLN)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))='GLN'
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=FOH*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C        -- 3. FIND BACKBONE INTERACTION --
C
         FBKB=1.20
         DIS1=2.00
         DIS2=3.50
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(XH1-XBKO).LT.DIS2.AND.
     $         ABS(YH1-YBKO).LT.DIS2.AND.
     $         ABS(ZH1-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH2-XBKO).LT.DIS2.AND.
     $         ABS(YH2-YBKO).LT.DIS2.AND.
     $         ABS(ZH2-ZBKO).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XBKO)**2+(YH1-YBKO)**2+(ZH1-ZBKO)**2)
             DISH2=SQRT((XH2-XBKO)**2+(YH2-YBKO)**2+(ZH2-ZBKO)**2)
             VECNM=SQRT((XBKC-XBKO)**2+(YBKC-YBKO)**2+(ZBKC-ZBKO)**2)
             XVCO=(XBKO-XBKC)/VECNM
             YVCO=(YBKO-YBKC)/VECNM
             ZVCO=(ZBKO-ZBKC)/VECNM
             XVOH1=(XH1-XBKO)/DISH1
             YVOH1=(YH1-YBKO)/DISH1
             ZVOH1=(ZH1-ZBKO)/DISH1
             XVOH2=(XH2-XBKO)/DISH2
             YVOH2=(YH2-YBKO)/DISH2
             ZVOH2=(ZH2-ZBKO)/DISH2
             AGOH1=XVCO*XVOH1+YVCO*YVOH1+ZVCO*ZVOH1
             AGOH2=XVCO*XVOH2+YVCO*YVOH2+ZVCO*ZVOH2
             DIS=MIN(DISH1,DISH2)
             AGOH=MAX(AGOH1,AGOH2)
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NBKHIS(IHIS)=NBKHIS(IHIS)+1
               NAMBKB(2,IHIS,NBKHIS(IHIS))=NAMPRT(I-1)
               NUMBKB(2,IHIS,NBKHIS(IHIS))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(2,IHIS,NBKHIS(IHIS))=FBKB*AGOH*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALBKB(2,IHIS,NBKHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(2,IHIS)+NMASS(5,ILYS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF((ABS(XCT-XN).LT.DIS2.AND.
     $         ABS(YCT-YN).LT.DIS2.AND.
     $         ABS(ZCT-ZN).LT.DIS2    )     ) THEN
             DIS=SQRT((XCT-XN)**2+(YCT-YN)**2+(ZCT-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='LYS'
               IF(ILYS.EQ.1)NAMCOL(2,IHIS,NCLHIS(IHIS))='N+ '
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(2,IHIS)+NMASS(6,IARG).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XCT-XCZ).LT.DIS2.AND.
     $        ABS(YCT-YCZ).LT.DIS2.AND.
     $        ABS(ZCT-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XCT-XCZ)**2+(YCT-YCZ)**2+(ZCT-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='ARG'
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END DO
C
       END DO
C
C
C
C      **********************
C      STEP 6. CYS
C      **********************
C
C
       DO ICYS=1, NCYS
       IF(TYPCYS(ICYS).NE.'BONDED')THEN
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
C
C
C
C        -- 2. SIDECHAIN INTERACTION --
C
C
C           - FIND SER -
C
         FOH=-1.60
         DIS1=3.50
         DIS2=4.50
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XSG-XSER).LT.DIS2 .AND.
     $         ABS(YSG-YSER).LT.DIS2 .AND.
     $         ABS(ZSG-ZSER).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XSER)**2+(YSG-YSER)**2+(ZSG-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='SER'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FOH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=-1.60
         DIS1=3.50
         DIS2=4.50
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XSG-XTHR).LT.DIS2 .AND.
     $         ABS(YSG-YTHR).LT.DIS2 .AND.
     $         ABS(ZSG-ZTHR).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XTHR)**2+(YSG-YTHR)**2+(ZSG-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='THR'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FOH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FNH=-1.60
         DIS1=2.50
         DIS2=3.50
         DO IGLN=1,NGLN
           XGLN1=XGLNP1(IGLN)
           YGLN1=YGLNP1(IGLN)
           ZGLN1=ZGLNP1(IGLN)
           XGLN2=XGLNP2(IGLN)
           YGLN2=YGLNP2(IGLN)
           ZGLN2=ZGLNP2(IGLN)
           IF((ABS(XSG-XGLN1).LT.DIS2.AND.
     $         ABS(YSG-YGLN1).LT.DIS2.AND.
     $         ABS(ZSG-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XSG-XGLN2).LT.DIS2.AND.
     $         ABS(YSG-YGLN2).LT.DIS2.AND.
     $         ABS(ZSG-ZGLN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XSG-XGLN1)**2+(YSG-YGLN1)**2+(ZSG-ZGLN1)**2)
             DIS12=SQRT((XSG-XGLN2)**2+(YSG-YGLN2)**2+(ZSG-ZGLN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='GLN'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C           - FIND ASN -
C
         FNH=-1.60
         DIS1=2.50
         DIS2=3.50
         DO IASN=1,NASN
           XASN1=XASNP1(IASN)
           YASN1=YASNP1(IASN)
           ZASN1=ZASNP1(IASN)
           XASN2=XASNP2(IASN)
           YASN2=YASNP2(IASN)
           ZASN2=ZASNP2(IASN)
           IF((ABS(XSG-XASN1).LT.DIS2.AND.
     $         ABS(YSG-YASN1).LT.DIS2.AND.
     $         ABS(ZSG-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XSG-XASN2).LT.DIS2.AND.
     $         ABS(YSG-YASN2).LT.DIS2.AND.
     $         ABS(ZSG-ZASN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XSG-XASN1)**2+(YSG-YASN1)**2+(ZSG-ZASN1)**2)
             DIS12=SQRT((XSG-XASN2)**2+(YSG-YASN2)**2+(ZSG-ZASN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='ASN'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C           - FIND TRP -
C
         FNH=-1.60
         DIS1=2.50
         DIS2=3.50
         DO ITRP=1,NTRP
           XTRP1=XTRPP1(ITRP)
           YTRP1=YTRPP1(ITRP)
           ZTRP1=ZTRPP1(ITRP)
           IF((ABS(XSG-XTRP1).LT.DIS2.AND.
     $         ABS(YSG-YTRP1).LT.DIS2.AND.
     $         ABS(ZSG-ZTRP1).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XTRP1)**2+(YSG-YTRP1)**2+(ZSG-ZTRP1)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='TRP'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LTRPRS(ITRP)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C        -- 3. FIND BACKBONE INTERACTION --
C
         FBKB=-2.40
         DIS1=3.5
         DIS2=4.5
         DO I=2,NPRTON
           XP=XPRTON(I)
           YP=YPRTON(I)
           ZP=ZPRTON(I)
           XN=XNITRN(I)
           YN=YNITRN(I)
           ZN=ZNITRN(I)
           IF((ABS(XSG-XP).LT.DIS2.AND.
     $         ABS(YSG-YP).LT.DIS2.AND.
     $         ABS(ZSG-ZP).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XP)**2+(YSG-YP)**2+(ZSG-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPS=-(XP-XSG)/DIS
             YVPS=-(YP-YSG)/DIS
             ZVPS=-(ZP-ZSG)/DIS
             AGPS=XVNP*XVPS + YVNP*YVPS + ZVNP*ZVPS
             IF(DIS.LT.DIS2 .AND. AGPS.GT.0.001)THEN
               NBKCYS(ICYS)=NBKCYS(ICYS)+1
               NAMBKB(3,ICYS,NBKCYS(ICYS))=NAMPRT(I)
               NUMBKB(3,ICYS,NBKCYS(ICYS))=NUMPRT(I)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALBKB(3,ICYS,NBKCYS(ICYS))=FBKB*VALUE*AGPS
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALBKB(3,ICYS,NBKCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C        -- COULOMBIC INTERACTIONS --
C
C
C
C           - FIND CYS PAIRS -
C
         FSS=-1.6
         DIS1=3.00
         DIS2=5.00
         DO JCYS=1,NCYS
           IF(LCYSRS(ICYS).EQ.LCYSRS(JCYS)-3)THEN
             XS=X(LCYSSG(JCYS))
             YS=Y(LCYSSG(JCYS))
             ZS=Z(LCYSSG(JCYS))
             DIS=SQRT((XSG-XS)**2+(YSG-YS)**2+(ZSG-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='CYS'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LCYSRS(JCYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FSS*MIN(1.0,VALUE)
C
               NSDCYS(JCYS)=NSDCYS(JCYS)+1
               NAMSDC(3,JCYS,NSDCYS(JCYS))='CYS'
               NUMSDC(3,JCYS,NSDCYS(JCYS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,JCYS,NSDCYS(JCYS))=-FSS*MIN(1.0,VALUE)
C
C              STRONG H-BONDING IF BURIED
               IF(NMASS(3,ICYS)+NMASS(3,JCYS).GT.900 .OR.
     $    (NMASS(3,ICYS).GT.400.AND.NMASS(3,JCYS).GT.400))THEN
                 VALSDC(3,ICYS,NSDCYS(ICYS))=-3.60
                 VALSDC(3,JCYS,NSDCYS(JCYS))=+3.60
               END IF
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
               PK1CYS(JCYS)=PK1CYS(JCYS)+VALSDC(3,JCYS,NSDCYS(JCYS))
             END IF
           END IF
         END DO
C
C
C           - FIND TYR-OH -
C
         FOH=-0.80
         DIS1=3.50
         DIS2=4.50
         DO ITYR=1,NTYR
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF(ABS(XSG-XOH).LT.DIS2.AND.
     $        ABS(YSG-YOH).LT.DIS2.AND.
     $        ABS(ZSG-ZOH).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XOH)**2+(YSG-YOH)**2+(ZSG-ZOH)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='TYR'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FOH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
C
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='CYS'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=-FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C           - FIND TYR(-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO ITYR=1,NTYR
           IF(NMASS(3,ICYS)+NMASS(4,ITYR).GT.900 .OR.
     $     (NMASS(3,ICYS).GT.400.AND.NMASS(4,ITYR).GT.400))THEN
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF(ABS(XSG-XOH).LT.DIS2.AND.
     $        ABS(YSG-YOH).LT.DIS2.AND.
     $        ABS(ZSG-ZOH).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XOH)**2+(YSG-YOH)**2+(ZSG-ZOH)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))='CYS'
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND LYS H-BONDING
C
         FSN=-1.60
         DIS1=3.00
         DIS2=4.00
         DO ILYS=1,NLYS
           FSN=-1.60
           DIS2=4.00
           IF(ILYS.EQ.1)FSN=-2.40
           IF(ILYS.EQ.1)DIS2=4.50
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XSG-XN).LT.DIS2.AND.
     $        ABS(YSG-YN).LT.DIS2.AND.
     $        ABS(ZSG-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='LYS'
               IF(ILYS.EQ.1)NAMSDC(3,ICYS,NSDCYS(ICYS))='N+ '
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FSN*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(3,ICYS)+NMASS(5,ILYS).GT.900 .OR.
     $     (NMASS(3,ICYS).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XSG-XN).LT.DIS2.AND.
     $        ABS(YSG-YN).LT.DIS2.AND.
     $        ABS(ZSG-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='LYS'
               IF(ILYS.EQ.1)NAMCOL(3,ICYS,NCLCYS(ICYS))='N+ '
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
C
               NCLLYS(ILYS)=NCLLYS(ILYS)+1
               NAMCOL(5,ILYS,NCLLYS(ILYS))='CYS'
               NUMCOL(5,ILYS,NCLLYS(ILYS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(5,ILYS,NCLLYS(ILYS))=-FCOUL*MIN(1.0,VALUE)
               PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND ARG(+) H-BONDING
C
         FSN=-1.60
         DIS1=2.50
         DIS2=4.00
         DO IARG=1,NARG
           XH1=XARGP1(IARG)
           YH1=YARGP1(IARG)
           ZH1=ZARGP1(IARG)
           XH2=XARGP2(IARG)
           YH2=YARGP2(IARG)
           ZH2=ZARGP2(IARG)
           XH3=XARGP3(IARG)
           YH3=YARGP3(IARG)
           ZH3=ZARGP3(IARG)
           XH4=XARGP4(IARG)
           YH4=YARGP4(IARG)
           ZH4=ZARGP4(IARG)
           XH5=XARGP5(IARG)
           YH5=YARGP5(IARG)
           ZH5=ZARGP5(IARG)
           IF(
     $        (ABS(XSG-XH1).LT.DIS2.AND.
     $         ABS(YSG-YH1).LT.DIS2.AND.
     $         ABS(ZSG-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH2).LT.DIS2.AND.
     $         ABS(YSG-YH2).LT.DIS2.AND.
     $         ABS(ZSG-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH3).LT.DIS2.AND.
     $         ABS(YSG-YH3).LT.DIS2.AND.
     $         ABS(ZSG-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH4).LT.DIS2.AND.
     $         ABS(YSG-YH4).LT.DIS2.AND.
     $         ABS(ZSG-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH5).LT.DIS2.AND.
     $         ABS(YSG-YH5).LT.DIS2.AND.
     $         ABS(ZSG-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XSG-XH1)**2+(YSG-YH1)**2+(ZSG-ZH1)**2)
             DIS12=SQRT((XSG-XH2)**2+(YSG-YH2)**2+(ZSG-ZH2)**2)
             DIS13=SQRT((XSG-XH3)**2+(YSG-YH3)**2+(ZSG-ZH3)**2)
             DIS14=SQRT((XSG-XH4)**2+(YSG-YH4)**2+(ZSG-ZH4)**2)
             DIS15=SQRT((XSG-XH5)**2+(YSG-YH5)**2+(ZSG-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='ARG'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FSN*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(3,ICYS)+NMASS(6,IARG).GT.900 .OR.
     $     (NMASS(3,ICYS).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XSG-XCZ).LT.DIS2.AND.
     $        ABS(YSG-YCZ).LT.DIS2.AND.
     $        ABS(ZSG-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XSG-XCZ)**2+(YSG-YCZ)**2+(ZSG-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='ARG'
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
C
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))='CYS'
               NUMCOL(6,IARG,NCLARG(IARG))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=-FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END DO
C
C
       END IF
       END DO
C
C
C
C      **********************
C      STEP 7. TYR
C      **********************
C
C
       DO ITYR=1, NTYR
         XOH=X(LTYROH(ITYR))
         YOH=Y(LTYROH(ITYR))
         ZOH=Z(LTYROH(ITYR))
C
C
C
C        -- 2. SIDECHAIN INTERACTION --
C
C
C           - FIND SER -
C
         FOH=-0.80
         DIS1=3.50
         DIS2=4.50
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XOH-XSER).LT.DIS2 .AND.
     $         ABS(YOH-YSER).LT.DIS2 .AND.
     $         ABS(ZOH-ZSER).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XSER)**2+(YOH-YSER)**2+(ZOH-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='SER'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND THR -
C
         FOH=-0.80
         DIS1=3.50
         DIS2=4.50
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XOH-XTHR).LT.DIS2 .AND.
     $         ABS(YOH-YTHR).LT.DIS2 .AND.
     $         ABS(ZOH-ZTHR).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XTHR)**2+(YOH-YTHR)**2+(ZOH-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='THR'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND ASN -
C
         FNH=-0.80
         DIS1=2.50
         DIS2=3.50
         DO IASN=1,NASN
           XASN1=XASNP1(IASN)
           YASN1=YASNP1(IASN)
           ZASN1=ZASNP1(IASN)
           XASN2=XASNP2(IASN)
           YASN2=YASNP2(IASN)
           ZASN2=ZASNP2(IASN)
           IF((ABS(XOH-XASN1).LT.DIS2.AND.
     $         ABS(YOH-YASN1).LT.DIS2.AND.
     $         ABS(ZOH-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XOH-XASN2).LT.DIS2.AND.
     $         ABS(YOH-YASN2).LT.DIS2.AND.
     $         ABS(ZOH-ZASN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XOH-XASN1)**2+(YOH-YASN1)**2+(ZOH-ZASN1)**2)
             DIS12=SQRT((XOH-XASN2)**2+(YOH-YASN2)**2+(ZOH-ZASN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='ASN'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FNH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FNH=-0.80
         DIS1=2.50
         DIS2=3.50
         DO IGLN=1,NGLN
           XGLN1=XGLNP1(IGLN)
           YGLN1=YGLNP1(IGLN)
           ZGLN1=ZGLNP1(IGLN)
           XGLN2=XGLNP2(IGLN)
           YGLN2=YGLNP2(IGLN)
           ZGLN2=ZGLNP2(IGLN)
           IF((ABS(XOH-XGLN1).LT.DIS2.AND.
     $         ABS(YOH-YGLN1).LT.DIS2.AND.
     $         ABS(ZOH-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XOH-XGLN2).LT.DIS2.AND.
     $         ABS(YOH-YGLN2).LT.DIS2.AND.
     $         ABS(ZOH-ZGLN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XOH-XGLN1)**2+(YOH-YGLN1)**2+(ZOH-ZGLN1)**2)
             DIS12=SQRT((XOH-XGLN2)**2+(YOH-YGLN2)**2+(ZOH-ZGLN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='GLN'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FNH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C        -- 3. FIND BACKBONE INTERACTION --
C
         FBKB=-1.20
         DIS1=3.5
         DIS2=4.5
         DO I=2,NPRTON
           XP=XPRTON(I)
           YP=YPRTON(I)
           ZP=ZPRTON(I)
           XN=XNITRN(I)
           YN=YNITRN(I)
           ZN=ZNITRN(I)
           IF((ABS(XOH-XP).LT.DIS2.AND.
     $         ABS(YOH-YP).LT.DIS2.AND.
     $         ABS(ZOH-ZP).LT.DIS2    )     ) THEN
             DIS=SQRT((XOH-XP)**2+(YOH-YP)**2+(ZOH-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPS=-(XP-XOH)/DIS
             YVPS=-(YP-YOH)/DIS
             ZVPS=-(ZP-ZOH)/DIS
             AGPS=XVNP*XVPS + YVNP*YVPS + ZVNP*ZVPS
             IF(DIS.LT.DIS2 .AND. AGPS.GT.0.001)THEN
               NBKTYR(ITYR)=NBKTYR(ITYR)+1
               NAMBKB(4,ITYR,NBKTYR(ITYR))=NAMPRT(I)
               NUMBKB(4,ITYR,NBKTYR(ITYR))=NUMPRT(I)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALBKB(4,ITYR,NBKTYR(ITYR))=FBKB*VALUE*AGPS
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALBKB(4,ITYR,NBKTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C        -- COULOMBIC INTERACTIONS --
C
C
C           - FIND ARG H-BONDING
C
         FOH=-0.80
         DIS1=2.50
         DIS2=4.00
         DO IARG=1,NARG
           XH1=XARGP1(IARG)
           YH1=YARGP1(IARG)
           ZH1=ZARGP1(IARG)
           XH2=XARGP2(IARG)
           YH2=YARGP2(IARG)
           ZH2=ZARGP2(IARG)
           XH3=XARGP3(IARG)
           YH3=YARGP3(IARG)
           ZH3=ZARGP3(IARG)
           XH4=XARGP4(IARG)
           YH4=YARGP4(IARG)
           ZH4=ZARGP4(IARG)
           XH5=XARGP5(IARG)
           YH5=YARGP5(IARG)
           ZH5=ZARGP5(IARG)
           IF(
     $        (ABS(XOH-XH1).LT.DIS2.AND.
     $         ABS(YOH-YH1).LT.DIS2.AND.
     $         ABS(ZOH-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH2).LT.DIS2.AND.
     $         ABS(YOH-YH2).LT.DIS2.AND.
     $         ABS(ZOH-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH3).LT.DIS2.AND.
     $         ABS(YOH-YH3).LT.DIS2.AND.
     $         ABS(ZOH-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH4).LT.DIS2.AND.
     $         ABS(YOH-YH4).LT.DIS2.AND.
     $         ABS(ZOH-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH5).LT.DIS2.AND.
     $         ABS(YOH-YH5).LT.DIS2.AND.
     $         ABS(ZOH-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XOH-XH1)**2+(YOH-YH1)**2+(ZOH-ZH1)**2)
             DIS12=SQRT((XOH-XH2)**2+(YOH-YH2)**2+(ZOH-ZH2)**2)
             DIS13=SQRT((XOH-XH3)**2+(YOH-YH3)**2+(ZOH-ZH3)**2)
             DIS14=SQRT((XOH-XH4)**2+(YOH-YH4)**2+(ZOH-ZH4)**2)
             DIS15=SQRT((XOH-XH5)**2+(YOH-YH5)**2+(ZOH-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='ARG'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(4,ITYR)+NMASS(6,IARG).GT.900 .OR.
     $     (NMASS(4,ITYR).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XOH-XCZ).LT.DIS2.AND.
     $        ABS(YOH-YCZ).LT.DIS2.AND.
     $        ABS(ZOH-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XOH-XCZ)**2+(YOH-YCZ)**2+(ZOH-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))='ARG'
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
C
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))='TYR'
               NUMCOL(6,IARG,NCLARG(IARG))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=-FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END DO
C
       END DO
C
C
C
C
C
       DO I=1, 1000
         PK2CAR(I)=PK1CAR(I)
         PK2HIS(I)=PK1HIS(I)
         PK2CYS(I)=PK1CYS(I)
         PK2TYR(I)=PK1TYR(I)
         PK2LYS(I)=PK1LYS(I)
         PK2ARG(I)=PK1ARG(I)
         PKACAR(I)=PK1CAR(I)
         PKAHIS(I)=PK1HIS(I)
         PKACYS(I)=PK1CYS(I)
         PKATYR(I)=PK1TYR(I)
         PKALYS(I)=PK1LYS(I)
         PKAARG(I)=PK1ARG(I)
       END DO
       DO I=501,1000
         NSDCAR(I)=NSDCAR(I-500)
         NBKCAR(I)=NBKCAR(I-500)
         NCLCAR(I)=NCLCAR(I-500)
         NSDHIS(I)=NSDHIS(I-500)
         NBKHIS(I)=NBKHIS(I-500)
         NCLHIS(I)=NCLHIS(I-500)
         NSDCYS(I)=NSDCYS(I-500)
         NBKCYS(I)=NBKCYS(I-500)
         NCLCYS(I)=NCLCYS(I-500)
         NSDTYR(I)=NSDTYR(I-500)
         NBKTYR(I)=NBKTYR(I-500)
         NCLTYR(I)=NCLTYR(I-500)
         NSDLYS(I)=NSDLYS(I-500)
         NBKLYS(I)=NBKLYS(I-500)
         NCLLYS(I)=NCLLYS(I-500)
         NSDARG(I)=NSDARG(I-500)
         NBKARG(I)=NBKARG(I-500)
         NCLARG(I)=NCLARG(I-500)
       END DO
C
C
C      **********************
C      STEP 8. PKA ITERATION
C      **********************
C         - USUALLY 1 TO 3 ITERATIONS ARE REQUIRED
C
       DO 500 ITER=1,10
C
       DO I=1, 1000
         PK1CAR(I)=PK2CAR(I)
         PK1HIS(I)=PK2HIS(I)
         PK1CYS(I)=PK2CYS(I)
         PK1TYR(I)=PK2TYR(I)
         PK1LYS(I)=PK2LYS(I)
         PK1ARG(I)=PK2ARG(I)
       END DO
       DO I=1,500
         NSDCAR(I)=NSDCAR(I+500)
         NBKCAR(I)=NBKCAR(I+500)
         NCLCAR(I)=NCLCAR(I+500)
         NSDHIS(I)=NSDHIS(I+500)
         NBKHIS(I)=NBKHIS(I+500)
         NCLHIS(I)=NCLHIS(I+500)
         NSDCYS(I)=NSDCYS(I+500)
         NBKCYS(I)=NBKCYS(I+500)
         NCLCYS(I)=NCLCYS(I+500)
         NSDTYR(I)=NSDTYR(I+500)
         NBKTYR(I)=NBKTYR(I+500)
         NCLTYR(I)=NCLTYR(I+500)
         NSDLYS(I)=NSDLYS(I+500)
         NBKLYS(I)=NBKLYS(I+500)
         NCLLYS(I)=NCLLYS(I+500)
         NSDARG(I)=NSDARG(I+500)
         NBKARG(I)=NBKARG(I+500)
         NCLARG(I)=NCLARG(I+500)
       END DO
C
C      **********************
C      ASP/GLU 
C      **********************
C
C
       DO ICAR=1, NCAR
         K=NSDCAR(ICAR)
         DO J=K+1,30
           NAMSDC(1,ICAR,J)='000'
           NUMSDC(1,ICAR,J)=0
           VALSDC(1,ICAR,J)=0.0
         END DO
C
         K=NBKCAR(ICAR)
         DO J=K+1,30
           NAMBKB(1,ICAR,J)='000'
           NUMBKB(1,ICAR,J)=0
           VALBKB(1,ICAR,J)=0.0
         END DO
C
         K=NCLCAR(ICAR)
         DO J=K+1,30
           NAMCOL(1,ICAR,J)='000'
           NUMCOL(1,ICAR,J)=0
           VALCOL(1,ICAR,J)=0.0
         END DO
C
         XO1=X(LCARO1(ICAR))
         YO1=Y(LCARO1(ICAR))
         ZO1=Z(LCARO1(ICAR))
         XO2=X(LCARO2(ICAR))
         YO2=Y(LCARO2(ICAR))
         ZO2=Z(LCARO2(ICAR))
         XO=(XO1+XO2)/2.0
         YO=(YO1+YO2)/2.0
         ZO=(ZO1+ZO2)/2.0
C
C           - FIND HIS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IHIS=1,NHIS
         IF(PKACAR(ICAR).LT.PKAHIS(IHIS))THEN
           IF(NMASS(1,ICAR)+NMASS(2,IHIS).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(2,IHIS).GT.400))THEN
           XCG=X(LHISCG(IHIS))
           YCG=Y(LHISCG(IHIS))
           ZCG=Z(LHISCG(IHIS))
           XND=X(LHISND(IHIS))
           YND=Y(LHISND(IHIS))
           ZND=Z(LHISND(IHIS))
           XCE=X(LHISCE(IHIS))
           YCE=Y(LHISCE(IHIS))
           ZCE=Z(LHISCE(IHIS))
           XNE=X(LHISNE(IHIS))
           YNE=Y(LHISNE(IHIS))
           ZNE=Z(LHISNE(IHIS))
           XCD=X(LHISCD(IHIS))
           YCD=Y(LHISCD(IHIS))
           ZCD=Z(LHISCD(IHIS))
           XCT=(XCG+XND+XCE+XNE+XCD)/5.0
           YCT=(YCG+YND+YCE+YNE+YCD)/5.0
           ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
           IF(ABS(XO-XCT).LT.DIS2.AND.
     $        ABS(YO-YCT).LT.DIS2.AND.
     $        ABS(ZO-ZCT).LT.DIS2    )  THEN
             DIS=SQRT((XO-XCT)**2+(YO-YCT)**2+(ZO-ZCT)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCAR(ICAR)=NCLCAR(ICAR)+1
               NAMCOL(1,ICAR,NCLCAR(ICAR))='HIS'
               NUMCOL(1,ICAR,NCLCAR(ICAR))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
C
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))=NAMCAR(ICAR)
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=-FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND ASP/GLU
C
         FOO=-0.80
         DIS1=2.50
         DIS2=3.50
         IF(ICAR.LT.NCAR)THEN
         DO JCAR=ICAR+1,NCAR
             XJO1=X(LCARO1(JCAR))
             YJO1=Y(LCARO1(JCAR))
             ZJO1=Z(LCARO1(JCAR))
             XJO2=X(LCARO2(JCAR))
             YJO2=Y(LCARO2(JCAR))
             ZJO2=Z(LCARO2(JCAR))
             IF((ABS(XO1-XJO1).LT.DIS2.AND.
     $         ABS(YO1-YJO1).LT.DIS2.AND.
     $         ABS(ZO1-ZJO1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XJO1).LT.DIS2.AND.
     $         ABS(YO2-YJO1).LT.DIS2.AND.
     $         ABS(ZO2-ZJO1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XJO2).LT.DIS2.AND.
     $         ABS(YO1-YJO2).LT.DIS2.AND.
     $         ABS(ZO1-ZJO2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XJO2).LT.DIS2.AND.
     $         ABS(YO2-YJO2).LT.DIS2.AND.
     $         ABS(ZO2-ZJO2).LT.DIS2    )     ) THEN
               DIS11=SQRT((XO1-XJO1)**2
     $                   +(YO1-YJO1)**2
     $                   +(ZO1-ZJO1)**2)
               DIS12=SQRT((XO1-XJO2)**2
     $                   +(YO1-YJO2)**2
     $                   +(ZO1-ZJO2)**2)
               DIS21=SQRT((XO2-XJO1)**2
     $                   +(YO2-YJO1)**2
     $                   +(ZO2-ZJO1)**2)
               DIS22=SQRT((XO2-XJO2)**2
     $                   +(YO2-YJO2)**2
     $                   +(ZO2-ZJO2)**2)
               DIS=MIN(DIS11,DIS12)
               DIS=MIN(DIS,DIS21)
               DIS=MIN(DIS,DIS22)
               IF(DIS.LT.DIS2)THEN
               IF(PKACAR(ICAR).LT.PKACAR(JCAR))THEN
                 NSDCAR(ICAR)=NSDCAR(ICAR)+1
                 NAMSDC(1,ICAR,NSDCAR(ICAR))=NAMCAR(JCAR)
                 NUMSDC(1,ICAR,NSDCAR(ICAR))=LCARRS(JCAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,ICAR,NSDCAR(ICAR))=FOO*MIN(1.0,VALUE)
C
                 NSDCAR(JCAR)=NSDCAR(JCAR)+1
                 NAMSDC(1,JCAR,NSDCAR(JCAR))=NAMCAR(ICAR)
                 NUMSDC(1,JCAR,NSDCAR(JCAR))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,JCAR,NSDCAR(JCAR))=-FOO*MIN(1.0,VALUE)
C
C                STRONG H-BONDING IF BURIED
C
                 IF(NMASS(1,ICAR)+NMASS(1,JCAR).GT.900 .OR.
     $           (NMASS(1,ICAR).GT.400.AND.NMASS(1,JCAR).GT.400))THEN
                   VALSDC(1,ICAR,NSDCAR(ICAR))=-1.60
                   VALSDC(1,JCAR,NSDCAR(JCAR))=+1.60
                 END IF
               ELSE
                 NSDCAR(ICAR)=NSDCAR(ICAR)+1
                 NAMSDC(1,ICAR,NSDCAR(ICAR))=NAMCAR(JCAR)
                 NUMSDC(1,ICAR,NSDCAR(ICAR))=LCARRS(JCAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,ICAR,NSDCAR(ICAR))=-FOO*MIN(1.0,VALUE)
C
                 NSDCAR(JCAR)=NSDCAR(JCAR)+1
                 NAMSDC(1,JCAR,NSDCAR(JCAR))=NAMCAR(ICAR)
                 NUMSDC(1,JCAR,NSDCAR(JCAR))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,JCAR,NSDCAR(JCAR))=FOO*MIN(1.0,VALUE)
C
C                STRONG H-BONDING IF BURIED
C
                 IF(NMASS(1,ICAR)+NMASS(1,JCAR).GT.900 .OR.
     $           (NMASS(1,ICAR).GT.400.AND.NMASS(1,JCAR).GT.400))THEN
                   VALSDC(1,ICAR,NSDCAR(ICAR))=+1.60
                   VALSDC(1,JCAR,NSDCAR(JCAR))=-1.60
                 END IF
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
               PK1CAR(JCAR)=PK1CAR(JCAR)+VALSDC(1,JCAR,NSDCAR(JCAR))
C
               END IF
             END IF
         END DO
         END IF
C
C
C           - FIND ASP/GLU (-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO JCAR=1,NCAR
           IF(JCAR.NE.ICAR .AND. PKACAR(ICAR).GT.PKACAR(JCAR))THEN
             IF(NMASS(1,ICAR)+NMASS(1,JCAR).GT.900 .OR.
     $         (NMASS(1,ICAR).GT.400.AND.NMASS(1,JCAR).GT.400))THEN
             XJO1=X(LCARO1(JCAR))
             YJO1=Y(LCARO1(JCAR))
             ZJO1=Z(LCARO1(JCAR))
             XJO2=X(LCARO2(JCAR))
             YJO2=Y(LCARO2(JCAR))
             ZJO2=Z(LCARO2(JCAR))
             XJO=(XJO1+XJO2)/2.0
             YJO=(YJO1+YJO2)/2.0
             ZJO=(ZJO1+ZJO2)/2.0
             IF(ABS(XO-XJO).LT.DIS2.AND.
     $          ABS(YO-YJO).LT.DIS2.AND.
     $          ABS(ZO-ZJO).LT.DIS2    ) THEN
               DIS=SQRT((XO-XJO)**2
     $                 +(YO-YJO)**2
     $                 +(ZO-ZJO)**2)
               IF(DIS.LT.DIS2)THEN
                 NCLCAR(ICAR)=NCLCAR(ICAR)+1
                 NAMCOL(1,ICAR,NCLCAR(ICAR))=NAMCAR(JCAR)
                 NUMCOL(1,ICAR,NCLCAR(ICAR))=LCARRS(JCAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
                 PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
               END IF
             END IF
             END IF
           END IF
         END DO
C
       END DO
C
C
C
C      **********************
C      HIS
C      **********************
C
C
       DO IHIS=1, NHIS
         K=NSDHIS(IHIS)
         DO J=K+1,30
           NAMSDC(2,IHIS,J)='000'
           NUMSDC(2,IHIS,J)=0
           VALSDC(2,IHIS,J)=0.0
         END DO
C
         K=NBKHIS(IHIS)
         DO J=K+1,30
           NAMBKB(2,IHIS,J)='000'
           NUMBKB(2,IHIS,J)=0
           VALBKB(2,IHIS,J)=0.0
         END DO
C
         K=NCLHIS(IHIS)
         DO J=K+1,30
           NAMCOL(2,IHIS,J)='000'
           NUMCOL(2,IHIS,J)=0
           VALCOL(2,IHIS,J)=0.0
         END DO
C
         XH1=XHISP1(IHIS)
         YH1=YHISP1(IHIS)
         ZH1=ZHISP1(IHIS)
         XH2=XHISP2(IHIS)
         YH2=YHISP2(IHIS)
         ZH2=ZHISP2(IHIS)
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XCT=(XCG+XND+XCE+XNE+XCD)/5.0
         YCT=(YCG+YND+YCE+YNE+YCD)/5.0
         ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
C
C
C
C           - FIND HIS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO JHIS=1,NHIS
         IF(IHIS.NE.JHIS.AND.PKAHIS(IHIS).LT.PKAHIS(JHIS))THEN
           IF(NMASS(2,IHIS)+NMASS(2,JHIS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(2,JHIS).GT.400))THEN
           XCG=X(LHISCG(JHIS))
           YCG=Y(LHISCG(JHIS))
           ZCG=Z(LHISCG(JHIS))
           XND=X(LHISND(JHIS))
           YND=Y(LHISND(JHIS))
           ZND=Z(LHISND(JHIS))
           XCE=X(LHISCE(JHIS))
           YCE=Y(LHISCE(JHIS))
           ZCE=Z(LHISCE(JHIS))
           XNE=X(LHISNE(JHIS))
           YNE=Y(LHISNE(JHIS))
           ZNE=Z(LHISNE(JHIS))
           XCD=X(LHISCD(JHIS))
           YCD=Y(LHISCD(JHIS))
           ZCD=Z(LHISCD(JHIS))
           XCTJ=(XCG+XND+XCE+XNE+XCD)/5.0
           YCTJ=(YCG+YND+YCE+YNE+YCD)/5.0
           ZCTJ=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
           IF(ABS(XCT-XCTJ).LT.DIS2.AND.
     $        ABS(YCT-YCTJ).LT.DIS2.AND.
     $        ABS(ZCT-ZCTJ).LT.DIS2    ) THEN
             DIS=SQRT((XCT-XCTJ)**2+(YCT-YCTJ)**2+(ZCT-ZCTJ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='HIS'
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LHISRS(JHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND CYS(-) COULOMBICS
C
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED'.AND.
     $                   PKAHIS(IHIS).GT.PKACYS(ICYS))THEN
           IF(NMASS(2,IHIS)+NMASS(3,ICYS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XCT-XS).LT.DIS2.AND.
     $         ABS(YCT-YS).LT.DIS2.AND.
     $         ABS(ZCT-ZS).LT.DIS2    )     ) THEN
             DIS=SQRT((XCT-XS)**2+(YCT-YS)**2+(ZCT-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='CYS'
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
C
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='HIS'
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=-FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
       END DO
C
C
C
C      **********************
C      CYS
C      **********************
C
C
       DO ICYS=1, NCYS
         K=NSDCYS(ICYS)
         DO J=K+1,30
           NAMSDC(3,ICYS,J)='000'
           NUMSDC(3,ICYS,J)=0
           VALSDC(3,ICYS,J)=0.0
         END DO
C
         K=NBKCYS(ICYS)
         DO J=K+1,30
           NAMBKB(3,ICYS,J)='000'
           NUMBKB(3,ICYS,J)=0
           VALBKB(3,ICYS,J)=0.0
         END DO
C
         K=NCLCYS(ICYS)
         DO J=K+1,30
           NAMCOL(3,ICYS,J)='000'
           NUMCOL(3,ICYS,J)=0
           VALCOL(3,ICYS,J)=0.0
         END DO
C
       IF(TYPCYS(ICYS).NE.'BONDED')THEN
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
C
C
C           - FIND CYS(-)
C
         FCOUL=+2.40
         DIS1=4.00
         DIS2=7.00
         DO JCYS=1,NCYS
         IF(TYPCYS(JCYS).NE.'BONDED'.AND.
     $            PKACYS(ICYS).GT.PKACYS(JCYS))THEN
           IF(NMASS(3,ICYS)+NMASS(3,JCYS).GT.900 .OR.
     $    (NMASS(3,ICYS).GT.400.AND.NMASS(3,JCYS).GT.400))THEN
             XS=X(LCYSSG(JCYS))
             YS=Y(LCYSSG(JCYS))
             ZS=Z(LCYSSG(JCYS))
             DIS=SQRT((XSG-XS)**2+(YSG-YS)**2+(ZSG-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='CYS'
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LCYSRS(JCYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
             END IF
           END IF
         END IF
         END DO
C
       END IF
       END DO
C
C
C
C      **********************
C      TYR
C      **********************
C
C
       DO ITYR=1, NTYR
         K=NSDTYR(ITYR)
         DO J=K+1,30
           NAMSDC(4,ITYR,J)='000'
           NUMSDC(4,ITYR,J)=0
           VALSDC(4,ITYR,J)=0.0
         END DO
C
         K=NBKTYR(ITYR)
         DO J=K+1,30
           NAMBKB(4,ITYR,J)='000'
           NUMBKB(4,ITYR,J)=0
           VALBKB(4,ITYR,J)=0.0
         END DO
C
         K=NCLTYR(ITYR)
         DO J=K+1,30
           NAMCOL(4,ITYR,J)='000'
           NUMCOL(4,ITYR,J)=0
           VALCOL(4,ITYR,J)=0.0
         END DO
C
         XOH=X(LTYROH(ITYR))
         YOH=Y(LTYROH(ITYR))
         ZOH=Z(LTYROH(ITYR))
C
C           - FIND TYR-OH
C
         FOO=0.80
         DIS1=3.50
         DIS2=4.50
         IF(ITYR.LT.NTYR)THEN
         DO JTYR=ITYR+1,NTYR
           XO=X(LTYROH(JTYR))
           YO=Y(LTYROH(JTYR))
           ZO=Z(LTYROH(JTYR))
           IF(ABS(XOH-XO).LT.DIS2.AND.
     $        ABS(YOH-YO).LT.DIS2.AND.
     $        ABS(ZOH-ZO).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XO)**2+(YOH-YO)**2+(ZOH-ZO)**2)
             IF(DIS.LT.DIS2)THEN
             IF(PKATYR(ITYR).GT.PKATYR(JTYR))THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='TYR'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LTYRRS(JTYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOO*MIN(1.0,VALUE)
C
               NSDTYR(JTYR)=NSDTYR(JTYR)+1
               NAMSDC(4,JTYR,NSDTYR(JTYR))='TYR'
               NUMSDC(4,JTYR,NSDTYR(JTYR))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,JTYR,NSDTYR(JTYR))=-FOO*MIN(1.0,VALUE)
             ELSE
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='TYR'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LTYRRS(JTYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=-FOO*MIN(1.0,VALUE)
C
               NSDTYR(JTYR)=NSDTYR(JTYR)+1
               NAMSDC(4,JTYR,NSDTYR(JTYR))='TYR'
               NUMSDC(4,JTYR,NSDTYR(JTYR))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,JTYR,NSDTYR(JTYR))=FOO*MIN(1.0,VALUE)
             END IF
             PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             PK1TYR(JTYR)=PK1TYR(JTYR)+VALSDC(4,JTYR,NSDTYR(JTYR))
             END IF
           END IF
         END DO
         END IF
C
C           - FIND TYR(-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO JTYR=1,NTYR
         IF(ITYR.NE.JTYR .AND. PKATYR(ITYR).GT.PKATYR(JTYR))THEN
           IF(NMASS(4,ITYR)+NMASS(4,JTYR).GT.900 .OR.
     $       (NMASS(4,ITYR).GT.400.AND.NMASS(4,JTYR).GT.400))THEN
           XOJ=X(LTYROH(JTYR))
           YOJ=Y(LTYROH(JTYR))
           ZOJ=Z(LTYROH(JTYR))
           IF(ABS(XOH-XOJ).LT.DIS2.AND.
     $        ABS(YOH-YOJ).LT.DIS2.AND.
     $        ABS(ZOH-ZOJ).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XOJ)**2+(YOH-YOJ)**2+(ZOH-ZOJ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))='TYR'
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LTYRRS(JTYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND LYS H-BONDING
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ILYS=1,NLYS
           FOH=-0.80
           DIS2=4.00
           IF(ILYS.EQ.1)FOH=-1.20
           IF(ILYS.EQ.1)DIS2=4.50
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XOH-XN).LT.DIS2.AND.
     $        ABS(YOH-YN).LT.DIS2.AND.
     $        ABS(ZOH-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XN)**2+(YOH-YN)**2+(ZOH-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               IF(PKATYR(ITYR).LT.PKALYS(ILYS))THEN
                 NSDTYR(ITYR)=NSDTYR(ITYR)+1
                 NAMSDC(4,ITYR,NSDTYR(ITYR))='LYS'
                 IF(ILYS.EQ.1)NAMSDC(4,ITYR,NSDTYR(ITYR))='N+ '
                 NUMSDC(4,ITYR,NSDTYR(ITYR))=LLYSRS(ILYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
                 PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
               END IF
             END IF
           END IF
         END DO
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(4,ITYR)+NMASS(5,ILYS).GT.900 .OR.
     $        (NMASS(4,ITYR).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XOH-XN).LT.DIS2.AND.
     $        ABS(YOH-YN).LT.DIS2.AND.
     $        ABS(ZOH-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XN)**2+(YOH-YN)**2+(ZOH-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               IF(PKATYR(ITYR).LT.PKALYS(ILYS))THEN
                 NCLTYR(ITYR)=NCLTYR(ITYR)+1
                 NAMCOL(4,ITYR,NCLTYR(ITYR))='LYS'
                 IF(ILYS.EQ.1)NAMCOL(4,ITYR,NCLTYR(ITYR))='N+ '
                 NUMCOL(4,ITYR,NCLTYR(ITYR))=LLYSRS(ILYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
                 PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
C
                 NCLLYS(ILYS)=NCLLYS(ILYS)+1
                 NAMCOL(5,ILYS,NCLLYS(ILYS))='TYR'
                 NUMCOL(5,ILYS,NCLLYS(ILYS))=LTYRRS(ITYR)
                 VALCOL(5,ILYS,NCLLYS(ILYS))=-FCOUL*MIN(1.0,VALUE)
                 PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
               END IF
             END IF
           END IF
           END IF
         END DO
C
       END DO
C
C
C      **********************
C      LYS
C      **********************
C
C
       DO ILYS=1, NLYS
         K=NSDLYS(ILYS)
         DO J=K+1,30
           NAMSDC(5,ILYS,J)='000'
           NUMSDC(5,ILYS,J)=0
           VALSDC(5,ILYS,J)=0.0
         END DO
C
         K=NBKLYS(ILYS)
         DO J=K+1,30
           NAMBKB(5,ILYS,J)='000'
           NUMBKB(5,ILYS,J)=0
           VALBKB(5,ILYS,J)=0.0
         END DO
C
         K=NCLLYS(ILYS)
         DO J=K+1,30
           NAMCOL(5,ILYS,J)='000'
           NUMCOL(5,ILYS,J)=0
           VALCOL(5,ILYS,J)=0.0
         END DO
C
         XNZ=X(LLYSNZ(ILYS))
         YNZ=Y(LLYSNZ(ILYS))
         ZNZ=Z(LLYSNZ(ILYS))
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO JLYS=1,NLYS
         IF(ILYS.NE.JLYS .AND. PKALYS(ILYS).LT.PKALYS(JLYS))THEN
           IF(NMASS(5,ILYS)+NMASS(5,JLYS).GT.900 .OR.
     $     (NMASS(5,ILYS).GT.400.AND.NMASS(5,JLYS).GT.400))THEN
           XN=X(LLYSNZ(JLYS))
           YN=Y(LLYSNZ(JLYS))
           ZN=Z(LLYSNZ(JLYS))
           IF(ABS(XNZ-XN).LT.DIS2.AND.
     $        ABS(YNZ-YN).LT.DIS2.AND.
     $        ABS(ZNZ-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XNZ-XN)**2+(YNZ-YN)**2+(ZNZ-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLYS(ILYS)=NCLLYS(ILYS)+1
               NAMCOL(5,ILYS,NCLLYS(ILYS))='LYS'
               IF(JLYS.EQ.1)NAMCOL(5,ILYS,NCLLYS(ILYS))='N+ '
               NUMCOL(5,ILYS,NCLLYS(ILYS))=LLYSRS(JLYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(5,ILYS,NCLLYS(ILYS))=FCOUL*MIN(1.0,VALUE)
               PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
         IF(NMASS(5,ILYS)+NMASS(6,IARG).GT.900 .OR.
     $   (NMASS(5,ILYS).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XNZ-XCZ).LT.DIS2.AND.
     $        ABS(YNZ-YCZ).LT.DIS2.AND.
     $        ABS(ZNZ-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XNZ-XCZ)**2+(YNZ-YCZ)**2+(ZNZ-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               IF(PKALYS(ILYS).LT.PKAARG(IARG))THEN
                 NCLLYS(ILYS)=NCLLYS(ILYS)+1
                 NAMCOL(5,ILYS,NCLLYS(ILYS))='ARG'
                 NUMCOL(5,ILYS,NCLLYS(ILYS))=LARGRS(IARG)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(5,ILYS,NCLLYS(ILYS))=FCOUL*MIN(1.0,VALUE)
                 PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
               ELSE
                 NCLARG(IARG)=NCLARG(IARG)+1
                 NAMCOL(6,IARG,NCLARG(IARG))='LYS'
                 IF(ILYS.EQ.1)NAMCOL(6,IARG,NCLARG(IARG))='N+ '
                 NUMCOL(6,IARG,NCLARG(IARG))=LLYSRS(ILYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(6,IARG,NCLARG(IARG))=FCOUL*MIN(1.0,VALUE)
                 PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
               END IF
             END IF
           END IF
         END IF
         END DO
C
       END DO
C
C
C      **********************
C      ARG
C      **********************
C
C
       DO IARG=1, NARG
         K=NSDARG(IARG)
         DO J=K+1,30
           NAMSDC(6,IARG,J)='000'
           NUMSDC(6,IARG,J)=0
           VALSDC(6,IARG,J)=0.0
         END DO
C
         K=NBKARG(IARG)
         DO J=K+1,30
           NAMBKB(6,IARG,J)='000'
           NUMBKB(6,IARG,J)=0
           VALBKB(6,IARG,J)=0.0
         END DO
C
         K=NCLARG(IARG)
         DO J=K+1,30
           NAMCOL(6,IARG,J)='000'
           NUMCOL(6,IARG,J)=0
           VALCOL(6,IARG,J)=0.0
         END DO
C
         X1=X(LARGN1(IARG))
         Y1=Y(LARGN1(IARG))
         Z1=Z(LARGN1(IARG))
         X2=X(LARGN2(IARG))
         Y2=Y(LARGN2(IARG))
         Z2=Z(LARGN2(IARG))
         X3=X(LARGN3(IARG))
         Y3=Y(LARGN3(IARG))
         Z3=Z(LARGN3(IARG))
         XCZ=X(LARGCZ(IARG))
         YCZ=Y(LARGCZ(IARG))
         ZCZ=Z(LARGCZ(IARG))
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO JARG=1,NARG
         IF(IARG.NE.JARG.AND.PKAARG(IARG).LT.PKAARG(JARG))THEN
           IF(NMASS(6,IARG)+NMASS(6,JARG).GT.900 .OR.
     $     (NMASS(6,IARG).GT.400.AND.NMASS(6,JARG).GT.400))THEN
           XCZJ=X(LARGCZ(JARG))
           YCZJ=Y(LARGCZ(JARG))
           ZCZJ=Z(LARGCZ(JARG))
           IF(ABS(XCZ-XCZJ).LT.DIS2.AND.
     $        ABS(YCZ-YCZJ).LT.DIS2.AND.
     $        ABS(ZCZ-ZCZJ).LT.DIS2    )  THEN
             DIS=SQRT((XCZ-XCZJ)**2+(YCZ-YCZJ)**2+(ZCZ-ZCZJ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))='ARG'
               NUMCOL(6,IARG,NCLARG(IARG))=LARGRS(JARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END IF
         END DO
C
       END DO
C
C
C
       CONV=.TRUE.
       DO I=1,1000
C        IF(PK1CAR(I).NE.PKACAR(I))
C    $ WRITE(11,*)NAMCAR(I),LCARRS(I),PK1CAR(I),PKACAR(I)
C        IF(PK1HIS(I).NE.PKAHIS(I))
C    $ WRITE(11,*)'HIS',LHISRS(I),PK1HIS(I),PKAHIS(I)
C        IF(PK1CYS(I).NE.PKACYS(I))
C    $ WRITE(11,*)'CYS',LCYSRS(I),PK1CYS(I),PKACYS(I)
C        IF(PK1TYR(I).NE.PKATYR(I))
C    $ WRITE(11,*)'TYR',LTYRRS(I),PK1TYR(I),PKATYR(I)
C        IF(PK1LYS(I).NE.PKALYS(I))
C    $ WRITE(11,*)'LYS',LLYSRS(I),PK1LYS(I),PKALYS(I)
C        IF(PK1ARG(I).NE.PKAARG(I))
C    $ WRITE(11,*)'ARG',LARGRS(I),PK1ARG(I),PKAARG(I)
         IF(PK1CAR(I).NE.PKACAR(I))CONV=.FALSE.
         IF(PK1HIS(I).NE.PKAHIS(I))CONV=.FALSE.
         IF(PK1CYS(I).NE.PKACYS(I))CONV=.FALSE.
         IF(PK1TYR(I).NE.PKATYR(I))CONV=.FALSE.
         IF(PK1LYS(I).NE.PKALYS(I))CONV=.FALSE.
         IF(PK1ARG(I).NE.PKAARG(I))CONV=.FALSE.
         PKACAR(I)=PK1CAR(I)
         PKAHIS(I)=PK1HIS(I)
         PKACYS(I)=PK1CYS(I)
         PKATYR(I)=PK1TYR(I)
         PKALYS(I)=PK1LYS(I)
         PKAARG(I)=PK1ARG(I)
       END DO
C
       IF(CONV)GOTO 510
 500   CONTINUE
 510   CONTINUE
C
C
C      **********************
C      STEP 9. PRINT RESULTS
C      **********************
C
       DO ICAR=1, NCAR
          PKACAR(ICAR)=PKACAR(ICAR)+6.00*NHBCAR(ICAR)
       END DO
       DO IHIS=1, NHIS
          PKAHIS(IHIS)=PKAHIS(IHIS)-6.00*NHBHIS(IHIS)
       END DO
C
C
C
       WRITE(11,*)' '
       WRITE(11,'(95(1H-))')
       WRITE(11,'(2(1H-),30X,A,30X,2(1H-))')
     $          'PROPKA: A PROTEIN PKA PREDICTOR'
       WRITE(11,'(2(1H-),28X,A,28X,2(1H-))')
     $          'VERSION 1.0,  04/25/2004, IOWA CITY'
       WRITE(11,'(2(1H-),41X,A,41X,2(1H-))')
     $          'BY HUI LI'
       WRITE(11,'(2(1H-),25X,41X,25X,2(1H-))')
       WRITE(11,'(2(1H-),22X,A,22X,2(1H-))')
     $          'PROPKA PREDICTS PROTEIN PKA VALUES ACCORDING TO'
       WRITE(11,'(2(1H-),30X,A,30X,2(1H-))')
     $          'THE EMPIRICAL RULES PROPOSED BY'
       WRITE(11,'(2(1H-),23X,A,23X,2(1H-))')
     $          'HUI LI, ANDREW D. ROBERTSON AND JAN H. JENSEN'
       WRITE(11,'(95(1H-))')
       WRITE(11,*)' '
       WRITE(11,'(7(1H-),3X,5(1H-),3X,6(1H-),3X,9(1H-),
     $        11(1H-),3(3X,13(1H-)))')
       WRITE(11,'(15X,12X,A,3(3X,A13))')
     $   'DESOLVATION  EFFECTS',
     $   '    SIDECHAIN',
     $   '     BACKBONE',
     $   '    COULOMBIC'
       WRITE(11,'(A15,6(3X,A))')
     $   'RESIDUE     pKa ',
     $   'LOCATE',
     $   '  MASSIVE',
     $   '   LOCAL',
     $   'HYDROGEN BOND',
     $   'HYDROGEN BOND',
     $   '  INTERACTION'
       WRITE(11,'(7(1H-),3X,5(1H-),3X,6(1H-),3X,9(1H-),
     $        3X,8(1H-),3(3X,13(1H-)))')
       WRITE(11,*)' '
C
C
       DO ICAR=1,NCAR
         WRITE(11,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   NAMCAR(ICAR), LCARRS(ICAR), PKACAR(ICAR),
     $   TYPCAR(ICAR),TOLMAS(1,ICAR),NMASS(1,ICAR),
     $   TOLLOC(1,ICAR),NLOCAL(1,ICAR),
     $   VALSDC(1,ICAR,1),NAMSDC(1,ICAR,1),NUMSDC(1,ICAR,1),
     $   VALBKB(1,ICAR,1),NAMBKB(1,ICAR,1),NUMBKB(1,ICAR,1),
     $   VALCOL(1,ICAR,1),NAMCOL(1,ICAR,1),NUMCOL(1,ICAR,1)
         DO J=2,30
           IF(J.LE.NSDCAR(ICAR) .OR.
     $        J.LE.NBKCAR(ICAR) .OR.
     $        J.LE.NCLCAR(ICAR)    )THEN
             WRITE(11,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       NAMCAR(ICAR), LCARRS(ICAR),
     $       VALSDC(1,ICAR,J),NAMSDC(1,ICAR,J),NUMSDC(1,ICAR,J),
     $       VALBKB(1,ICAR,J),NAMBKB(1,ICAR,J),NUMBKB(1,ICAR,J),
     $       VALCOL(1,ICAR,J),NAMCOL(1,ICAR,J),NUMCOL(1,ICAR,J)
           END IF
         END DO
         WRITE(11,*)' '
       END DO
C
C
       DO IHIS=1,NHIS
         WRITE(11,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'HIS', LHISRS(IHIS), PKAHIS(IHIS),
     $   TYPHIS(IHIS),TOLMAS(2,IHIS),NMASS(2,IHIS),
     $   TOLLOC(2,IHIS),NLOCAL(2,IHIS),
     $   VALSDC(2,IHIS,1),NAMSDC(2,IHIS,1),NUMSDC(2,IHIS,1),
     $   VALBKB(2,IHIS,1),NAMBKB(2,IHIS,1),NUMBKB(2,IHIS,1),
     $   VALCOL(2,IHIS,1),NAMCOL(2,IHIS,1),NUMCOL(2,IHIS,1)
         DO J=2,30
           IF(J.LE.NSDHIS(IHIS) .OR.
     $        J.LE.NBKHIS(IHIS) .OR.
     $        J.LE.NCLHIS(IHIS)    )THEN
             WRITE(11,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'HIS', LHISRS(IHIS),
     $       VALSDC(2,IHIS,J),NAMSDC(2,IHIS,J),NUMSDC(2,IHIS,J),
     $       VALBKB(2,IHIS,J),NAMBKB(2,IHIS,J),NUMBKB(2,IHIS,J),
     $       VALCOL(2,IHIS,J),NAMCOL(2,IHIS,J),NUMCOL(2,IHIS,J)
           END IF
         END DO
         WRITE(11,*)' '
       END DO
C
C
       DO ICYS=1,NCYS
         WRITE(11,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'CYS', LCYSRS(ICYS), PKACYS(ICYS),
     $   TYPCYS(ICYS),TOLMAS(3,ICYS),NMASS(3,ICYS),
     $   TOLLOC(3,ICYS),NLOCAL(3,ICYS),
     $   VALSDC(3,ICYS,1),NAMSDC(3,ICYS,1),NUMSDC(3,ICYS,1),
     $   VALBKB(3,ICYS,1),NAMBKB(3,ICYS,1),NUMBKB(3,ICYS,1),
     $   VALCOL(3,ICYS,1),NAMCOL(3,ICYS,1),NUMCOL(3,ICYS,1)
         DO J=2,30
           IF(J.LE.NSDCYS(ICYS) .OR.
     $        J.LE.NBKCYS(ICYS) .OR.
     $        J.LE.NCLCYS(ICYS)    )THEN
             WRITE(11,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'CYS', LCYSRS(ICYS),
     $       VALSDC(3,ICYS,J),NAMSDC(3,ICYS,J),NUMSDC(3,ICYS,J),
     $       VALBKB(3,ICYS,J),NAMBKB(3,ICYS,J),NUMBKB(3,ICYS,J),
     $       VALCOL(3,ICYS,J),NAMCOL(3,ICYS,J),NUMCOL(3,ICYS,J)
           END IF
         END DO
         WRITE(11,*)' '
       END DO
C
C
       DO ITYR=1,NTYR
         WRITE(11,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'TYR', LTYRRS(ITYR), PKATYR(ITYR),
     $   TYPTYR(ITYR),TOLMAS(4,ITYR),NMASS(4,ITYR),
     $   TOLLOC(4,ITYR),NLOCAL(4,ITYR),
     $   VALSDC(4,ITYR,1),NAMSDC(4,ITYR,1),NUMSDC(4,ITYR,1),
     $   VALBKB(4,ITYR,1),NAMBKB(4,ITYR,1),NUMBKB(4,ITYR,1),
     $   VALCOL(4,ITYR,1),NAMCOL(4,ITYR,1),NUMCOL(4,ITYR,1)
         DO J=2,30
           IF(J.LE.NSDTYR(ITYR) .OR.
     $        J.LE.NBKTYR(ITYR) .OR.
     $        J.LE.NCLTYR(ITYR)    )THEN
             WRITE(11,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'TYR', LTYRRS(ITYR),
     $       VALSDC(4,ITYR,J),NAMSDC(4,ITYR,J),NUMSDC(4,ITYR,J),
     $       VALBKB(4,ITYR,J),NAMBKB(4,ITYR,J),NUMBKB(4,ITYR,J),
     $       VALCOL(4,ITYR,J),NAMCOL(4,ITYR,J),NUMCOL(4,ITYR,J)
           END IF
         END DO
         WRITE(11,*)' '
       END DO
C
C
       DO ILYS=1,NLYS
       IF(ILYS.EQ.1)THEN
         WRITE(11,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'N+ ', LLYSRS(ILYS), PKALYS(ILYS),
     $   TYPLYS(ILYS),TOLMAS(5,ILYS),NMASS(5,ILYS),
     $   TOLLOC(5,ILYS),NLOCAL(5,ILYS),
     $   VALSDC(5,ILYS,1),NAMSDC(5,ILYS,1),NUMSDC(5,ILYS,1),
     $   VALBKB(5,ILYS,1),NAMBKB(5,ILYS,1),NUMBKB(5,ILYS,1),
     $   VALCOL(5,ILYS,1),NAMCOL(5,ILYS,1),NUMCOL(5,ILYS,1)
         DO J=2,30
           IF(J.LE.NSDLYS(ILYS) .OR.
     $        J.LE.NBKLYS(ILYS) .OR.
     $        J.LE.NCLLYS(ILYS)    )THEN
             WRITE(11,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'N+ ', LLYSRS(ILYS),
     $       VALSDC(5,ILYS,J),NAMSDC(5,ILYS,J),NUMSDC(5,ILYS,J),
     $       VALBKB(5,ILYS,J),NAMBKB(5,ILYS,J),NUMBKB(5,ILYS,J),
     $       VALCOL(5,ILYS,J),NAMCOL(5,ILYS,J),NUMCOL(5,ILYS,J)
           END IF
         END DO
         WRITE(11,*)' '
       ELSE
         WRITE(11,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'LYS', LLYSRS(ILYS), PKALYS(ILYS),
     $   TYPLYS(ILYS),TOLMAS(5,ILYS),NMASS(5,ILYS),
     $   TOLLOC(5,ILYS),NLOCAL(5,ILYS),
     $   VALSDC(5,ILYS,1),NAMSDC(5,ILYS,1),NUMSDC(5,ILYS,1),
     $   VALBKB(5,ILYS,1),NAMBKB(5,ILYS,1),NUMBKB(5,ILYS,1),
     $   VALCOL(5,ILYS,1),NAMCOL(5,ILYS,1),NUMCOL(5,ILYS,1)
         DO J=2,30
           IF(J.LE.NSDLYS(ILYS) .OR.
     $        J.LE.NBKLYS(ILYS) .OR.
     $        J.LE.NCLLYS(ILYS)    )THEN
             WRITE(11,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'LYS', LLYSRS(ILYS),
     $       VALSDC(5,ILYS,J),NAMSDC(5,ILYS,J),NUMSDC(5,ILYS,J),
     $       VALBKB(5,ILYS,J),NAMBKB(5,ILYS,J),NUMBKB(5,ILYS,J),
     $       VALCOL(5,ILYS,J),NAMCOL(5,ILYS,J),NUMCOL(5,ILYS,J)
           END IF
         END DO
         WRITE(11,*)' '
       END IF
       END DO
C
C
       DO IARG=1,NARG
         WRITE(11,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'ARG', LARGRS(IARG), PKAARG(IARG),
     $   TYPARG(IARG),TOLMAS(6,IARG),NMASS(6,IARG),
     $   TOLLOC(6,IARG),NLOCAL(6,IARG),
     $   VALSDC(6,IARG,1),NAMSDC(6,IARG,1),NUMSDC(6,IARG,1),
     $   VALBKB(6,IARG,1),NAMBKB(6,IARG,1),NUMBKB(6,IARG,1),
     $   VALCOL(6,IARG,1),NAMCOL(6,IARG,1),NUMCOL(6,IARG,1)
         DO J=2,30
           IF(J.LE.NSDARG(IARG) .OR.
     $        J.LE.NBKARG(IARG) .OR.
     $        J.LE.NCLARG(IARG)    )THEN
             WRITE(11,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'ARG', LARGRS(IARG),
     $       VALSDC(6,IARG,J),NAMSDC(6,IARG,J),NUMSDC(6,IARG,J),
     $       VALBKB(6,IARG,J),NAMBKB(6,IARG,J),NUMBKB(6,IARG,J),
     $       VALCOL(6,IARG,J),NAMCOL(6,IARG,J),NUMCOL(6,IARG,J)
           END IF
         END DO
         WRITE(11,*)' '
       END DO
C
C
       WRITE(11,'(95(1H-))')
       WRITE(11,'(A)')'SUMMARY OF THIS PREDICTION'
       DO ICAR=1,NCAR
         IF(NAMCAR(ICAR).EQ.'ASP')THEN
           WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $     NAMCAR(ICAR),LCARRS(ICAR),SPACE1,LCARCH(ICAR),PKACAR(ICAR)
         END IF
       END DO
       DO ICAR=1,NCAR
         IF(NAMCAR(ICAR).EQ.'GLU')THEN
           WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $     NAMCAR(ICAR), LCARRS(ICAR),SPACE1,LCARCH(ICAR),PKACAR(ICAR)
         END IF
       END DO
       DO ICAR=1,NCAR
         IF(NAMCAR(ICAR).EQ.'C- ')THEN
           WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $     NAMCAR(ICAR), LCARRS(ICAR),SPACE1,LCARCH(ICAR),PKACAR(ICAR)
         END IF
       END DO
       DO IHIS=1,NHIS
         WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $   'HIS', LHISRS(IHIS), SPACE1, LHISCH(IHIS), PKAHIS(IHIS)
       END DO
       DO ICYS=1,NCYS
         WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $   'CYS', LCYSRS(ICYS), SPACE1, LCYSCH(ICYS), PKACYS(ICYS)
       END DO
       DO ITYR=1,NTYR
         WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $   'TYR', LTYRRS(ITYR), SPACE1, LTYRCH(ITYR), PKATYR(ITYR)
       END DO
       WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $   'N+ ', LLYSRS(1), SPACE1, LLYSCH(1), PKALYS(1)
       DO ILYS=2,NLYS
         WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $   'LYS', LLYSRS(ILYS), SPACE1, LLYSCH(ILYS), PKALYS(ILYS)
       END DO
       DO IARG=1,NARG
         WRITE(11,'(3X,A3,I5,A1,A1,F8.2)')
     $   'ARG', LARGRS(IARG), SPACE1, LARGCH(IARG), PKAARG(IARG)
       END DO
       WRITE(11,'(95(1H-))')
C
       CLOSE(11)
C
C
       END
