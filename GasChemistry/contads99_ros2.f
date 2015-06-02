ccc g77 -fno-silent contads99_ros2.f adjros2.f saprcnov_stem_model.f 
ccc saprcnov_stem_hessian.f saprcnov_stem_linalg.f 
ccc saprcnov_stem_util.f


C                                                                  
C **************************************************************** 



      BLOCK DATA GLOBAL_DATA

      INCLUDE 'saprcnov_stem.h'


      data CFACTOR / 2.4476e+13 / 
      
      data LOOKAT /
     *  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
     * 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
     * 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
     * 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
     * 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
     * 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
     * 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
     * 85, 86, 87, 88, 89, 90, 91, 92, 93 / 


      data SLOOKAT /
     *'H2SO4','BC','OC','SSF','SSC',
     *'PM10','PM25','DMS','DST1','DST2',
     *'DST3','CO2','HCOOH','CCO_OH','RCO_OH',
     *'CCO_OOH','RCO_OOH','XN','XC','SO2',
     *'O1D','C2H6','BACL','PAN','PAN2',
     *'PBZN','MA_PAN','H2O2','C3H8','ETOH',
     *'N2O5','HONO','ALK3','TBU_O','ALK5',
     *'ARO2','HNO4','COOH','HOCOO','BZNO2_O',
     *'MEOH','ALK4','ARO1','DCB3','DCB2',
     *'CRES','C2H2','DCB1','BALD','ROOH',
     *'NPHE','PHEN','CO','MGLY','ACET',
     *'HNO3','ETHENE','C3H6','GLY','BZ_O',
     *'ISOPRENE','R2O2','TERP','METHACRO','OLE1',
     *'ISOPROD','OLE2','MVK','CCHO','HCHO',
     *'RNO3','O3P','RCHO','MEK','PROD2',
     *'O3','BZCO_O2','HO2','OH','MA_RCO3',
     *'NO','NO3','C_O2','CCO_O2','NO2',
     *'RO2_N','RO2_R','RCO_O2','AIR','O2',
     *'H2O','H2','CH4' / 


      data MONITOR /
     * 76 / 


      data SMONITOR /
     *'O3' / 





      data( SEQN(i), i = 1, 24 ) /
     *'                       NO2 --> O3P + NO         ',
     *'            O3P + AIR + O2 --> O3               ',
     *'                  O3P + O3 --> 2O2              ',
     *'            O3P + NO + AIR --> NO2              ',
     *'                 O3P + NO2 --> NO               ',
     *'                 O3P + NO2 --> NO3              ',
     *'                   O3 + NO --> NO2              ',
     *'                  O3 + NO2 --> NO3              ',
     *'                  NO + NO3 --> 2NO2             ',
     *'                  2NO + O2 --> 2NO2             ',
     *'                 NO3 + NO2 --> N2O5             ',
     *'                      N2O5 --> NO3 + NO2        ',
     *'                N2O5 + H2O --> 2HNO3            ',
     *'                 NO3 + NO2 --> NO + NO2         ',
     *'                       NO3 --> NO               ',
     *'                       NO3 --> O3P + NO2        ',
     *'                        O3 --> O3P              ',
     *'                        O3 --> O1D              ',
     *'                 O1D + H2O --> 2OH              ',
     *'                 O1D + AIR --> O3P              ',
     *'                   OH + NO --> HONO             ',
     *'                      HONO --> OH + NO          ',
     *'                      HONO --> HO2 + NO2        ',
     *'                 HONO + OH --> NO2              ' / 

      data( SEQN(i), i = 25, 48 ) /
     *'                  OH + NO2 --> HNO3             ',
     *'                  OH + NO3 --> HO2 + NO2        ',
     *'                 HNO3 + OH --> NO3              ',
     *'                      HNO3 --> OH + NO2         ',
     *'                   CO + OH --> HO2              ',
     *'                   O3 + OH --> HO2              ',
     *'                  HO2 + NO --> OH + NO2         ',
     *'                 HO2 + NO2 --> HNO4             ',
     *'                      HNO4 --> HO2 + NO2        ',
     *'               HNO4 --> 1HO2 + 0OH + 0NO3 + 1NO2',
     *'                 HNO4 + OH --> NO2              ',
     *'                  O3 + HO2 --> OH               ',
     *'                      2HO2 --> H2O2             ',
     *'                2HO2 + H2O --> H2O2             ',
     *'                HO2 + NO3 --> 0HNO3 + 1OH + 1NO2',
     *'                      2NO3 --> 2NO2             ',
     *'                      H2O2 --> 2OH              ',
     *'                 H2O2 + OH --> HO2              ',
     *'                  HO2 + OH --> O2 + H2O         ',
     *'                  SO2 + OH --> H2SO4 + HO2      ',
     *'                   OH + H2 --> HO2              ',
     *'                 NO + C_O2 --> HCHO + HO2 + NO2 ',
     *'                HO2 + C_O2 --> COOH             ',
     *'                NO3 + C_O2 --> HCHO + HO2 + NO2 ' / 

      data( SEQN(i), i = 49, 72 ) /
     *'                     2C_O2 --> MEOH + HCHO      ',
     *'                     2C_O2 --> 2HCHO + 2HO2     ',
     *'                NO + RO2_R --> HO2 + NO2        ',
     *'               HO2 + RO2_R --> ROOH             ',
     *'               NO3 + RO2_R --> HO2 + NO2        ',
     *'            C_O2 + RO2_R --> 0MEOH + 1HCHO + HO2',
     *'                    2RO2_R --> HO2              ',
     *'                 R2O2 + NO --> NO2              ',
     *'                R2O2 + HO2 --> HO2              ',
     *'                R2O2 + NO3 --> NO2              ',
     *'               R2O2 + C_O2 --> C_O2             ',
     *'              R2O2 + RO2_R --> RO2_R            ',
     *'                     2R2O2 --> 2R2O2            ',
     *'                NO + RO2_N --> RNO3             ',
     *'               HO2 + RO2_N --> ROOH             ',
     *'     C_O2 + RO2_N --> 0MEOH + 1HCHO + 0MEK + 0PROD2 + ',
     *'               NO3 + RO2_N --> MEK + HO2 + NO2  ',
     *'           RO2_N + RO2_R --> 0MEK + 0PROD2 + HO2',
     *'              R2O2 + RO2_N --> RO2_N            ',
     *'                    2RO2_N --> MEK + PROD2 + HO2',
     *'              CCO_O2 + NO2 --> PAN              ',
     *'                       PAN --> CCO_O2 + NO2     ',
     *'               NO + CCO_O2 --> C_O2 + NO2       ',
     *'       HO2 + CCO_O2 --> 0CCO_OH + 1CCO_OOH + 0O3' / 

      data( SEQN(i), i = 73, 96 ) /
     *'              NO3 + CCO_O2 --> C_O2 + NO2       ',
     *'             C_O2 + CCO_O2 --> CCO_OH + HCHO    ',
     *'            CCO_O2 + RO2_R --> CCO_OH           ',
     *'             R2O2 + CCO_O2 --> CCO_O2           ',
     *'            CCO_O2 + RO2_N --> CCO_OH + PROD2   ',
     *'                   2CCO_O2 --> 2C_O2            ',
     *'              NO2 + RCO_O2 --> PAN2             ',
     *'                      PAN2 --> NO2 + RCO_O2     ',
     *'              NO + RCO_O2 --> CCHO + NO2 + RO2_R',
     *'       HO2 + RCO_O2 --> 0RCO_OH + 1RCO_OOH + 0O3',
     *'             NO3 + RCO_O2 --> CCHO + NO2 + RO2_R',
     *'             C_O2 + RCO_O2 --> RCO_OH + HCHO    ',
     *'            RO2_R + RCO_O2 --> RCO_OH           ',
     *'             R2O2 + RCO_O2 --> RCO_O2           ',
     *'            RO2_N + RCO_O2 --> RCO_OH + PROD2   ',
     *'         CCO_O2 + RCO_O2 --> CCHO + C_O2 + RO2_R',
     *'                   2RCO_O2 --> 2CCHO + 2RO2_R   ',
     *'             BZCO_O2 + NO2 --> PBZN             ',
     *'                      PBZN --> BZCO_O2 + NO2    ',
     *'              BZCO_O2 + NO --> BZ_O + R2O2 + NO2',
     *'      BZCO_O2 + HO2 --> 0RCO_OH + 1RCO_OOH + 0O3',
     *'             BZCO_O2 + NO3 --> BZ_O + R2O2 + NO2',
     *'            BZCO_O2 + C_O2 --> RCO_OH + HCHO    ',
     *'           BZCO_O2 + RO2_R --> RCO_OH           ' / 

      data( SEQN(i), i = 97, 120 ) /
     *'            R2O2 + BZCO_O2 --> BZCO_O2          ',
     *'           BZCO_O2 + RO2_N --> RCO_OH + PROD2   ',
     *'         BZCO_O2 + CCO_O2 --> BZ_O + R2O2 + C_O2',
     *' BZCO_O2 + RCO_O2 --> BZ_O + R2O2 + CCHO + RO2_R',
     *'                  2BZCO_O2 --> 2BZ_O + 2R2O2    ',
     *'             MA_RCO3 + NO2 --> MA_PAN           ',
     *'                    MA_PAN --> MA_RCO3 + NO2    ',
     *'            MA_RCO3 + NO --> HCHO + CCO_O2 + NO2',
     *'      HO2 + MA_RCO3 --> 0RCO_OH + 1RCO_OOH + 0O3',
     *'           MA_RCO3 + NO3 --> HCHO + CCO_O2 + NO2',
     *'            MA_RCO3 + C_O2 --> RCO_OH + HCHO    ',
     *'           MA_RCO3 + RO2_R --> RCO_OH           ',
     *'            R2O2 + MA_RCO3 --> MA_RCO3          ',
     *'           MA_RCO3 + RO2_N --> 2RCO_OH          ',
     *'       MA_RCO3 + CCO_O2 --> HCHO + C_O2 + CCO_O2',
     *' MA_RCO3 + RCO_O2 --> CCHO + HCHO + CCO_O2 + RO2_R',
     *'BZCO_O2 + MA_RCO3 --> BZ_O + R2O2 + HCHO + CCO_O2',
     *'                  2MA_RCO3 --> 2HCHO + 2CCO_O2  ',
     *'               TBU_O + NO2 --> RNO3             ',
     *'                     TBU_O --> ACET + C_O2      ',
     *'                BZ_O + NO2 --> NPHE             ',
     *'                BZ_O + HO2 --> PHEN             ',
     *'                      BZ_O --> PHEN             ',
     *'             BZNO2_O + NO2 --> 2XN + 6XC        ' / 

      data( SEQN(i), i = 121, 144 ) /
     *'             BZNO2_O + HO2 --> NPHE             ',
     *'                   BZNO2_O --> NPHE             ',
     *'                      HCHO --> CO + 2HO2        ',
     *'                      HCHO --> CO               ',
     *'                 HCHO + OH --> CO + HO2         ',
     *'                HCHO + HO2 --> HOCOO            ',
     *'                     HOCOO --> HCHO + HO2       ',
     *'                HOCOO + NO --> HCOOH + HO2 + NO2',
     *'                HCHO + NO3 --> CO + HNO3 + HO2  ',
     *'                 CCHO + OH --> CCO_O2           ',
     *'                      CCHO --> CO + HO2 + C_O2  ',
     *'                CCHO + NO3 --> HNO3 + CCO_O2    ',
     *'        RCHO + OH --> 0CO + 0CCHO + 0RO2_N + 0RO2_R + ',
     *'                RCHO --> CO + CCHO + HO2 + RO2_R',
     *'                RCHO + NO3 --> HNO3 + RCO_O2    ',
     *'              ACET + OH --> R2O2 + HCHO + CCO_O2',
     *'                      ACET --> C_O2 + CCO_O2    ',
     *'         MEK + OH --> 1R2O2 + 0CCHO + 0HCHO + 0RCHO + ',
     *'                   MEK --> CCHO + CCO_O2 + RO2_R',
     *'                 MEOH + OH --> HCHO + HO2       ',
     *'        ETOH + OH --> 1CCHO + 0HCHO + 1HO2 + 0RO2_R',
     *'               COOH + OH --> 0HCHO + 0OH + 1C_O2',
     *'                      COOH --> HCHO + HO2 + OH  ',
     *'               ROOH + OH --> RCHO + 1OH + 0RO2_R' / 

      data( SEQN(i), i = 145, 168 ) /
     *'                      ROOH --> RCHO + HO2 + OH  ',
     *'                       GLY --> 2CO + 2HO2       ',
     *'                       GLY --> CO + HCHO        ',
     *'               GLY + OH --> 1CO + 1HO2 + 0RCO_O2',
     *'        GLY + NO3 --> 1CO + HNO3 + 1HO2 + 0RCO_O2',
     *'                      MGLY --> CO + HO2 + CCO_O2',
     *'                 MGLY + OH --> CO + CCO_O2      ',
     *'               MGLY + NO3 --> CO + HNO3 + CCO_O2',
     *'                      BACL --> 2CCO_O2          ',
     *'             PHEN + OH --> 0GLY + 0BZ_O + 1RO2_R',
     *'                PHEN + NO3 --> HNO3 + BZ_O      ',
     *'            CRES + OH --> 0MGLY + 0BZ_O + 1RO2_R',
     *'                CRES + NO3 --> HNO3 + BZ_O      ',
     *'                NPHE + NO3 --> BZNO2_O + HNO3   ',
     *'                 BALD + OH --> BZCO_O2          ',
     *'                      BALD --> 7XC              ',
     *'                BALD + NO3 --> HNO3 + BZCO_O2   ',
     *'    METHACRO + OH --> 0CO + 0MGLY + 0HCHO + 0MEK + 0MA',
     *'    METHACRO + O3 --> 0HCOOH + 0CO + 1MGLY + 0HCHO + 0',
     *'   METHACRO + NO3 --> 0CO + 0HNO3 + 0MA_RCO3 + 0RO2_R',
     *'            METHACRO + O3P --> RCHO             ',
     *'         METHACRO --> 1CO + 1HCHO + 0HO2 + 0OH + 0MA_R',
     *'         MVK + OH --> 0MGLY + 1R2O2 + 0HCHO + 1RCHO + ',
     *'         MVK + O3 --> 0HCOOH + 0CO + 1MGLY + 0HCHO + 0' / 

      data( SEQN(i), i = 169, 192 ) /
     *'                 MVK + O3P --> 0RCHO + 1MEK     ',
     *'              MVK --> 1CO + 1PROD2 + 0MA_RCO3 + 0C_O2',
     *'     ISOPROD + OH --> 0CO + 0MGLY + 0GLY + 0CCHO + 0HC',
     *'     ISOPROD + O3 --> 0HCOOH + 0RCO_OH + 0CO + 1MGLY +',
     *'    ISOPROD + NO3 --> 1CO + 0MGLY + 0HNO3 + 0HCHO + 1R',
     *'          ISOPROD --> 1CO + 0CCHO + 0HCHO + 0MEK + 1HO',
     *'       PROD2 + OH --> 0CCHO + 0HCHO + 1RCHO + 0MEK + 0',
     *'            PROD2 --> 1R2O2 + 0CCHO + 1HCHO + 1RCHO + ',
     *'        RNO3 + OH --> 0ACET + 1R2O2 + 0CCHO + 0HCHO + ',
     *'             RNO3 --> 0ACET + 0R2O2 + 0CCHO + 0HCHO + ',
     *'                 DCB1 + OH --> CO + RCHO + RO2_R',
     *'            DCB1 + O3 --> 2CO + GLY + 2HO2 + 0OH',
     *'              DCB2 + OH --> R2O2 + RCHO + CCO_O2',
     *'             DCB2 --> CO + 0MGLY + 0GLY + R2O2 + 0HO2 ',
     *'              DCB3 + OH --> R2O2 + RCHO + CCO_O2',
     *'             DCB3 --> CO + 0MGLY + 0GLY + R2O2 + 0HO2 ',
     *'                  OH + CH4 --> C_O2             ',
     *'           ETHENE + OH --> 0CCHO + 2HCHO + RO2_R',
     *'      ETHENE + O3 --> 0HCOOH + 0CO + HCHO + 0HO2 + 0OH',
     *'              ETHENE + NO3 --> RCHO + RO2_R     ',
     *'     ETHENE + O3P --> 0CO + 0GLY + 0CCHO + 0HCHO + 0HO',
     *'    ISOPRENE + OH --> 0R2O2 + 0METHACRO + 0ISOPROD + 0',
     *'    ISOPRENE + O3 --> 0HCOOH + 0RCO_OH + 0CO + 0R2O2 +',
     *'   ISOPRENE + NO3 --> 0R2O2 + 1ISOPROD + 0NO2 + 0RO2_N' / 

      data( SEQN(i), i = 193, 216 ) /
     *'   ISOPRENE + O3P --> 0R2O2 + 0HCHO + 1PROD2 + 0MA_RCO',
     *'        TERP + OH --> 0R2O2 + 0HCHO + 0RCHO + 0PROD2 +',
     *' 0RO2_R + 0RCO_O2 --> 0HCOOH + 0RCO_OH + 0BACL + 0CO +',
     *'       TERP + NO3 --> 1R2O2 + 0RNO3 + 0RCHO + 0NO2 + 0',
     *'                TERP + O3P --> 0RCHO + 1PROD2   ',
     *'                 C2H6 + OH --> CCHO + RO2_R     ',
     *'        C3H8 + OH --> 1ACET + 0RCHO + 0RO2_N + 1RO2_R',
     *'        C2H2 + OH --> 0HCOOH + 0CO + 1GLY + 0HCHO + 0H',
     *'        ALK3 + OH --> 0TBU_O + 0ACET + 1R2O2 + 0CCHO +',
     *'        ALK4 + OH --> 0CO + 0ACET + 1R2O2 + 0CCHO + 0H',
     *'        ALK5 + OH --> 0ACET + 1R2O2 + 0CCHO + 0HCHO + ',
     *'        ARO1 + OH --> 0DCB3 + 0DCB2 + 0CRES + 0DCB1 + ',
     *'        ARO2 + OH --> 0BACL + 0DCB3 + 0DCB2 + 0CRES + ',
     *'        OLE1 + OH --> 0ACET + 0R2O2 + 0CCHO + 1HCHO + ',
     *'                  --> 0HCOOH + 0CCO_OH + 0RCO_OH + 0CO',
     *'       OLE1 + NO3 --> 0ACET + 0R2O2 + 0CCHO + 1RNO3 + ',
     *'            OLE1 + O3P --> 0RCHO + 0MEK + 0PROD2',
     *'        OLE2 + OH --> 0BALD + 0ACET + 0R2O2 + 0METHACR',
     *'OD2 + 0HO2 + 0OH + 0C_O2 + 0CCO_O2 + 0RO2_N + 0RO2_R +',
     *'       OLE2 + NO3 --> 0BALD + 0ACET + 1R2O2 + 0MVK + 1',
     *'       OLE2 + O3P --> 0CO + 0METHACRO + 0RCHO + 1MEK +',
     *'           C2H2 + O3 --> 0CO2 + 2CO + 2HO2 + 0OH',
     *'        C3H6 + OH --> 0XC + 1CCHO + 1HCHO + 0RO2_N + 1',
     *'        C3H6 + O3 --> 0CO2 + 0HCOOH + 0CCO_OH + 0XC + ' / 

      data( SEQN(i), i = 217, 235 ) /
     *'       C3H6 + NO3 --> XN + 3XC + 0RO2_N + 1RO2_R',
     *'               C3H6 + O3P --> 1XC + 0RCHO + 1MEK',
     *'                       SO2 --> H2SO4            ',
     *'                       HO2 -->                  ',
     *'                       SO2 -->                  ',
     *'                     H2SO4 -->                  ',
     *'                      HNO3 -->                  ',
     *'                      H2O2 -->                  ',
     *'                        BC -->                  ',
     *'                        OC -->                  ',
     *'                       SSF -->                  ',
     *'                       SSC -->                  ',
     *'                      PM10 -->                  ',
     *'                      PM25 -->                  ',
     *'                      DST1 -->                  ',
     *'                      DST2 -->                  ',
     *'                      DST3 -->                  ',
     *'                       DMS -->                  ',
     *'                       CO2 -->                  '
     * / 


C User defined global variables                                    

C End user defined global variables                                

      END


C **************************************************************** 
C                                                                  
C INITVAL - function to update initial concentrations              
C   Parameters :                                                   
C                                                                  
C **************************************************************** 

      SUBROUTINE INITVAL ( )

      INCLUDE 'saprcnov_stem.h'


      INTEGER i
      REAL*8 x

      x = (0.0E0)*CFACTOR
      do i = 1, NVAR
        VAR(i) = x
      end do

      x = (0.0E0)*CFACTOR
      do i = 1, NRAD
        RAD(i) = x
      end do

      x = (0.0E0)*CFACTOR
      do i = 1, NFIX
        FIX(i) = x
      end do

      VAR(13) = (6.77E-4)*CFACTOR
      VAR(14) = (1.16E-3)*CFACTOR
      VAR(15) = (3.92E-4)*CFACTOR
      VAR(20) = (5.E-2)*CFACTOR
      VAR(32) = (1.E-3)*CFACTOR
      VAR(33) = (4.69E-2)*CFACTOR
      VAR(35) = (3.06E-2)*CFACTOR
      VAR(36) = (8.74E-3)*CFACTOR
      VAR(41) = (5.89E-3)*CFACTOR
      VAR(42) = (4.17E-2)*CFACTOR
      VAR(43) = (1.18E-2)*CFACTOR
      VAR(46) = (5.60E-4)*CFACTOR
      VAR(49) = (7.51E-5)*CFACTOR
      VAR(52) = (6.06E-4)*CFACTOR
      VAR(54) = (8.37E-5)*CFACTOR
      VAR(55) = (5.07E-3)*CFACTOR
      VAR(57) = (1.89E-2)*CFACTOR
      VAR(59) = (1.21E-4)*CFACTOR
      VAR(61) = (4.33E-4)*CFACTOR
      VAR(63) = (8.20E-4)*CFACTOR
      VAR(64) = (1.30E-3)*CFACTOR
      VAR(65) = (1.04E-2)*CFACTOR
      VAR(66) = (8.93E-5)*CFACTOR
      VAR(67) = (7.97E-3)*CFACTOR
      VAR(69) = (2.32E-3)*CFACTOR
      VAR(70) = (1.12E-2)*CFACTOR
      VAR(73) = (1.72E-3)*CFACTOR
      VAR(74) = (3.26E-3)*CFACTOR
      VAR(75) = (1.93E-3)*CFACTOR
      VAR(81) = (1.0E-1)*CFACTOR
      VAR(85) = (5.0E-2)*CFACTOR
      FIX(1) = (1.0E6)*CFACTOR
      FIX(2) = (2.09E+5)*CFACTOR
      FIX(3) = (2.0E+04)*CFACTOR
      FIX(4) = 0.0E0      
      Fix(5) = 1.0E0*CFACTOR

      do i = 1, NSPEC
        C_DEFAULT(i) = C(i)
      end do
             
C User defined initialisations                                     

        STEPMIN=1.0d-3
        STEPMAX= 3600.0d0
        ISNOTAUTONOM=0
        STEPSTART=STEPMIN


	TSTART = 0.0D0
	TEND = TSTART + 24.0d0*3600.0d0
	DT = 3600.D0
        TEMP = 300.0D0


C End user defined initialisations                                 

      RETURN
      END

C End of INITVAL function                                          
C **************************************************************** 


      SUBROUTINE Update_SUN()

      INCLUDE 'saprcnov_stem.h'

      REAL*8 SunRise, SunSet
      REAL*8 Thour, Tlocal, Ttmp 
   
      SunRise = 4.5
      SunSet  = 19.5
      Thour = TIME/3600.
      Tlocal = Thour - (INT(Thour)/24)*24

      IF ((Tlocal.GE.SunRise).AND.(Tlocal.LE.SunSet)) THEN
        Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise)
        IF (Ttmp.GT.0) THEN
          Ttmp =  Ttmp*Ttmp
        ELSE
          Ttmp = -Ttmp*Ttmp
        END IF
        SUN = ( 1.0 + COS(PI*Ttmp) )/2.0 
      ELSE
        SUN = 0.0
      END IF

      RETURN
      END

C End of Update_SUN function                                       
C **************************************************************** 


C **************************************************************** 
C                                                                  
C INTEGRATE - Integrator routine                                   
C   Parameters :                                                   
C      TIN       - Start Time for Integration                      
C      TOUT      - End Time for Integration                        
C                                                                  
C **************************************************************** 

      SUBROUTINE INTEGRATE( TIN, TOUT )
 
      INCLUDE 'saprcnov_stem.h'

C TIN - Start Time
      REAL*8 TIN
C TOUT - End Time
      REAL*8 TOUT

      INTEGER    INFO(5)

      EXTERNAL FUN, JAC_SP

      INFO(1) = 1 - ISNOTAUTONOM
 
      call ROS2(NVAR,TIN,TOUT,STEPMIN,STEPMAX,
     +                   STEPSTART,VAR,ATOL,RTOL,
     +                   Info,FUN,JAC_SP)
 

      RETURN
      END
  

 
 
      SUBROUTINE ROS2(N,T,Tnext,Hmin,Hmax,Hstart,
     +                   y,AbsTol,RelTol,
     +                   Info,FUN,JAC_SP)
     
      INCLUDE 'saprcnov_stem_s.h'                                                                                                  
      INCLUDE 'storage.h'
      
C  INPUT ARGUMENTS:
C     y = Vector of (NVAR) concentrations, contains the
C         initial values on input
C     [T, Tnext] = the integration interval
C     Hmin, Hmax = lower and upper bounds for the selected step-size.
C          Note that for Step = Hmin the current computed
C          solution is unconditionally accepted by the error
C          control mechanism.
C     AbsTol, RelTol = (NVAR) dimensional vectors of
C          componentwise absolute and relative tolerances.
C     FUN = name of routine of derivatives. KPP syntax.
C          See the header below.
C     JAC_SP = name of routine that computes the Jacobian, in
C          sparse format. KPP syntax. See the header below.
C     Info(1) = 1  for  autonomous   system
C             = 0  for nonautonomous system
C     Info(2) = 1  for third order embedded formula
C             = 0  for first order embedded formula 
C
C   Note: Stage 3 used to build strongly A-stable order 3 formula for error control
C            Embed3 = (Info(2).EQ.1)
C         if Embed3 = .true. then the third order embedded formula is used
C                     .false. then a first order embedded  formula is used
C
C
C  OUTPUT ARGUMENTS:
C     y = the values of concentrations at Tend.
C     T = equals Tend on output.
C     Info(2) = # of FUN calls.
C     Info(3) = # of JAC_SP calls.
C     Info(4) = # of accepted steps.
C     Info(5) = # of rejected steps.
 
      REAL*8     K1(NVAR), K2(NVAR), K3(NVAR)
      REAL*8     F1(NVAR), JAC(LU_NONZERO_V)
      REAL*8     DFDT(NVAR)
      REAL*8     Hmin,Hmax,Hstart,ghinv,uround      
      REAL*8     y(NVAR), ynew(NVAR)
      REAL*8     AbsTol(NVAR), RelTol(NVAR)
      REAL*8     T, Tnext, H, Hold, Tplus
      REAL*8     ERR, factor, facmax
      INTEGER    n,nfcn,njac,Naccept,Nreject,i,j
      INTEGER    Info(5)
      LOGICAL    IsReject, Autonom, Embed3
      EXTERNAL    FUN, JAC_SP                                                                                                

      REAL*8 gamma, m1, m2, alpha, beta1, beta2, delta, w, e
 
c     Initialization of counters, etc.
      Autonom = Info(1) .EQ. 1
      Embed3  = Info(2) .EQ. 1
      uround  = 1.d-15
      dround  = dsqrt(uround)
      H = DMAX1(1.d-8, Hmin)
      Tplus = T
      IsReject = .false.
      Naccept  = 0
      Nreject  = 0
      Nfcn     = 0
      Njac     = 0
      
C     Method Parameters      
      gamma = 1.d0 + 1.d0/sqrt(2.d0)
      a21   = - 1.d0/gamma
      m1 = -3.d0/(2.d0*gamma)
      m2 = -1.d0/(2.d0*gamma)
      c31 = -1.0D0/gamma**2*(1.0D0-7.0D0*gamma+9.0D0*gamma**2)
     &         /(-1.0D0+2.0D0*gamma)
      c32 = -1.0D0/gamma**2*(1.0D0-6.0D0*gamma+6.0D0*gamma**2)
     &          /(-1.0D0+2.0D0*gamma)/2
      gamma3 = 0.5D0 - 2*gamma
      d1 = ((-9.0D0*gamma+8.0D0*gamma**2+2.0D0)/gamma**2/
     &          (-1.0D0+2*gamma))/6.0D0
      d2 = ((-1.0D0+3.0D0*gam)/gamma**2/
     &          (-1.0D0+2.0D0*gamma))/6.0D0
      d3 = -1.0D0/(3.0D0*gamma)
 
C === Starting the time loop ===
 10    CONTINUE
       Tplus = T + H
       if ( Tplus .gt. Tnext ) then
          H = Tnext - T
          Tplus = Tnext
       end if
 
       call JAC_SP(NVAR, T, y, JAC)
 
       Njac = Njac+1
       ghinv = -1.0d0/(gamma*H)
       DO 20 j=1,NVAR
         JAC(LU_DIAG_V(j)) = JAC(LU_DIAG_V(j)) + ghinv
 20    CONTINUE
       call KppDecomp (NVAR, JAC, ier)
 
       if (ier.ne.0) then
         if ( H.gt.Hmin) then
            H = 5.0d-1*H
            go to 10
         else
            print *,'IER <> 0, H=',H
            stop
         end if
       end if
 
       call FUN(NVAR, T, y, F1)
 
C ====== NONAUTONOMOUS CASE ===============
       IF (.not. Autonom) THEN
         tau = dsign(dround*dmax1( 1.0d-6, dabs(T) ), T)
         call FUN(NVAR, T+tau, y, K2)
         nfcn=nfcn+1
         DO 30 j = 1,NVAR
           DFDT(j) = ( K2(j)-F1(j) )/tau
 30      CONTINUE
       END IF ! .NOT.Autonom
 
C ----- STAGE 1 -----
       delta = gamma*H
       DO 40 j = 1,NVAR
          K1(j) =  F1(j) 
 40    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 45 j = 1,NVAR
           K1(j) = K1(j) + delta*DFDT(j)
 45      CONTINUE       
       END IF ! .NOT.Autonom
       call KppSolve (JAC, K1)
 
C ----- STAGE 2  -----
       DO 50 j = 1,NVAR
         ynew(j) = y(j) + a21*K1(j)
 50    CONTINUE
       call FUN(NVAR, T+H, ynew, F1)
       nfcn=nfcn+1
       beta = 2.d0/(gamma*H)
       delta = -gamma*H
       DO 55 j = 1,NVAR
         K2(j) = F1(j) + beta*K1(j) 
 55    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 56 j = 1,NVAR
           K2(j) = K2(j) + delta*DFDT(j)
 56      CONTINUE       
       END IF ! .NOT.Autonom
       call KppSolve (JAC, K2)
 
C ----- STAGE 3  -----
       IF (Embed3) THEN
       beta1 = -c31/H
       beta2 = -c32/H
       delta = gamma3*H
       DO 57 j = 1,NVAR
         K3(j) = F1(j) + beta1*K1(j) + beta2*K2(j) 
 57    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 58 j = 1,NVAR
           K3(j) = K3(j) + delta*DFDT(j)
 58      CONTINUE       
       END IF ! .NOT.Autonom
       CALL KppSolve (JAC, K3)
       END IF ! Embed3
 


C ---- The Solution ---
       DO 120 j = 1,NVAR
         ynew(j) = y(j) + m1*K1(j) + m2*K2(j)
 120   CONTINUE
 
 
C ====== Error estimation ========
 
        ERR=0.d0
        DO 130 i=1,NVAR
           w = AbsTol(i) + RelTol(i)*DMAX1(DABS(y(i)),DABS(ynew(i)))
	   IF ( Embed3 ) THEN
	     e = d1*K1(i) + d2*K2(i) + d3*K3(i)
	   ELSE
             e = 1.d0/(2.d0*gamma)*(K1(i)+K2(i))
	   END IF  ! Embed3 
           ERR = ERR + ( e/w )**2
 130    CONTINUE
        ERR = DMAX1( uround, DSQRT( ERR/NVAR ) )
	
C ======= Choose the stepsize ===============================
 
        IF ( Embed3 ) THEN
            elo = 3.0D0 ! estimator local order
	ELSE
	    elo = 2.0D0
	END IF    
        factor = DMAX1(2.0D-1,DMIN1(6.0D0,ERR**(1.0D0/elo)/.9D0))
        Hnew   = DMIN1(Hmax,DMAX1(Hmin, H/factor))
        Hold = H
         
C ======= Rejected/Accepted Step ============================
 
        IF ( (ERR.gt.1).and.(H.gt.Hmin) ) THEN
          IsReject = .true.
	  H = DMIN1(H/10,Hnew)
          Nreject  = Nreject+1
        ELSE
          DO 140 i=1,NVAR
             y(i)  = ynew(i)
 140      CONTINUE
          T = Tplus
	  IF (.NOT.IsReject) THEN
	      H = Hnew   ! Do not increase stepsize if previous step was rejected
	  END IF    
          IsReject = .false.
          Naccept = Naccept+1
	  CALL TRAJISTORE(y,hold)    
        end if        
C ======= End of the time loop ===============================
      if ( T .lt. Tnext ) go to 10 
 
C ======= Output Information =================================
      Info(2) = Nfcn
      Info(3) = Njac
      Info(4) = Naccept
      Info(5) = Nreject      
     
      RETURN
      END
 

 
        SUBROUTINE FUN(N, T, Y, P)
        INCLUDE 'saprcnov_stem.h'
        REAL*8   T, Told
        REAL*8   Y(NVAR), P(NVAR)
        Told = TIME
        TIME = T
cc        CALL UPDATE_SUN()
cc        CALL UPDATE_RCONST()
        CALL FunVar( Y, RAD, FIX, RCONST, P )	
        TIME = Told
        RETURN
        END
 
        SUBROUTINE JAC_SP(N, T, Y, J)
        INCLUDE 'saprcnov_stem.h'
        REAL*8   Told, T
        REAL*8   Y(NVAR), J(LU_NONZERO_V)
        Told = TIME
        TIME = T
cc        CALL UPDATE_SUN()
cc        CALL UPDATE_RCONST()
        CALL JacVar_SP( Y, RAD, FIX, RCONST, J )
        TIME = Told
        RETURN
        END                                                                                                                 





                                                                                                                            
C End of INTEGRATE function                                        
C **************************************************************** 

      PROGRAM driver

      INCLUDE 'saprcnov_stem.h'
      INCLUDE 'storage.h'
      
      real*8 advar(nvar),fdvar(nvar),cadvar(nvar)
      real*8 t, pert,h
      real*8 tendeq
      EXTERNAL   FUN, JAC_SP                                                                                                
      INTEGER    Info(5)
      
      do i=1,NVAR
        RTOL(i) = 1.0d-3
        ATOL(i) = 1.0d0
      end do 
                 
      CALL INITVAL()  
      
ccc run for equilibrium 
      
      time = tstart 
      tendeq = tstart + 24.0d0*3600.0d0*1.0d0    
      istore = 0
      do while (time .lt. TENDEQ)              
        CALL UPDATE_SUN() 
        CALL UPDATE_RCONST()	    
        CALL INTEGRATE( TIME, TIME+DT )        	
      end do   
ccc                

 100  FORMAT(F32.16)                     
    
        
      CALL InitSaveData ()        
      time = TSTART                        	      
      
      iperiod = 1 ! counter for integration period
      istore = 0  ! don't store trajectory inside dt period
      
      call TRAJPSTORE() ! save initial state

      do while (time .lt. TEND)
        iperiod = iperiod + 1
	istep = 1 ! counter for nr. of steps inside iperiod
        CALL SaveData ()  
        CALL UPDATE_SUN() 
        CALL UPDATE_RCONST()	    
        CALL INTEGRATE( TIME, TIME+DT )
        CALL TRAJPSTORE()	
      end do
            
      iperiodsave = iperiod
      
      CALL SaveData ()   
      CALL CloseSaveData ()
      CALL GenerateMatlab (' ')
               
      write(6,*)'end forward integration'
    	
ccc   open file for nr. of steps statistics - optional    
      open(44, file = 'istep.dat')        
ccc continuous adjoint       
      
      time = tend
      call INITADJ(nvar,advar) ! initialize adjoint variables
ccc   open file for adjoint trajectory storage - optional             
      open(12, file = 'cadjtraj_ros2.dat')
      write(12,999)(time-tstart)/3600.0d0,(advar(i),i=1,nvar)
            
      do while (iperiod .gt. 1)
        iperiod = iperiod - 1
	time = time - dt
		
c restore fwd concentration at iperiod      
        do i=1,nvar
          var(i) = fwdtrajp(iperiod,i)	 
        enddo

c forward trajectory recomputation and storage
        CALL UPDATE_SUN() 
        CALL UPDATE_RCONST()	
        istore = 1  ! store trajectory inside dt period		
	istep = 1
	do i=1,nvar
	   fwdtraji(istep,i) = fwdtrajp(iperiod,i)
	enddo				
        CALL INTEGRATE( TIME, TIME+DT )	
	write(44,902)istep-1	     
c continuous adjoint backward integration	
	do while (istep .gt. 1)	            
c set the time step of the backward integration       
          h= stepvect(istep-1)	  
c restore forward variables
          do i=1,nvar
            var(i) = fwdtraji(istep,i)	 
          enddo	  
c backward integration 
          call CADINTEGRATE(advar,h,time)
	  time = time - h
          istep = istep - 1
	  write(12,999)(time-tstart)/3600.0d0,(advar(i),i=1,nvar)
        end do  ! istep loop  
ccc        write(12,999)(time-tstart)/3600.0d0,(advar(i),i=1,nvar) 
      end do ! iperiod loop, end backward integration
      close(12)
      close(44)
cccc end continuous adjoint
      write(6,*)'end adjoint integration'
      
ccc  finite differences check 

!!!      go to 1000

cccc begin finite differences check 
      pert = 1.d-6
      open(15, file = 'cadjsensit_ros2.dat') 
      write(6,*)'finite differences check, eps = ',pert,' ppm' 
      write(6,900)       
      istore = 0     
      do k=1,nvar
       do i=1,nvar
          var(i) = fwdtrajp(1,i)
       enddo     
       var(k) = var(k) + pert*cfactor
       time = TSTART                
       do while (time .lt. TEND)
         CALL UPDATE_SUN() 
         CALL UPDATE_RCONST()       
         CALL INTEGRATE( TIME, TIME+DT )		
       end do 
       fdvar(k)=(var(i_o3)-fwdtrajp(iperiodsave,i_o3))/(pert*cfactor)	        
       write(6,901) k, SLOOKAT(k), fdvar(k), advar(k),
     & dabs(advar(k)-fdvar(k))
       write(15,*)fdvar(k), advar(k), dabs(advar(k)-fdvar(k))
      enddo  
      close(15)         
cccc end finite differences check            
      
1000  continue                 	       
      
900   FORMAT('KPPNo.', 2X, 'chem spec',7X,'fd sensit',14X, 
     &       'dadj sensit', 10X, 'abs(adj-fd)')
901   FORMAT(I3,5X, A8, 6X, E16.8, 7X, E16.8, 5X, E16.8)  
902   FORMAT(I3)   
999   FORMAT(E24.16,100(1X,E24.16))
      stop
      end

      
      SUBROUTINE INITADJ(n,advar)
      INCLUDE 'saprcnov_stem.h'      
      
      integer n
      real*8 advar(n)
      
      do i=1,n
        advar(i) = 0.0d0
      enddo
      
      advar(I_O3) = 1.0d0
      
      return
      end
      
      SUBROUTINE TRAJPSTORE()
      INCLUDE 'saprcnov_stem.h'      
      INCLUDE 'storage.h'      
      
      integer i
      
      do i=1,nvar
       fwdtrajp(iperiod,i) = var(i) 
      enddo 
      
      RETURN
      END
      
      SUBROUTINE TRAJISTORE(y, h)
      INCLUDE 'saprcnov_stem.h'      
      INCLUDE 'storage.h'      
      
      real*8 h,y(nvar)
      integer  i
      
      if (istore .ne. 1) return	
          
      stepvect(istep) = h
      istep = istep+1
      if(istep .gt. nsmax) then
	 write(6,*) 'max nr of steps istep =',istep, 
     &              ' at time = ',time
         stop
      end if	          
      do i=1,nvar
         fwdtraji(istep,i) = y(i) 
      enddo       
      
      RETURN
      END      
                   
      SUBROUTINE CADINTEGRATE(ady,h,t)
      	
      INCLUDE 'saprcnov_stem.h' 
      REAL*8 ady(nvar)
      REAL*8 h,t
      
      call cadj_ros2(var,ady,h,t) 
      
      return
      end
       
