!  Copyright(c) 2005 Argonne National Laboratory
!  This header file is used to control the data set assignments in the DIF3D code
!     8 bit INTEGER -128 to 128
!    16 bit INTEGER -32,768 to 32,767
!    32 bit INTEGER -2,147,483,648 to 2,147,483,647
!    64 bit INTEGER -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807
!   128 bit INTEGER why?
!    32 bit REAL 1.17549435E-38 to 3.40282347E38
!    64 bit REAL 2.2250738585072013D-308 to 1.7976931348623158D308
!   128 bit REAL 6.4751751194380251109244389582276465524996Q-4966 to 1.189731495357231765085759326628007016196477Q4932
!
!
! DIF3D wide
! --------------
!#define DIF3D_Debug
#define DIF3D_CallTree
#define DIF3D_US_NonErrorSystem   -2147483648
#define DIF3D_NonErrorSystem      0
#define DIF3D_System8bitcheck
#define DIF3D_R32_Tolerance       1.0d-6
#define DIF3D_R64_Tolerance       1.0d-13
#define DIF3D_R64_LowerCutoff     8.232d-11
#define DIF3D_R64_HigherCutoff    1.304d19
#define DIF3D_PeakToAverage       1.0D-7
#define DIF3D_Definition_PI       3.141592653589793d0
#define DIF3D_Definition_Avogadro 0.602214179d0
#define DIF3D_Convert_To_MB       9.53674316D-7
#define DIF3D_Tchebychev_DR_low   0.1d0
#define DIF3D_Tchebychev_DR_high  1.0d0
#define DIF3D_MaxFileUnits        9999
#define DIF3D_Definition_Yes      '      Yes       '
#define DIF3D_Definition_No       '      No        '
#define DIF3D_BC_Reflected        'REFLECTIVE      '
#define DIF3D_BC_Vacuum           'VOID            '
#define DIF3D_BC_GenericInterface 'GENERICINTERFACE'
#define DIF3D_ULO_Error_Eigenvalue 1.0D-3
#define DIF3D_LLO_Error_Eigenvalue 1.0D-10
#define DIF3D_DO_Error_Eigenvalue  1.0D-6
#define DIF3D_ULO_Error_Fission    1.0D-3
#define DIF3D_LLO_Error_Fission    1.0D-10
#define DIF3D_DO_Error_Fission     5.0D-6
#define DIF3D_ULO_Error_Flux       1.0D-3
#define DIF3D_LLO_Error_Flux       1.0D-10
#define DIF3D_DO_Error_Flux        5.0D-7
!  Variable definitions
#define DIF3D_Int_08bit     INTEGER(KIND=1)
#define DIF3D_Int_16bit     INTEGER(KIND=2)
#define DIF3D_Int_32bit     INTEGER(KIND=4)
#define DIF3D_Int_64bit     INTEGER(KIND=8)
#define DIF3D_Int           INTEGER(KIND=4)
#define DIF3D_Real_32bit    REAL(KIND=4)
#define DIF3D_Real_64bit    REAL(KIND=8)
#define DIF3D_Real          REAL(KIND=8)
#define DIF3D_Log           LOGICAL(KIND=4)
#define DIF3D_FileNames     CHARACTER*128
#define DIF3D_Char          CHARACTER*8
#define DIF3D_Char16        CHARACTER*16
#define DIF3D_Char32        CHARACTER*32
#define DIF3D_Undefined     'UNKNOWN'

! NUBOW

! Processor Angle-Space decomposition style
! #define DIF3D_OLDPROCSETUP

!  The matching byte size definitions
#define DIF3D_Int_08bit_Size  1
#define DIF3D_Int_16bit_Size  2
#define DIF3D_Int_32bit_Size  4
#define DIF3D_Int_64bit_Size  8
#define DIF3D_Int_Size        4
#define DIF3D_Real_32bit_Size 4
#define DIF3D_Real_64bit_Size 8
#define DIF3D_Real_Size       8
#define DIF3D_Log_Size        4
#define DIF3D_FileNames_Size  128
#define DIF3D_Char_Size       8
#define DIF3D_Char32_Size     32
#define DIF3D_MESSAGE_LENGTH  120

! BUILD ZPR Model
#define ZPR_Literal_Void 'VOIDDUDE'

! ISOTXS
! --------------
! #define ISOTXS_Debug
!  Fixed constants
#define ISOTXS_MaxLegendre     11
#define ISOTXS_MaxGroups      1000
!  Variable definitions
#define ISOTXS_Char           CHARACTER*8
#define ISOTXS_Int            INTEGER(KIND=4)
#define ISOTXS_Real           REAL(KIND=4)
!  The matching byte size definitions
#define ISOTXS_Char_Size      8
#define ISOTXS_Int_Size       4
#define ISOTXS_Real_Size      4

! DLAYXS
! --------------
! #define DLAYXS_Debug
!  Variable definitions
#define DLAYXS_Char           CHARACTER*8
#define DLAYXS_Int            INTEGER(KIND=4)
#define DLAYXS_Real           REAL(KIND=4)
!  The matching byte size definitions
#define DLAYXS_Char_Size      8
#define DLAYXS_Int_Size       4
#define DLAYXS_Real_Size      4

! COMPXS
! --------------
! #define COMPXS_Debug
!  Variable definitions
#define COMPXS_Char           CHARACTER*8
#define COMPXS_Int            INTEGER(KIND=4)
#define COMPXS_Real           REAL(KIND=8)
!  The matching byte size definitions
#define COMPXS_Char_Size      8
#define COMPXS_Int_Size       4
#define COMPXS_Real_Size      8

! GEODST
! --------------
! #define GEODST_Debug
!  Variable definitions
#define GEODST_Char           CHARACTER*8
#define GEODST_Int            INTEGER(KIND=4)
#define GEODST_Real           REAL(KIND=8)
#define GEODST_Real_5D        REAL(KIND=4)
!  The matching byte size definitions
#define GEODST_Char_Size      8
#define GEODST_Int_Size       4
#define GEODST_Real_Size      8
#define GEODST_Real_5D_Size   4

! LABELS
! --------------
! #define LABELS_Debug
!  Variable definitions
#define LABELS_Char           CHARACTER*8
#define LABELS_Int            INTEGER(KIND=4)
#define LABELS_Real           REAL(KIND=8)
!  The matching byte size definitions
#define LABELS_Char_Size      8
#define LABELS_Int_Size       4
#define LABELS_Real_Size      8

! FIXSRC
! --------------
! #define FIXSRC_Debug
!  Variable definitions
#define FIXSRC_Char           CHARACTER*8
#define FIXSRC_Int            INTEGER(KIND=4)
#define FIXSRC_Real           REAL(KIND=8)
!  The matching byte size definitions
#define FIXSRC_Char_Size      8
#define FIXSRC_Int_Size       4
#define FIXSRC_Real_Size      8

! NDXSRF
! --------------
! #define NDXSRF_Debug
!  Variable definitions
#define NDXSRF_Char           CHARACTER*8
#define NDXSRF_Int            INTEGER(KIND=4)
#define NDXSRF_Real           REAL(KIND=4)
!  The matching byte size definitions
#define NDXSRF_Char_Size      8
#define NDXSRF_Int_Size       4
#define NDXSRF_Real_Size      4

! PWDINT
! --------------
! #define PWDINT_Debug
!  Variable definitions
#define PWDINT_Char           CHARACTER*8
#define PWDINT_Int            INTEGER(KIND=4)
#define PWDINT_Real           REAL(KIND=4)
!  The matching byte size definitions
#define PWDINT_Char_Size      8
#define PWDINT_Int_Size       4
#define PWDINT_Real_Size      4

! RTFLUX
! --------------
! #define RTFLUX_Debug
!  Variable definitions
#define RTFLUX_Char           CHARACTER*8
#define RTFLUX_Int            INTEGER(KIND=4)
#define RTFLUX_Real           REAL(KIND=8)
!  The matching byte size definitions
#define RTFLUX_Char_Size      8
#define RTFLUX_Int_Size       4
#define RTFLUX_Real_Size      8

! RMFLUX
! --------------
! #define RTFLUX_Debug
!  Variable definitions
#define RMFLUX_Char           CHARACTER*8
#define RMFLUX_Int            INTEGER(KIND=4)
#define RMFLUX_Real           REAL(KIND=8)
!  The matching byte size definitions
#define RMFLUX_Char_Size      8
#define RMFLUX_Int_Size       4
#define RMFLUX_Real_Size      8

! RZMFLX
! --------------
! #define RZMFLX_Debug
!  Variable definitions
#define RZMFLX_Char           CHARACTER*8
#define RZMFLX_Int            INTEGER(KIND=4)
#define RZMFLX_Real           REAL(KIND=8)
!  The matching byte size definitions
#define RZMFLX_Char_Size      8
#define RZMFLX_Int_Size       4
#define RZMFLX_Real_Size      8

! NHFLUX
! --------------
! #define NHFLUX_Debug
!  Variable definitions
#define NHFLUX_Char           CHARACTER*8
#define NHFLUX_Int            INTEGER(KIND=4)
#define NHFLUX_Real4          REAL(KIND=4)
#define NHFLUX_Real           REAL(KIND=8)
!  The matching byte size definitions
#define NHFLUX_Char_Size      8
#define NHFLUX_Int_Size       4
#define NHFLUX_Real4_Size     4
#define NHFLUX_Real_Size      8

! RCTDEN
! --------------
! #define RCTDEN_Debug
!  Variable definitions
#define RCTDEN_Char           CHARACTER*8
#define RCTDEN_Int            INTEGER(KIND=4)
#define RCTDEN_Real4          REAL(KIND=4)
#define RCTDEN_Real           REAL(KIND=8)
!  The matching byte size definitions
#define RCTDEN_Char_Size      8
#define RCTDEN_Int_Size       4
#define RCTDEN_Real4_Size     4
#define RCTDEN_Real_Size      8

! RCTFLX
! --------------
! #define RCTFLX_Debug
!  Variable definitions
#define RCTFLX_Char           CHARACTER*8
#define RCTFLX_Int            INTEGER(KIND=4)
#define RCTFLX_Real4          REAL(KIND=4)
#define RCTFLX_Real           REAL(KIND=8)
!  The matching byte size definitions
#define RCTFLX_Char_Size      8
#define RCTFLX_Int_Size       4
#define RCTFLX_Real4_Size     4
#define RCTFLX_Real_Size      8

! RCTPWD
! --------------
! #define RCTPWD_Debug
!  Variable definitions
#define RCTPWD_Char           CHARACTER*8
#define RCTPWD_Int            INTEGER(KIND=4)
#define RCTPWD_Real4          REAL(KIND=4)
#define RCTPWD_Real           REAL(KIND=8)
!  The matching byte size definitions
#define RCTPWD_Char_Size      8
#define RCTPWD_Int_Size       4
#define RCTPWD_Real4_Size     4
#define RCTPWD_Real_Size      8

! ZNATDN
! --------------
! #define ZNATDN_Debug
!  Variable definitions
#define ZNATDN_Char           CHARACTER*8
#define ZNATDN_Int            INTEGER(KIND=4)
#define ZNATDN_Real           REAL(KIND=4)
!  The matching byte size definitions
#define ZNATDN_Char_Size      8
#define ZNATDN_Int_Size       4
#define ZNATDN_Real_Size      4

! NE_FreeForm
! --------------
! #define NE_FreeForm_Debug
!  Fixed constants
#define NE_FreeForm_MaxInputLength  1024
#define NE_FreeForm_MaxParseLength  128
#define NE_FreeForm_MaxDataPerLine  1024
#define NE_FreeForm_Error_Int       -123456789
#define NE_FreeForm_Error_Real      -1.2345678d9
#define NE_FreeForm_Error_Char      '-FAIL'
!  Variable definitions
#define NE_FreeForm_Log            LOGICAL(KIND=1)
#define NE_FreeForm_InputString    CHARACTER*1024
#define NE_FreeForm_Char           CHARACTER*128
#define NE_FreeForm_Int            INTEGER(KIND=4)
#define NE_FreeForm_Real           REAL(KIND=8)
!  The matching byte size definitions
#define NE_FreeForm_Log_Size           1
#define NE_FreeForm_InputString_Size  16
#define NE_FreeForm_Char_Size          4
#define NE_FreeForm_Int_Size           4
#define NE_FreeForm_Real_Size          8

! Large Memory Array (LMA)
! --------------
!#define NE_LMA_Debug
#define NE_LMA_MaxPartitions      5
#define NE_LMA_GuessPieces        20
#define NE_LMA_MaxPieces          100
#define NE_LMA_MaxRecordSize      1048576
#define NE_LMA_MinRecordSize      65536
#define NE_LMA_MaxAllocationSteps 25
#define NE_LMA_Max_LMA_MB         500000
#define NE_LMA_Char               CHARACTER*16
#define NE_LMA_Char_Size          16
#define NE_LMA_RecordBytesPerWord 2
! The maximum number of 8 byte words in a LMA file
! Note that 1 TB or 1024*1024 is the typical disk size as of 9/2012
#define NE_LMA_MaxFileSize_MB     1024
#define NE_LMA_MaxFiles           1024

! PERSENT
! --------------
#define PERSENT_Invalid_Material -1.23454321d0
#define PERSENT_NumNumerators  7
#define PERSENT_NumPrincipleRx 6
#define PERSENT_MaxModel 8192
#define PERSENT_DefaultSens 1.01d0
#define PERSENT_FatalEigError 0.00010d0
#define PERSENT_WarnEigError  0.000005d0
#define PERSENT_PertMaxToMin  1.0d-9
#define PERSENT_SensThreshold 1.0d-9
#define PERSENT_M01  'FIRST_ORDER_PT'
#define PERSENT_M02  'GENERALIZED_PT'
#define PERSENT_M03  'NS_FIRST_ORDER'
#define PERSENT_T01  'PERTURB_DENSITY'
#define PERSENT_T02  'PERTURB_ZONE'
#define PERSENT_T03  'PERTURB_XS'
#define PERSENT_T04  'LAMBDA_BETA'
#define PERSENT_R01  'REACTION_RATE'
#define PERSENT_R02  'REACTION_RATIO'
#define PERSENT_R03  'REACTION_WORTH'
#define PERSENT_R04  'POWER_FRACTION'
#define PERSENT_R05  'SENS_BETA'
#define PERSENT_R06  'SENS_LAMBDA'
#define PERSENT_R07  'EGPT_WORTH'
#define PERSENT_R08  'SENS_EIGENVALUE'
#define PERSENT_R09  'BILINEAR'
#define PERSENT_R10  'SENS_FILE'
#define PERSENT_Edit01 'PRINT_BY_MESH'
#define PERSENT_Edit02 'PRINT_BY_REGION'
#define PERSENT_Edit03 'PRINT_BY_AREA'
#define PERSENT_Edit04 'PRINT_BALANCE'
#define PERSENT_Edit05 'PRINT_BY_ISOTOPE'
#define PERSENT_Edit06 'PRINT_BY_GROUP'
#define PERSENT_Edit07 'EXPORT_VTK'
#define PERSENT_Edit08 'PRINT_PERTURBATION'
#define PERSENT_Edit09 'PRINT_BY_MASS'
#define PERSENT_Edit10 'PRINT_BY_UNIQUE'
#define PERSENT_Edit11 'PRINT_BY_FAMILY'

#define PERSENT_S_Principles 23
#define PERSENT_xs01 'TOTAL'
#define PERSENT_xs02 'NU'
#define PERSENT_xs03 'NUFISSION'
#define PERSENT_xs04 'FISSION'
#define PERSENT_xs05 'CAPTURE'
#define PERSENT_xs06 'GAMMA'
#define PERSENT_xs07 'ALPHA'
#define PERSENT_xs08 'PROTON'
#define PERSENT_xs09 'TRITIUM'
#define PERSENT_xs10 'DEUTERIUM'
#define PERSENT_xs11 'SCATTER'
#define PERSENT_xs12 'ELASTIC'
#define PERSENT_xs13 'INELASTIC'
#define PERSENT_xs14 'N2N'
#define PERSENT_xs15 'P1SCATTER'
#define PERSENT_xs16 'P1ELASTIC'
#define PERSENT_xs17 'P1INELASTIC'
#define PERSENT_xs18 'P1N2N'
#define PERSENT_xs19 'POWER'
#define PERSENT_xs20 'CHI'
#define PERSENT_xs21 'CHI_FD'
#define PERSENT_xs22 'STANDARDSET'
#define PERSENT_xs23 'EVERYTHING'

! PERSENT
! --------------
#define PINDETAIL_PINMAPTYPE1  'CARTESIAN'
#define PINDETAIL_PINMAPTYPE2  'HEXAGONAL'

! COMMARA
! --------------- Compare with the above
#define COMMARA_xs01 'CHI' 
#define COMMARA_xs02 'P1ELASTIC'
#define COMMARA_xs03 'NU'         
#define COMMARA_xs04 'FISSION'    
#define COMMARA_xs05 'GAMMA'     
#define COMMARA_xs06 'CAPTURE'    
#define COMMARA_xs07 'ELASTIC'    
#define COMMARA_xs08 'INELASTIC'  
#define COMMARA_xs09 'NXN'        
#define COMMARA_xs10 'DISCARD'

! MovingFuel code
#define MovingFuel_XS_PrecNuFis 1.0d-6
#define MovingFuel_XS_Scaling   1.0d-24
#define MovingFuel_DensityLimit 1000.0d0
