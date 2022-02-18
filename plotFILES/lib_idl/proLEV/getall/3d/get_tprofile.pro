PRO get_tprofile, dir, FILE_EXIST, NDXMAX, NDYMAX, NDZMAX, DIM_MU, DIM_PHI, NXOBS, INDX_XT1, INDX_YT1, INDX_ZT1, PR_ANGDEP, PRBAR_ANGDEP
;
;INPUT: NDXMAX, NDYMAX, NDZMAX, DIM_MU, DIM_PHI, NXOBS
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/tprof_in1.dat'
fname2=dir+'/tprof_t1.dat'
fname3=dir+'/tprof_t2.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname1)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname1
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname2)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname2
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname3)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname3
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
HEADER=''
;
OPENR, 1, fname1
   READF, 1, HEADER
   READF, 1, INDX_XT1, INDX_YT1, INDX_ZT1
CLOSE, 1
;
;
;
IF(INDX_XT1 EQ 0) THEN BEGIN
   PR_ANGDEP=FLTARR(NDXMAX,NXOBS, DIM_MU,DIM_PHI)*1.D0
ENDIF
;
IF(INDX_YT1 EQ 0) THEN BEGIN
   PR_ANGDEP=FLTARR(NDYMAX,NXOBS, DIM_MU,DIM_PHI)*1.D0
ENDIF
;
IF(INDX_ZT1 EQ 0) THEN BEGIN
   PR_ANGDEP=FLTARR(NDZMAX,NXOBS, DIM_MU,DIM_PHI)*1.D0
ENDIF

OPENR, 1, fname2, /F77_UNFORMATTED
   READU, 1, PR_ANGDEP
CLOSE, 1
;
;
;
PRBAR_ANGDEP=FLTARR(NDXMAX,NDYMAX,NDZMAX,DIM_MU,DIM_PHI)*1.D0
OPENR, 1, fname3, /F77_UNFORMATTED
   READU, 1, PRBAR_ANGDEP
CLOSE, 1
;
RETURN
;
END
