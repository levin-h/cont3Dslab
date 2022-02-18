PRO getOPAC3D, dir, FILE_EXIST, NDXMAX, NDYMAX, NDZMAX, OPATH3D
;
fname=dir+'/MODEL_OPATH3D.dat'
;
OPATH3D=FLTARR(NDXMAX, NDYMAX, NDZMAX)*1.D0
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, OPATH3D
CLOSE, 1
;
RETURN
;
END