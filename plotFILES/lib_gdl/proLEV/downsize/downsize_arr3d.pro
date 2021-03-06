pro downsize_arr3d, x, y, z, arr3d, ndxmax, ndymax, ndzmax, arr2=arr2, arr3=arr3, arr4=arr4, arr5=arr5, level=level, help=print_help
;
;+
; NAME:
;	downsize_arr3d
;
; PURPOSE:
;	This procedure calculates a new 3d-array out of a given 3d-array
;       but with smaller dimensions (defined by level)
;
; CALLING SEQUENCE:
;	downsize_arr3d, x, y, z, arr3d, ndxmax, ndymax, ndzmax
;
; INPUTS:
;	arr3d:	input 3d-array
;       x:      x-coordinate of the input array
;       y:      y-coordinate of the input array
;       z:      z-coordinate of the input array
;       ndxmax:   dimension of x
;       ndymax:   dimension of y
;       ndzmax:   dimension of z
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this
;               procedure
;       arr2:   Set this keyword to an additional 3d array that shall be
;               converted the same way
;       arr3:   Set this keyword to an additional 3d array that shall be
;               converted the same way
;       arr4:   Set this keyword to an additional 3d array that shall be
;               converted the same way
;       arr5:   Set this keyword to an additional 3d array that shall be
;               converted the same way
;       level:  Set this keyword to an integer value, to specified that each 
;               level-th grid-point shall be taken into account (default: 2)
;
; OUTPUTS:
;	arr3d:	new 3d-array
;       x:      x-coordinate of the new array
;       y:      y-coordinate of the new array
;       z:      z-coordinate of the new array
;       ndxmax:   new dimension of x
;       ndymax:   new dimension of y
;       ndzmax:   new dimension of z

;
; EXAMPLE:
;	downsize_arr3d, x, y, z, arr3d, ndxmax, ndymax, ndzmax
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'downsize_arr3d'
   return
endif
;
if(not keyword_set(level)) then level=2
;
;-------------------------calculate new x array-------------------------
;
;mark those grid points which are included in new array
flagx=indgen(ndxmax)*0
;
;include outer two boundary points
flagx(0)=1
flagx(1)=1
flagx(ndxmax-2)=1
flagx(ndxmax-1)=1
ndxmax_new=4
;
for i=2, ndxmax-3 do begin
   if(i mod level eq 0) then begin
;include every second point
      flagx(i)=1
      ndxmax_new=ndxmax_new+1
   endif else begin
      if(abs(x(i)-1.d0) lt 1.d-8) then begin
;include points 1 and -1, if not already included in every second point
         flagx(i)=1
         ndxmax_new=ndxmax_new+1
      endif else begin
         if(abs(x(i)) lt 1.d-8) then begin
;include origin, if not already included in every second point
            flagx(i)=1
            ndxmax_new=ndxmax_new+1
         endif
      endelse
   endelse
endfor
;
;-------------------------calculate new y array-------------------------
;
;mark those grid points which are included in new array
flagy=indgen(ndymax)*0
;
;include outer two boundary points
flagy(0)=1
flagy(1)=1
flagy(ndymax-2)=1
flagy(ndymax-1)=1
ndymax_new=4
;
for i=2, ndymax-3 do begin
   if(i mod level eq 0) then begin
;include every second point
      flagy(i)=1
      ndymax_new=ndymax_new+1
   endif else begin
      if(abs(y(i)-1.d0) lt 1.d-8) then begin
;include points 1 and -1, if not already included in every second point
         flagy(i)=1
         ndymax_new=ndymax_new+1
      endif else begin
         if(abs(y(i)) lt 1.d-8) then begin
;include origin, if not already included in every second point
            flagy(i)=1
            ndymax_new=ndymax_new+1
         endif
      endelse
   endelse
endfor
;
;-------------------------calculate new z array-------------------------
;
;mark those grid points which are included in new array
flagz=indgen(ndzmax)*0
;
;include outer two boundary points
flagz(0)=1
flagz(1)=1
flagz(ndzmax-2)=1
flagz(ndzmax-1)=1
ndzmax_new=4
;
for i=2, ndzmax-3 do begin
   if(i mod level eq 0) then begin
;include every second point
      flagz(i)=1
      ndzmax_new=ndzmax_new+1
   endif else begin
      if(abs(z(i)-1.d0) lt 1.d-8) then begin
;include points 1 and -1, if not already included in every second point
         flagz(i)=1
         ndzmax_new=ndzmax_new+1
      endif else begin
         if(abs(z(i)) lt 1.d-8) then begin
;include origin, if not already included in every second point
            flagz(i)=1
            ndzmax_new=ndzmax_new+1
         endif
      endelse
   endelse
endfor
;
;-----------------------------store new grid and data-------------------
;
;create new 3d arrays
arr3d_new=fltarr(ndxmax_new, ndymax_new, ndzmax_new)
if(keyword_set(arr2)) then begin
   arr2_new=fltarr(ndxmax_new, ndymax_new, ndzmax_new)
endif
if(keyword_set(arr3)) then begin
   arr3_new=fltarr(ndxmax_new,ndymax_new,ndzmax_new)
endif
if(keyword_set(arr4)) then begin
   arr4_new=fltarr(ndxmax_new,ndymax_new,ndzmax_new)
endif
if(keyword_set(arr5)) then begin
   arr5_new=fltarr(ndxmax_new,ndymax_new,ndzmax_new)
endif
;
;create new axes
x_new=x(where(flagx eq 1))
y_new=y(where(flagy eq 1))
z_new=z(where(flagz eq 1))
;
indx_x=-1
for i=0, ndxmax-1 do begin
   if(flagx(i) eq 1) then begin
      indx_x=indx_x+1
      indx_y=-1
      for j=0, ndymax-1 do begin
         if(flagy(j) eq 1) then begin
            indx_y=indx_y+1
            indx_z=-1
            for k=0, ndzmax-1 do begin
               if(flagz(k) eq 1) then begin
                  indx_z=indx_z+1 
                  arr3d_new(indx_x, indx_y, indx_z) = arr3d(i,j,k)
                  if(keyword_set(arr2)) then begin
                     arr2_new(indx_x, indx_y, indx_z) = arr2(i,j,k)
                  endif
                  if(keyword_set(arr3)) then begin
                     arr3_new(indx_x, indx_y, indx_z) = arr3(i,j,k)
                  endif
                  if(keyword_set(arr4)) then begin
                     arr4_new(indx_x, indx_y, indx_z) = arr4(i,j,k)
                  endif
                  if(keyword_set(arr5)) then begin
                     arr5_new(indx_x, indx_y, indx_z) = arr5(i,j,k)
                  endif
               endif
            endfor
         endif
      endfor
   endif
endfor
;
;overwrite old grids
ndxmax=ndxmax_new
ndymax=ndymax_new
ndzmax=ndzmax_new
x=x_new
y=y_new
z=z_new
arr3d=arr3d_new
if(keyword_set(arr2)) then arr2=arr2_new
if(keyword_set(arr3)) then arr3=arr3_new
if(keyword_set(arr4)) then arr4=arr4_new
if(keyword_set(arr5)) then arr5=arr5_new

;
end
