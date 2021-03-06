pro ng_expol1d, viit, nd
;
;-------same as ng_expol_3d, but: input and output in vector form-------
;
;print, '-------------------------calculating new iterate (ng)--------------------------'
;print, ''
;
a1=0.d0
a2=0.d0
b1=0.d0
b2=0.d0
c1=0.d0
c2=0.d0
;
for i=0L, nd-1 do begin
   a1 = a1+(viit(3,i)-2.d0*viit(2,i)+viit(1,i))^2.d0
   b1 = b1+(viit(3,i)-viit(2,i)-viit(1,i)+viit(0,i))*(viit(3,i)-2.d0*viit(2,i)+viit(1,i))
   c1 = c1+(viit(3,i)-viit(2,i))*(viit(3,i)-2.d0*viit(2,i)+viit(1,i))
   a2 = a2+(viit(3,i)-viit(2,i)-viit(1,i)+viit(0,i))*(viit(3,i)-2.d0*viit(2,i)+viit(1,i))
   b2 = b2+(viit(3,i)-viit(2,i)-viit(1,i)+viit(0,i))^2.d0
   c2 = c2+(viit(3,i)-viit(2,i))*(viit(3,i)-viit(2,i)-viit(1,i)+viit(0,i))
endfor
;
ap=(b2*c1-b1*c2)/(a1*b2-a2*b1)
bp=(a1*c2-a2*c1)/(a1*b2-a2*b1)
;
;extrapolation can yield to negative source-functions
;        if so, then skip extrapolation
negative=0

for i=0L, nd-1 do begin
   cp=(1.d0-ap-bp)*viit(3,i) + ap*viit(2,i)+bp*viit(1,i)
   if(cp lt 0.d0) then begin
;      print, 'ng-extrapolation yields negative source-functions => skip'
      negative=1
      break
    endif
endfor
;
if(negative eq 0) then begin
   viit(0,*)=(1.d0-ap-bp)*viit(3,*) + ap*viit(2,*)+bp*viit(1,*)
endif
;
;
;
end
