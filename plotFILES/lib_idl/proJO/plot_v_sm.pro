PRO PLOT_V_SM,X,Y,VRAD,N,x1,y1
C=1.-VRAD/2.9979E5
x1=x*c
y1=smooth(y,n)
PLOT,x1,y1,THI=2
RETURN
END