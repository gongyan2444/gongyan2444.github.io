;;;; an example of dust SED fitting assuming a single modified-blackbody spectrum

FUNCTION PLANCK, X, P
  RETURN, 1.53392e-08*2*6.6260755d-27*(X*1d9)^(3)/(2.99792458d10)^(2.)*1d/(exp(6.6260755d-27*(X*1d9)/1.380658d-16/P[0])-1d)*(1.-exp(-P[1]*1d22*2.8*1.6726231d-24*0.009*(X/230.)^(P[2])))/1d-23
END

;Solid angle
omega = 1.53392e-08
;7.7654585e-08
;1.5978312d-08
print, PLANCK(345, [20, 1.6, 2.0])
;FUNCTION PLANCK, X, P
;  RETURN, 1.5978312d-08*2*6.6260755d-27*(230*1d9)^(3)/(2.99792458d10)^(2.)*1d/(exp(6.6260755d-27*(X*1d9)/1.380658d-16/P[0])-1d)*(-P[1]*1d22*2.8*1.6726231d-24*0.09*(X/230.)^(P[2]+3))/1d-23
;END


;/(2.99792458d10)^(2.)1d/(exp(6.6260755d-27*(X*1d9)/1.380658d-16/P[0])-1d)
;*(1.-exp(-P[1]*2.8*1.6726231d-24*0.09*(X/230.)^(P[2])))

;h = 6.6260755d-27
;c = 2.99792458d10
;k = 1.380658d-16
;mp = 1.6726231d-24
;maj = 30d 
;min = 20d
;omega = 1.133*(maj/3600d*!dtor)*(min/3600*!dtor)
;tau = nh2 * 2.8*mp*0.09*(nu/230.)^(beta)    ;;; nu in units of GHz
;bv = 2*h*nu^(3)/c^(2)*1./(exp(h*nu/k/td)-1) ;;; nu in units of Hz
;sv = omega*bv*(1-exp(-tau))


;;X: nu in units of GHz; P[0]: Td;  P[1]: nh2; P[2]: beta; 

;FUNCTION PLANCK, X, P
;  RETURN, omega*2d*h*(X*1d9)^(3)/c^(2d)*1d/(exp(h*(X*1d9)/k/P[0])-1d)*(1d-exp(-P[1]*2.8*mp*0.09*(X/230.)^(P[2])))
;END

;FUNCTION PLANCK, X, P
;  RETURN, P[0]+GAUSS1(X, P[1:3])
;END

;;; sed data
;freq = [552.10397d, 582.68699d, 614.95889d, 686.02393d, 1108.2900d]
;sv   = [7.7d, 9.5d, 13.1d, 17.4d, 57.1d]

;freq = [250, 268, 345]
;sv   = [1.39d, 1.782, 2.96d]
;sverr = [0.209, 0.367, 0.444]

;freq = [250, 345]
;sv   = [1.39d, 2.96d]

freqall = [250, 268, 345, 599.584916, 856.54988, 1199.169832, 1873.7028625000003, 4282.749400000001, 1873.7029, 2141.3747, 3331.0273, 4612.1917] ;552.10397d, 582.68699d, 614.95889d, 686.02393d];, 1108.2900d]
svall   = [1.39d, 1.782, 2.96d, 23.142, 44.916, 127.186, 231.152, 359.565, 275.100, 302.900, 176.300, 313.10] ;7.7d, 9.5d, 13.1d, 17.4d];, 57.1d]


;freqall = [250, 268, 345, 599.584916, 856.54988, 1199.169832, 1873.7028625000003, 4282.749400000001, 1873.7029, 2141.3747, 3331.0273, 4612.1917, 552.10397d, 582.68699d, 614.95889d, 686.02393d, 1108.2900d]
;svall   = [1.39d, 1.782, 2.96d, 23.142, 44.916, 127.186, 231.152, 359.565, 275.100, 302.900, 176.300, 313.10, 7.7d, 9.5d, 13.1d, 17.4d, 57.1d]

;freq = [250, 268, 345, 599.584916, 856.54988, 1199.169832]
;sv   = [1.39d, 1.782, 2.96d, 23.142, 44.916, 127.186]

freq = [250, 268, 345, 599.584916, 856.54988, 1199.169832, 1873.7028625000003, 1873.7029, 2141.3747, 3331.0273] ;552.10397d, 582.68699d, 614.95889d, 686.02393d];, 1108.2900d]
sv   = [1.39d, 1.782, 2.96d, 23.142, 44.916, 127.186, 231.152, 275.100, 302.900, 176.300] ;7.7d, 9.5d, 13.1d, 17.4d];, 57.1d]


;freq = [250, 268, 345, 599.584916, 856.54988, 1199.169832, 1873.7028625000003, 1873.7029, 2141.3747, 3331.0273, 552.10397d, 582.68699d, 614.95889d, 686.02393d, 1108.2900d]
;sv   = [1.39d, 1.782, 2.96d, 23.142, 44.916, 127.186, 231.152, 275.100, 302.900, 176.300, 7.7d, 9.5d, 13.1d, 17.4d, 57.1d]


start = [20d, 0.55, 1.9]
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
;parinfo[0].fixed = 1
;parinfo[0].fixed = 1
;parinfo[0].limited[0]=1
;parinfo[0].limits[0]=10
;parinfo[0].limits[1]=200
;parinfo[1].limited[0]=1
;parinfo[1].limits[0]=0.1
;parinfo[1].limits[1]=20
;parinfo[2].limited[0]=1
;parinfo[2].limits[0]=0.5
;parinfo[2].limits[1]=3.0

result = MPFITFUN('PLANCK', freq, sv, sv*0.2, start,parinfo=parinfo, PERROR=myerr)

print, result
print, myerr

mass  = omega*(3.8*1000d*3.086d18)^(2.)*result[1]*1d22*2.8*1.67d-24/1.989d33
errmass = omega*(3.8*1000d*3.086d18)^(2.)*myerr[1]*1d22*2.8*1.67d-24/1.989d33
print, "the total mass is ", mass,"+-",errmass

;aa = findgen(4000)
;plot, aa, PLANCK(aa, [20.,0.5, 2.])
;print,PLANCK(250, [20, 5.5d21, 0.77])
myfreq = findgen(10000)
ccc = 299792458
mywave = ccc/myfreq/1d9*1d6
mysv  = PLANCK(myfreq, result)




entry_device=!d.name
set_plot,'ps'
ysize=20
xsize=25
device,filename='sed.eps',encapsulated=0,/color,/portrait,xsize=xsize,ysize=ysize,yoffset=0.5

tvlct,0,0,0,0
tvlct,255,0,0,1
tvlct,0,255,0,2
tvlct,0,0,255,3
tvlct,128,128,128,4
tvlct,255,255,255,5



;print,ccc/freq/1d9*1d6
;ploterror, ccc/freqall/1d9*1d6, svall, svall*0.15,
;psym=1,xrange=[10,5000],yrange=[0.1,500],/xlog,/ylog,/xst,/yst
plotsym, 0 , 1,/fill, thick =3
ploterror, freqall, svall, svall*0.2, psym=8,xrange=[100,10000],yrange=[0.1,500],/xlog,/ylog,xstyle=8, ystyle=1, xtitle='Frequency (GHz)', ytitle='S!d'+textoidl('\nu')+'!n (Jy)',color=0,xthick=4,ythick=4, charthick=4, charsize=2,position=[0.15, 0.15, 0.95, 0.9],ERRTHICK=4,ERRCOLOR=0

;oploterror,freqall, svall, svall*0.2, psym=8
oplot, myfreq, mysv, linestyle =2 ,color=4,thick=4
;oplot, mywave, mysv, linestyle =1,color='xx0000'xl

axis,xaxis=1,xrange=[ccc/100/1d9*1d6,ccc/10000/1d9*1d6],/xst,/xlog,xthick=4,charthick=4,charsize=2,xtitle='Wavelength ('+textoidl('\mu')+'m)',color=0
;axis,yaxis=1,yrange=[0.1, 500],/yst,/ylog,ythick=4, charthick=4, charsize=2,color=0

xyouts, 0.2, 0.82, "T!dd!n=25K",charsize=2,charthick=4,/normal,ALIGNMENT=0
xyouts, 0.2, 0.75, "N!dH2!n=4.4e22!ncm!u-2!n",charsize=2,charthick=4,/normal,ALIGNMENT=0
xyouts, 0.2, 0.68, textoidl('\beta')+"=1.7",charsize=2,charthick=4,/normal,ALIGNMENT=0

device,/close_file
set_plot,entry_device






end 
