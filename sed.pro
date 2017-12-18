;;;; an example of dust SED fitting assuming a single modified-blackbody spectrum
FUNCTION PLANCK, X, P
  RETURN, 7.7654585e-08*2*6.6260755d-27*(X*1d9)^(3)/(2.99792458d10)^(2.)*1d/(exp(6.6260755d-27*(X*1d9)/1.380658d-16/P[0])-1d)*(1.-exp(-P[1]*1d22*2.8*1.6726231d-24*0.009*(X/230.)^(P[2])))/1d-23
END
;Solid angle
;7.7654585e-08
;1.5978312d-08
print, PLANCK(345, [20, 0.55, 0.77])
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

freq = [250, 268, 345] ;552.10397d, 582.68699d, 614.95889d, 686.02393d];, 1108.2900d]
sv   = [1.39d, 1.782, 2.96d] ;7.7d, 9.5d, 13.1d, 17.4d];, 57.1d]


start = [15d, 0.55, 0.77]
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
parinfo[0].fixed = 1
;parinfo[1].fixed = 1
;parinfo[0].limited[0]=1
;parinfo[0].limits[0]=10
;parinfo[0].limits[1]=200
;parinfo[1].limited[0]=1
;parinfo[1].limits[0]=0.1
;parinfo[1].limits[1]=20
;parinfo[2].limited[0]=1
;parinfo[2].limits[0]=0.5
;parinfo[2].limits[1]=3.0

result = MPFITFUN('PLANCK', freq, sv, sv*0.15, start,parinfo=parinfo, PERROR=myerr)

print, result
print, myerr
;aa = findgen(4000)
;plot, aa, PLANCK(aa, [20.,0.5, 2.])
;print,PLANCK(250, [20, 5.5d21, 0.77])


end 
