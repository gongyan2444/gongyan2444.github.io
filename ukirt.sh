#convert, convert to sdf files
> fits2ndf in="w${date}_${num}_sf_st.fit" out="temp" container=true
#kappa, extract individual files
> ndfcopy in='temp.HDU_1' out='tempw' reset 
> ndfcopy in='temp.HDU_2' out='tempx' reset 
> ndfcopy in='temp.HDU_3' out='tempy' reset 
> ndfcopy in='temp.HDU_4' out='tempz' reset 
#compadd, compress
#kappa, align wcs
> wcsalign "in=*_bin ref=z101_bin lbnd=! acc=0.2 method=bilinear params=[0,2] out=*_wcs reset"
#ccdpack, mosaic
> makemos in='*_wcs' out=dr21_mosaic method=broad zero reset
