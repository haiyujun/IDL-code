loadct,3
nxx=512.0
ntt=639.0;原数据640差分后
fits2map,'/Users/lihaiyu/Desktop/IDL code/c_corr/ccr_asymmetry/differentdirection_asymmetry/dataset_quiet/Oct03_3_V.fits',map   
data0 = map.data

print,'part_0 finished'

data=fltarr(nxx,nxx,ntt)
for i = 0, n_elements(data[0,0,*])-1 do begin
    data[*,*,i]=data0[*,*,i+1]-data0[*,*,i]
endfor  ;进行差分运算，滤掉低频对流


;kmin=1/(12.15*0.03*nxx),Mm^-1
;vmin=1000/(ntt*45),mHz

k=fltarr(nxx^2)
a=fltarr(ceil(sqrt(2)*nxx),ntt)
ftb=reform((abs(fft(data)))^2)

for i = 0, nxx-1 do begin
    for j = 0, nxx-1 do begin
        k[nxx*i+j]=round(sqrt((i*1.0)^2+(j*1.0)^2))  
        ;for x = 0, 639 do a[k[500*i+j],x]=a[k[500*i+j],x]+ftb[i,j,x]
        a[k[nxx*i+j],*]=a[k[nxx*i+j],*]+ftb[i,j,*]
    endfor
endfor
;v=findgen(640)*1000/(640*45)    ;频率序列
;kp=findgen(142)/(12.15*0.03*100)    ;波数序列
;contour,a(*,0:319),kp,v(0:319),nlevels = 60,/fill,max_value = 20


tvim,a[0:nxx/2,0:ntt/2],xrange=findgen((nxx+1)/2)/(12.15*0.03*nxx),yrange=findgen((ntt+1)/2)*1000.0/(ntt*45)
;tvim,a[0:nxx/2,0:ntt/2],xrange=findgen((nxx+1)/2)/(12.15*0.03*nxx)*4373.0972,yrange=findgen((ntt+1)/2)*1000.0/(ntt*45)
end