fits2map,'/Users/lihaiyu/Desktop/IDL code/c_corr/ccr_asymmetry/differentdirection_asymmetry/dataset_spot/Oct10_2_V_1.fits',map   
data0 = map.data

print,'part_0 finished'

centerx=235
centery=256
spotr=48;黑子的中心x，y坐标以及半径

data=fltarr(512,512,639)
for i = 0, n_elements(data[0,0,*])-1 do begin
    data[*,*,i]=data0[*,*,i+1]-data0[*,*,i]
endfor  ;进行差分运算，滤掉低频对流

data_spot=fltarr(512,512,639)
data_out=fltarr(512,512,639)
for i = centerx-spotr-1, centerx+spotr+1 do begin
    for j = centery-spotr-1, centery+spotr+1 do begin
        if (i-centerx)^2+(j-centery)^2 le (spotr*spotr) then data_spot[i,j,*]=data[i,j,*]
        ;if (((i-centerx)^2+(j-centery)^2 ge --) and ((i-centerx)^2+(j-centery)^2 le --)) then data0[i,j,*]=0
    endfor
endfor
data_out=data-data_spot    ;将滤波后数据分为黑子范围和外部范围，所需范围之外全为0

ncircle=fltarr(151)
circlex=fltarr(151,round(2*!pi*150))
circley=fltarr(151,round(2*!pi*150))
;+
;dbg=fltarr(301,301)
;-
for r = 20, 150 do begin
    ncircle[r]=max([round(2*!pi*r),1])  ;取半径为r的圆上的点数;这一行的round可以去掉试试。
    for i = 1, ncircle[r] do begin
        circlex[r,i-1]=round(r*cos(2*!pi*i/ncircle[r]))
        circley[r,i-1]=round(r*sin(2*!pi*i/ncircle[r]))  ;半径为r的圆上点相对坐标
    endfor
endfor 
print,'part_1 finished'
;_____________________________________________________________________________________________________________
;m=fltarr(301,311,151)
;n=0
tmpalpha0=0
tmpalpha1=0
tmpalpha2=0
temp=fltarr(639)
ccr_1=fltarr(12,151,301)
for l = 0, 5 do begin;l=0～11代表12个方向，为了进行并行运算把程序分成4份分别计算l=0～2，l=3～5...
    for i = centerx-spotr-1, centerx+spotr+1 do begin
        for j = centery-spotr-1, centery+spotr+1 do begin
            tmpalpha0=180/!PI*atan(1.0*(j-centery)/(i-centerx))
            tmpalpha1=180/!PI*atan(j-centery,i-centerx)
            tmpalpha2=180/!PI*acos((i-centerx)/sqrt((i-centerx)*(i-centerx)+(j-centery)*(j-centery)))
            for r = 20, 150 do begin
                temp=fltarr(639)
                ;n=0
                for k = 0, ncircle[r]-1 do begin
                    if((i+circlex[r,k] ge 0) and (j+circley[r,k] ge 0) and (i+circlex[r,k] le 511) and (j+circley[r,k] le 511) $
                    and (((tmpalpha2 le 150) and (tmpalpha1-180/!PI*atan(circley[r,k],circlex[r,k]) ge -30) and (tmpalpha1-180/!PI*atan(circley[r,k],circlex[r,k]) le 30)) or ((tmpalpha2 gt 150) and (180/!PI*acos(circlex[r,k]/sqrt(circlex[r,k]*circlex[r,k]+circley[r,k]*circley[r,k])) gt 120) and (tmpalpha0-180/!PI*atan(circley[r,k]/circlex[r,k]) ge -30) and (tmpalpha0-180/!PI*atan(circley[r,k]/circlex[r,k]) le 30))))$
                    ;and (tmpalpha2-180/!PI*acos(circlex[r,k]/sqrt(circlex[r,k]*circlex[r,k]+circley[r,k]*circley[r,k])) ge -15) and (tmpalpha2-180/!PI*acos(circlex[r,k]/sqrt(circlex[r,k]*circlex[r,k]+circley[r,k]*circley[r,k])) le 15)) $
                    then temp=temp+data_out[i+circlex[r,((k-l*ncircle[r]/12+ncircle[r]) mod ncircle[r])],j+circley[r,((k-l*ncircle[r]/12+ncircle[r]) mod ncircle[r])],*]
                        ;n=n+1
                        ;dbg[150+circlex[r,((k-l*ncircle[r]/12) mod ncircle[r])],150+circley[r,((k-l*ncircle[r]/12) mod ncircle[r])]]=1               
                endfor
                ;m[i,j,r]=n/ncircle[r]
                if total(temp) ne 0 and total(data_spot[i,j,*]) ne 0 then $
                   ccr_1[l,r,*]=ccr_1[l,r,*]+c_correlate(data_spot[i,j,*],temp,findgen(301)-150)
            endfor
        endfor
        print,i
    endfor
    print,'l=',l
endfor




end