fits2map,'/Users/lihaiyu/Desktop/IDL code/c_corr/ccr_asymmetry/differentdirection_asymmetry/dataset_quiet/Oct02_1_V.fits',map   
data0 = map.data

print,'part_0 finished'

data=fltarr(512,512,639)
for i = 0, n_elements(data[0,0,*])-1 do begin
    data[*,*,i]=data0[*,*,i+1]-data0[*,*,i]
endfor  ;进行差分运算，滤掉低频对流

data_spot=fltarr(512,512,639)
data_out=fltarr(512,512,639)
for i = 199, 300 do begin
    for j = 199, 350 do begin
        if((i-250.0)^2+(j-257.0)^2 le 2500.0) then data_spot[i,j,*]=data[i,j,*]
    endfor
endfor
data_out=data-data_spot    ;将滤波后数据分为黑子范围和外部范围，所需范围之外全为0（对于宁静区这块代码没用）

ncircle=fltarr(151)
circlex=fltarr(151,round(2*!pi*150))
circley=fltarr(151,round(2*!pi*150))
for r = 20, 150 do begin
    ncircle[r]=max([round(2*!pi*r),1])  ;取半径为r的圆上的点数;这一行的round可以去掉试试。
    for i = 1, ncircle[r] do begin
        circlex[r,i-1]=round(r*cos(2*!pi*i/ncircle[r]))
        circley[r,i-1]=round(r*sin(2*!pi*i/ncircle[r]))  ;半径为r的圆上点相对坐标
    endfor
endfor 
print,'part_1 finished'

temp=fltarr(639)
ccr=fltarr(151,301)
lag=findgen(301)-150
for r = 20, 150 do begin
    for i = 0, 511 do begin
        for j = 0, 511 do begin
            temp=fltarr(639)
            for k = 0, ncircle[r]-1 do begin
                if((i+circlex[r,k] ge 0) and (j+circley[r,k] ge 0) and (i+circlex[r,k] le 511) and (j+circley[r,k] le 511)) then temp=temp+data[i+circlex[r,k],j+circley[r,k],*]
            endfor
            if((total(temp) ne 0) and (total(data[i,j,*]) ne 0)) then ccr[r,*]=ccr[r,*]+c_correlate(data[i,j,*],temp,lag)
        endfor
        print,i
    endfor
    print,r
endfor

end