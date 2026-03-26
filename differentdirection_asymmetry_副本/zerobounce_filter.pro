;cosine bell
restore,'/Users/lihaiyu/Downloads/t_d_curv_4_use.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/c_corr/ccr_asymmetry/differentdirection_asymmetry/result/ccrangle/ccrangle_total.sav'
;tt[1]=19*45/60.0
y=interpol(tt,1440*2+1)*60/45
;=interpol(dd,1450*2)
ccr1=fltarr(12,151,301)
tmp0=fltarr(151,151)
tmp=fltarr(151,301)
;y[20:40]=findgen(21)*12/20+9

for i = 20, 150 do begin
   tmp0[i,round(y[i])-5:round(y[i])+14]=1
   tmp0[i,round(y[i])-14:round(y[i])-5]=reverse(0.5*(cos(findgen(10)/9.*!pi)+1))
   tmp0[i,round(y[i])+14:round(y[i])+23]=0.5*(cos(findgen(10)/9.*!pi)+1)
   ;print,round(y[i])
endfor
tmp[*,150:300]=tmp0
tmp[*,0:150]=reverse(tmp0,2)

for l = 0, 11 do begin
    ccr1[l,*,*]=ccrangle_total[l,*,*]*tmp
endfor
;ccr1=ccr0*tmp0
;ccr2=ccr[*,0:150]*tmp2


end