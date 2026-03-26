restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/cs.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/H_1.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/N2_1.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/va.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/vax.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/vaz.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/wc_1.sav'
n_x=511.0
n_z=512.0
centerxy=255.0

D3_cs=dblarr(n_x,n_x,n_z)
D3_H=D3_cs
D3_N2=D3_cs
D3_wc=D3_cs
D3_vaz=D3_cs
D3_varr=D3_cs
D3_vaxx=D3_cs
D3_vayy=D3_cs
D3_va=D3_cs

for i = 0, n_x-1 do begin
    tic
    for j = 0, n_x-1 do begin
        if ((i-centerxy)^2+(j-centerxy)^2 LE centerxy^2) AND ((i-centerxy)^2+(j-centerxy)^2 GT 0) then begin
            ;标量
            rij=round(sqrt((i-centerxy)^2 + (j-centerxy)^2))
            D3_cs[i,j,*]=cs[centerxy+rij,*]
            D3_H[i,j,*]=H[centerxy+rij,*]
            D3_N2[i,j,*]=N2[centerxy+rij,*]
            D3_wc[i,j,*]=wc[centerxy+rij,*]
            D3_vaz[i,j,*]=vaz[centerxy+rij,*]
            D3_va[i,j,*]=va[centerxy+rij,*]
            ;矢量
            D3_varr[i,j,*]=vax[centerxy+rij,*]
            D3_vaxx[i,j,*]=D3_varr[i,j,*]*(i-centerxy)/sqrt((i-centerxy)^2 + (j-centerxy)^2)
            D3_vayy[i,j,*]=D3_varr[i,j,*]*(j-centerxy)/sqrt((i-centerxy)^2 + (j-centerxy)^2)
        endif
    endfor
    print,i
    toc
endfor



end