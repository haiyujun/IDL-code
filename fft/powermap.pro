data = mrdfits('/Users/lihaiyu/Downloads/Oct10_2_V.fits')

f1=complexarr(640)
d=fltarr(6,512,512)

for j = 2, 7 do begin
    low=ceil(j*0.64*45)
    up=floor((j+1)*0.64*45)
    for x = 0, 511 do begin
        for y = 0, 511 do begin
            for i = low, up do begin
                f_data=fft(data[x,y,*])
                f1[i]=f_data[i]
                f1[640-i]=f_data[640-i]
            endfor
            d1=real_part(fft(f1,1))
            d[j-2,x,y]=total(d1^2)
        endfor
        print,x
    endfor
    print,j
endfor
end