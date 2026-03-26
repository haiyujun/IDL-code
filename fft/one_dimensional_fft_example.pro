time=FINDGEN(100)
series=10.0*SIN(2*!PI*time+!PI/4.0)
f_series=FFT(series)
print,f_series
plot, f_series.real() ,position=[0.1,0.1,0.45,0.9]
plot, f_series.imaginary() ,position=[0.55,0.1,0.9,0.9],/noerase
end