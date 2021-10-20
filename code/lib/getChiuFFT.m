function [mX,pX]=getChiuFFT(pattern)

N = length(pattern); 
X = fft(pattern)/N; 
mX = abs(X); 
pX = angle(X); 



