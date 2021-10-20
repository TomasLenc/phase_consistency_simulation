function x=getPinkNoise(N,fs)

% tau = 200; 

T = N/fs; 
hz = [0:N-1]/T; 

if mod(N,2)==0
    hN          = N/2-1; 
%     magnitudes  = rand(1,hN) .* exp(-[1:hN]/tau); 
    magnitudes  = rand(1,hN) ./ hz(2:hN+1); 
    phases      = exp(1j*2*pi*rand(1,hN)); 
    X           = [0, magnitudes, 0, magnitudes(end:-1:1)] .* [0, phases, 0, conj(phases(end:-1:1))]; 
else
    hN          = floor(N/2); 
%     magnitudes  = rand(1,hN) .* exp(-[1:hN]/tau); 
    magnitudes  = rand(1,hN) ./ hz(2:hN+1); 
    phases      = exp(1j*2*pi*rand(1,hN)); 
    X           = [0, magnitudes, magnitudes(end:-1:1)] .* [0, phases, conj(phases(end:-1:1))]; 
end

x = ifft(X); 



