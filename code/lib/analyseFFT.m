function [itpc,aXbeat,meanVector,meanAngle,mXavgFreq,mXavgTime] = analyseFFT(eeg,N,fs,beatFrex,maxFreqLim,snrmin,snrmax)

maxFreqIdx = round(maxFreqLim/fs*N)+1; 
beatFrexIdx = round(beatFrex/fs*N)+1; 


%% PHASE 

% get complex-valued FFT
X = fft(eeg,[],2) / N * 2; 
X = X(:,1:maxFreqIdx); 

% get angles
aX = angle(X); 

% get angles at beat frequency
aXbeat = aX(:,beatFrexIdx); 

% get mean vector, it's angle and length (aka itpc)
meanVector  = mean(exp(1j*aX),1); 
meanAngle   = angle(meanVector); 
itpc        = abs(meanVector); 


%% MAGNITUDE

% average trials in frequency domain (i.e. average magnitudes)
mXavgFreq = mean(abs(X),1); 
% set DC to 0
mXavgFreq(1) = 0; 
% subtract surrounding bins
mXavgFreq = subtractSNR(mXavgFreq,snrmin,snrmax); 


% average trials in the time domain (i.e. average complex numbers)
mXavgTime = abs(fft(mean(eeg,1)) / N * 2); 
mXavgTime = mXavgTime(:,1:maxFreqIdx); 
% set DC to 0
mXavgTime(1) = 0; 
% subtract surrounding bins
mXavgTime = subtractSNR(mXavgTime,snrmin,snrmax); 
