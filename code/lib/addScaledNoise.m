function eeg = addScaledNoise(eeg,SNR,N,fs)
% INPUT: 
%     eeg             matrix, size m (trials) x n (timepoints) 


% get rms of EEG 
eeg_rms = rms(eeg,2); 

for triali=1:size(eeg,1)
    % generate pink noise
    noise = getPinkNoise(N,fs); 
    % scale it to achieve requested SNR
    noise = (noise / rms(noise)) * (eeg_rms(triali)/SNR); 
    % add noise to EEG signal
    eeg(triali,:) = eeg(triali,:) + noise;     
end

