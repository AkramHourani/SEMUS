% Match the interference signal with the orignal SOI template signal
% Take FFT of the interference signal to convert to frequency domain
S_InfR = fft1d2(sInfR);                                % FFT for the time domain signal (FFT along each eta row)
% Apply the matched filter
S_I  = repmat(G,size(S_InfR,1),1).*S_InfR;               
