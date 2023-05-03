% Match the interference signal with the orignal SOI template signal
% Take FFT of the interference signal to convert to frequency domain
S_InfR = fft1d2(sInfR);                                % FFT for the time domain signal (FFT along each eta row)
% Apply the matched filter
S_I  = repmat(G,size(S_InfR,1),1).*S_InfR;
% Get the leaked interference 
Interf = Src - S_I;
interf = ifft1d2(Interf);
% imagesc(real(interf)
% Interferece power calculation
Pi = 10*log10(sum((abs(interf)).^2,'all')/size(interf,1)); 
% Signal power calculation
Pr = 10*log10(sum((abs(sqd_ref)).^2,'all')/size(sqd_ref,1)); 
% Signal to interference ratio SIR
SIR =10*log10(Pr/Pi)