% Generate binary data stream
data=randi([0 1],QPSK.NofBits,1);                               % Input binary data stream (0's and 1's) to modulate
L = 2 * QPSK.OF;                                                % Samples in each symbol (QPSK has 2 bits in each symbol)
ak = 2 * data-1;                                                % NRZ encoding 0-> -1, 1->+1
I = ak(1:2:end);                                                % I: Odd bit streams - baseband I channel waveform (no carrier)
Q = ak(2:2:end);                                                % Q: Even bit streams - baseband Q channel waveform (no carrier)
I = repmat(I,1,L).'; Q=repmat(Q,1,L).';                         % Even/odd streams at 1/2Tb baud
I = I(:).'; Q = Q(:).';                                         % Serialize
t = 0:1/RadPar.fs:(length(I)-1)/RadPar.fs;                      % Time base vector for the carrier modulated signal
% iChannel = I.*cos(2*pi*QPSK.fc*t);                              % In-phase channel
% qChannel = -Q.*sin(2*pi*QPSK.fc*t);                             % Quadrature-phase channel
% signal_IQ = iChannel + 1i* qChannel;                            % QPSK modulated signal with carrier
signal_IQ = iChannel + 1i* qChannel;                            % QPSK modulated signal with carrier
signal_IQ = signal_IQ.';
%% Define QPSK Timing
TxTime = Param.tg;
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all"));

QPSKTime = length(signal_IQ)/RadPar.fs;
Discard_Q = round(length(signal_IQ)* TxTime/QPSKTime);         % Number of samples on LORA data to be removed
Keep_Q = round(length(signal_IQ)* RxTime/QPSKTime);            % Number of samples on LORA data to be considered
%% Creat QPSK Interference Matrix
% First define Number of section in QPSK data
sections = round(length(signal_IQ) / (Discard_Q + Keep_Q));
% Creat Matrix
for i = 1: sections-1
  sqpsk(i,:) = signal_IQ(i+(i*Discard_Q)+((i-1)*Keep_Q):i+(i*Discard_Q)+(i*Keep_Q));
end
% Adjust size of QPSK signal to match the SAR raw data
sQPSK = zeros(size(sqd,1),size(sqd,2));
for i = 1 : QPSK.NumberofQPSK
scale(i) = 0.25 + (1 - 0.25) * rand(1,1);                       % 0.25 is the min limit of the scale and 1 is the max limit of the scale
QPSKi{i}(:,:) = padarray(sqpsk,[(round(size(sqd,1)*scale(i))-size(sqpsk,1)), (size(sqd,2)-size(sqpsk,2))],0,'pre');
sQPSKi = padarray(QPSKi{i}(:,:),[(size(sqd,1)-size(QPSKi{i}(:,:),1)), (size(sqd,2)-size(QPSKi{i}(:,:),2))],0,'post');
sQPSK = sQPSKi +sQPSK;
end
% imagesc(real(sQPSK))
NI04a_QPSKGain
sQPSK = sQPSK .* sqrt(PQPSK);
