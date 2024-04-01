% add folders and sub-folders in the current working directory path
currentDir = pwd;
addpath(genpath(currentDir));

% LDPC parameters
N = 3000;  % Codeword length
K = 7040;  % Message length
R = K / N; % Code rate

% Generate random message bits
messageBits = randi([0 1], K, 1);

% LDPC encoding
% H = dvbs2ldpc(R); % LDPC parity-check matrix
% encoder = comm.LDPCEncoder(H);
% encodedBits = encoder(messageBits);
encodedBits = LDPCEncode(messageBits,1);

% QPSK modulation
qpskSymbols = qpskmod(encodedBits);

% Add AWGN noise
EbNo = 5;  % Energy per bit to noise power spectral density ratio (dB)
SNR = EbNo + 10*log10(R); % Signal-to-noise ratio (dB)
noiseVar = 10^(-SNR/10);
receivedSymbols = qpskSymbols + sqrt(noiseVar/2)*(randn(size(qpskSymbols)) + 1i*randn(size(qpskSymbols)));

% QPSK demodulation
demodulatedBits = qpskdemod(receivedSymbols);

% LDPC decoding
decoder = comm.LDPCDecoder(H);
decodedBits = decoder(demodulatedBits);

% Check if the decoded bits match the original message bits
bitErrors = sum(decodedBits ~= messageBits);
bitErrorRate = bitErrors / K;

disp(['Bit Error Rate: ' num2str(bitErrorRate)]);
