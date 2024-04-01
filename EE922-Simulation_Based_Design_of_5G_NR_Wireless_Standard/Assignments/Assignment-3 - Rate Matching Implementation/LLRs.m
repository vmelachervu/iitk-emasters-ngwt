% Received QPSK symbols
receivedSymbols = [-1-1i, 1+1i, -1+1i, 1-1i]; % Example received symbols

% QPSK demodulation
demodulatedBits = QPSK.qpsk_modulation(receivedSymbols);

% Log-likelihood ratios (LLRs) for '1' and '0'
LLR_1 = log((1 + exp(-abs(receivedSymbols - 1i*pi/4)).^2) ./ (1 + exp(-abs(receivedSymbols - 3i*pi/4)).^2));
LLR_0 = log((1 + exp(-abs(receivedSymbols + 1i*pi/4)).^2) ./ (1 + exp(-abs(receivedSymbols - 1i*pi/4)).^2));

% Display the LLRs
disp('Log-Likelihood Ratios (LLRs):');
disp(['LLR for 1: ' num2str(LLR_1)]);
disp(['LLR for 0: ' num2str(LLR_0)]);
