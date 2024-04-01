% QPSK Modulation and Demodulation Example

% Parameters
nBits = 10;     % Number of bits to transmit
EbNo_dB = 10;     % Eb/No in dB

% Generate random bits
bits = randi([0, 1], 1, nBits);

% QPSK modulation
modulatedSymbols = qpskModulation(bits);

% Add AWGN noise to the channel
receivedSymbols = addAWGN(modulatedSymbols, EbNo_dB);

% QPSK demodulation
demodulatedBits = qpskDemodulation(receivedSymbols);

% Calculate bit error rate (BER)
ber = calculateBER(bits, demodulatedBits);

% Display results
fprintf('Bit Error Rate (BER): %.4f\n', ber);

% QPSK modulation function
function symbols = qpskModulation(bits)
    % Convert bits to symbols
    symbols = sqrt(0.5) * (2 * bits(1:2:end) - 1 + 1i * (2 * bits(2:2:end) - 1));
end

% Add AWGN noise function
function receivedSymbols = addAWGN(symbols, EbNo_dB)
    % Convert Eb/No from dB to linear scale
    EbNo = 10^(EbNo_dB / 10);
    
    % Calculate symbol energy
    Es = mean(abs(symbols).^2);
    
    % Calculate noise variance
    N0 = Es / (2 * EbNo);
    
    % Generate complex Gaussian noise
    noise = sqrt(N0/2) * (randn(size(symbols)) + 1i * randn(size(symbols)));
    
    % Add noise to symbols
    receivedSymbols = symbols + noise;
end

% QPSK demodulation function
function demodulatedBits = qpskDemodulation(symbols)
    % Calculate decision threshold
    threshold = 0;
    
    % Perform symbol to bit conversion
    realBits = real(symbols) > threshold;
    imagBits = imag(symbols) > threshold;
    
    % Combine real and imaginary bits
    demodulatedBits = reshape([realBits; imagBits], 1, []);
end

% Calculate bit error rate (BER) function
function ber = calculateBER(bits, demodulatedBits)
    % Count bit errors
    nErrors = sum(bits ~= demodulatedBits);
    
    % Calculate bit error rate
    ber = nErrors / numel(bits);
end
