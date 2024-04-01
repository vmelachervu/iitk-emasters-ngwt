%{
    Library name : Library for QPSK Modulation and Demodulation for AWGN channel with Plots
    eMasters - Communication Systems - Simulation-based Design of 5G NR Wireless Standard - EE922    
    Roll number : 23156022
    Student Name : Venkateswar Reddy Melachervu    
    email : vmela23@iitk.ac.in

    History:
    V1.0.0  -   Initial complete solution - 16-06-2023        
    (C) Venkateswar Reddy Melachervu. 2023-2024.
%}


% QPSK test code

% Noise parameters
nBits = 10;     % Number of source bits to transmit
EbNo_dB = 0.5;     % Eb/No in dB

% Generate random source bits
% bits = randi([0, 1], 1, nBits);
bits = [0 1 0 1 1 1 0 0 1 1 0 1 0 1 0 1 0 0 1 1 1 1 1 0 0 0]; % information
num_of_bits = length(bits);
source_bits_string = sprintf("%d", bits);
disp(['Source bits data: ' source_bits_string]);
disp(['Length of source bits: ' num_of_bits]);

% Let's plot the source bits data
figure(1)
stem(bits, 'linewidth',3), grid on;
title('Source Data to QPSK Modulator');
% axis([ 0 11 0 1.5]);
axis([ 0 num_of_bits 0 1.5]);

% qpsk modulation
modulated_symbols = qpsk_modulation(bits);
mod_symbols_string = sprintf("%s", modulated_symbols);
disp(['Modulated data: ' mod_symbols_string]);
disp(['Length of modulated data: ' num2str(length(modulated_symbols))]);

% let's add some AWGN noise
disp(['Adding a noise value of Eb/N0 in dB:' num2str(EbNo_dB)])
received_symbols = add_AWGN(modulated_symbols, EbNo_dB);
received_symbols_string = sprintf("%s", received_symbols);
disp(['Received noisy symbols are: ' received_symbols_string]);

% qpsk - demodulation
demodulated_bits = qpsk_demodulation(received_symbols);
demodulated_bits_string = sprintf("%d", demodulated_bits);
disp(['Demodulated bits are: ' demodulated_bits_string]);

% Calculate bit error rate (BER)
ber = calculate_BER(bits, demodulated_bits);
disp(['BER is: ' num2str(ber)]);


% let's plot demodulated data
figure(3)
stem(demodulated_bits,'linewidth',3) 
title('Demodulated Data');
% axis([ 0 11 0 1.5]), grid on;
axis([ 0 num_of_bits 0 1.5]), grid on;


% QPSK modulation function
function symbols = qpsk_modulation(bits)
    % Convert bits to symbols
    symbols = sqrt(0.5) * (2 * bits(1:2:end) - 1 + 1i * (2 * bits(2:2:end) - 1));
end

% Add some AWGN noise function
function noisy_symbols = add_AWGN(symbols, EbNo_dB)
    % Convert Eb/No from dB to linear scale
    EbNo = 10^(EbNo_dB / 10);
    
    % Calculate symbol energy
    Es = mean(abs(symbols).^2);
    
    % Calculate noise variance
    N0 = Es / (2 * EbNo);
    
    % Generate complex Gaussian noise
    noise = sqrt(N0/2) * (randn(size(symbols)) + 1i * randn(size(symbols)));
    
    % Add noise to symbols
    noisy_symbols = symbols + noise;
end

% QPSK demodulation function
function demodulated_bits = qpsk_demodulation(symbols)
    % Calculate decision threshold
    threshold = 0;
    
    % Perform symbol to bit conversion
    realBits = real(symbols) > threshold;
    imagBits = imag(symbols) > threshold;
    
    % Combine real and imaginary bits
    demodulated_bits = reshape([realBits; imagBits], 1, []);
end

% Calculate bit error rate (BER) function
function ber = calculate_BER(bits, demodulatedBits) 
    % Count bit errors
    nErrors = sum(bits ~= demodulatedBits);
    
    % Calculate bit error rate
    ber = nErrors / numel(bits);
end
    

