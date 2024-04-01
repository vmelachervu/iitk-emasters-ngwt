%{
    Library name : Library for QPSK
    eMasters - Communication Systems - Simulation-based Design of 5G NR Wireless Standard - EE922    
    Roll number : 23156022
    Student Name : Venkateswar Reddy Melachervu    
    email : vmela23@iitk.ac.in

    Description:
    Modulation, demodulation, adding AWGN noise are defined as static
    methods of the class QPSK. After adding this file path, these QPSK methods can be called as 
    below from other matlab files/functions:
        - Calling modulation function
            - symbols = QPSK.qpsk_modulation(source_data_bits)
        - Calling demodulation function
            - demodulated_bits = QPSK.qpsk_demodulation(source_data_bits)
        - Calling function to add AWGn noise
            - noisy_symbols = QPSK.add_AWGN(symbols, snr_in_db)

    History:
    V1.0.0  -   Initial complete solution - 16-06-2023        
    (C) Venkateswar Reddy Melachervu. 2023-2024.
%}

classdef QPSK
    methods(Static)
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
            % realBits = real(symbols) > threshold;
            % imagBits = imag(symbols) > threshold;
            realBits = real(symbols);
            imagBits = imag(symbols);
            
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
    end
end     

