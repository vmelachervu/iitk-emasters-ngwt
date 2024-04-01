%{
    eMasters - Communication Systems
    Simulation-based Design of 5G NR Wireless Standard - EE910
    Assignment 1 - CRC Implementation
    Roll number : 23156022
    Student Name : Venkateswar Reddy Melachervu    
%}

% Define globals
generatorPoly = [];
transport_block_bits_array = [];
bitsArrayAsString = '';

% Prompt user if he wants to provide transport block
disp("Would you like to type in tranport block binary data for CRC generation and validation?");
disp("1. Yes");
disp("2. No");
disp("3. Exit");

% Prompt the user for input
choice = input("Please enter your choice (1-3): ");
% Handle the user's choice
switch choice
    case 1
        disp("You chose to type in the transport block binary data");
        disp("Please be aware that this program assumes you are typing binary number and does not validate the typed in data!");
        % Prompt the user to enter the transport block binary data
        transport_block_bits_array = read_binary_data();
        if isempty(bits)
            return;
        else
            bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
            disp('You have entered :');                        
            disp(bitsArrayAsString);
        end
    case 2
        disp("You chose to let the program generate random transport block binary data");       
    case 3
        disp("You chose to exit the program. Exiting...");
        return;        
    otherwise
        disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 3.");
        return;
end

% Prompt user to choose the generator polynomial for CRC
disp('Please select the 5G NR cyclic generator polynomial you would like to use');
disp("1.gCRC24A(D) - 24-bit CRC - [D^24+D^23+D^18+D^17+D^14+D^11+D^10+D^7+D^6+D^5+D^4+D^3+D+1]");
disp("2.gCRC24B(D) - 24-bit CRC - [D^24+D^23+D^6+D^5+D+1]");
disp("3.gCRC24C(D) - 24-bit CRC - [D^24+D^23+D^21++D^20+D^17+D^15+D^13+D^12+D^8+D^4+D^2+D+1]");
disp("4.gCRC16(D) - 16-bit CRC - [D^16+D^12+D^5+1]");
disp("5.gCRC11(D) - 11-bit CRC - [D^11+D^10+D^9+D^5+1]");
disp("6.gCRC6(D) - 6-bit CRC - [D^6+D^5+1]");
disp("7.Exit");

choice = input("Please enter your choice (1-7): ");
% Handle the user's choice
switch choice
    case 1
        disp("You chose option 1: gCRC24A - 24-bit cyclic generator polynomial");        
        generatorPoly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];   
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 2
        disp("You chose option 2: gCRC24B - 24-bit cyclic generator polynomial");
    case 3
        disp("You chose option 3: gCRC24C - 24-bit cyclic generator polynomial");
    case 4
        disp("You chose option 4: gCRC16 - 16-bit cyclic generator polynomial");
    case 5
        disp("You chose option 5: gCRC11 - 11-bit cyclic generator polynomial");
    case 6
        disp("You chose option 6: gCRC6 - 6-bit cyclic generator polynomial");
    case 7
        disp("You chose to exit the program. Exiting...");
        return;        
    otherwise
        disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 7.");
        return;
end

% Append padding bits to the data for CRC calculation
disp("Appending zeros to the data for CRC calculation...");
numCRCBits = length(generatorPoly) - 1;
paddedData = [transport_block_bits_array,zeros(1,numCRCBits)];
padded_data = sprintf('%d ', paddedData);
disp(['Padded data is: ' padded_data]);

% Calculating CRC natively
disp("Calculating CRC natively without any library...");
remainder = crcDivision(paddedData, generatorPoly);
numCRCBits = length(generatorPoly) - 1;
% Extract the CRC bits
size = length(remainder);
nativeCRCBits1 = remainder(size - numCRCBits + 1:size);

% Display calculated CRC bits
crcBitsString1 = sprintf('%d ', nativeCRCBits1);
disp("Calculated CRC bits are:");
disp(crcBitsString1);

disp(['Input data is: ' bitsArrayAsString]);

dataWithCRC = [transport_block_bits_array nativeCRCBits1];
dataWithCRCString = sprintf('%d', dataWithCRC);
disp("CRC appended data is:");
disp(dataWithCRCString);

% Validate the calculated CRC
crc_5g_nr_validate(dataWithCRC,generatorPoly);



% If user selects program generation, generate a random transport block with user provided block length

% Prompt user for selecting the generator polynomial for CRC generation from the standard

% Display the transmit block and validate the CRC

% Error case - flip the bits in the transmit block, validate and display
% check-sum failure

% Function for reading transport block binary data from user 
function bit_array = read_binary_data()
    binary_string = input('Enter binary data: ', 's');
    bit_array = [];
    for i = 1:length(binary_string)
        if binary_string(i) == '0'
            bit_array = [bit_array, 0];
        elseif binary_string(i) == '1'
            bit_array = [bit_array, 1];
        else
            fprintf('Invalid binary data! Exiting the program. Please retry with valid data.\n');
            bit_array = [];
            return;
        end
    end
end


% Function to calculate 5G NR CRC
function crc = calculate_5g_nr_crc(transport_block_bits_array, generatorPoly)    
    numCRCBits = length(generatorPoly) - 1;
    % Append padding bits to the data for CRC calculation
    disp("Appending zeros to the data for CRC calculation...");
   
    paddedData = [transport_block_bits_array,zeros(1,numCRCBits)];
    padded_data = sprintf('%d ', paddedData);
    disp(['Padded data is: ' padded_data]);
    
    remainder = crcDivision(paddedData, generatorPoly);    
    % Extract the CRC bits
    crc = remainder((end - numCRCBits + 1):end);
end

%Function to validate 5G NR CRC
function is_5g_nr_crc_valid = crc_5g_nr_validate(receivedMessage, generatorPolynomial)
    % Perform CRC division on the received message
    remainder = crcDivision(receivedMessage, generatorPolynomial);
    
    % Check if the remainder is all zeros (indicating no error)
    if all(remainder == zeros(1, length(remainder)))
        disp('CRC is valid. No error detected.');
    else
        disp('CRC is invalid. Error detected.');
    end
end

% crcDivision fuction
function remainder = crcDivision(message, generatorPolynomial)
    % Perform CRC division using bitwise XOR
    remainder = message;
    divisor = generatorPolynomial;
    
    for i = 1 : length(message) - length(generatorPolynomial) + 1
        if remainder(i)
            remainder(i:i + length(generatorPolynomial) - 1) = bitxor(remainder(i:i + length(generatorPolynomial) - 1), divisor);
        end
    end
    
    % Trim leading zeros from the remainder
    % [~, startIndex] = max(remainder);
    % remainder = remainder(startIndex:end);
end