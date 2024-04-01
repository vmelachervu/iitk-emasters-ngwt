%{
    eMasters - Communication Systems
    Simulation-based Design of 5G NR Wireless Standard - EE910
    Assignment 1 - CRC Implementation
    Professor : Rohit Budhiraja
    Roll number : 23156022
    Student Name : Venkateswar Reddy Melachervu    
    email : vmela23@iitk.ac.in
    (C) Venkateswar Reddy Melachervu.
%}

% Globals
generatorPoly = [];
transport_block_bits_array = [];
bitsArrayAsString = '';

% Prompt user if he wants to provide transport block
fprintf('\n');
disp("Welcome to interactive 5G NR CRC calculator as per 3GPP Standard TS38.212-f20");
fprintf('\n');
disp("Would you like to type in input tranport block binary data for CRC generation and validation?");
disp("1. Yes");
disp("2. No - I want program to auto-generate input transport data block");
disp("3. I am already feeling lucky!");
% Prompt the user for input
choice = input("Please enter your choice (1-3): ");
% Handle the user's choice
switch choice
    case 1
        disp("You chose to type in the transport block binary data");                
        transport_block_bits_array = read_binary_data();
        if isempty(transport_block_bits_array)
            return;
        else
            bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
            disp('The input transport block binary data you have typed in is:');                        
            disp(bitsArrayAsString);
        end
    case 2
        disp("You chose to let the program generate random transport block binary data");  
        disp("Generating a random transport block of random size between 24 and 319784 bits (inclusive) length...");  
        % let's set the lower and upper bounds of the transport block length as per 5G NR
        lowerBound = 24;
        upperBound = 319784;   
        size = randi([lowerBound, upperBound]);   
        disp("Randomly selected size for the input transmit block is:");
        fprintf('%d\n', round(size));
        transport_block_bits_array = randi([0, 1], 1, size);
        bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
        disp('Randomly generated transport block is :');                        
        disp(bitsArrayAsString);
        prompt = input("Press any key to continue...");
    case 3
        disp("You are feeling lucky already! Exiting the program...");
        return;        
    otherwise
        disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 3.");
        return;
end
fprintf('\n');

% Prompt user to choose the 5G NR generator polynomial for CRC
disp('Please select the 5G NR cyclic generator polynomial you would like to use:');
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
        generatorPoly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];   
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 3
        disp("You chose option 3: gCRC24C - 24-bit cyclic generator polynomial");
        generatorPoly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];   
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 4
        disp("You chose option 4: gCRC16 - 16-bit cyclic generator polynomial");
        generatorPoly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];   
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 5
        disp("You chose option 5: gCRC11 - 11-bit cyclic generator polynomial");
        generatorPoly = [1 1 1 0 0 0 1 0 0 0 0 1];   
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 6
        disp("You chose option 6: gCRC6 - 6-bit cyclic generator polynomial");
        generatorPoly = [1 1 0 0 0 0 1];   
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 7
        disp("You chose to exit the program. Exiting...");
        return;        
    otherwise
        disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 7.");
        return;
end
fprintf('\n');

% Append padding bits to the data for CRC calculation
disp("Appending zeros to the data for CRC calculation...");
numCRCBits = length(generatorPoly) - 1;
paddedData = [transport_block_bits_array,zeros(1,numCRCBits)];
padded_data = sprintf('%d ', paddedData);
disp(['Padded data is: ' padded_data]);
fprintf('\n');

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
fprintf('\n');

% Validate the calculated CRC
prompt = input("Press any key to continue to validate the calculated CRC...");
disp(['Input data is: ' bitsArrayAsString]);
dataWithCRC = [transport_block_bits_array nativeCRCBits1];
dataWithCRCString = sprintf('%d', dataWithCRC);
disp(['CRC appended input data is:' dataWithCRCString]);
disp("Now validating the CRC...");
is_crc_valid = crc_5g_nr_validate(dataWithCRC,generatorPoly);
if is_crc_valid
    disp('CRC is valid. No error detected.');
else
    disp('CRC is invalid. Error detected.');
end
fprintf('\n');

% Error case validation
prompt = input("Press any key to continue to test error case...");
disp('For testing the error case, input data bits are flipped, calculated CRC of non-flipped input is appended and validated.');
prompt = input("Press any key to continue...");
disp('Flipping bits in the input data (without CRC bits)...')
flippedBits = ~transport_block_bits_array;
flippedBitString = sprintf('%d', flippedBits);
disp(['Flipped input data is:' flippedBitString]);
disp("Calculated CRC bits of un-flipped original input binary data:");
disp(crcBitsString1);
disp('Appending calculated CRC bits (of non-flipped input) to flipped input data...');
flippedInputWithCRC = [flippedBits nativeCRCBits1];
flippedInputWithCRCString = sprintf('%d', flippedInputWithCRC);
disp(['Flipped input data with CRC is:' flippedInputWithCRCString]);
disp('Validating flipped input data with CRC...');
is_crc_valid = crc_5g_nr_validate(flippedInputWithCRC,generatorPoly);
if is_crc_valid
    disp('CRC is valid. No error detected.');
else
    disp('CRC is invalid. Error detected.');
end
fprintf('\n');
disp('Thanks for testing and using this program - Venkateswar Reddy Melachervu. More about me at http://www.linkedin.com/in/vmelachervu');
% Error case validation


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


% CRC modulo-2 division fuction
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

% Function to validate 5G NR CRC
function crc_valid = crc_5g_nr_validate(receivedMessage, generatorPolynomial)
    % Perform CRC division on the received message
    remainder = crcDivision(receivedMessage, generatorPolynomial);    
    % Check if the remainder is all zeros (indicating no error)
    if all(remainder == zeros(1, length(remainder)))
        crc_valid = 1;
        % disp('CRC is valid. No error detected.');
    else
        crc_valid = 0;
        % disp('CRC is invalid. Error detected.');        
    end
end