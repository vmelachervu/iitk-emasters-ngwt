% ------------------------------------------------------------------------------------------------------%
%{
    Program name : 5G NR transceiver implementation with CRC and LDPC as per 3GPP TS38.212-f20 standard
    eMasters - Communication Systems - Simulation-based Design of 5G NR Wireless Standard - EE910
    Assignment 2 - 5G Transceiver Implementation - I
    Roll number : 23156022
    Student Name : Venkateswar Reddy Melachervu    
    email : vmela23@iitk.ac.in
    History:
    V1.0.0 - Initial complete solution - 05-06-2023    

    (C) Venkateswar Reddy Melachervu.
%}
% ------------------------------------------------------------------------------------------------------%

% 3GPP TS38.212-f20 standard constants 
    % transport block length bounds as per 5G NR
        lowerBound = 24;
        upperBound = 319784;   
    % base graphs
    BASE_GRAPH_1 = 1;
    BASE_GRAPH_2 = 2;

    % code block lengths
        K_cb_base_graph_1_cb_length = 8448;
        K_cb_base_graph_2_cb_length = 3840;
    % Lengths
        B = 0;
    % gCRC24A - 24-bit cyclic generator polynomial for TB segmentation
        gCRC24A_generatorPoly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];  
    % gCRC24B - 24-bit cyclic generator polynomial for TB segmentation    
        gCRC24B_generatorPoly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];   
    % gCRC24C - 24-bit cyclic generator polynomial for TB segmentation             
        gCRC24C_generatorPoly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];   
    % gCRC16 - 16-bit cyclic generator polynomial for TB segmentation       
        gCRC16_generatorPoly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1]; 
    
% variables 
generator_poly = [];
tb_crc_size = 0;
transport_block_bits_array = [];
bitsArrayAsString = '';


% TB length with CRC
B = 0;
% CB CRC length
L = 24;
% total number of segmented code blocks
C = 0;
% base graph matrix to use
base_graph_to_use = 0;
% code block length to use
K_cb = 0;
% net payload
B_prime = 0;
% user selected crc size for CRC computation
selected_crc_size= 0;
% selected generator polynomial for segmentation
selected_segmentation_generatorPoly = gCRC24C_generatorPoly;

% clear the command window
clc;

% Prompt user if he wants to provide transport block
transport_block_bits_array = user_prompt_menu_for_tb_data();

disp("Welcome to interactive 5G NR transceiver implementation (Matlab Asisgnment - 2 of Cohort-3 EE922) as per 3GPP Standard TS38.212-f20");
disp("Would you like to type in input tranport block binary data for CRC generation and validation?");
disp("1. Yes");
disp("2. No - I want the program to auto-generate input transport data block");
disp("3. I am already feeling lucky - that is it and that is all ;-)");
% Prompt the user for TB data input
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
        size = randi([lowerBound, upperBound]);   
        disp(['Randomly selected size for the input transmit block is:' num2str(size)]);       
        transport_block_bits_array = randi([0, 1], 1, size);
        bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
        disp(['Randomly generated transport block data is :']);                        
        disp(bitsArrayAsString);
    case 3
        disp("You are feeling lucky already! That is it an that is all. Exiting the program...");
        return;        
    otherwise
        disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 3.");
        return;
end
fprintf('\n');

% Prompt user to choose the 5G NR generator polynomial for CRC
[generator_poly, tb_crc_size] = read_user_selection_of_gp_for_crc_calculation();

disp('Please select the 5G NR cyclic generator polynomial - CGP - you would like to use for the CRC calculation');
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
        disp("You chose option 1: gCRC24A - 24-bit cyclic generator polynomial for CRC");  
        disp("As per the 3GPP standard this is an invalaid choice of CRC CGP for TB segmentation. ");
        disp("Please re-run the program and select a valid choice for the CGP - one of {gCRC24B, gCRC24C, gCRC16}.");
        disp("Exiting the program.");
        return;
        % generatorPoly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];   
        % selected_crc_size = size(generatorPoly);
        % generatorPolyString = sprintf('%d ', generatorPoly);
        % disp(['Generator polynomial is: ' generatorPolyString]);
    case 2
        disp("You chose option 2: gCRC24B - 24-bit cyclic generator polynomial");
        generatorPoly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];   
        L = length(generatorPoly) - 1;
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 3
        disp("You chose option 3: gCRC24C - 24-bit cyclic generator polynomial");
        generatorPoly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];   
        L = length(generatorPoly) - 1;
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 4
        disp("You chose option 4: gCRC16 - 16-bit cyclic generator polynomial");
        generatorPoly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];   
        L = length(generatorPoly) - 1;
        generatorPolyString = sprintf('%d ', generatorPoly);
        disp(['Generator polynomial is: ' generatorPolyString]);
    case 5
        disp("You chose option 5: gCRC11 - 11-bit cyclic generator polynomial");
        disp("As per the 3GPP standard this is an invalaid choice of CRC CGP for TB segmentation. ");
        disp("Please re-run the program and select a valid choice for the CGP - one of {gCRC24B, gCRC24C, gCRC16}.");
        disp("Exiting the program.");
        return;
        % generatorPoly = [1 1 1 0 0 0 1 0 0 0 0 1];
        % selected_crc_size = size(generatorPoly);
        % generatorPolyString = sprintf('%d ', generatorPoly);
        % disp(['Generator polynomial is: ' generatorPolyString]);
    case 6
        disp("You chose option 6: gCRC6 - 6-bit cyclic generator polynomial");
        disp("As per the 3GPP standard this is an invalaid choice of CRC CGP for TB segmentation. ");
        disp("Please re-run the program and select a valid choice for the CGP - one of {gCRC24B, gCRC24C, gCRC16}.");
        disp("Exiting the program.");
        return;
        % generatorPoly = [1 1 0 0 0 0 1];   
        % selected_crc_size = size(generatorPoly);
        % generatorPolyString = sprintf('%d ', generatorPoly);
        % disp(['Generator polynomial is: ' generatorPolyString]);
    case 7
        disp("You chose to exit the program. Exiting...");
        return;        
    otherwise
        disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 7.");
        return;
end
% compte B - TB length plus selected CRC length
B = length(transport_block_bits_array) + L;

% setting the selected base graph and code block length global variables per the base graph selection
base_graph_to_use = BASE_GRAPH_1;
K_cb = set_segmentation_code_block_length(base_graph_to_use);

% compute total number of code block segments and netpayload
[C, B_prime] = needed_segmentation_code_blocks_and_net_payload(B, L, K_cb);


fprintf('\n');
% ----------------------- end of the main program --------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------------------------------
% utility and computation functions
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


%Function to validate 5G NR CRC
function crc_5g_nr_validate(receivedMessage, generatorPolynomial)
    % Perform CRC division on the received message
    remainder = crcDivision(receivedMessage, generatorPolynomial);
    
    % Check if the remainder is all zeros (indicating no error)
    if all(remainder == zeros(1, length(remainder)))
        disp('CRC is valid. No error detected.');
    else
        disp('CRC is invalid. Error detected.');
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
end


% Calculate CRC with padding and return CRC as per 3GPP TS38.212-f20 standard
function calculated5GNRCRC = calculate5GNRCRC(message, generatorPoly)
    % Append padding bits to the data for CRC calculation
    disp("Appending zeros to the data for CRC calculation...");
    numCRCBits = length(generatorPoly) - 1;
    paddedData = [message,zeros(1,numCRCBits)];
    padded_data = sprintf('%d ', paddedData);
    disp(['Padded data is: ' padded_data]);    
    fprintf('\n');

    % Calculating CRC natively
    disp("Calculating CRC natively without any library...");
    remainder = crcDivision(paddedData, generatorPoly);    
    % Extract the CRC bits
    size = length(remainder);
    calculated5GNRCRC = remainder(size - numCRCBits + 1:size);

    % Display calculated CRC bits
    calculated5GNRCRCString = sprintf('%d ', calculated5GNRCRC);
    disp("Calculated CRC bits are:");
    fprintf(calculated5GNRCRCString);
    fprintf('\n');
end

% Selection of base graph and setting the segmented code block length
function [seg_code_block_length] = set_segmentation_code_block_length(selected_base_graph)
    % 5G NR standard values
    BASE_GRAPH_1 = 1;
    BASE_GRAPH_2 = 2;
    K_cb_base_graph_1_cb_length = 8448;
    K_cb_base_graph_2_cb_length = 3840;

    switch selected_base_graph
        case BASE_GRAPH_1           
            % code block length to use
            seg_code_block_length = K_cb_base_graph_1_cb_length;
            disp('Selected base graph is : Base Graph 1');
            disp(['Segmentation code block length to be used is: ' num2str(seg_code_block_length)]);
        case BASE_GRAPH_2           
            % code block length to use
            seg_code_block_length = K_cb_base_graph_2_cb_length;
            disp('Selected base graph is : Base Graph 2');
            disp(['Segmentation code block length to be used is: ' num2str(seg_code_block_length)]);
    end
end


% Compute total number of segmentation code blocks from the avaliable user selection
function [needed_code_blocks, net_pay_load] = needed_segmentation_code_blocks_and_net_payload(b, l, k_cb)   
    if (b > k_cb)
        needed_code_blocks = b/(code_block_length_to_use - l);
        net_pay_load = b + (needed_code_blocks * l);
    else
        needed_code_blocks = 1;
        net_pay_load = b;
        disp(['Segmentation is NOT needed - as TB + L=' num2str(b) ' is less than K_cb=' num2str(k_cb)]);
    end
    disp(['Total number of code blocks ' num2str(needed_code_blocks)]);    
end

% Function to prompt user for TB data
function transport_block_bits_array = user_prompt_menu_for_tb_data()
    % Prompt user if he wants to provide transport block
    disp("Welcome to interactive 5G NR transceiver implementation (Matlab Asisgnment - 2 of Cohort-3 EE922) as per 3GPP Standard TS38.212-f20");
    disp("Would you like to type in input tranport block binary data for CRC generation and validation?");
    disp("1. Yes");
    disp("2. No - I want the program to auto-generate input transport data block");
    disp("3. I am already feeling lucky - that is it and that is all ;-)");
    % Prompt the user for TB data input
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
            size = randi([lowerBound, upperBound]);   
            disp(['Randomly selected size for the input transmit block is:' num2str(size)]);       
            transport_block_bits_array = randi([0, 1], 1, size);
            bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
            disp(['Randomly generated transport block data is :']);                        
            disp(bitsArrayAsString);
        case 3
            disp("You are feeling lucky already! That is it an that is all. Exiting the program...");
            return;        
        otherwise
            disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 3.");
            return;
    end
end

% Function to prompt user to choose the 5G NR generator polynomial for CRC
function [generatorPoly, selected_crc_size] = read_user_selection_of_gp_for_crc_calculation()
    disp('Please select the 5G NR cyclic generator polynomial - CGP - you would like to use for the CRC calculation');
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
            disp("You chose option 1: gCRC24A - 24-bit cyclic generator polynomial for CRC");                         
            generatorPoly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 2
            disp("You chose option 2: gCRC24B - 24-bit cyclic generator polynomial");
            generatorPoly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 3
            disp("You chose option 3: gCRC24C - 24-bit cyclic generator polynomial");
            generatorPoly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 4
            disp("You chose option 4: gCRC16 - 16-bit cyclic generator polynomial");
            generatorPoly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 5
            disp("You chose option 5: gCRC11 - 11-bit cyclic generator polynomial");            
            generatorPoly = [1 1 1 0 0 0 1 0 0 0 0 1];
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 6
            disp("You chose option 6: gCRC6 - 6-bit cyclic generator polynomial");            
            generatorPoly = [1 1 0 0 0 0 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 7
            disp("You chose to exit the program. Exiting...");
            return;        
        otherwise
            disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 7.");
            return;
    end
end 