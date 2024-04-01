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
    V1.2.0 - User prompt menus added - 05-06-2023
    V1.2.2 - In-memory table implementation for base graphs, Zc table,finding the min Zc etc. - 05-06-2023
    V1.3.0 - Asking user for TB length  - inclusion in the user menu prompt - 06-06-2023

    (C) Venkateswar Reddy Melachervu. 2023-2024.
%}
% ------------------------------------------------------------------------------------------------------%

% 3GPP TS38.212-f20 standard constants 
    % transport block length bounds as per 5G NR
    lowerBound = 24;
    upperBound = 319784;    
    BASE_GRAPH_1_NUM = 1;
    BASE_GRAPH_2_NUM = 2;

    
    % 5G NR base graph key-value structures
    % base graph 1  
    BASE_GRAPH_1 = containers.Map();
    BASE_GRAPH_1('R_min_code_rate') = '1/3';
    BASE_GRAPH_1('base_matrix_size') = '46x68';
    BASE_GRAPH_1('K_b_systematic_columns') = 22;
    BASE_GRAPH_1('K_cb_max_info_block_size') = 8448;
    BASE_GRAPH_1('num_non_zero_elements') = 316;
    
    % base graph 2    
    BASE_GRAPH_2 = containers.Map();
    BASE_GRAPH_2('R_min_code_rate') = '1/5';
    BASE_GRAPH_2('base_matrix_size') = '42x62';
    BASE_GRAPH_2('K_b_systematic_columns') = 10;
    BASE_GRAPH_2('K_cb_max_info_block_size') = 3840;
    BASE_GRAPH_2('num_non_zero_elements') = 197;

    % lifting sizes table  
    LIFTING_SIZES_TABLE = containers.Map();
    LIFTING_SIZES_TABLE('0') = [2 4 8 16 32 64 128 256];
    LIFTING_SIZES_TABLE('1') = [3 6 12 24 48 96 192 384];
    LIFTING_SIZES_TABLE('2') = [5 10 20 40 80 160 320];
    LIFTING_SIZES_TABLE('3') = [7 14 28 56 112 224];
    LIFTING_SIZES_TABLE('4') = [9 18 36 72 144 288];
    LIFTING_SIZES_TABLE('5') = [11 22 44 88 176 352];
    LIFTING_SIZES_TABLE('6') = [13 26 52 104 208];
    LIFTING_SIZES_TABLE('7') = [15 30 60 120 240];    
    
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
% if tb input data empty, exit
if isempty(transport_block_bits_array)
    disp('Thanks for testing and using this interactive program written by Venkateswar Reddy Melachervu - Roll number:23156022');
    disp('More about me at http://www.linkedin.com/in/vmelachervu');
    return;
end
fprintf('\n');

% Prompt user to choose the 5G NR generator polynomial for TB CRC
[generator_poly, tb_crc_size] = read_user_selection_of_gp_for_crc_calculation();
% if generator poly is empty, exit
if isempty(generator_poly)
    disp('Thanks for testing and using this interactive program written by Venkateswar Reddy Melachervu - Roll number:23156022');
    disp('More about me at http://www.linkedin.com/in/vmelachervu');
    return;
end
fprintf('\n');

% Calculate CRC as per 5G NR standard
disp("Now computing the CRC as per 3GPP Standard TS38.212-f2...");
calculated5GNRCRCBits = calculate5GNRCRC(transport_block_bits_array,generator_poly);
fprintf('\n');
promptCRC = input("Press any key to continue...");

% Prompt user to choose the 5G NR generator polynomial for segmentation CRC GP
[seg_crc_generator_poly, seg_crc_length] = get_user_choice_of_segmentation_gp_crc_with_length();
if seg_crc_length == 0
    disp('Thanks for testing and using this interactive program written by Venkateswar Reddy Melachervu - Roll number:23156022');
    disp('More about me at http://www.linkedin.com/in/vmelachervu');
    return;
end
fprintf('\n');

% Compute B - TB length plus selected CRC length
B = length(transport_block_bits_array) + L;

% setting the selected base graph and code block length global variables per the base graph selection
%{
    Base Graph 1 : Optimized for large information block sizes and high code rates - when the SNR is good - K_bc = 8448
    Base Graph 2 : Optimized for smaller information block sizes and lower code rates - when the SNR is lower - K_bc = 3840
%}
base_graph_to_use = user_choice_of_base_graph();
if base_graph_to_use == 0
    disp('Thanks for testing and using this interactive program written by Venkateswar Reddy Melachervu - Roll number:23156022');
    disp('More about me at http://www.linkedin.com/in/vmelachervu');
    return;
end
K_cb = set_segmentation_code_block_length(base_graph_to_use);
fprintf('\n');

% compute total number of code block segments and effective payload
[C, B_prime] = needed_segmentation_code_blocks_and_effective_payload(B, L, K_cb);
fprintf('\n');

% effective paylod bits per code block - K' - K_prime
K_prime = B_prime/C;
Z_c_min = 0;
i_LS = '';
filler_bits = 0;
K_b = 0;

% find Z_c_min
switch base_graph_to_use
    case BASE_GRAPH_1_NUM    
        K_b = BASE_GRAPH_1('K_b_systematic_columns');
        disp(['The systematic/message columns value for the selected base graph is ' num2str(K_b) ' ...']);
        [Z_c_min, i_LS, filler_bits] = find_min_lifting_size_Z_c_and_filler_bits(LIFTING_SIZES_TABLE, K_b, K_prime);
    case BASE_GRAPH_2_NUM       
        K_b = BASE_GRAPH_2('K_b_systematic_columns');
        disp(['The systematic/message columns value for the selected base graph is ' num2str(K_b) ' ...']);
        [Z_c_min, i_LS, filler_bits] = find_min_lifting_size_Z_c_and_filler_bits(LIFTING_SIZES_TABLE, K_b, K_prime);
end 
if Z_c_min == 0
    disp('Unable to find a valid Z_c in the selected base graph! Exiting the system.');
    disp('Thanks for testing and using this interactive program written by Venkateswar Reddy Melachervu - Roll number:23156022');
    disp('More about me at http://www.linkedin.com/in/vmelachervu');
    return;
end
disp(['The Z_c minimum value found is : ' num2str(Z_c_min)]);
disp(['The lifting index i_ls is : ' i_LS]);
disp(['The K (K_b*Z_c_min) is : ' num2str(K_b*Z_c_min)]);
disp(['The K_prime (B_prime/C) is : ' num2str(B_prime/C)]);
disp(['The number of filler bits (K-K_prime) needed - prior to ceiis : ' num2str(filler_bits)]);
fprintf('\n');

% base graph matrix u needs to be transformed into parity check matrix H using lifting factor Z_c



disp('Thanks for testing and using this interactive program written by Venkateswar Reddy Melachervu - Roll number:23156022');
disp('More about me at http://www.linkedin.com/in/vmelachervu');

% ---------------------------------- End of the main program ---------------------------------------------------------------



% ---------------------------------- Begin of functions and utilities ------------------------------------------------------
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
        case BASE_GRAPH_2           
            % code block length to use
            seg_code_block_length = K_cb_base_graph_2_cb_length;
            disp('Selected base graph is : Base Graph 2');            
    end
    disp(['Segmentation code block length to be used, if needed : ' num2str(seg_code_block_length)]);
end


% Compute total number of segmentation code blocks from the avaliable user selection
function [needed_code_blocks, effective_pay_load] = needed_segmentation_code_blocks_and_effective_payload(b, l, k_cb)   
    disp(['Transport block length inclusive of the bits for the selected CRC (B) is : ' num2str(b)]);    
    if (b > k_cb)
        needed_code_blocks = ceil(b/(k_cb - l));
        effective_pay_load = b + (needed_code_blocks * l);
        disp(['Segmentation is needed - as B (TB + L) = ' num2str(b) ' is greater than K_cb = ' num2str(k_cb)]);
    else
        needed_code_blocks = 1;
        effective_pay_load = b;
        disp(['Segmentation is NOT needed - as B (TB + L) = ' num2str(b) ' is not greater than K_cb = ' num2str(k_cb)]);
    end
    disp(['Total number of code blocks needed (C) : ' num2str(needed_code_blocks)]);    
    disp(['Effective payload length (B_prime) is : ' num2str(effective_pay_load)]);    
    disp(['K_prime (B_prime/C) is : ' num2str(effective_pay_load/needed_code_blocks)]);    
end

% Function to prompt user for TB data
function transport_block_bits_array = user_prompt_menu_for_tb_data()
    transport_block_bits_array = [];
    % 3GPP standard constants
    lowerBound = 24;
    upperBound = 319784;   
    % Prompt user if he wants to provide transport block
    disp("Welcome to interactive 5G NR transceiver implementation (Matlab Asisgnment - 2 of Cohort-3 EE922) as per 3GPP Standard TS38.212-f20");
    disp("Would you like to type in input tranport block binary data for CRC generation and validation?");
    disp("1. Yes");
    disp("2. No - I want the program to auto generate input transport block data of a specific length");
    disp("3. No - I want the program to auto generate input transport block data between (inclusive) 24 and 3,19,784 bits length");
    disp("4. I am already feeling lucky - that is it and that is all ;-)");
    % Prompt the user for TB data input
    choice = input("Please enter your choice (1-3) : ");
    % Handle the user's choice
    switch choice
        case 1
            disp("You chose to type in the transport block binary data");                
            transport_block_bits_array = read_binary_data();
            if isempty(transport_block_bits_array)
                return;
            else
                bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
                disp('The input transport block binary data you have typed in is :');                        
                disp(bitsArrayAsString);
            end
        case 2
            disp("You chose to let the program auto generate input transport block data of a specific length");  
            tb_length = input("Please enter the TB bit length for auto generating TB data - between 24 and 3,19,784 : ");
            if ((tb_length < lowerBound) || (tb_length > upperBound))
                disp(['You have entered a length value which is illigitimate as per the 3GPP standard - ' num2str(tb_length)]);
                disp('Exiting the program. Please retry with legitimate value.');
                transport_block_bits_array = [];                
                return;
            end
            disp(['Auto generating transport block data of size ' num2str(tb_length) ' bits...']);                                              
            transport_block_bits_array = randi([0, 1], 1, tb_length);
            bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
            disp(['Auto generated transport block data is :']);                        
            disp(bitsArrayAsString);
        case 3
            disp("You chose to let the program auto generate random transport block binary data of random size.");  
            disp("Generating random transport block data of random size between 24 and 3,19,784 bits (inclusive) length...");          
            size = randi([lowerBound, upperBound]);               
            disp(['Randomly selected size for the input transmit block is:' num2str(size)]);       
            transport_block_bits_array = randi([0, 1], 1, size);
            bitsArrayAsString = sprintf('%d ', transport_block_bits_array);
            disp('Auto generated transport block data is :');   
            disp(bitsArrayAsString);
        case 4
            disp("You are feeling lucky already! That is it and that is all. Exiting the program...");
            return;        
        otherwise
            disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 3.");
            return;
    end
end

% Function to prompt user to choose the 5G NR generator polynomial for CRC
function [generatorPoly, selected_crc_size] = read_user_selection_of_gp_for_crc_calculation()
    selected_crc_size = 0;
    generatorPoly = [];
    disp('Please select the 5G NR generator polynomial you would like to use for the transport block CRC calculation:');
    disp("1.gCRC24A - 24-bit CRC - [D^24+D^23+D^18+D^17+D^14+D^11+D^10+D^7+D^6+D^5+D^4+D^3+D+1]");
    disp("2.gCRC24B - 24-bit CRC - [D^24+D^23+D^6+D^5+D+1]");
    disp("3.gCRC24C - 24-bit CRC - [D^24+D^23+D^21++D^20+D^17+D^15+D^13+D^12+D^8+D^4+D^2+D+1]");
    disp("4.gCRC16 - 16-bit CRC - [D^16+D^12+D^5+1]");
    disp("5.gCRC11 - 11-bit CRC - [D^11+D^10+D^9+D^5+1]");
    disp("6.gCRC6 - 6-bit CRC - [D^6+D^5+1]");
    disp("7.Exit");
    
    choice = input("Please enter your choice (1-7) : ");
    % Handle the user's choice
    switch choice
        case 1
            disp("You chose option 1 : gCRC24A - 24-bit cyclic generator polynomial for CRC");                         
            generatorPoly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 2
            disp("You chose option 2 : gCRC24B - 24-bit cyclic generator polynomial");
            generatorPoly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 3
            disp("You chose option 3 : gCRC24C - 24-bit cyclic generator polynomial");
            generatorPoly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 4
            disp("You chose option 4 : gCRC16 - 16-bit cyclic generator polynomial");
            generatorPoly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];   
            selected_crc_size = size(generatorPoly) - 1;
            generatorPolyString = sprintf('%d ', generatorPoly);
            disp(['Generator polynomial is: ' generatorPolyString]);
        case 5
            disp("You chose option 5 : gCRC11 - 11-bit cyclic generator polynomial");            
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

% Function to prompt user for selecting segmentation CRC polynomial
function [seg_crc_generator_poly, seg_crc_length] = get_user_choice_of_segmentation_gp_crc_with_length()
    seg_crc_length = 0;

    % Prompt user to choose the 5G NR segementation generator polynomial
    disp('Please select 5G NR generator polynomial you would like to use for segemented code blocks:');
    disp("1.gCRC24A - 24-bit CRC - [D^24+D^23+D^18+D^17+D^14+D^11+D^10+D^7+D^6+D^5+D^4+D^3+D+1]");
    disp("2.gCRC24B - 24-bit CRC - [D^24+D^23+D^6+D^5+D+1]");
    disp("3.gCRC24C - 24-bit CRC - [D^24+D^23+D^21++D^20+D^17+D^15+D^13+D^12+D^8+D^4+D^2+D+1]");
    disp("4.gCRC16 - 16-bit CRC - [D^16+D^12+D^5+1]");
    disp("5.gCRC11 - 11-bit CRC - [D^11+D^10+D^9+D^5+1]");
    disp("6.gCRC6 - 6-bit CRC - [D^6+D^5+1]");
    disp("7.Exit");
    
    choice = input("Please enter your choice (1-7): ");    
    % Handle the user's choice
    switch choice
        case 1
            disp("You chose option 1: gCRC24A - 24-bit cyclic generator polynomial for CRC");  
            disp("As per the 3GPP standard, this is an invalid choice of CRC CGP for TB segmentation. ");
            disp("Please re-run the program and select a valid choice for the CGP - one of {gCRC24B, gCRC24C, gCRC16}.");
            disp("Exiting the program.");
            return;           
        case 2
            disp("You chose option 2: gCRC24B - 24-bit cyclic generator polynomial");
            seg_crc_generator_poly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];   
            seg_crc_length = length(seg_crc_generator_poly) - 1;
            seg_crc_gp_poly_string = sprintf('%d ', seg_crc_generator_poly);
            disp(['Generator polynomial is: ' seg_crc_gp_poly_string]);
        case 3
            disp("You chose option 3: gCRC24C - 24-bit cyclic generator polynomial");
            seg_crc_generator_poly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];   
            seg_crc_length = length(seg_crc_generator_poly) - 1;
            seg_crc_gp_poly_string = sprintf('%d ', seg_crc_generator_poly);
            disp(['Generator polynomial is : ' seg_crc_gp_poly_string]);
        case 4
            disp("You chose option 4: gCRC16 - 16-bit cyclic generator polynomial");
            seg_crc_generator_poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];   
            seg_crc_length = length(seg_crc_generator_poly) - 1;
            seg_crc_gp_poly_string = sprintf('%d ', seg_crc_generator_poly);
            disp(['Generator polynomial is : ' seg_crc_gp_poly_string]);
        case 5
            disp("You chose option 5: gCRC11 - 11-bit cyclic generator polynomial");
            disp("As per the 3GPP standard, this is an invalid choice of CRC CGP for TB segmentation. ");
            disp("Please re-run the program and select a valid choice for the CGP - one of {gCRC24B, gCRC24C, gCRC16}.");
            disp("Exiting the program.");
            return;           
        case 6
            disp("You chose option 6: gCRC6 - 6-bit cyclic generator polynomial");
            disp("As per the 3GPP standard, this is an invalid choice of CRC CGP for TB segmentation. ");
            disp("Please re-run the program and select a valid choice for the CGP - one of {gCRC24B, gCRC24C, gCRC16}.");
            disp("Exiting the program.");
            return;           
        case 7
            disp("You chose to exit the program. Exiting...");
            return;        
        otherwise
            disp("Invalid choice. Exiting the program. Please re-run the program and type a number between 1 and 7.");
            return;
    end
end

% Function to read user choice for the base graph to be used for the segmentation process
function chosen_base_graph = user_choice_of_base_graph()    
    BASE_GRAPH_1 = 1;
    BASE_GRAPH_2 = 2;    
    chosen_base_graph = BASE_GRAPH_1;
    if chosen_base_graph > BASE_GRAPH_2
        disp('Invalid base graph chosen! Base graph needs to be between 1 and 2. Exiting the system. Please retry with legitimate value.')
        chosen_base_graph = 0;        
        return;
    elseif chosen_base_graph < BASE_GRAPH_1
        disp(['Invalid base graph chosen! Base graph needs to be between 1 and 2. Exiting the system. Please retry with legitimate value.'])
        chosen_base_graph = 0;                    
        return;
    end
end

% Function to find Z_c - minimum lifting size for the parity check matrix needed for LDPC
function [Z_c_min, i_LS, filler_bits] = find_min_lifting_size_Z_c_and_filler_bits(ls_table, K_b, K_prime)             
    % iterate over the map using for loop to find min Z_c
    disp('Iterating over Lifting size table to find minimum Z_c...');
    i_ls_indices = ls_table.keys;
    Z_c_rows = ls_table.values;
    Z_c_min = 0;    
    for i = 1:numel(i_ls_indices)        
        i_ls_current = i_ls_indices{i};
        disp(['Currently searching for legitimate Z_c_min in the row with the index ' num2str(i_ls_current) ' ...']);
        Z_c_current_row = Z_c_rows{i};
        for j = 1:numel(Z_c_current_row)
            Z_c_current = Z_c_current_row(j);             
            if (K_b * Z_c_current) >= K_prime
                disp(['Found a Z_c=' num2str(Z_c_current) ' where ' num2str(K_b) '*' num2str(Z_c_current) '=' num2str(K_b*Z_c_current) ' is greater than or equals to K_prime=' num2str(K_prime)]);
                % initialize Z_c_min with first hit
                if Z_c_min == 0                
                    Z_c_min = Z_c_current;            
                    i_LS = i_ls_current;
                    disp(['The initial Z_c_min hit value is: ' num2str(Z_c_min)]);
                end
                % check if this legitimate Z_c_current is lesser than previous hit
                if (Z_c_current < Z_c_min)
                    disp(['Found this Z_c=' num2str(Z_c_current) ' is less than previously found Z_c=' num2str(Z_c_min)]);
                    Z_c_min = Z_c_current;
                    i_LS = i_ls_current;
                end 
            end 
        end
    end         
    filler_bits = (K_b*Z_c_min) - K_prime;
end