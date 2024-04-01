% ------------------------------------------------------------------------------------------------------%
%{
    Program name : Writing to a circular buffer
    eMasters - Communication Systems - Simulation-based Design of 5G NR Wireless Standard - EE910
    Assignment 3 - 5G Transceiver Implementation - II
    Roll number : 23156022
    Student Name : Venkateswar Reddy Melachervu    
    email : vmela23@iitk.ac.in
    History:
    V1.0.0  -   Initial complete solution - 18-06-2023        
      
    (C) Venkateswar Reddy Melachervu. 2023-2024.
%}
% ------------------------------------------------------------------------------------------------------%

clc;
buff_size = 10;
demod_cb = [1 0 1 0 1 0 1 0 1 1 0 0 0 0 0 ];
rv_id = 1;
disp(['Circular buffer size is:' num2str(buff_size)]);
disp('Demod bits are :');
disp(demod_cb);
disp(['Demod bits size is:' num2str(length(demod_cb))]);
disp(['RV offset value is: ' num2str(rv_id)]);


rate_recovered_cb = rate_recovery(buff_size, demod_cb, rv_id);
disp('Rate recovered bits are:');
disp(rate_recovered_cb);

function [rate_recovered_cb] = rate_recovery(circ_buff_size, demod_code_block, rv_id)
    % Number of code blocks    
    neutral_info = 0.5;

    % length of the decoded code block segment
    length_demod_code_block = length(demod_code_block);
    
    k = 0; % running index for rate recovery
    j = 0; % running index for computing the circular index    
    if (rv_id == -1 )
        rv_id = rv0;
    end
    k_0 = rv_id; % redundancy version offset

    % Initialization
    rate_recovered_cb = [];
    while k < circ_buff_size
        index = mod((k_0 + j), circ_buff_size);
        if k < length_demod_code_block
            rate_recovered_cb(index + 1) = demod_code_block( k + 1); %#ok<*AGROW>            
        else
            rate_recovered_cb(index + 1) = neutral_info; 
        end
        k = k + 1;
        j = j + 1;
    end
end