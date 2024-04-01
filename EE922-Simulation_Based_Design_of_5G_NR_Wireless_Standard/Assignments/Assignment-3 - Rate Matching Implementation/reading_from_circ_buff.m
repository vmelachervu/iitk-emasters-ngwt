E = 12;
ldpc_code_block = [1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 -1 -1 -1 -1 -1];
rv_id = 19;

rate_matched_cb = rate_matcher_v2(E, ldpc_code_block,rv_id);
E_string = sprintf('%d', rate_matched_cb);
disp(['Rate matched CB is:' E_string]);

% function rate matcher
function [rate_matched_cb] = rate_matcher_v2(E, ldpc_code_block, rv_id)      
    rate_matched_cb = [];
    % length of the ldpc encoded block segment
    length_ldpc_code_block = length(ldpc_code_block);
    
    % for skipping copying of filler bit
    NULL = -1;
                 
    k_0 = rv_id; % redundancy version offset
    k = 0; % running index for e_k - the rate matched outout array - rate_matched_cb
    j = 0; % running index for d - the ldpc code block and input to rate matching - ldpc_code_block

    % populate rate matched code block
    % If transmittable bits length is less than ldpc encode output - need
    % to drop/puncture ldpc encode output to match the transmittable
    % rate matched block length - E with offset RV0 = 0
    % If transmittable bits length is more than ldpc encoded output - need
    % to add redundant bits with offset RV0=0
    while k < E
        index = mod((k_0 + j), length_ldpc_code_block);       
        rate_matched_cb(k + 1) = ldpc_code_block(index + 1); 
        k = k + 1;        
        j = j + 1;
    end
end