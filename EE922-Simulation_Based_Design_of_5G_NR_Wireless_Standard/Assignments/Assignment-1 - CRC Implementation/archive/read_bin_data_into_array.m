% Example usage
bits = read_binary_data();
disp('Bit array:');
disp(bits);
function bit_array = read_binary_data()
    binary_string = input('Enter binary data: ', 's');
    bit_array = [];

    for i = 1:length(binary_string)
        if binary_string(i) == '0'
            bit_array = [bit_array, 0];
        elseif binary_string(i) == '1'
            bit_array = [bit_array, 1];
        else
            fprintf('Invalid binary data. Exiting.\n');
            bit_array = [];
            return;
        end
    end
end

