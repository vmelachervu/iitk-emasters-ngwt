binaryArray = [1, 0, 1, 1, 0, 1];

% Randomly select the index of the bit to flip
indexToFlip = randi(numel(binaryArray));

% Flip the selected bit
binaryArray(indexToFlip) = ~binaryArray(indexToFlip);
disp(['flip bit index is:']);
disp(indexToFlip);

% Display the modified binary array
disp(binaryArray);

n = 10;  % Upper bound of the range
p = randperm(n);  % Generate a random permutation

disp(p);  % Display the random permutation

% Size of the bit array
n = 5;   % Number of bits to flip
m = 10;  % Number of additional bits

% Generate the bit array
bitArray = randi([0, 1], 1, 15);

% Set of indices to flip
%indicesToFlip = randperm(15, 5);
indicesToFlip = randperm(15, 15);

% Flip the bits at the specified indices
bitArray(indicesToFlip) = ~bitArray(indicesToFlip);

% Display the modified bit array
disp(bitArray);

n = 10;  % Upper bound of the range
fprintf('Flipping %d bits\n', n);
p = randperm(n);  % Generate a random permutation

disp(p);  % Display the random permutation