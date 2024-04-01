function info = getCBSInfo(B,bgn)
% getCBSInfo provides the information of code block segments
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   INFO = nr5g.internal.getCBSInfo(B,BGN) provides the information of the
%   code block segments based on the inputs, block length B and base graph
%   number BGN
%   1) Block length (B) should be a scalar nonnegative integer
%   2) Base graph number (BGN) should be either 1 or 2
%
%   INFO contains the following fields:
%   C   - Number of code block segments
%   CBZ - Number of bits in each code block (excluding CB-CRC bits and
%         filler bits)
%   Lcb - Number of parity bits in each code block
%   F   - Number of filler bits in each code block
%   K   - Number of bits in each code block (including CB-CRC bits and 
%         filler bits)
%   Zc  - Selected lifting size
%   Z   - Full lifting size set
%
%   Example:
%   % Get the information of code block segments for block length 3850
%   % and base graph number 2
%
%   info = nr5g.internal.getCBSInfo(3850,2)

%   Copyright 2018 The MathWorks, Inc.

%#codegen

    % Cast B to double, to make all the output fields have same data type
    B = cast(B,'double');

    % Get the maximum code block size
    if bgn == 1
      Kcb = 8448;
    else
      Kcb = 3840;
    end

    % Get number of code blocks and length of CB-CRC coded block
    if B <= Kcb
      L = 0;
      C = 1;
      Bd = B;
    else
      L = 24; % Length of the CRC bits attached to each code block
      C = ceil(B/(Kcb-L));
      Bd = B+C*L;
    end

    % Obtain the number of bits per code block (excluding CB-CRC bits)
    cbz = ceil(B/C);

    % Get number of bits in each code block (excluding filler bits)
    Kd = ceil(Bd/C);

    % Find the minimum value of Z in all sets of lifting sizes in 38.212
    % Table 5.3.2-1, denoted as Zc, such that Kb*Zc>=Kd
    if bgn == 1
      Kb = 22;
    else
      if B > 640
        Kb = 10;
      elseif B > 560
        Kb = 9;
      elseif B > 192
        Kb = 8;
      else
        Kb = 6;
      end
    end
    Zlist = [2:16 18:2:32 36:4:64 72:8:128 144:16:256 288:32:384];
    Zc =  min(Zlist(Kb*Zlist >= Kd));

    % Get number of bits (including <NULL> filler bits) to be input to the LDPC
    % encoder
    if bgn == 1
      K = 22*Zc;
    else
      K = 10*Zc;
    end

    info.C   = C;       % Number of code block segments
    info.CBZ = cbz;     % Number of bits in each code block (excluding CB-CRC bits and filler bits)
    info.Lcb = L;       % Number of parity bits in each code block
    info.F   = K-Kd;    % Number of filler bits in each code block
    info.K   = K;       % Number of bits in each code block (including CB-CRC bits and filler bits)
    info.Zc  = Zc;      % Selected lifting size
    info.Z   = Zlist;   % Full lifting size set

end