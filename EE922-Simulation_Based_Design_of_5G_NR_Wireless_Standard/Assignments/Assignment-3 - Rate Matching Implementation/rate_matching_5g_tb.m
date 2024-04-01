%% 5G LDPC Block Error Rate Simulation Using the Cloud or a Cluster
% This example shows how to use the cloud or a cluster for block error rate
% (BLER) simulation of low-density parity-check (LDPC) coding for the 5G NR
% downlink shared transport channel (DL-SCH).
%
% For ultra-reliable low-latency communication in 5G systems, some use cases
% require a BLER as low as 1e-6 <#12 [1]>. In this low BLER region, obtaining
% accurate results requires simulating many millions of blocks. On a single
% desktop computer, this simulation can take days. You can use the cloud or
% a cluster to reduce simulation time. For example, a 256-worker cluster on
% AWS is about 42 times as fast as a 6-core desktop in one test scenario.
% For more details, see <#7 Sample Results>.
%
% In this example, you first calculate one point on a BLER curve by using a
% single desktop computer. You then use MATLAB Parallel Server in the cloud
% or on a cluster on your local network to calculate the BLER curve across a
% range of signal-to-noise ratios.

% Copyright 2022 The MathWorks, Inc.

%% DL-SCH with LDPC Coding
% First, simulate one transport block for the 5G NR DL-SCH with LDPC
% coding. This code is the basis of the |helperLDPCBLERSim| function which
% uses |parfor| to simulate many transport blocks in parallel.

% Set up DL-SCH coding parameters
TBS = 3816;            % Transport block size, a positive integer
codeRate = 308/1024;   % Target code rate, a real number between 0 and 1
rv = 0;                % Redundancy version, 0-3
modulation = 'QPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1;           % Number of layers, 1-4 for a transport block
cbsInfo = nrDLSCHInfo(TBS,codeRate);
disp('DL-SCH coding parameters')
disp(cbsInfo)

switch modulation
    case 'QPSK'
        bitsPerSymbol = 2;
    case '16QAM'
        bitsPerSymbol = 4;
    case '64QAM'
        bitsPerSymbol = 6;
    case '256QAM'
        bitsPerSymbol = 8;
end

% Set up AWGN channel
EbNo = 1.25; % in dB
outlen = ceil(TBS/codeRate);
snrdB = convertSNR(EbNo,"ebno",...
    BitsPerSymbol=bitsPerSymbol,CodingRate=TBS/outlen);

% Random transport block data generation
in = randi([0 1],TBS,1,'int8');
% Transport block CRC attachment
tbIn = nrCRCEncode(in,cbsInfo.CRC);
% Code block segmentation and CRC attachment
cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
% LDPC encoding
enc = nrLDPCEncode(cbsIn,cbsInfo.BGN);
% Rate matching and code block concatenation
chIn = nrRateMatchLDPC(enc,outlen,rv,modulation,nlayers);
% Symbol mapping
symOut = nrSymbolModulate(chIn,modulation);
% AWGN channel
[rxSig, noiseVar] = awgn(symOut,snrdB);
% Symbol demapping
rxllr = nrSymbolDemodulate(rxSig,modulation,noiseVar);
% Rate recovery
raterec = nrRateRecoverLDPC(rxllr,TBS,codeRate,rv,modulation,nlayers);
% LDPC decoding, with early termination and at most 12 iterations
decBits = nrLDPCDecode(raterec,cbsInfo.BGN,12);
% Code block desegmentation and CRC decoding
[blk,~] = nrCodeBlockDesegmentLDPC(decBits,cbsInfo.BGN,TBS+cbsInfo.L);
% Transport block CRC decoding
[out,~] = nrCRCDecode(blk,cbsInfo.CRC);
% Compare
blockError = any(out~=in)