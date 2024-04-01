% ------------------------------------------------------------------------------------------------------%
%{
    Library name : Library for QPSK Modulator and Demodulator with Plots
    eMasters - Communication Systems - Simulation-based Design of 5G NR
    Wireless Standard - EE92    
    Roll number : 23156022
    Student Name : Venkateswar Reddy Melachervu    
    email : vmela23@iitk.ac.in

    History:
    V1.0.0  -   Initial complete solution - 16-06-2023        
    (C) Venkateswar Reddy Melachervu. 2023-2024.
%}
% ------------------------------------------------------------------------------------------------------%

% test code
% data = randi([0 1], 1, 1024);
data=[0 1 0 1 1 1 0 0 1 1 0 1 0 1 0 1 0 0 1 1 1 1 1 0 0 0]; % information
num_of_bits = length(data);

figure(1)
stem(data, 'linewidth',3), grid on;
title('Input Data to QPSK Modulator');
% axis([ 0 11 0 1.5]);
axis([ 0 num_of_bits 0 1.5]);

% qpsk modulation
[Tx_sig, tt, y_in, y_qd, f, t, T] = qpsk_modulate(data);
figure(2)

subplot(3,1,1);
plot(tt,y_in,'linewidth',3), grid on;
title(' QPSK Modulated Waveform - Inphase Waveform Component');
xlabel('Time - Sec');
ylabel('Amplitude');

subplot(3,1,2);
plot(tt,y_qd,'linewidth',3), grid on;
title(' QPSK Modulated Waveform - Quadrature Waveform Component');
xlabel('Time - Sec');
ylabel('Amplitude');

subplot(3,1,3);
plot(tt,Tx_sig,'r','linewidth',3), grid on;
title('QPSK Modulated Signal - Sum of Inphase and Quadrature Phase Signals)');
xlabel('Time - Sec');
ylabel('Apmplitude');

% qpsk - demodulation
[Rx_data] = qpsk_demod(Tx_sig, data, f, t, T);

figure(3)
stem(Rx_data,'linewidth',3) 
title('Demodulated Data');
% axis([ 0 11 0 1.5]), grid on;
axis([ 0 num_of_bits 0 1.5]), grid on;


% function - qpsk modulation
function [Tx_sig, tt, y_in, y_qd, f, t, T] = qpsk_modulate(data)    
    % initialization
    data_NZR=2*data-1; % Data Represented at NZR form for QPSK modulation
    s_p_data=reshape(data_NZR,2,length(data)/2);  % S/P convertion of data       
    br=10.^6; %transmission bit rate  1000000
    f=br; % minimum carrier frequency
    T=1/br; % bit duration
    t=T/99:T/99:T; % Time vector for one bit information     
    
    % qpsk modulation
    y=[];
    y_in=[];
    y_qd=[];
    for(i=1:length(data)/2)
        y1=s_p_data(1,i)*cos(2*pi*f*t); % inphase component
        y2=s_p_data(2,i)*sin(2*pi*f*t) ;% Quadrature component
        y_in=[y_in y1]; % inphase signal vector
        y_qd=[y_qd y2]; %quadrature signal vector
        y=[y y1+y2]; % modulated signal vector
    end
    Tx_sig=y; % transmitting signal after modulation
    tt=T/99:T/99:(T*length(data))/2;
end 

% qpsk demodulation
function [Rx_data] = qpsk_demod(Tx_sig, data, f, t, T)    
    Rx_data=[];
    disp(['Source data length is:' num2str(length(data))]);
    disp(['Tx_sig data length is:' num2str(length(Tx_sig)/2)]);
    Rx_sig=Tx_sig; % Received signal
    for(i=1:1:length(data)/2)    
        % inphase coherent dector
        Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)).*cos(2*pi*f*t); 
        % above line indicate multiplication of received & inphase carred signal        
        Z_in_intg=(trapz(t,Z_in))*(2/T);% integration using trapizodial rull

        if(Z_in_intg>0) % Decession maker
            Rx_in_data=1;
        else
           Rx_in_data=0; 
        end
        
        % quadrature coherent dector
        Z_qd=Rx_sig((i-1)*length(t)+1:i*length(t)).*sin(2*pi*f*t);
        %above line indicat multiplication ofreceived & Quadphase carred signal
        
        Z_qd_intg=(trapz(t,Z_qd))*(2/T); %integration using trapizodial rull
        if (Z_qd_intg>0) % Decession Maker
            Rx_qd_data=1;
        else
            Rx_qd_data=0; 
        end                       
        Rx_data=[Rx_data  Rx_in_data  Rx_qd_data]; % Received Data vector
    end
end
    

