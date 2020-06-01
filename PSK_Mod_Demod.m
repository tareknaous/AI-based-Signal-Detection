%>>>>>>>>> Training Data Generation of PSK Mod/Demod with errors caused by a noisy channel>>>>>>>%

clc;
clear all;
close all;

%Specify the number of samples the dataset will be composed of
SizeOfDataset = 5;

%Specify the value of the target SNR
SNR = -5;

%Initialize the matrix that will contain the transmitted (T) and received
%symbols (R)
T = [0 0 0 0 0 0 0];
R = [0 0 0 0 0 0 0];

%Begin iterating to generate data samples

for i=1:1:SizeOfDataset
    
%Randomize the bit sequence transmitted in each iteration    
    
t1 = round(rand(1,1));
t2 = round(rand(1,1));
t3 = round(rand(1,1));
t4 = round(rand(1,1));
t5 = round(rand(1,1));
t6 = round(rand(1,1));
t7 = round(rand(1,1));

%Transmitted Sequence
 
x=[ t1 t2 t3 t4 t5 t6 t7];                                    
T=[T;x];

% Specify the value of the bit period
bp=.000001;  

%Print the transmitted bit squence (optional)
%disp('Transmitted Bit Stream :')
%disp(x)

%>>>Representing the transmitted sequence as a digital signal (pulses)<<<%
bit=[]; 
for n=1:1:length(x)
    if x(n)==1;
       se=ones(1,100);
    else x(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];

end

%Visualize the transmitted digital signal
t1=bp/100:bp/100:100*length(x)*(bp/100);
subplot(3,1,1);
plot(t1,bit,'lineWidth',2.5);grid on;
axis([ 0 bp*length(x) -.5 1.5]);
ylabel('Amplitude (V)');
xlabel('Time (sec)');
title('Transmitted Digital Signal');



%>>> Binary-PSK modulation <<<%
A=5;                                          % Amplitude of carrier signal 
br=1/bp;                                                         %Bit rate
f=br*2;                                                 %Carrier frequency 
t2=bp/99:bp/99:bp;                 
ss=length(t2);
m=[];

for (i=1:1:length(x))
    if (x(i)==1)
        y=A*cos(2*pi*f*t2);
        
    else
        y=A*cos(2*pi*f*t2+pi);  
        
    end
    m=[m y];
end

%Obtained modulated waveform 
t3=bp/99:bp/99:bp*length(x);

%Addition of AWGN to the modulated waveform
t3_snr = awgn(t3,SNR,'measured');

%Visualization of the noisy waveform
subplot(3,1,2);
plot(t3_snr,m);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('Binary PSK modulated waveform');




%>>> Binary PSK demodulation <<<%
mn=[];
for n=ss:ss:length(m)
  t=bp/99:bp/99:bp;
  y=cos(2*pi*f*t);                                        % carrier signal 
  mm=y.*m((n-(ss-1)):n);
  
  %Demodulated waveform
  t4=bp/99:bp/99:bp;
  
  %Demodulated waveform with noise
  t4_snr=awgn(t4,SNR,'measured');
  z=trapz(t4_snr,mm);                                             
  zz=round((2*z/bp));                                     
  if(zz>0)                                       
                       
    a=1;
  else
    a=0;
  end
  mn=[mn a];
end

%Print received sequence (optional)
%disp(' Received Bit Stream :');
%disp(mn);

R=[R;mn];


%>>>Representing the received sequence as a digital signal (pulses)<<<%
bit=[];
for n=1:length(mn);
    if mn(n)==1;
       se=ones(1,100);
    else mn(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];

end

t4=bp/100:bp/100:100*length(mn)*(bp/100);
subplot(3,1,3);
plot(t4,bit,'LineWidth',2.5);grid on;
axis([ 0 bp*length(mn) -.5 1.5]);
ylabel('Amplitude (V)');
xlabel(' Time (sec)');
title('Received Signal after PSK Demodulation');



end


%Append Transmitted and Received sequence in one matrix that represents the
%dataset ready to be exported to excel
D = [ T R ];


%>>> END <<<%
