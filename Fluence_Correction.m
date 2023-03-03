%%%% loading image prebeamformed

%% load images and model fluence 
%{
fnameBase = {'bld10ml+gel20ml+int5ml Tran1 T2.iq';};
fnameBase = char(fnameBase(1));
fname = strcat(fnameBase, '.pamode');
fnameXml = strcat(fnameBase, '.xml');
% This determines how many quadrants are used for full-width
% Apperture
AS = 0; % zero-based
AE = 255; % zero-based
% MS250/LZ250 
a = 0.25e-3; % m
pitch = 90e-6; %m
iframe = 3;
% Speed of sound: ct=transverse, cl=longitudinal
ct = 1540; %m/s
cl = 2340; %m/s
% Set the region of the dataset to beamform (for speed)
% Defined as fractional percentage of the total number of lines or samples
SL = 0.0; 
EL = 1.0;          
SS = 0.01;     
ES = 1.0;       

[param, IntFac, fs_int, f_rf, IdataRo, QdataRo] = PA_Open_IQ_File_EH(fname, fnameXml, pitch, AS, AE, iframe);
%%%%%%%%%%%%%%%%%%%%% RECONSTRUCT RF DATA AND BEAMFORM %%%%%%%%%%%%%%%%%%%%
[BfData_eno, RfData] = PA_BF_EH(IdataRo, QdataRo, param, IntFac, fs_int, f_rf, pitch, a, SS, ES, SL, EL, ct);

%sig1 = real(BfData_eno(X(i):(X(i)+P-1),190:215,3));
signal = real(BfData_eno(400:end,:));
%signal = RfData(300:end,:);
N_max = length(signal);
fs = 672e6;
signal1 = signal.*exp(-0.*(1:N_max)./fs.*1540000)';
[BfData_eno, RfData] = PA_BF_EH_A(IdataRo, QdataRo, param, IntFac, fs_int, f_rf, pitch, a, SS, ES, SL, EL, ct);
signal = real(BfData_eno(400:end,:));
signal2 = signal;
%}

%% load the images 
signal1 = real(BfData_eno(1:end,1:255,8));
signal2 = real(BfData_eno(1:end,1:255,18));

% load('BF_phantom7002.mat')
% signal1 = RF(:,1:255);
% load('BF_phantom9002.mat')
% signal2 = RF(:,1:255);
N_max = length(signal1);

S1 = max(abs(real(signal1(300:1500,:))));
S2 = max(abs(real(signal2(300:1500,:))));
signal1 = signal1./S1;
signal2 = signal2./S2;

window_size = [700];
clear SO SOC SOC1 SOC2
SO = zeros(234,1);
fs = 672e6;

%% calculate and compute the fluece difference

h = waitbar(0,'Initializing waitbar...');
for fil = 1: length(window_size)        %% testing for different window sizes
    clear Sl 
    n_window = window_size(fil);
    %overlap_win = round(n_window/45);
    overlap_win = 20;
    X = 1:overlap_win:(N_max - n_window);
    for i = 1:length(X)
        waitbar(i./length(X),h,sprintf('Fluence Correcting %d%% along...',round(100.*i./length(X))))
        sig1 = real(signal1(X(i):(X(i)+n_window-1),:));
        sig2 = real(signal2(X(i):(X(i)+n_window-1),:));

        N = 5000;        k = 0:N-1;                %create a vector from 0 to N-1
        T = N/fs;                 %get the frequency interval
        freq = k/T;          fsh1 = fft(sig1,N);        fsh2 = fft(sig2,N);
        f1 = fft(sig1,N);        f2 = fft(sig2,N);
        w = 20;
        S = find(freq<=10e6,1,'last');
        E = find(freq>=30e6,1,'first');
        for ii = 1:234
            x = ((freq(S:E)')\(mean(20*log10(abs(fsh2(S:E,(ii):(ii+w)))),2)-mean(20*log10(abs(fsh1(S:E,(ii):(ii+w)))),2)));
%             x = polyfit(freq(S:E)',mean(20*log10(abs(fsh2(S:E,(ii):(ii+w)))),2)-mean(20*log10(abs(fsh1(S:E,(ii):(ii+w)))),2),1);
            SO(ii,i+1) = x(1);
        end
        SOC1 = (cumsum(SO(:,1:(end-1)),2).*(overlap_win/n_window));
        SOC(:,i) = (SOC1(:,end) + SO(:,end)).*(1 + 1i);
        filter2 = 10.^((freq(S:E))'*SOC(:,i)'./20);
        %filter2 = filter./filter(1,:).*(1 + 1i);  
        LTF2 = [repmat(filter2(1,:),S-1,1); filter2; repmat(filter2(end,:),N-(E*2-1),1); flip(filter2); repmat(filter2(1,:),S-2,1)] ;

        sig1c = ifft(f1(:,11:244).*LTF2,N);
%         sig1O(X(i):(X(i)+n_window-1),:,i) = sig1(1:n_window,:);
        sig1C(X(i):(X(i)+n_window-1),:,i) = sig1c(1:n_window,:);
    end
%     RefO = (sum(sum(reshape((20*log10(abs(sig1O)) > -500),[],size(sig1O,2),i),3),2)./(size(sig1O,2)-1));
    RefC = (sum(sum(reshape((20*log10(abs(sig1C)) > -500),[],size(sig1C,2),i),3),2)./(size(sig1C,2)));
%     signal1O = sum(sig1O,3)./RefC;
    signal1C = sum(sig1C,3)./RefC;
end
close(h)