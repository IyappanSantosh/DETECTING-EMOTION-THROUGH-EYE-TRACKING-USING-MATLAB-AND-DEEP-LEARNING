
clc;
close all;
clear all;
%------------ Image Reading ------------------------------------------
delete Fdata.mat;
FilePathStorage=pwd;
FilePathname1=strcat(FilePathStorage,'/Fear/');
FileNameListArray=dir([FilePathname1 '*.jpg']);
FileNameListArray1=dir([FilePathname1 '*.mat']);
Totalfiles=size(FileNameListArray1,1);
Train_dataF=zeros(Totalfiles,18);
Train_LabF=zeros(2,Totalfiles);
Train_LabF1=zeros(Totalfiles,1);
for PPP=1:Totalfiles
TempArray=sprintf('Fear/%d.mat',PPP);
TempArray1=sprintf('Fear/%d.jpg',PPP);
load(TempArray);
DataImg = imread(TempArray1);
DataArray=imresize(DataImg,[150,200]);
%-------------Calculating Height and Width of picture----------------
[ImH,ImW,C]=size(DataArray);
disp('Image Heigt');
disp(ImH);
disp('Image Width');
disp(ImW);
disp('Color Channel');
disp(C);
rminiris=30;
rmaxiris=100;
if(C==3)
DataArray=rgb2gray(DataArray);
end
[ci,cp,o]=thresh(DataArray,rminiris,rmaxiris);
% figure,
% imshow(o)
% title('Iris and Pupil Boundary');
circleImage = false(ImW,ImH);
[x, y] = meshgrid(1:ImW, 1:ImH); 
circleImage((x - ci(2)).^2 + (y - ci(1)).^2 <= ci(3).^2) = true; 
% circleImage((x - cp(2)).^2 + (y - cp(1)).^2 <= cp(3).^2) = true; 
maskedImage = DataArray;
maskedImage(~circleImage) = 255;
circleImage1 = false(ImW,ImH);
circleImage1((x - cp(2)).^2 + (y - cp(1)).^2 <= cp(3).^2) = true; 
maskedImage1 = maskedImage;
maskedImage1(~circleImage1) = 255;
segment=zeros(ImH,ImW);
for ii=1:ImH
    for jj=1:ImW
        if(maskedImage1(ii,jj)<100)
            segment(ii,jj)=maskedImage(ii,jj);
        else
            segment(ii,jj)=255;
        end
    end
end
% figure,
% imshow(uint8(segment))
% title('Pupil Detection');
% 
Diameter=uint8(segment);
Diameter=double(reshape(Diameter,size(Diameter,1)*size(Diameter,2),1));
mean_P=mean(Diameter);
std_P=std(Diameter);
%DE Features
time=1:length(Diameter);
signal = Diameter;                                 % Signal data in Time-Domain
N = length(signal);                         % Number Of Samples
fs = [0.2 0.4 0.6 1];
DE=zeros(1,length(fs));
for i=1:length(fs)
freq =fs(i)*(-N/2:(N-1)/2);
FFT_Data = fftshift(fft(signal,N));
DE(1,i)=entropy(real(FFT_Data)).*fs(i);
end
Adata = disperse(Diameter);
Disp_mean=mean(Adata);
Disp_Std=std(Adata);
E_Feat=[mean_P std_P DE Disp_mean Disp_Std];
disp('Eye Feature');
disp(E_Feat)
ls = length(EEGdata);
% figure,
% plot(EEGdata);
% title('Input EEG Signal'),grid on
% pause(0.5);
%  DCT basis vector based 3-level multirate filterbank structure
% Level-1
fs=120;
Tw = 25;          
Ts = 10;
alpha = 0.97;
R = [1 60];  % frequency range to consider
M = 2;% number of filterbank channels 
N=2;
C = 5; % number of  coefficients
L = 10;
Nw = 600;
Ns = 1;
nfft = 2^nextpow2( Nw );     % length of FFT analysis 
K = (nfft/2)+1;               % length of the unique part of the FFT 
hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));
hz2mel = @( hz )( 1127*log(1+hz/700) );     % Hertz to mel warping function
mel2hz = @( mel )( 700*exp(mel/1127)-700 ); % mel to Hertz warping function
dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M) ...
                                       .* repmat(pi*([1:M]-0.5)/M,N,1) ) );
ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );
H = trifbank( M, K, R, fs, hz2mel, mel2hz ); 
EEG_data = filter( [1 -alpha], 1, EEGdata );
frames = vect2array( EEG_data, Nw, Ns, 'cols',hamming,false);
MAG = abs(fft(frames,nfft,2));
FBE = H*MAG(1:K,:); 
DCT = dctm(N,M);
CC=  DCT * log(FBE);
lifter = ceplifter( N, L );
D_Filt = diag( lifter ) * CC;
Fdata1=D_Filt(1,:);
Fdata2=D_Filt(2,:);
% figure,
% plot(Fdata1);
% title('0-30Hz Signal'),grid on
% pause(0.5);
% figure,
% plot(Fdata2);
% title('Gamma rhythm-30-60Hz'),grid on
% pause(0.5);
% Level-2
M = 5;% number of filterbank channels 
N=5;
Nw = 600;
Ns = 1;
nfft = 2^nextpow2( Nw );     % length of FFT analysis 
K = (nfft/2)+1;               % length of the unique part of the FFT 
hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));
hz2mel = @( hz )( 1127*log(1+hz/700) );     % Hertz to mel warping function
mel2hz = @( mel )( 700*exp(mel/1127)-700 ); % mel to Hertz warping function
dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M) ...
                                       .* repmat(pi*([1:M]-0.5)/M,N,1) ) );
ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );
H = trifbank( M, K, R, fs, hz2mel, mel2hz ); 
EEG_data = filter( [1 -alpha], 1, Fdata1);
frames = vect2array(EEG_data, Nw, Ns, 'cols',hamming,false);
MAG = abs(fft(frames,nfft,2));
FBE = H*MAG(1:K,:); 
DCT = dctm(N,M);
CC=  DCT * log(FBE);
lifter = ceplifter( N, L );
D_Filt = diag( lifter ) * CC;
Fdata3=D_Filt(1,:)+D_Filt(2,:);
Fdata4=D_Filt(3,:)+D_Filt(4,:)+D_Filt(5,:);
% figure,
% plot(Fdata3);
% title('0-12Hz Signal'),grid on
% pause(0.5);
% figure,
% plot(Fdata4);
% title('Beta rhythm-12-30Hz'),grid on
% pause(0.5);
% Level-3
M = 3;% number of filterbank channels 
N=3;
Nw = 600;
Ns = 1;
nfft = 2^nextpow2( Nw );     % length of FFT analysis 
K = (nfft/2)+1;               % length of the unique part of the FFT 
hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));
hz2mel = @( hz )( 1127*log(1+hz/700) );     % Hertz to mel warping function
mel2hz = @( mel )( 700*exp(mel/1127)-700 ); % mel to Hertz warping function
dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M) ...
                                       .* repmat(pi*([1:M]-0.5)/M,N,1) ) );
ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );
H = trifbank( M, K, R, fs, hz2mel, mel2hz ); 
EEG_data = filter( [1 -alpha], 1, Fdata3);
frames = vect2array(EEG_data, Nw, Ns, 'cols',hamming,false);
MAG = abs(fft(frames,nfft,2));
FBE = H*MAG(1:K,:); 
DCT = dctm(N,M);
CC=  DCT * log(FBE);
lifter = ceplifter( N, L );
D_Filt = diag( lifter ) * CC;
Fdata5=D_Filt(1,:);
Fdata6=D_Filt(2,:);
Fdata7=D_Filt(3,:);
% figure,
% plot(Fdata5);
% title('delta rhythm 0-4Hz'),grid on
% pause(0.5);
% figure,
% plot(Fdata6);
% title('theta rhythm 4-8Hz'),grid on
% pause(0.5);
% figure,
% plot(Fdata7);
% title('alpha rhythm 8-12Hz'),grid on
% pause(0.5);
% Feature Extraction
% PSD
delta=Fdata5;
theta=Fdata6;
alpha=Fdata7;
beta=Fdata4;
gamma=Fdata2;
Fs=1000;
[psdx_delta,freq_delta]=PSD(delta,Fs);
[psdx_theta,freq_theta]=PSD(theta,Fs);
[psdx_alpha,freq_alpha]=PSD(alpha,Fs);
[psdx_beta,freq_beta]=PSD(beta,Fs);
[psdx_gamma,freq_gamma]=PSD(gamma,Fs);
% figure,
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% subplot(2,3,1)
% plot(freq_delta,10*log10(psdx_delta))
% grid on
% title('Periodogram Using FFT(delta)')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% subplot(2,3,2)
% plot(freq_theta,10*log10(psdx_theta))
% grid on
% title('Periodogram Using FFT(theta)')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% subplot(2,3,3)
% plot(freq_alpha,10*log10(psdx_alpha))
% grid on
% title('Periodogram Using FFT(alpha)')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% subplot(2,3,4)
% plot(freq_beta,10*log10(psdx_beta))
% grid on
% title('Periodogram Using FFT(beta)')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% subplot(2,3,5)
% plot(freq_gamma,10*log10(psdx_gamma))
% grid on
% title('Periodogram Using FFT(beta)')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% hold on
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 1,'\bf PSD feature','HorizontalAlignment','center','VerticalAlignment', 'top');
% pause(0.5);
psdx_delta=mean(psdx_delta);
psdx_theta=mean(psdx_theta);
psdx_alpha=mean(psdx_alpha);
psdx_beta=mean(psdx_beta);
psdx_gamma=mean(psdx_gamma);
PSD_Feat=[psdx_delta psdx_theta psdx_alpha psdx_beta psdx_gamma];
disp('PSD Feature')
disp(PSD_Feat)
% DE Feature
%DE Features
signal=zeros(1,length(delta));
FFT_Data=zeros(1,length(delta));
signal(1,:)=delta;
signal(2,:)=theta;
signal(3,:)=alpha;
signal(4,:)=beta;
signal(5,:)=gamma;
fs=100;
time=1:length(delta);        
N = length(delta);                        
DE=zeros(1,5);
for i=1:5
freq =fs*(-N/2:(N-1)/2);
FFT_Data(i,:) = fftshift(fft(signal(i,:),N));
DE(1,i)=entropy(real(FFT_Data(i,:)));
end
% figure,
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% subplot(2,3,1)
% plot(real(FFT_Data(1,:)));
% grid on
% title('STFT(delta)')
% xlabel('Frequency (Hz)')
% ylabel('STFT')
% subplot(2,3,2)
% plot(real(FFT_Data(2,:)));
% grid on
% title('STFT(theta)')
% xlabel('Frequency (Hz)')
% ylabel('STFT')
% subplot(2,3,3)
% plot(real(FFT_Data(3,:)));
% grid on
% title('STFT(alpha)')
% xlabel('Frequency (Hz)')
% ylabel('STFT')
% subplot(2,3,4)
% plot(real(FFT_Data(4,:)));
% grid on
% title('STFT(beta)')
% xlabel('Frequency (Hz)')
% ylabel('STFT')
% subplot(2,3,5)
% plot(real(FFT_Data(5,:)));
% grid on
% title('STFT(gamma)')
% xlabel('Frequency (Hz)')
% hold on
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 1,'\bf Short-term Fourier transforms(STFT)','HorizontalAlignment','center','VerticalAlignment', 'top');
% ylabel('STFT')
disp('DE Features')
disp(DE)
EEG_Feat=[PSD_Feat DE];
disp('EEG Features')
disp(EEG_Feat)
% 
T_feat=[EEG_Feat E_Feat];
disp('Final Feature')
disp(T_feat)
Train_dataF(PPP,:)=T_feat;
Train_LabF1(PPP,1)=4;
Train_LabF(1,PPP)=0;
Train_LabF(2,PPP)=0;
Train_LabF(3,PPP)=0;
Train_LabF(4,PPP)=1;
end
save Fdata.mat Train_dataF Train_LabF Train_LabF1