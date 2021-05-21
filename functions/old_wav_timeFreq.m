function [wav_power]= wav_timeFreq(data,fs,Min_freq,Max_freq)
% definitions, selections...
% chan2use = 'fcz';

min_freq =  Min_freq;
max_freq = Max_freq;
% min_freq =  1;
% max_freq = 100;
num_frex = round(max_freq-min_freq);

pnts= size(data,1);
trials= size(data,2);

len_time= size(data,1);
times= 1/fs:1/fs:len_time/fs;

% define wavelet parameters
time = -1:1/fs:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% define convolution parameters
n_wavelet            = length(time);
n_data               = pnts*trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
datfft = fft(reshape(data,1,pnts*trials),n_conv_pow2);

% initialize
datpower = zeros(num_frex,pnts); % frequencies X time X trials

baseidx = dsearchn(times',[1/fs 500/fs]');

% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    datconv = ifft(wavelet.*datfft);
    datconv = datconv(1:n_convolution);
    datconv = datconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(datconv,pnts,trials)).^2,2);
    datpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end
wav_power= datpower;

% % Checking visually
% figure
% subplot(121)
% contourf(times,frex,datpower,40,'linecolor','none')
% set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
% title('Logarithmic frequency scaling')
% 
% 
% subplot(122)
% contourf(times,frex,datpower,20,'linecolor','none')
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% colorbar


