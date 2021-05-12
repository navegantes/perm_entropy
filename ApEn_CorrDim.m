%% Using Unidimensional time series with varied Complexity Levels (time-dependent complexity) 
close all
clc
clear
tot_time= 3; % In seconds 

% Sampling parameters and total time
fs= 1000;
len_time= fs.*tot_time;

% Creating two sets of data, each with 500 trials 
data_1= wgn(3000,500,0); % generate random data points for data_1
data_2= wgn(3000,500,0); % generate random data points for data_2
% For each set a non linear data initiates at approx different narrow time
% window
% Setting data_1: nonlinearity starts earlier than data_2...
...between sample 1001 and 1020
for id1= 1:size(data_1,2)
    for i= ((len_time/2)-randi([1 20],1,1)):(0.85*len_time) % generate change of data complexity
        data_1(i,id1) = 0.05 + 0.6*data_1(i-1,id1) + 0.2*data_1(i-2,id1)-0.1*data_1(i-3,id1)+wgn(1,1,0);
    end
end
% Setting data_1: nonlinearity starts latter than data_1...
...between sample 1480 and 1499
 % generate random data points
for id2= 1:size(data_2,2)
    for i= ((len_time/3)+randi([1 20],1,1)):(0.85*len_time) % generate change of data complexity
        data_2(i,id2) = 0.05 + 0.6*data_2(i-1,id2) + 0.2*data_2(i-2,id2)-0.1*data_2(i-3,id2)+wgn(1,1,0);
    end
end

%% Checking time and frequency domain
% Time domain
rand_data1= data_1(:,randi([1 500],1,4));
rand_data2= data_2(:,randi([1 500],1,4));
time= 1/fs:1/fs:len_time/fs;
% Plotting time domain examples
figure(1)
for i=1:4
    subplot(4,1,i)
    plot(time,rand_data1(:,i));
    xlabel('Time (s)');
end

figure(2)
for i=1:4
    subplot(4,1,i)
    plot(time,rand_data2(:,i));
    xlabel('Time (s)');
end

% Frequency domain; Complex Morlet Wavelet Convolution
% frequency range parameters
min_freq= 1; max_freq= 100;
num_frex= round(max_freq-min_freq);
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
% calculating Morlet Wavelets time-frequency 
dat1_spec= wav_timeFreq(data_1,1000,min_freq,max_freq);
dat2_spec= wav_timeFreq(data_2,1000,min_freq,max_freq);

% Plotting frequency domain
figure(3)
subplot(1,2,1)
contourf(time,frex,dat1_spec,20,'linecolor','none')
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar
subplot(1,2,2)
contourf(time,frex,dat2_spec,20,'linecolor','none')
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar
%% Calculating permutation entropy for each data trial...
    ...not for mean data across trials
% Calculating permutatin entropy for signal and surrogates
% Permutation Entropy Parameters
delay = 1; % delay 1 between points in ordinal patterns (successive points)
order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
windowSize = 0.1*fs; % 100 miliseconds sliding window

% Calculating permutation entropy of data_1 for each trial and...
...averaging across trials
avPE_data1= cell(1,size(data_1,2)); 
avPE_data2= cell(1,size(data_1,2));
for i=1:size(data_1,2)
    avPE_data1{1,i} = PE(data_1(:,i),delay, order, windowSize)';
    avPE_data2{1,i} = PE(data_2(:,i),delay, order, windowSize)';
end
% mean(cell2mat(avPE_data1),2);
avPE_data1= cell2mat(avPE_data1);
avPE_data2= cell2mat(avPE_data2);

% Plotting both perm entropy data sets
figure(4)
plot(length(data_1)-length(mean(avPE_data1,2))+1:length(data_1), mean(avPE_data1,2), 'r', 'LineWidth', 0.2 );
hold on;
plot(length(data_2)-length(mean(avPE_data2,2))+1:length(data_2), mean(avPE_data2,2), 'b', 'LineWidth', 0.2 );
hold off;
legend('Data1','Data2');
title('Permutation entropy');

figure(5)
avPE_data1= avPE_data1';
avPE_data2= avPE_data2';
new_time= length(data_1)-length(mean(avPE_data1,1))+1:length(data_1);
subplot(1,2,1)
contourf(new_time,1:500,avPE_data1,20,'linecolor','none')
xlabel('Time (s)');
ylabel('Trials');
colorbar
subplot(1,2,2)
contourf(new_time,1:500,avPE_data2,20,'linecolor','none')
xlabel('Time (s)');
ylabel('Trials');
colorbar
%% Now calculating AppEn considering trials as dimesions
% 50 and 100 ms will be used
close all;
slwin1= 0.1*fs;
slwin2= 0.05*fs;
ApEn_TrialDim1= zeros(length(data_1)/slwin1,2);
ApEn_TrialDim2= zeros(length(data_1)/slwin2,2);
% For 100 ms slwin
count_ind1= [1 slwin1];
for i=1:length(ApEn_TrialDim1)
    ApEn_TrialDim1(i,1)= approximateEntropy(data_1(count_ind1(1):count_ind1(2),:));
    ApEn_TrialDim1(i,2)= approximateEntropy(data_2(count_ind1(1):count_ind1(2),:));
    count_ind1= count_ind1 + slwin1; 
end

count_ind2= [1 slwin2];
for i=1:length(ApEn_TrialDim2)
    ApEn_TrialDim2(i,1)= approximateEntropy(data_1(count_ind2(1):count_ind2(2),:));
    ApEn_TrialDim2(i,2)= approximateEntropy(data_2(count_ind2(1):count_ind2(2),:));
    count_ind2= count_ind2 + slwin2; 
end

% Plotting
% 100 ms slwin
figure(6)
plot(slwin1:slwin1:length(data_1), ApEn_TrialDim1(:,1),'b','linewidth',3); hold on;
plot(slwin1:slwin1:length(data_1), ApEn_TrialDim1(:,2),'r','linewidth',3); hold off;
legend('Data1', 'Data2');
title('Approximate Entropy for 100 ms slwin')

figure(7)
plot(slwin2:slwin2:length(data_1), ApEn_TrialDim2(:,1),'b','linewidth',3); hold on;
plot(slwin2:slwin2:length(data_1), ApEn_TrialDim2(:,2),'r','linewidth',3); hold off;
legend('Data1', 'Data2');
title('Approximate Entropy for 50 ms slwin')

% Doing the same but for one surrogate created for each trial
% Calculating Surrogates
surr1= zeros(size(data_1,1),size(data_1,2));
surr2= zeros(size(data_2,1),size(data_2,2));
for i=1:size(data_1,2)
    surr1(:,i) = IAAFTsur(data_1(:,i),1);
    surr2(:,i) = IAAFTsur(data_2(:,i),1);
end

ApEn_TrialDimSurr1= zeros(length(surr1)/slwin1,2);
ApEn_TrialDimSurr2= zeros(length(surr1)/slwin2,2);
% For 100 ms slwin
count_ind1= [1 slwin1];
for i=1:length(ApEn_TrialDimSurr1)
    ApEn_TrialDimSurr1(i,1)= approximateEntropy(surr1(count_ind1(1):count_ind1(2),:));
    ApEn_TrialDimSurr1(i,2)= approximateEntropy(surr2(count_ind1(1):count_ind1(2),:));
    count_ind1= count_ind1 + slwin1; 
end

count_ind2= [1 slwin2];
for i=1:length(ApEn_TrialDimSurr2)
    ApEn_TrialDimSurr2(i,1)= approximateEntropy(surr1(count_ind2(1):count_ind2(2),:));
    ApEn_TrialDimSurr2(i,2)= approximateEntropy(surr2(count_ind2(1):count_ind2(2),:));
    count_ind2= count_ind2 + slwin2; 
end

% Plotting
% 100 ms slwin
figure(8)
plot(slwin1:slwin1:length(data_1), ApEn_TrialDimSurr1(:,1),'b','linewidth',3); hold on;
plot(slwin1:slwin1:length(data_1), ApEn_TrialDimSurr1(:,2),'r','linewidth',3); hold off;
legend('Data1', 'Data2');
title('Approximate Entropy for 100 ms slwin: Surrogates')

figure(9)
plot(slwin2:slwin2:length(data_1), ApEn_TrialDimSurr2(:,1),'b','linewidth',3); hold on;
plot(slwin2:slwin2:length(data_1), ApEn_TrialDimSurr2(:,2),'r','linewidth',3); hold off;
legend('Data1', 'Data2');
title('Approximate Entropy for 50 ms slwin: : Surrogates')

%% Calculating the Correlation Dimension
% dim= size(data_1,2);
clc
close all
fs= 1000;
slwin= 0.05*fs;
% corDim1= zeros(length(data_1)/slwin,1);
% corDim2= zeros(length(data_2)/slwin,1);
corDim1= zeros(length(data_1)/slwin,size(data_1,2));
corDim2= zeros(length(data_2)/slwin,size(data_2,2));

keep_msg= cell((length(data_1)/slwin)*size(data_1,2),1);
count= 1;
for i= 1:size(data_1,2)
    count_ind= [1 slwin];
    while count_ind(1) < length(data_1)
        for ii= 1:length(data_1)/slwin
            [~,lag1,dim1] = phaseSpaceReconstruction(data_1(count_ind(1):count_ind(2),i));
            corDim1(ii,i)= correlationDimension(data_1(count_ind(1):count_ind(2),i)',lag1,dim1,'MaxRadius',2.5,'NumPoints',200);
            [~,lag2,dim2] = phaseSpaceReconstruction(data_2(count_ind(1):count_ind(2),i));
            corDim2(ii,i)= correlationDimension(data_2(count_ind(1):count_ind(2),i)',lag2,dim2,'MaxRadius',2.5,'NumPoints',200);
            keep_msg{count,1} = lastwarn;
            count= count+1;
            count_ind= count_ind + slwin;
        end
        disp(string(count_ind(1)) +" "+ string(length(data_1)));
    end
end

%%
%Plotting
figure(10)
plot(slwin:slwin:length(data_1),mean(corDim2,2));

contourf(slwin:slwin:length(data_1),1:500,corDim1',20,'linecolor','none')
xlabel('Time (s)');
ylabel('Trials');
colorbar

[XR,lag1,dim1] = phaseSpaceReconstruction(data_1(:,1));
teste= correlationDimension(data_1(count_ind(1):count_ind(2),i)',lag1,dim1,'MaxRadius',2.5,'NumPoints',200);
