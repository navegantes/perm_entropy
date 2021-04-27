
chdir 'D:\Users\Sources\Entropy\perm_entropy';
%% -----------------------------------------------------------------

close all;
clear;
clc;

addpath('PE');
addpath('functions')

tmin = 10;
tmax = 1930;
RWD_label = 'RWD-250';
ev_duration = '.250<=2'; %(s)

% ------------------------------------------------------------------------
% Read from csv file
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
csv_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.csv';
edf_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.edf';

file = csvread(csv_file, 2, 0);
file = file';

rawCSV = pop_importdata( 'setname', 'JLC-210120_S2-csv', ...
                         'data', file, ...
                         'nbchan', 13, ...
                         'dataformat', 'array', ...
                         'srate', 256);
                     
EEGcsv = pop_select( rawCSV, 'time',[tmin+1 tmax+1] );

%Read from edf file
EEG = pop_biosig(edf_file, 'importevent','off');
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', 'JLC-210120_S2-edf', 'gui','off');
EEG = pop_select( EEG,'time',[tmin tmax] );
% EEG = eeg_checkset( EEG );
% eeglab redraw;

% muda os labels dos canais
labels = {'EEG','ECG'};
for i=1:2
    EEG.chanlocs(i).labels = labels{i};
end

% Fitragem passa-faixa.
EEG = pop_eegfiltnew(EEG, 1,100,900,0,[],0);

% ------------------------------------------------------------------------
% Adiciona os canais (11, 2, 5, 8) RWD e Amplitudes theta SMR hibeta, do csv.
labels = {'RWD' 'Theta:Amp' 'SMR:Amp' 'Hibeta:Amp' 'Theta:Power'...
          'SMR:Power' 'Hibeta:Power' 'Alpha:Power'};
chanlist = [11 2 5 8];
for chan=1:length(chanlist)
    EEG.data(end+1,:) = EEGcsv.data(chanlist(chan), :);
    EEG.nbchan = EEG.nbchan + 1;
    EEG.chanlocs(chan+2).labels = labels{chan};
end

% ------------------------------------------------------------------------
% Adiciona potencias instantaneas theta, SMR, hbeta, alpha

fminmax = {[4 7] [12 15] [20 30] [8 12]};
lenfmM = length(fminmax);
for f=1:length(fminmax)
    freq = fminmax{f};
    EEGfilt = pop_eegfiltnew(EEG, freq(1),freq(2),900,0,[],0);
    hilb = hilbert(EEGfilt.data(1, :));
    EEG.data(end+1,:) = abs(hilb).^2;
    EEG.nbchan = EEG.nbchan + 1;
    EEG.chanlocs(f+6).labels = labels{f+lenfmM};
end

%% -----------------------------------------------------------------
% Calcula entropia permutação e gera surrogate, inclui no objeto EEG
[Hdatas, EEG] = fntools.hsurdata(EEG);

%% -----------------------------------------------------------------
% Cria os eventos e extrai trials        % EEGnfb = pop_select( EEG,'time',[60 1860] );
ev_range = [-.3 .5];
rwdchan = 3;
latency = '60<=1860';
EEGev = fntools.create_events(EEG, rwdchan, latency, RWD_label, ev_duration, ev_range);

%% -----------------------------------------------------------------
% -------------------------------------------------------------------------
% VISUALIZACAO
% -------------------------------------------------------------------------
% Mostra o espectro (confirmar resultado do feedback)
f_range = [1 50];
data = permute(EEGev.data(1,:,:), [2 3 1]);
[time, frex] = fntools.show_evspec(data, EEGev.srate, ev_range, f_range);

%% -----------------------------------------------------------------
[~] = fntools.show_erp_events(EEGev, RWD_label);

%% -----------------------------------------------------------------
Hperm = fntools.evpermentropy(EEGev);

%% -----------------------------------------------------------------
wess_dist = zeros(1,size(Hperm,1));
for i=1:size(Hperm,1)
    wess_dist(i) = ws_distance(EEGev.data(12,:,i), EEGev.data(13,:,i),1);
end

%% -----------------------------------------------------------------
figure;
hold on;
plot(Hdatas.permH(1:end-5*256-2));
plot(Hdatas.surpermH(1:end-5*256-2));
hold off;

%% -----------------------------------------------------------------
figure;
hold on;
for i=1:20
    plot(Hperm(i,1:size(Hperm,2)-round(256/12)-2));
end
hold off;

%% -----------------------------------------------------------------
figure;
hold on;
[~, len, num] = size(EEGev.data);
for i=1:20
    plot(EEGev.data(13,:, i))
end
hold off;
%% -----------------------------------------------------------------
delay = 1; % delay 1 between points in ordinal patterns (successive points)
order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
windowSize = 4*EEG.srate;
H_perm = zeros(1, size(EEG.data(1,:), 2));
h = PE(EEG.data(1,:)', delay, order, windowSize);
H_perm(1:length(h)) = h;
%Computing Entropy
EEG.data(end+1,:) = H_perm; 

EEG.nbchan = size(EEG.data,1);
% [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% eeglab redraw;
pop_eegplot( EEG, 1, 1, 1);

%% -----------------------------------------------------------------
% define os eventos  disp([newline 'Def events']);

% EEG = pop_chanevent( EEG, 3, 'edge', 'leading', 'duration', 'on',...
%                      'typename', 'RWD', 'delchan', 'off', 'edgelen', 1);
% 
% % Seleciona eventos duracao maior que 250 ms
% EEG = pop_selectevent( EEG, 'latency','60<=1860','duration',ev_duration,...
%                        'renametype',RWD_label,'oldtypefield','RWD','deleteevents','on');
% EEG = pop_epoch( EEG, { RWD_label }, [-0.3 0.5], ...
%                  'newname', 'JLC-210120_S2-epochs', 'epochinfo', 'yes');
             
%% -----------------------------------------------------------------
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
% EEG = eeg_checkset( EEG );
% eeglab redraw;
% 
% sz = size(EEG.data);
% hdatas   = cell(sz(3));
% delay = 1; % delay 1 between points in ordinal patterns (successive points)
% order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
% windowSize = 4*EEG.srate;

% for i=1:
% hdatas = PE(sbj_dt(EEG.data(1,:, :)', delay, order, windowSize);

%% -----------------------------------------------------------------

% EEG = fn_tools.show_erp_events(EEG, RWD_label);

% function EEG = show_erp_events(EEG)
% 
%     RWD_label = 'RWD-250';
%     pop_eegplot( EEG, 1, 1, 1);
% 
%     smooth = 5;
%     chans = [1 4 5 6 7 8 9 10];
%     band_labels = {'EEG' 'Theta Amp' 'SMR Amp' 'Hibeta Amp' 'Theta Power'...
%                    'SMR Power' 'Hibeta Power' 'Alpha Power'};
%     ylims = {[-1 1] [2 8] [3 5] [1 7] [4 21] [5 9] [1 10] [14 25]};
% 
%     for i=1:length(chans)
%     figure;
%     pop_erpimage(EEG,1, chans(i),[],band_labels{i},smooth,1,{RWD_label},[],'epoch' ,...
%         'yerplabel','\muV','erp','on','cbar','on', 'limits',[NaN NaN ylims{i} NaN NaN NaN NaN] );
%     end
% end



% 
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% EEG = pop_biosig('D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.edf', 'importevent','off');
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
% pop_saveh( EEG.history, 'eeglabhist.m', 'D:\Users\Sources\Entropy\ENTROPY FOR RAPHAEL\');
% eeglab redraw;
