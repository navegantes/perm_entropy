% EEGLAB history file generated on the 05-Jan-2021
% ------------------------------------------------
close all;
clear;
clc;

addpath('PE');



tmin = 10;
tmax = 1930;
RWD_label = 'RWD-250';
ev_duration = '.250<=2'; %(s)

% Read from csv file
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
fl = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.csv';
file = csvread(fl, 2, 0);
file = file';
EEG = pop_importdata( 'setname', 'JLC-210120_S2-csv', 'data', file, 'nbchan', 13,...
                      'dataformat', 'array', 'srate', 256);
EEGcsv = pop_select( EEG,'time',[tmin tmax] );

%Read from edf file
fl = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.edf';
EEG = pop_biosig('D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.edf', 'importevent','off');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', 'JLC-210120_S2-edf', 'gui','off');
EEG = pop_select( EEG,'time',[tmin-1 tmax-1] );
EEG = eeg_checkset( EEG );
eeglab redraw;

% Fitragem passa-faixa.
EEG    = pop_eegfiltnew(EEG, 1,100,900,0,[],0);
% RWD and Amp bands, theta smr hibeta, from csv file
chanlist = [11 2 5 8];
for chan=1:length(chanlist)
    EEG.data(end+1,:) = EEGcsv.data(chanlist(chan), :);
end

% EEGcsv = pop_eegfiltnew(EEGcsv, 1,100,900,0,[],0);

fminmax = {[4 7] [12 15] [20 30] [8 12]};
for f=1:length(fminmax)
    freq = fminmax{f};
    EEGfilt = pop_eegfiltnew(EEG, freq(1),freq(2),900,0,[],0);
    hilb = hilbert(EEGfilt.data(1, :));
    EEG.data(end+1,:) = abs(hilb).^2;
end

delay = 1; % delay 1 between points in ordinal patterns (successive points)
order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
windowSize = 4*EEG.srate;
H_perm = zeros(1, size(EEG.data(1,:), 2));
h = PE(EEG.data(1,:)', delay, order, windowSize);
H_perm(1:length(h)) = h;
%Computing Entropy
EEG.data(end+1,:) = H_perm; 

EEG.nbchan = size(EEG.data,1);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
% pop_eegplot( EEG, 1, 1, 1);


% define os eventos  disp([newline 'Def events']);
EEG = pop_chanevent( EEG, 3, 'edge', 'leading', 'duration', 'on',...
                     'typename', 'RWD', 'delchan', 'off', 'edgelen', 1);
% disp([newline 'Selecting events']);

% Seleciona eventos duracao maior que 250 ms
EEG = pop_selectevent( EEG, 'latency','60<=1860','duration',ev_duration,...
                       'renametype',RWD_label,'oldtypefield','RWD','deleteevents','on');
% disp([newline 'Epochs']);
EEG = pop_epoch( EEG, { RWD_label }, [-0.3 0.5], ...
                 'newname', 'JLC-210120_S2-epochs', 'epochinfo', 'yes');
%% --------------------------------------------------------------------------
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
EEG = eeg_checkset( EEG );
eeglab redraw;

hdatas   = cell(sz(3));
sz = size(EEG.data);
delay = 1; % delay 1 between points in ordinal patterns (successive points)
order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
windowSize = 4*EEG.srate;

for i=1:
hdatas = PE(sbj_dt(EEG.data(1,:, :)', delay, order, windowSize);



show_erp_events(EEG);

function EEG = show_erp_events(EEG)
    pop_eegplot( EEG, 1, 1, 1);

    smooth = 5;
    chans = [1 4 5 6 7 8 9 10];
    band_labels = {'EEG' 'Theta Amp' 'SMR Amp' 'Hibeta Amp' 'Theta Power'...
                   'SMR Power' 'Hibeta Power' 'Alpha Power'};
    ylims = {[-1 1] [2 8] [3 5] [1 7] [4 21] [5 9] [1 10] [14 25]};

    for i=1:length(chans)
    figure;
    pop_erpimage(EEG,1, [chans(i)],[[]],band_labels{i},smooth,1,{RWD_label},[],'epoch' ,...
        'yerplabel','\muV','erp','on','cbar','on', 'limits',[NaN NaN ylims{i} NaN NaN NaN NaN] );
    end
end

% 
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% EEG = pop_biosig('D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.edf', 'importevent','off');
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
% pop_saveh( EEG.history, 'eeglabhist.m', 'D:\Users\Sources\Entropy\ENTROPY FOR RAPHAEL\');
% eeglab redraw;
