
chdir 'D:\Users\Sources\Entropy\perm_entropy';
%% -----------------------------------------------------------------

close all;
clear;
clc;

rootpath = 'D:\Users\Sources\Entropy\perm_entropy';
chdir(rootpath)

addpath('PE');
addpath('functions')

tmin = 10;
tmax = 1939;

% -------------------------------------------------------------------------
% Read from csv file
% csv_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-200120\JLC-200120_S1.csv';
% edf_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-200120\JLC-200120_S1.edf';
% csv_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.csv';
% edf_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-210120\JLC-210120_S2.edf';
% csv_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-220120\JLC-220120_S3.csv';
% edf_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-220120\JLC-220120_S3.edf';
% csv_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-230120\JLC-230120_S4.csv';
% edf_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-230120\JLC-230120_S4.edf';
% csv_file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-240120\JLC-240120_S5.csv';
file = 'D:\Users\NFB\Pacientes\JLC\NFB\nfb-240120\JLC-240120_S5'; %.edf';
% -------------------------------------------------------------------------
% Init_process()
csv_file = join([file, '.csv']);
edf_file = join([file, '.edf']);

spltpath = strsplit(edf_file, {'\', '.'});
filename = spltpath(8);
sbjName = spltpath(5);
session = strsplit(filename{1}, '_');
session = session(2);

file = csvread(csv_file, 2, 0);
file = file';

rawCSV = pop_importdata( 'setname', 'JLC-210120_S2-csv', ...
                         'data', file, ...
                         'nbchan', 13, ...
                         'dataformat', 'array', ...
                         'srate', 256);
                     
EEGcsv = pop_select( rawCSV, 'time',[tmin+1 tmax+1] );

% -------------------------------------------------------------------------
%Read from edf file
EEG = pop_biosig(edf_file, 'importevent','off');
EEG = pop_select( EEG,'time',[tmin tmax] );

% Muda as labels dos canais
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
    EEG.chanlocs(f+lenfmM+2).labels = labels{f+lenfmM};
end

% clear EEGfilt;

disp('Done...');
% Init_process() Fim

%% -----------------------------------------------------------------
% Gera surrogate e calcula entropia permutação, inclui no objeto EEG
% [PEdata, EEG] = fntools.pesurdata(EEG, filename{1});
[EEG] = fntools.pesurdata(EEG, filename{1});

%% -----------------------------------------------------------------
% |-------|-------------------|----------|
% 0  bs1  60       nfb       1860  bs2  1920
% Separa trechos baseline
% EEG = pop_select( EEG,'time',[0 1920] );
NFB = pop_select( EEG,'time',[60 1860] );
BS1 = pop_select( EEG,'time',[0 60] );
BS2 = pop_select( EEG,'time',[1860 1920] );

%% -------------------------------------------------------------------------
% Cria os eventos e extrai trials        % EEGnfb = pop_select( EEG,'time',[60 1860] );
ev_range = [-.3 .5]; % ms
rwdchan = 3;
latency = '0<=1800';
RWD_label = 'RWD-250';
ev_duration = '.250<=2'; %(s)
EEGev = fntools.create_events(NFB, rwdchan, latency, RWD_label, ev_duration, ev_range);

disp('Eventos criados...');

%% -------------------------------------------------------------------------
% Calcula entropia permutação por trials   % 1:dado, 11:surrogate, 12:PermH, %13:surpermH
PEdata_bytrials = fntools.pe_bytrials(EEGev, 1);
PESurr_bytrials = fntools.pe_bytrials(EEGev, 11);

disp('Entropia e Surrgate terminado...');

%% -----------------------------------------------------------------
% -------------------------------------------------------------------------
% VISUALIZACAO
% -------------------------------------------------------------------------
% Mostra o espectro (confirmar resultado do feedback)
% figDir = join([currpath, "fig", sbjName], "\");
show = true;

f_range = [1 50];
chans = [1 11];
[time, frex, currfig] = vis.show_erpspectrum(EEGev, EEGev.srate, ev_range, f_range, chans);
suptitle(session);
% -------------------------------------------------------------------------
% Saving fig
fntools.savefigure(currfig, rootpath, filename{1}, "speccom");

%% -----------------------------------------------------------------
% Mostra espectro potencia PSD
currfig = figure;
pop_spectopo(EEGev, 1, [], 'EEG', 'freqrange', [1 45], 'plotchans', [1 11]);
legend('C3 channel','Surrogate');
fntools.savefigure(currfig, rootpath, filename{1}, "spectopo");

% % -----------------------------------------------------------------
% Visualiza PE dado e surrogate por trials;
pefig = vis.show_evtrials(EEGev, ev_range);
fntools.savefigure(pefig, rootpath, filename{1}, "petrials");

%% -----------------------------------------------------------------
% Comparação PE entre baseline Pre e Pos
% Correlação
basefig = vis.show_pebaselines(EEG);
fntools.savefigure(basefig, rootpath, filename{1}, "pebases");

% % -----------------------------------------------------------------
% % Checagem das condicoes de reconpensa NFB
% % pop_eegplot( EEGev, 1, 1, 1);
% [~] = vis.show_erp_events(EEGev, RWD_label);


%% -----------------------------------------------------------------
wess_dist = zeros(1,size(Hperm,1));
for i=1:size(Hperm,1)
    wess_dist(i) = ws_distance(EEGev.data(12,:,i), EEGev.data(13,:,i),1);
end

%% -----------------------------------------------------------------
figure;
hold on;
plot(PEdata.permH(1:end-5*256-2));
plot(PEdata.surpermH(1:end-5*256-2));
hold off;

%% -----------------------------------------------------------------
figure;
[ntrials, lentime] = size(PEdata_bytrials);
hold on;
for trial=1:20
    plot(time(1:139), PEdata_bytrials(trial, 1:139));
end
hold off;
xlabel('Time (ms)');
ylabel('Permutation Entropy');

%% ----
figure;
hold on;
for trial=1:20
    plot(time(1:139), PESurr_bytrials(trial, 1:139));
end
hold off;
xlabel('Time (ms)');
ylabel('Permutation Entropy');

%% -----------------------------------------------------------------
% 11:surrogate, 12:PermH, %13:surpermH
figure;
hold on;
[~, len, num] = size(EEGev.data);
for i=1:20
    plot(time(1:139), EEGev.data(12,1:139, i));
end
hold off;

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
