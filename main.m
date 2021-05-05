close all;
clear;
clc;

rootpath = 'D:\Users\Sources\Entropy\perm_entropy';
datapath = 'D:\Users\NFB\Pacientes\';

chdir(rootpath);

addpath('PE');
addpath('functions');

subj_list = ["JLC"]; % ["EYK", "JLC", "JRJ", "SAJ"];
SBJ_DT = struct();
TASKS = struct();

SBJ_DT = fntools.gendatastruct(datapath, SBJ_DT, subj_list );
% % -----------------------------------------------------------------
for suj=1:length(subj_list)
    filepath = SBJ_DT(suj).filespath;
end

for indx=1:length(filepath)
    tminmax = [10 1939];
    EEG(indx) = init_process(filepath{indx}, tminmax);

    % % -----------------------------------------------------------------
    % Gera surrogate e calcula entropia permutação, inclui no objeto EEG
    splitpath = strsplit(filepath{indx}, {'\'});
    filename = splitpath{end};

%     [PEdata, EEG] = fntools.pesurdata(EEG, filename);
    [EEG(indx)] = fntools.pesurdata(EEG(indx), filename);

    % % -----------------------------------------------------------------
    % |-------|-------------------|----------|------|
    % 0  bs1  60       nfb       1860  bs2  1920   1829
    % Separa trechos baseline  % EEG = pop_select( EEG,'time',[0 1920] );
    TASKS(indx).nfb = pop_select( EEG(indx),'time',[60 1860] );
    TASKS(indx).bs1 = pop_select( EEG(indx),'time',[0 60] );
    TASKS(indx).bs2 = pop_select( EEG(indx),'time',[1860 1920] );

    % % -----------------------------------------------------------------
    % Cria os eventos e extrai trials
    ev_range = [-.3 .5]; % ms
    rwdchan = 3;
    latency = '0<=1800';
    RWD_label = 'RWD-250';
    ev_duration = '.250<=2'; %(s)
    EEGev(indx) = fntools.create_events(TASKS(indx).nfb, rwdchan, latency, RWD_label, ev_duration, ev_range);

    disp('Eventos criados...');

end

% % -----------------------------------------------------------------
% % Calcula entropia permutação por trials   % 1:dado, 11:surrogate, 12:PermH, %13:surpermH
% PEdata_bytrials = fntools.pe_bytrials(EEGev, 1);
% PESurr_bytrials = fntools.pe_bytrials(EEGev, 11);
% disp('Entropia e Surrgate terminado...');

%% -----------------------------------------------------------------
% -------------------------------------------------------------------------
% VISUALIZACAO
% -------------------------------------------------------------------------
% Mostra o espectro (confirmar resultado do feedback)
% figDir = join([currpath, "fig", sbjName], "\");

f_range = [1 50];
chans = [1 11];
for indx=1:length(filepath)
    splitpath = strsplit(filepath{indx}, {'\'});
    filename = splitpath{end};
    [time, frex, currfig] = vis.show_erpspectrum(EEGev(indx), ...
                                                 EEGev(indx).srate, ...
                                                 ev_range, ...
                                                 f_range, ...
                                                 chans);
    supt = suptitle(filename);
    set(supt,'Interpreter','none');
    % ---------------------------------------------------------------------
    % Saving fig
    vis.savefigure(currfig, rootpath, filename, "speccom");
    close(currfig);
    % % -----------------------------------------------------------------
    % Mostra espectro potencia PSD
    specfig = figure('visible','off');
%     set(specfig, 'visible', 'off');
    pop_spectopo(EEGev(indx), 1, [], 'EEG', 'freqrange', [1 45], 'plotchans', [1 11]);
    legend('C3 channel','Surrogate');
    vis.savefigure(specfig, rootpath, filename, "psdspec");
    close(specfig);
    % % -----------------------------------------------------------------
    % Visualiza PE dado e surrogate por trials;
    pefig = vis.show_evtrials(EEGev(indx), ev_range);
    vis.savefigure(pefig, rootpath, filename, "htrials");
    % % -----------------------------------------------------------------
    % Comparação PE entre baseline Pre e Pos
    % Correlação
    basefig = vis.show_pebaselines(EEG(indx));
    vis.savefigure(basefig, rootpath, filename, "pebases");
    close(basefig);
end
%% -----------------------------------------------------------------









