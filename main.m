close all;
clear;
clc;

rootpath = 'D:\Users\Sources\Entropy\perm_entropy';
datapath = 'D:\Users\NFB\Pacientes\';
savesurpath = 'D:\Users\Sources\Entropy\perm_entropy\datas';

chdir(rootpath);

addpath('PE');
addpath('functions');

subj_list = ["JLC"]; % ["EYK", "JLC", "JRJ", "SAJ"];
SBJ_DT = struct();

SBJ_DT = fntools.gendatastruct(datapath, SBJ_DT, subj_list );

% % -----------------------------------------------------------------
TASKS = struct();
bs_data = struct();
metrics = struct();

t_slice = [10 1939];
frange = {[4 7] [12 15] [20 30] [8 12]};

% Parametros criaçao dos eventos
ev_range = [-.3 .5]; % ms
RWD_chan = 3;
RWD_label = 'RWD-250';
latency = '0<=1800';
ev_duration = '.250<=2'; %(s)
rejspecevent = true; % rejeitar epocas baseado no espectro
% Task slices em segundos
nfb_slice = [60 1860];
bs1_slice = [0 60];
bs2_slice = [1860 1920];

loadsurdata = true;
savesurdata = false;

% % -----------------------------------------------------------------
for suj=1:length(subj_list)
    filepath = SBJ_DT(suj).filespath;
    
    for sess=1:length(filepath)
        EEG(sess) = init_process(filepath{sess}, t_slice);

        % % -----------------------------------------------------------------
        % Gera surrogate e calcula entropia permutação, inclui no objeto EEG
        splitpath = strsplit(filepath{sess}, {'\'});
        filename = splitpath{end};  % XXX-DDMMAA_Sn
        
        % istherefile = false; calc surrogate
        % savesurdata
        [EEG(sess)] = fntools.pesurdata(EEG(sess), filename, savesurpath, loadsurdata, savesurdata);
        
        % -----------------------------------------------------------------
        % |-------|-------------------|----------|------|
        % 0  bs1  60       nfb       1860  bs2  1920   1829
        % Separa trechos baseline  % EEG = pop_select( EEG,'time',[0 1920] );
        TASKS(sess).nfb = pop_select( EEG(sess),'time', nfb_slice );
        TASKS(sess).bs1 = pop_select( EEG(sess),'time', bs1_slice );
        TASKS(sess).bs2 = pop_select( EEG(sess),'time', bs2_slice );
        % 12:dados 13:surrogates
        bs_data(sess).bs = { TASKS(sess).nfb.data(12,:,:), TASKS(sess).nfb.data(13,:,:);
                             TASKS(sess).bs1.data(12,:,:), TASKS(sess).bs1.data(13,:,:);
                             TASKS(sess).bs2.data(12,:,:), TASKS(sess).bs2.data(13,:,:) };    
    %     disp(">> CORRELATION WSDIST....");

        for bs=1:3
            metrics(sess).corr(bs) = { corrcoef( bs_data(sess).bs{bs,1}, ...
                                                 bs_data(sess).bs{bs, 2} )};
            metrics(sess).wass(bs) = ws_distance(bs_data(sess).bs{bs,1}, ...
                                                 bs_data(sess).bs{bs, 2});
        end

        % Cria os eventos e extrai trials por sessao
        SBJ_DT(suj).events(sess) = fntools.create_events( TASKS(sess).nfb, ...
                                                          RWD_chan, latency, ...
                                                          RWD_label, ev_duration, ...
                                                          ev_range, ...
                                                          filename, rejspecevent); % true - rejspec
        disp("Eventos criados. Sessao: " + string(sess));
        
%         bandsPWR = fntools.calc_bandpower(SBJ_DT(suj).events(sess), frange);
        SBJ_DT(suj).bandspower{sess,1} = fntools.calc_bandpower(SBJ_DT(suj).events(sess), frange);
    end

    SBJ_DT(suj).tasks = TASKS;
%     SBJ_DT(suj).events = EEGev;
    SBJ_DT(suj).metrics = metrics;
end

%% -----------------------------------------------------------------
% -------------------------------------------------------------------------
% VISUALIZACAO
% -------------------------------------------------------------------------

% BOXPLOT
bndPWR = SBJ_DT.bandspower;
g = [];
banDatas = cellfun(@(x) x(:,1), bndPWR,'UniformOutput',false);
bandsvec = vertcat(banDatas{:});

for sess=1:5
    s = split(SBJ_DT.filespath{sess}, "_");
    g = [g; repmat(sess,size(banDatas{sess},1),1) ];
end
boxchart(g, bandsvec); %,'notch','on');

%%
figure
for sess=1:5
    plot(bndPWR{sess}(:,2));
    hold on;
end
hold off;
%%
% figure
for i=1:5
    data = permute(SBJ_DT.events(i).data(12,:,:), [2 3 1]);
    dtmean = mean(data, 2);
    setname = SBJ_DT.events(i).setname;
    times = SBJ_DT.events(i).times;
    
    datasur = permute(SBJ_DT.events(i).data(13,:,:), [2 3 1]);
    dtmeansur = mean(datasur, 2);
    
    figure
    plot(times, dtmean,'linewidth',2); %'DisplayName', setname(12:13));
    hold on;
    plot(times, dtmeansur,'linewidth',2); %'DisplayName', setname(12:13));
    hold off;
    legend({'Data', 'Surrogate'});
end
title('Permutation Entropy  Data');
legend('show');
% hold off;
%%

figure
for i=1:5
    datasur = permute(SBJ_DT.events(i).data(13,:,:), [2 3 1]);
    dtmeansur = mean(datasur, 2);
    setname = SBJ_DT.events(i).setname;
    times = SBJ_DT.events(i).times;
    
    plot(times, dtmeansur,'linewidth',2,'DisplayName', setname(12:13));
    hold on;
end
title('Permutation Entropy Surrogate');
legend('show');
hold off;

%% -----------------------------------------------------------------

% Comparação PE entre baseline Pre, Pos, NFB
% Correlação e Ws distance
show_peBaselines(SBJ_DT, rootpath);

% % -----------------------------------------------------------------
% Mostra o espectro (confirmar resultado do feedback)
f_range = [1 45];
chans = [1 11];
show_erpSpec(SBJ_DT, rootpath, chans, ev_range, f_range);

% % -----------------------------------------------------------------
% Mostra espectro potencia PSD
show_PSDSpec(SBJ_DT, rootpath, chans, f_range);

% % -----------------------------------------------------------------
% Visualiza PE dado e surrogate por trials;
show_evPETrials(SBJ_DT, rootpath, ev_range);

% % -----------------------------------------------------------------

%% -----------------------------------------------------------------
% LOCAL VIS FUNCTIONS

function show_peBaselines(SBJ_DT, rootpath)

    subj_list = [ SBJ_DT.names ];

    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEG = SBJ_DT(suj).tasks;
        metrics = SBJ_DT(suj).metrics;

        for sess=1:length(filepath)
            splitpath = strsplit(filepath{sess}, {'\'});
            filename = splitpath{end};
            basefig = vis.show_pebaselines(EEG(sess), metrics(sess));
            vis.savefigure(basefig, rootpath, filename, "pebases");
            close(basefig);
        end
    end
end

function show_erpSpec(SBJ_DT, rootpath, chans, ev_range, f_range)

    subj_list = [ SBJ_DT.names ];

    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEGev = SBJ_DT(suj).events;

        for sess=1:length(SBJ_DT(suj).filespath)
            splitpath = strsplit(SBJ_DT(suj).filespath{sess}, {'\'});
            filename = splitpath{end};

            [~, ~, currfig] = vis.show_erpspectrum(EEGev(sess), ...
                                                   EEGev(sess).srate, ...
                                                   ev_range, ...
                                                   f_range, ...
                                                   chans);
%             supt = suptitle(filename);
%             set(supt,'Interpreter','none');
            vis.savefigure(currfig, rootpath, filename, "speccom");
            close(currfig);
        end
    end
end

function show_PSDSpec(SBJ_DT, rootpath, channels, f_range)

    subj_list = [ SBJ_DT.names ];
    
    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEGev = SBJ_DT(suj).events;

        for sess=1:length(filepath)
            splitpath = strsplit(filepath{sess}, {'\'});
            filename = splitpath{end};

            specfig = figure('visible','off');
            pop_spectopo(EEGev(sess), 1, [], 'EEG', ...
                         'freqrange', f_range, ...
                         'plotchans', channels);
            legend('C3 channel','Surrogate');
            vis.savefigure(specfig, rootpath, filename, "psdspec");
            close(specfig);
        end
    end
end

function show_evPETrials(SBJ_DT, rootpath, ev_range)

    subj_list = [ SBJ_DT.names ];

    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEGev = SBJ_DT(suj).events;

        for sess=1:length(filepath)
            splitpath = strsplit(filepath{sess}, {'\'});
            filename = splitpath{end};

            pefig = vis.show_evtrials(EEGev(sess), ev_range);
            vis.savefigure(pefig, rootpath, filename, "htrials");
            close(pefig);
        end
    end
end


% % -----------------------------------------------------------------
% % Calcula entropia permutação por trials   % 1:dado, 11:surrogate, 12:PermH, %13:surpermH
% PEdata_bytrials = fntools.pe_bytrials(EEGev, 1);
% PESurr_bytrials = fntools.pe_bytrials(EEGev, 11);
% disp('Entropia e Surrgate terminado...');

% % -----------------------------------------------------------------
% f_range = [1 50];
% chans = [1 11];
% for suj=1:length(subj_list)
%     filepath = SBJ_DT(suj).filespath;
%     EEGev = SBJ_DT(suj).events;
%     EEG = SBJ_DT(suj).tasks;
%     metrics = SBJ_DT(suj).metrics;
%     
%     for sess=1:length(filepath)
%         splitpath = strsplit(filepath{sess}, {'\'});
%         filename = splitpath{end};
%         
%         [time, frex, currfig] = vis.show_erpspectrum(EEGev(sess), ...
%                                                      EEGev(sess).srate, ...
%                                                      ev_range, ...
%                                                      f_range, ...
%                                                      chans);
%         supt = suptitle(filename);
%         set(supt,'Interpreter','none');
%         vis.savefigure(currfig, rootpath, filename, "speccom");
%         close(currfig);
% 
%         % % -----------------------------------------------------------------
%         % Mostra espectro potencia PSD
%         specfig = figure('visible','off');
%     %     set(specfig, 'visible', 'off');
%         pop_spectopo(EEGev(sess), 1, [], 'EEG', 'freqrange', [1 45], 'plotchans', [1 11]);
%         legend('C3 channel','Surrogate');
%         vis.savefigure(specfig, rootpath, filename, "psdspec");
%         close(specfig);
% 
%         % % -----------------------------------------------------------------
%         % Visualiza PE dado e surrogate por trials;
%         pefig = vis.show_evtrials(EEGev(sess), ev_range);
%         vis.savefigure(pefig, rootpath, filename, "htrials");

%         % % -----------------------------------------------------------------
%         % Comparação PE entre baseline Pre e Pos
%         % Correlação
%         basefig = vis.show_pebaselines(EEG(sess), metrics(sess));
%         vis.savefigure(basefig, rootpath, filename, "pebases");
%         close(basefig);
%     end
% end








