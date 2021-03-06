%
close all;
clear;
clc;

datapath = 'D:\Users\NFB\Pacientes\';
% rootpath = 'E:\Datas\nfb-data';
rootpath = 'D:\Users\Sources\Entropy\perm_entropy';
scriptspath = 'D:\Users\Sources\Entropy\perm_entropy';

% chdir(scriptspath);
chdir(rootpath);
addpath('PE');
addpath('functions');

subj_list = ["EYK"]; % ["EYK", "JLC", "JRJ", "SAJ"];
SBJ_DT = struct();

SBJ_DT = fntools.gendatastruct(datapath, SBJ_DT, subj_list );

% % -----------------------------------------------------------------

TASKS = struct();
bs_data = struct();
metrics = struct();

% Parametros cria�ao dos eventos
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
zscore_norm = false;
savesbjdata = false;

lang = 'ptbr'; %en';

% % -----------------------------------------------------------------

t_slice = [10 1939];
frange = {[4 7] [12 15] [20 30] [8 12]};

for suj=1:length(subj_list)
    filepath = SBJ_DT(suj).filespath;
    
    for sess=1:length(filepath)
        EEG(sess) = init_process(filepath{sess}, t_slice, zscore_norm);

        % % -----------------------------------------------------------------
        % Gera surrogate e calcula entropia permuta��o, inclui no objeto EEG
        splitpath = strsplit(filepath{sess}, {'\'});
        filename = splitpath{end};  % XXX-DDMMAA_Sn
        
        [EEG(sess)] = fntools.pesurdata(EEG(sess), filename, rootpath, loadsurdata, savesurdata);
        
        % -----------------------------------------------------------------
        % |-------|-------------------|----------|------|
        % 0  bs1  60       nfb       1860  bs2  1920   1829
        % Separa trechos baseline
        TASKS(sess).nfb = pop_select( EEG(sess),'time', nfb_slice );
        TASKS(sess).bs1 = pop_select( EEG(sess),'time', bs1_slice );
        TASKS(sess).bs2 = pop_select( EEG(sess),'time', bs2_slice );
        % 12:dados 13:surrogates
        bs_data(sess).bs = { TASKS(sess).nfb.data(12,:,:), TASKS(sess).nfb.data(13,:,:);
                             TASKS(sess).bs1.data(12,:,:), TASKS(sess).bs1.data(13,:,:);
                             TASKS(sess).bs2.data(12,:,:), TASKS(sess).bs2.data(13,:,:) };    

        for bs=1:3
            metrics(sess).corr(bs) = { corrcoef( bs_data(sess).bs{bs,1}, bs_data(sess).bs{bs, 2} )};
            metrics(sess).wass(bs) = ws_distance( bs_data(sess).bs{bs,1}, bs_data(sess).bs{bs, 2} );
        end
        
        % Cria os eventos e extrai trials por sessao
        SBJ_DT(suj).events(sess) = fntools.create_events( TASKS(sess).nfb, ...
                                                          RWD_chan, latency, ...
                                                          RWD_label, ev_duration, ...
                                                          ev_range, ...
                                                          filename, rejspecevent); % true - rejspec
        disp("..." +newline+ ...
             "Eventos criados. Sessao: " + string(sess));
        SBJ_DT(suj).bandspower{sess,1}    = fntools.calc_bandpower(SBJ_DT(suj).events(sess), 1, frange);
        SBJ_DT(suj).surbandspower{sess,1} = fntools.calc_bandpower(SBJ_DT(suj).events(sess), 11, frange);
    end

    SBJ_DT(suj).tasks = TASKS;
    SBJ_DT(suj).metrics = metrics;
end

%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% VISUALIZACAO
% -------------------------------------------------------------------------
%% ------------------------------------------------------------------------

% Comparação espectro nfb/surrogate (confirmar resultado do feedback)
f_range = [1 45];
chans = [1 11];
taskinfo = ["nfb", "speccom"];

time2plotnfb = show_ERPSpec__(SBJ_DT, rootpath, chans, ev_range, f_range, taskinfo, lang);
% % Comparação espectro baseline pre-pos
% bs_evrange = [0 1];
% taskinfo = ["bs1", "specbs1"];
% time2plotbs1 = show_ERPSpec__(SBJ_DT, rootpath, chans, bs_evrange, f_range, taskinfo, lang);
% taskinfo = ["bs2", "specbs2"];
% time2plotbs2 = show_ERPSpec__(SBJ_DT, rootpath, chans, bs_evrange, f_range, taskinfo, lang);

%% ------------------------------------------------------------------------
% Visualiza PE dado e surrogate por trials;
show_evPETrials__(SBJ_DT, rootpath, ev_range, lang);
%% ------------------------------------------------------------------------


show_Entropies__(SBJ_DT, rootpath, lang);
%% ------------------------------------------------------------------------

showCoefsBoxchart__(SBJ_DT, rootpath, lang);
%% ------------------------------------------------------------------------

if savesbjdata
    num = char(string(length(SBJ_DT)));
    save([rootpath, '\surdatas\sbjdata_', num], 'SBJ_DT');
end
%% ------------------------------------------------------------------------



% Comparação PE entre baseline Pre, Pos, NFB
% Correlação e Ws distance
% showPEBaselines__(SBJ_DT, rootpath, lang);

% % ------------------------------------------------------------------------
% % Mostra espectro potencia PSD
% show_PSDSpec__(SBJ_DT, rootpath, chans, f_range, lang);
% % % ------------------------------------------------------------------------
% taskinfo = ["bs1", "psdbs1"];
% show_PSDSpec__(SBJ_DT, rootpath, chans, f_range, taskinfo, lang);
% taskinfo = ["bs2", "psdbs2"];
% show_PSDSpec__(SBJ_DT, rootpath, chans, f_range, taskinfo, lang);

% % ------------------------------------------------------------------------
% % Boxchart potencias bandas
% disp("..." +newline+ "Generating boxplot...");
% showBandsBoxchart__(SBJ_DT, rootpath, lang);
% 
% disp("..." +newline+ "Visualization done...");

%% ------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                TESTES                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------------------------------------------------

% marca a posição dos eventos de feedback
showevents(SBJ_DT);
%% ------------------------------------------------------------------------

dt = SBJ_DT.events(1);
d = reshape(dt.data(12,:,:),[],1)';
chan = 12;
m=mean(dt.data(12, :,:),3);
mM = minmax(d);
% 0.8556 .8572
% 0.828  .832
pop_erpimage(dt,1, chan,[[]],dt.chanlocs(chan).labels,1,1,{},[],'' , ...
             'yerplabel','H','erp',1, 'renorm', 'no',...
             'limits',[NaN NaN 0.8556 .8572 NaN NaN NaN NaN], ...
             'cbar','on', 'caxis', [mM(1) 1], 'erpstd','on');
%% ------------------------------------------------------------------------
figure
bndPWR = SBJ_DT(1).bandspower;
for sess=1:5
    plot(bndPWR{sess}{2}(:,2), 'o-', 'linewidth', 1.);
    hold on;
end
hold off;
%% ------------------------------------------------------------------------
% figure data e surrogate por sessao
suj = 1;
for i=1:5
%     data = permute(SBJ_DT(suj).events(i).data(12,:,:), [2 3 1]);
    data = SBJ_DT(suj).events(i).data(12,:,:);
    dtmean = mean(data, 3);
%     setname = SBJ_DT(suj).events(i).setname;
    times = SBJ_DT(suj).events(i).times;
    
    datasur = permute(SBJ_DT(suj).events(i).data(13,:,:), [2 3 1]);
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
%% ------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                LOCAL VIS FUNCTIONS                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
chans = [12 13];
for c=1:2
    [SBJ_DT] = show_SlopesAVG__(SBJ_DT, chans(c), rootpath, lang);
% [SBJ_DT] = show_SlopesAVG__(SBJ_DT, 13, rootpath, lang);
end
%% ------------------------------------------------------------------------

% SBJ_DT = show_statsTrials__(SBJ_DT, rootpath);

SBJ_DT = show_crosStatsTrials__(SBJ_DT, rootpath, lang);
%% ------------------------------------------------------------------------

show_bandpower(SBJ_DT, rootpath);
%% ------------------------------------------------------------------------

function SBJ_DT = show_crosStatsTrials__(SBJ_DT, rootpath, lang) %, ev_range, lang)
    
    if (nargin < 3 || strcmp(lang, 'ptbr'))
        lang = 'ptbr'; %'en';
    else
        lang = 'en';
    end

    subj_list = [ SBJ_DT.names ];
    trange    = [0 250];
    
    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEGev    = SBJ_DT(suj).events;
        statsdt  = {cell(5,1), cell(5,1)};
        statsum  = {zeros(5,2), zeros(5,2)};
        ntrials  = zeros(5,1);
        
        for sess=1:length(filepath)
            tindx = get_evtimeindxs(EEGev(sess).times, trange);
            
            numtrials = EEGev(sess).trials;
            pe        = permute(EEGev(sess).data(12,:,:), [2 3 1]);
            pesu      = permute(EEGev(sess).data(13,:,:), [2 3 1]);
            se        = permute(EEGev(sess).data(14,:,:), [2 3 1]);
            sesu      = permute(EEGev(sess).data(15,:,:), [2 3 1]);
            entrData  = {pe, pesu, se, sesu};
            
%             numdata = length(entrData);
            pval = {zeros(numtrials, 2), zeros(numtrials, 2)};
            hval = {zeros(numtrials, 2), zeros(numtrials, 2)};
            Sval = {cell(numtrials, 2), cell(numtrials, 2)};
            
            test = {'ranksum', 'ttest'};
            for sts=1:2
                [pval{sts}, hval{sts}, rS] = stats.cross_stats(pe, pesu, tindx, test{sts});
                for rwd=1:2
                    for t=1:numtrials
                        Sval{sts}{t,rwd} = rS{rwd}{t,1};
                    end
                end
                
                statsdt{sts}{sess}.p = pval{sts};
                statsdt{sts}{sess}.h = hval{sts};
                statsdt{sts}{sess}.S = Sval{sts};
                
                statsum{sts}(sess, :) = sum(hval{sts});
                ntrials(sess,1)       = numtrials;
            end
        end
        
        Crosstats.ranksum = statsdt{1};
        Crossum.ranksum   = statsum{1};
        Crosstats.ttest   = statsdt{2};
        Crossum.ttest     = statsum{2};
        
        SBJ_DT(suj).Crosstats = Crosstats;
        SBJ_DT(suj).Crossum   = Crossum;
        
        splitpath = strsplit(filepath{sess}, {'\'});
        filename  = splitpath{end};
        
        info.title = '';
        info.legend = {'EP-EPsu: nRWD', 'EP-EPsu: RWD'};
        
        [fig, fign] = vis.show_numtrials(Crossum.ranksum, ntrials, info, lang);
%         vis.savefigure(fig,  rootpath, filename(1:3), "stats-ranksum");
%         vis.savefigure(fign, rootpath, filename(1:3), "stats-ranksum-n");
%         close(fig); close(fign);
%         
%         [fig, fign] = vis.show_numtrials(soma.ttest,   ntrials, test{2}, lang);
%         vis.savefigure(fig,  rootpath, filename(1:3), "stats-ttest");
%         vis.savefigure(fign, rootpath, filename(1:3), "stats-ttest-n");
%         close(fig); close(fign);
    end
end

function show_bandpower(SBJ_DT, rootpath)

    subj_list = [ SBJ_DT.names ];
    
    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        splitpath = strsplit(filepath{1}, {'\'});
        filename = splitpath{end};
        
        bandpwr    = SBJ_DT.bandspower;
        surbandpwr = SBJ_DT.surbandspower;

        bandnames = {'Teta', 'SMR', 'Hibeta', 'Alpha'};
        lenband = length(bandnames);

        bandmean    = {zeros(5,lenband), zeros(5,lenband)};
        bandstd     = {zeros(5,lenband), zeros(5,lenband)};
        surbandmean = {zeros(5,lenband), zeros(5,lenband)};
        surbandstd  = {zeros(5,lenband), zeros(5,lenband)};

        data    = bandpwr;
        surdata = surbandpwr;

        for rwd=1:2
            for i=1:5
                for band=1:lenband
                    bandmean{rwd}(i,band)    = mean(data{i}{rwd}(:,band));
                    bandstd{rwd}(i,band)     = std(data{i}{rwd}(:,band));
                    surbandmean{rwd}(i,band) = mean(surdata{i}{rwd}(:,band));
                    surbandstd{rwd}(i,band)  = std(surdata{i}{rwd}(:,band));
                end
            end
        end
        
        infos.rootpath  = rootpath;
        infos.filename  = filename(1:3);
        infos.legend    = {'Dados', 'Surr'};
        infos.sufixname = "-datasurr";
        show_bars(bandmean{2}, bandstd{2}, surbandmean{2}, surbandstd{2}, bandnames, infos);
        
        infos.legend    = {'Antes Reforço', 'Depois Reforço'};
        infos.sufixname = "-nrwdrwd";
        show_bars(bandmean{1}, bandstd{1}, bandmean{2}, bandstd{2}, bandnames, infos);
    end
end

function show_bars(data1, std1, data2, std2, bandnames, infos)
    
    fig = cell(1,4);
    lenband = length(bandnames); 
    hexcolors = {'#D95319', '#0072BD'};
    
    for bi=1:lenband
        fig{bi} = figure;
        meanbar = bar([data1(:,bi), data2(:,bi)]); hold on;
        meanbar(1).FaceColor = hexcolors{1};
        meanbar(2).FaceColor = hexcolors{2};
        
        er1 = errorbar(meanbar(1).XEndPoints, meanbar(1).YData, std1(:,bi));
        er1.Color = hexcolors{1}; %[0 0 0]; % meanbar(bi).FaceColor;
        er1.LineStyle = '--';
        er1.LineWidth = 1.2;
        er2 = errorbar(meanbar(2).XEndPoints, meanbar(2).YData, std2(:,bi));
        er2.Color = hexcolors{2}; %[0 0 0];
        er2.LineStyle = '-';
        er2.LineWidth = 1.2;
        
        m1 = max([data1(:,bi); data2(:,bi)]);
        m2 = max([std1(:,bi); std2(:,bi)]);
        
        ylim([-0.01 m1+m2+0.05]);
        
        title(bandnames{bi}, 'FontSize', 16);
        legend(infos.legend);
        ylabel('Potencia Relativa (%)', 'FontSize', 14);
        xlabel('Sessões', 'FontSize', 14);
        
        sfxname = bandnames{bi} + infos.sufixname ;
        vis.savefigure(fig{bi}, infos.rootpath, infos.filename, sfxname, 'bands' )
        close(fig{bi});
    end
%     for f=1:lenband
%         sfxname = bandnames{f} + infos.sufixname ;
%         disp(filename(1:3))
%         vis.savefigure(barfig{f}, rootpath, filename(1:3), sfxname, 'bands' )
%         close(fig{bi});
%     end
end

function showevents(SBJ_DT)

    sess = 2;
    suj  = 1;
    
    dt      = SBJ_DT(suj).tasks(sess).nfb.data(12,:,:);
    times   = SBJ_DT(suj).tasks(sess).nfb.times;
    EEGev   = SBJ_DT(suj).events(sess);
    urevent = EEGev.urevent;
%     pop_eegplot(EEGev,1,1,1);
    
    numevents = length(urevent);
    durations = [urevent.duration];
    lat       = [urevent.latency];
    
    figure;
    plot(times/1000, dt); hold on;
    
    for ev=1:numevents
        if durations(ev)>=64
            val = times(lat(ev))/1000;
            plot(val, dt(lat(ev)), 'og', ...
                'linewidth', .5, ...
                'MarkerFaceColor', 'g', ...
                'MarkerSize',2);
        end
    end
    hold off;
end

function show_Entropies__(SBJ_DT, rootpath, lang)
    
    if (nargin < 3 || strcmp(lang, 'ptbr'))
        lang = 'ptbr'; %'en';
    else
        lang = 'en';
    end
    
    subj_list = [ SBJ_DT.names ];
    
    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEG = SBJ_DT(suj); %.tasks;
        
        offset = 0.1726;
        boxpos = 0.8060;
        
        hfig = figure('Position', [282,132,800,840]);
        for sess=1:length(filepath)
            splitpath = strsplit(filepath{sess}, {'\'});
            filename = splitpath{end};
            
            subplot(5,1, sess,'Parent',hfig);
            ax = gca;
            annonbox = [0.7552 boxpos 0.1447 0.03214];
            vis.show_entropies(EEG, ax, sess, annonbox, lang);
            boxpos = boxpos - offset;
        end
        vis.savefigure(hfig, rootpath, filename(1:end-3), "specpermh");
%         close(hfig);
    end
end

function SBJ_DT = show_statsTrials__(SBJ_DT, rootpath, lang) %, ev_range, lang)
    
    if (nargin < 3 || strcmp(lang, 'ptbr'))
        lang = 'ptbr'; %'en';
    else
        lang = 'en';
    end

    subj_list = [ SBJ_DT.names ];
    trange    = [0 250];
    
    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEGev    = SBJ_DT(suj).events;
        statsdt  = {cell(5,1), cell(5,1)};
        statsum  = {zeros(5,4), zeros(5,4)};
        ntrials  = zeros(5,1);
        
        for sess=1:length(filepath)
            tindx = get_evtimeindxs(EEGev(sess).times, trange);
            
            numtrials = EEGev(sess).trials;
            pe        = permute(EEGev(sess).data(12,:,:), [2 3 1]);
            pesu      = permute(EEGev(sess).data(13,:,:), [2 3 1]);
            se        = permute(EEGev(sess).data(14,:,:), [2 3 1]);
            sesu      = permute(EEGev(sess).data(15,:,:), [2 3 1]);
            entrData  = {pe, pesu, se, sesu};
            
            numdata = length(entrData);
            pval = {zeros(numtrials, numdata), zeros(numtrials, numdata)};
            hval = {zeros(numtrials, numdata), zeros(numtrials, numdata)};
            Sval = {cell(numtrials, numdata), cell(numtrials, numdata)};
            
            test = {'ranksum', 'ttest'};
            for sts=1:2
                for dt=1:numdata
                    [rp, rh, rS] = stats.calc_stats(entrData{dt}, numtrials, tindx, test{sts});
                    pval{sts}(:,dt) = rp;
                    hval{sts}(:,dt) = rh;
                    for t=1:numtrials
                        Sval{sts}{t,dt} =  rS{t};
                    end
                end
                
                statsdt{sts}{sess}.p = pval{sts};
                statsdt{sts}{sess}.h = hval{sts};
                statsdt{sts}{sess}.S = Sval{sts};
                
                statsum{sts}(sess, :) = sum(hval{sts});
                ntrials(sess,1)       = numtrials;
            end
        end
        
        statsdata.ranksum = statsdt{1};
        soma.ranksum      = statsum{1};
        statsdata.ttest   = statsdt{2};
        soma.ttest        = statsum{2};
        
        SBJ_DT(suj).StatsData = statsdata;
        SBJ_DT(suj).StatSum   = soma;
        
        splitpath = strsplit(filepath{sess}, {'\'});
        filename  = splitpath{end};
        
        info.title = test{1};
        info.legend = {'EP', 'EPSu', 'EE', 'EESu'};
        [fig, fign] = vis.show_numtrials(soma.ranksum, ntrials, info, lang);
        vis.savefigure(fig,  rootpath, filename(1:3), "stats-ranksum");
        vis.savefigure(fign, rootpath, filename(1:3), "stats-ranksum-n");
        close(fig); close(fign);
        
        info.title = test{2};
        [fig, fign] = vis.show_numtrials(soma.ttest,   ntrials, info, lang);
        vis.savefigure(fig,  rootpath, filename(1:3), "stats-ttest");
        vis.savefigure(fign, rootpath, filename(1:3), "stats-ttest-n");
        close(fig); close(fign);
    end
end

function tindex = get_evtimeindxs(times, trange)
    t_zero  = dsearchn(times', trange(1));
    t_evend = dsearchn(times', trange(2));
    t_evdne = dsearchn(times', -1*trange(2));
    tindex = [t_evdne t_zero t_evend];
end

function [SBJ_DT] = show_SlopesAVG__(SBJ_DT, chan, rootpath, lang)

    if (nargin < 3 || strcmp(lang, 'ptbr'))
        lang = 'ptbr'; %'en';
    else
        lang = 'en';
    end
    
    % Slope das medias Permutation Entropy
    subj_list = [ SBJ_DT.names ];
    trange    = [0 250]; % em ms

    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        
        numsess      = length(SBJ_DT(suj).filespath);
        dtmeansur    = cell(numsess, 1);
        mpcoefs      = cell(numsess, 1);
        mpolyn       = cell(numsess, 1);
        trialslpcoef = cell(numsess, 1);
        prwdcoefs    = cell(numsess, 1);
        trialpoly    = cell(numsess, 1);
        prwdpolyn    = cell(numsess, 1);

        for sess=1:numsess
            EEGev          = SBJ_DT(suj).events(sess);
            chanlabels     = {EEGev.chanlocs.labels}; 
            datasur        = EEGev.data(chan,:,:);
            times          = SBJ_DT(suj).events(sess).times;
            [tslice, tind] = timeslice(times, trange);
            
            dtmeansur{sess} = mean(datasur, 3);%             [mpcoefs{sess}, mpolyn{sess}, tmslope] = polyslope__(dtmeansur{sess}, times);
            mpcoefs{sess}   = polyfit(tslice{2}, dtmeansur{sess,1}(tind(2):tind(3)), 1);
            mpolyn{sess}    = polyval(mpcoefs{sess}, tslice{2}); % + mpcoefs{sess}(2);
            
            [trialslpcoef{sess}, trialpoly{sess}, prwdcoefs{sess}, prwdpolyn{sess}] = calc_slopes(EEGev, chan, tind);
            
            SBJ_DT(suj).metrics(sess).slopecoefs.(chanlabels{chan})  = trialslpcoef{sess};
            SBJ_DT(suj).metrics(sess).pslopecoefs.(chanlabels{chan}) = prwdcoefs{sess};
        end
        
        coefs  = {prwdcoefs, mpcoefs};
        polyn  = {prwdpolyn, mpolyn, tind};
        slpfig = vis.show_slopesavg(SBJ_DT(suj), chan, coefs, polyn, tslice, lang);
        
        chanlabels = {SBJ_DT(1).events(1).chanlocs.labels};
        splitpath  = strsplit(filepath{sess}, {'\'});
        filename   = splitpath{end};
        vis.savefigure(slpfig, rootpath, filename(1:end-3), "slopesavg-" + chanlabels{chan});
        close(slpfig);
    end
end

function [tslice, tindex] = timeslice(times, trange)
    tzero  = dsearchn(times', trange(1));
    tevend = dsearchn(times', trange(2));
    tevdne = dsearchn(times', -1*trange(2));
    tindex = [tevdne tzero tevend];
    tslice = {times(tevdne:tzero), times(tzero:tevend)};
end

function [pcoefs, polyn, prwdcoefs, prwdpolyn] = calc_slopes(EEGev, chan, trange)
    
    petrials  = permute(EEGev.data(chan,:,:), [2 3 1]); % (lentime, numtrials)
%     surpetrials = permute(EEGev.data(13,:,:), [2 3 1])';
    times     = EEGev.times;
    tevdne    = trange(1);
    tzero     = trange(2);
    tevend    = trange(3);
    numtrials = size(petrials, 2);
    lentime   = (tevend-tzero)+1; % <<<<<<%%%
    tslice    = {times(tevdne:tzero), times(tzero:tevend)};
    
    pcoefs    = zeros(numtrials, 2);
    prwdcoefs = zeros(numtrials, 2);
    polyn     = zeros(lentime, numtrials); %cell(numtrials, 1);
    prwdpolyn = zeros(lentime, numtrials);
    
    for trial=1:numtrials
%         [pcoefs(trial,:), polyn(:,trial), ~] = polyslope__(petrials(:,trial), times);
        pcoefs(trial,:)    = polyfit(tslice{2}, petrials(tzero:tevend, trial), 1);
        polyn(:,trial)     = polyval(pcoefs(trial,:), tslice{2}); % + pcoefs(trial,2);
        prwdcoefs(trial,:) = polyfit(times(tevdne:tzero), petrials(tevdne:tzero, trial), 1);
        prwdpolyn(:,trial) = polyval(prwdcoefs(trial,:), times(tevdne:tzero));
    end

end

function time2plot = show_ERPSpec__(SBJ_DT, rootpath, chans, ev_range, f_range, taskinfo, lang)

    if nargin < 6
        taskinfo = ["nfb", "speccom"];
    end
    
    endlabel = taskinfo(2);
    subj_list = [ SBJ_DT.names ];
    infodata.taskinfo = taskinfo;
    
    if (nargin < 7 || strcmp(lang, 'ptbr'))
        infodata.figtitle = { "Dado " + "Original " + upper(taskinfo(1)), "Dado \it{Surrogates}"};
    elseif strcmp(lang, 'en')
        infodata.figtitle = { "Original " + upper(taskinfo(1)) + " Data", "Surrogates Data"};
    else
        infodata.figtitle = { "Dado " + "Original " + upper(taskinfo(1)), "Dado \it{Surrogates}"};
    end
    
    for suj=1:length(subj_list)
        
        EEGev = get_taskdata__( SBJ_DT(suj), ev_range, taskinfo(1) );
        
        for sess=1:length(SBJ_DT(suj).filespath)
            splitpath = strsplit(SBJ_DT(suj).filespath{sess}, {'\'});
            filename = splitpath{end};
            
            [time2plot, currfig] = vis.show_erpspectrum(EEGev(sess), ...
                                                        EEGev(sess).srate, ...
                                                        ev_range, ...
                                                        f_range, ...
                                                        chans, ...
                                                        infodata, lang);
            
            vis.savefigure(currfig, rootpath, filename, endlabel);
            close(currfig);
        end
    end
end

function EEGev = get_taskdata__(SUJDT, tslice, taskinfo)
    if isempty(taskinfo)
        taskinfo = "nfb";
    end
    
    switch taskinfo
        case "nfb"
            EEGev = SUJDT.events;
            return
        case "bs1"
            bsdt = [ SUJDT.tasks.bs1 ];
        case "bs2"
            bsdt = [ SUJDT.tasks.bs2 ];
        otherwise
            EEGev = SUJDT.events;
            return
    end
    % Cria segmentos, tslice(2) segundos, no baseline selecionado
    for s=1:length(bsdt)
        EEGepch = eeg_regepochs( bsdt(s), 'recurrence', tslice(2), ...
                                  'limits', tslice, ...
                                  'eventtype', 'bsevent');
        EEGepch = fntools.rejtrends(EEGepch);
        EEGev(s) = fntools.rejspec(EEGepch);
    end
end


function show_evPETrials__(SBJ_DT, rootpath, ev_range, lang)

    subj_list = [ SBJ_DT.names ];

    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEGev = SBJ_DT(suj).events;

        for sess=1:length(filepath)
            splitpath = strsplit(filepath{sess}, {'\'});
            filename = splitpath{end};

            pefig = vis.show_evtrials(EEGev(sess), ev_range, lang);
            vis.savefigure(pefig, rootpath, filename, "htrials");
            close(pefig);
        end
    end
end

function showBandsBoxchart__(SBJ_DT, rootpath, lang)
    
    subj_list = [ SBJ_DT.names ];
    
    if (nargin < 3 || strcmp(lang, 'ptbr'))
        bandlabels = {'Teta' 'RSM' 'Beta-Alto' 'Alfa'};
    elseif strcmp(lang, 'en')
        bandlabels = {'Theta' 'SMR' 'Hibeta' 'Alpha'};
    else
        bandlabels = {'Teta' 'RSM' 'Beta-Alto' 'Alfa'};
    end
%     labels = SBJ_DT.filespath;
    for suj=1:length(subj_list)
        filename = SBJ_DT(suj).names;
        
        boxfig = vis.gen_boxchart(SBJ_DT(suj), bandlabels, lang);
        
        for band=1:length(bandlabels)
            vis.savefigure(boxfig{band}, rootpath, filename, "box-"+bandlabels{band});
            close(boxfig);
        end
    end
end

function showCoefsBoxchart__(SBJ_DT, rootpath, lang)
    
    subj_list = [ SBJ_DT.names ];
%     bandlabels = {'Theta' 'SMR' 'Hibeta' 'Alpha'};
    
    for suj=1:length(subj_list)
        filename = SBJ_DT(suj).names;
%         aux      = cellfun(@(x) split(x, "_"), SBJ_DT(suj).filespath,'UniformOutput',false);
%         labels   = cellfun(@(x) x{2}, aux,'UniformOutput',false);
        
        boxfig = vis.gen_coefsboxchart(SBJ_DT(suj), lang);
        
%         for band=1:length(bandlabels)
        vis.savefigure(boxfig, rootpath, filename, "slopebox");
%         close(boxfig);
%         end
    end
end

function showPEBaselines__(SBJ_DT, rootpath, lang)

    subj_list = [ SBJ_DT.names ];

    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEG = SBJ_DT(suj).tasks;
        metrics = SBJ_DT(suj).metrics;

        for sess=1:length(filepath)
            splitpath = strsplit(filepath{sess}, {'\'});
            filename = splitpath{end};
            basefig = vis.show_pebaselines(EEG(sess), metrics(sess), lang);
            vis.savefigure(basefig, rootpath, filename, "pebases");
            close(basefig);
        end
    end
end

function show_PSDSpec__(SBJ_DT, rootpath, channels, f_range, taskinfo, lang)

    if nargin < 5
        taskinfo = ["nfb", "psdnfb"];
    end
    
    if (nargin < 6 || strcmp(lang, 'ptbr'))
        chan = 'Canal C3';
        surchan = 'Dado \it{Surrogate}';
    elseif strcmp(lang, 'en')
        chan = 'C3 channel';
        surchan = 'Surrogate';
    else
        chan = 'Canal C3';
        surchan = 'Dado \it{Surrogate}';
    end
    
    ev_range = [0 1];
    endlabel = taskinfo(2);
    subj_list = [ SBJ_DT.names ];
    
    for suj=1:length(subj_list)
        filepath = SBJ_DT(suj).filespath;
        EEGev = get_taskdata__( SBJ_DT(suj), ev_range, taskinfo(1) );

        for sess=1:length(filepath)
            splitpath = strsplit(filepath{sess}, {'\'});
            filename = splitpath{end};

            specfig = figure('visible','off');
            pop_spectopo(EEGev(sess), 1, [], 'EEG', ...
                         'freqrange', f_range, ...
                         'plotchans', channels);
            legend(chan, surchan);
            
            vis.savefigure(specfig, rootpath, filename, endlabel);
            close(specfig);
        end
    end
end

%% ------------------------------------------------------------------------
