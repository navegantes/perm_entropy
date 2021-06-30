classdef vis
methods (Static)
    
    function [fig, fign] = show_numtrials(data, ntrials, info, lang)
        
        ttitle = info.title;
        legs = info.legend;
    
        if (nargin < 4 || strcmp(lang, 'ptbr'))
            y_label  = 'Numero \it{Trials}';
            yn_label = 'Numero Relativo \it{Trials}';
            x_label  = 'Sessões';
        elseif strcmp(lang, 'en')
            y_label  = 'Number of Trials';
            yn_label = 'Relative Number of Trials';
            x_label  = 'Sessions';
        else
            y_label  = 'Numero \it{Trials}';
            yn_label = 'Numero Relativo \it{Trials}';
            x_label  = 'Sessões';
        end

        fig = figure;
        plot(data, 'o-', 'LineWidth', 1.7, 'MarkerFaceColor', 'auto', 'MarkerSize', 4);
    %     p.MarkerFaceColor = 'auto';
%         legs = {'EP', 'EPSu', 'EE', 'EESu'};
        legend(legs);
        xticks([1 2 3 4 5]);
        xticklabels({"S1 " + "("+ string(ntrials(1)) + ")", ...
                     "S2 " + "("+ string(ntrials(2)) + ")", ...
                     "S3 " + "("+ string(ntrials(3)) + ")", ...
                     "S4 " + "("+ string(ntrials(4)) + ")", ...
                     "S5 " + "("+ string(ntrials(5)) + ")"});
        xlim([0.5 7]);
        ylabel(y_label, 'FontSize', 14);
        xlabel(x_label, 'FontSize', 14);
    %     title(ttitle);

        fign = figure;
        plot(data./ntrials, 'o-', 'LineWidth', 1.7, 'MarkerFaceColor', 'auto', 'MarkerSize', 4);
%         legs = {'EP', 'EPSu', 'EE', 'EESu'};
        legend(legs);
        xticks([1 2 3 4 5]);
        xticklabels({"S1 ($ \frac{" + string(data(1)) + "}{" + string(ntrials(1)) + "}$)", ...
                     "S2 ($ \frac{" + string(data(2)) + "}{" + string(ntrials(2)) + "}$)", ...
                     "S3 ($ \frac{" + string(data(3)) + "}{" + string(ntrials(3)) + "}$)", ...
                     "S4 ($ \frac{" + string(data(4)) + "}{" + string(ntrials(4)) + "}$)", ...
                     "S5 ($ \frac{" + string(data(5)) + "}{" + string(ntrials(5)) + "}$)"} );
        xlim([0.5 7]);
        ylabel(yn_label, 'FontSize', 14);
        xlabel(x_label,  'FontSize', 14);
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
%     ax.FontSize = 8;
%     title(ttitle);
    end
% -------------------------------------------------------------------------
    function show_entropies(EEGsuj, fig, sess, annonbox, lang)

        EEG      = EEGsuj.tasks(sess);
        EEGev    = EEGsuj.events(sess);
%         pe       = EEG.nfb.data(12,:,:);
%         surpe    = EEG.nfb.data(13,:,:);
        times    = EEG.nfb.times;
%         [te, se] = fntools.calcSpecH(EEG);
        
%         lat = [EEGev.event.latency];
%         for ev=1:length(EEGev.event)
%             plot([times(lat(ev)) times(lat(ev))],[0 1], 'b--' ); hold on;
%         end
        
        if (nargin < 5 || strcmp(lang, 'ptbr'))
            y_label = "Entropia";
            subtitle = "Sessão ";
            x_label = 'Tempo (segundo)';
            legs = {'EP','EE', 'EPSu', 'EESu'};
        elseif strcmp(lang, 'en')
            y_label = "Entropy";
            subtitle = "Session ";
            x_label = 'Time (second)';
            legs = {'PE','SE', 'SuPE', 'SuSE'};
        end
        
        chanorder = [12, 14, 13, 15];
        for i=1:4
            chan = chanorder(i);
            plot(times/1000, EEG.nfb.data(chan,:,:), 'LineWidth',1.2); hold on;
            xlim([-50 2150]);
        end
        
        ylabel(y_label);
        title(subtitle + string(sess));
        
%         if sess~=5
%            set(fig,'xtick',[]);
%         else
        if sess==5
            xlabel(x_label);
        end
        
        if sess==1
            legend('show', legs, ... 
                   'Location', 'northeast', ...
                   'FontWeight','bold', ...
                   'FontSize', 8, ...
                   'EdgeColor','none', ...
                   'Color','none');
        end
        
%         intspecH = interp1(te,se, times, 'spline');
%         ccoef = corrcoef(pe, intspecH);
%         infodata = {sprintf("r: %7.4f", string(ccoef(1,2))), 0, annonbox};
%         setannon__(fig, infodata)
        hold off;
    end
% -------------------------------------------------------------------------
    function fig = show_slopesavg(SBJ_DT, chan, coefs, polyn, tmslope, lang)
        
        mpcoefs     = coefs{1};
        prwdpcoefs  = coefs{2};
        mpolyn      = polyn{2};
        prwdmpolyn  = polyn{1};
        trange      = polyn{3};
    %     subj_list = [ SBJ_DT.names ];
    %     for suj=1:length(subj_list)
        fig = figure('Position', [282,132,960,840]);
        numsess = length(SBJ_DT.filespath);
        mM = [Inf -Inf];
        
        if (nargin < 6 || strcmp(lang, 'ptbr'))
            sesstitlelabel = 'Sessão ';
            y_label = 'Entropia de Permutação';
            x_label = 'Tempo (segundo)';
            prepos = {'pre', 'pos'};
        elseif strcmp(lang, 'en')
            sesstitlelabel = 'Session ';
            y_label = 'Permutation Entropy';
            x_label = 'Time (ms)';
%             legs = {'PE','SE', 'Sur'};
            prepos = {'bf', 'af'};
        end
%         colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], ...
%             [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
        for sess=1:numsess
            EEGev = SBJ_DT.events(sess);
            datasur = EEGev.data(chan,:,:);
            times = SBJ_DT.events(sess).times;
            dtmeansur = mean(datasur, 3);
            mM = check_minmax(dtmeansur, mM);
            
%             sesslabel = SBJ_DT.events(sess).setname(12:13);
            mprwdpcoefs = mean(prwdpcoefs{sess}(:,1),1);
            
            lgnd = sprintf(['\n', sesstitlelabel, '%d'...
                            '\nb_1 ', prepos{1}, ' %4.2e\n', '     ',  prepos{2}, ' %4.2e\n'], ...
                            sess, mprwdpcoefs(1), mpcoefs{sess}(1));
            plot(times, dtmeansur,'linewidth',4, 'DisplayName', lgnd); %sesslabel);
            hold on;
        end
%         title('Slope of the Permutation Entropy Means','FontSize',14);
        legend('show', 'Location', 'bestoutside', ...
               'FontWeight','bold', ...
               'FontSize',11, ...
               'AutoUpdate','off');
        xlabel(x_label, 'FontSize',16);
        ylabel(y_label, 'FontSize',16);
        % plot slopes %
        X = "\bf{Slopes}";
        for sess=1:numsess
            mprwdpolyn = mean(prwdmpolyn{sess}, 2)';
%             X = [X, sprintf("S" +string(sess)+" %6.4g", mpcoefs{sess}(1))];
            plot(tmslope{2}, mpolyn{sess}, 'k-.','linewidth',3); %, 'DisplayName', X);
            hold on;
            plot(tmslope{1}, mprwdpolyn, '-.','linewidth',3, 'Color', [.64 .08 .18]); %, 'DisplayName', X);
            hold on;
        end
%         infodata = {X, 0};
%         setannon__(fig, infodata );
        color = [0 0.4470 0.7410];
        offset = .001;
        edges = [mM(1)-offset mM(1)-offset mM(2)+offset mM(2)+offset];
        ylim([mM(1)-offset mM(2)+offset]);
        xlim([-300 500]);
        p  = patch([0 250 250 0], edges, color, 'LineStyle', ':');
%         p2 = patch([-250 0 0 -250], edges, color, 'LineStyle', '-.');
        p.FaceAlpha  = .08;
%         p2.FaceAlpha = .08;
        plot([0 0], edges(2:3), 'm--','linewidth',3);
        hold off;
    %     end
    end
% -------------------------------------------------------------------------
    function [time2plot, fig] = show_erpspectrum(EEG, fs, t_range, f_range, chans, infodata, lang)
%         function [time, frex, fig] = show_erpspectrum(EEG, fs, t_range, f_range, chans, infodata)

        data = { permute(EEG.data(chans(1),:,:), [2 3 1]), ...
                 permute(EEG.data(chans(2),:,:), [2 3 1])};

        figtitle = infodata.figtitle;
%             figtitle = { 'Original NFB Data', 'Surrogates Data'};

        num_frex= round(f_range(2)-f_range(1));
        frex = logspace(log10(f_range(1)), log10(f_range(2)), num_frex);

        d_spec = cell(1,2);
        time2plot = cell(1,2);
        time = t_range(1):1/fs:t_range(2);
        t_zero = dsearchn(time', 0);
        for i=1:length(data)
           [d_spec{i}, time2plot{i}] = wav_timeFreq(data{i}, fs, f_range(1), f_range(2), t_zero); 
        end

        if infodata.taskinfo{1} ~= "nfb"
            len_time= size(data{1},1);
            time = time(1:len_time).*1000; % to ms
%                 time = linspace(t_range(1), t_range(2), fs);
        end

        mM = minmax(minmax(cell2mat(d_spec))');
        
        if (nargin < 7 || strcmp(lang, 'ptbr'))
            y_label = 'Frequência (Hz)';
            x_label = 'Tempo (ms)';
        elseif strcmp(lang, 'en')
            y_label = 'Frequency (Hz)';
            x_label = 'Time (ms)';
        end

        fig = figure('visible','off');
%             set(fig, 'Visible', 'off');
        fig.Position = [100 100 1080 720];
        for d=1:length(data)
%                 d_spec= wav_timeFreq(data{d}, fs, f_range(1), f_range(2));
            subplot(1,2,d)
            hold on;
%                 contourf(time, frex, d_spec{d}, 20, 'linecolor', 'none');
            if infodata.taskinfo{1} == "nfb"
                contourf(time2plot{d}, frex, d_spec{d}, 20, 'linecolor', 'none');
            else
                contourf(time, frex, d_spec{d}, 20, 'linecolor', 'none');
            end

            plot([0 0], [1 frex(end)], '--k', 'LineWidth',2);
            xlabel(x_label);
            ylabel(y_label);
%                 if ~isempty(figtitle{d})
            title(figtitle{d});
%                 end
            caxis([mM(1,1), mM(2,2)]);
            colorbar;
        end
        hold off;
    end
% -------------------------------------------------------------------------
    function [fig] = show_evtrials(EEGev, ev_range, lang)

        Hperm = permute(EEGev.data(12,:,:), [2 3 1]);
        HpermSur = permute(EEGev.data(13,:,:), [2 3 1]);

        mM = minmax(minmax([Hperm; HpermSur])');

        [len_time, trials] = size(Hperm);
        time = ev_range(1):1/EEGev.srate:ev_range(2);
        time = time(1:len_time)*1000;
        frex = 1:trials; %logspace(log10(1), log10(trials), trials);
        
        if (nargin < 3 || strcmp(lang, 'ptbr'))
            y_label = 'Épocas';
            x_label = 'Tempo (ms)';
            petitle  = 'Entropia Permutação ERP';
            sutitle = 'Entropia \it{Surrogate} ERP';
        elseif strcmp(lang, 'en')
            y_label = 'Trials';
            x_label = 'Time (ms)';
            petitle = 'Permutation Entropy ERP';
            sutitle = 'PE Surrogate ERP';
        end

        fig = figure('visible','off');

        fig.Position = [200 200 1080 720];
        subplot(1,2,1)
        contourf(time, frex, Hperm', 200, 'linecolor', 'none');
        xlabel(x_label);
        ylabel(y_label);
        title(petitle);
        caxis([mM(1,1), mM(2,2)]);
        colorbar;

        % figure;
        subplot(1,2,2)
        contourf(time, frex, HpermSur', 200, 'linecolor', 'none');
        caxis([mM(1,1), mM(2,2)]);
        colorbar;
        xlabel(x_label);
        ylabel(y_label);
        title(sutitle);
        % [time, frex] = fntools.show_evspec(permute(Hperm, [2 1]), EEGev.srate, ev_range, f_range, 'Reward Data');
    end
% -------------------------------------------------------------------------
    function EEG = show_erp_events(EEG, RWD_label)

        smooth = 5;
        chans = [1 4 5 6 7 8 9 10 11 12 13];
        band_labels = {'EEG' 'Theta Amp' 'SMR Amp' 'Hibeta Amp' 'Theta Power'...
                       'SMR Power' 'Hibeta Power' 'Alpha Power', ...
                       'Surrgate', 'EntrPerm', 'SurrEntrPerm'};
%         band_labels = {'EEG' 'Teta Amp' 'SMR Amp' 'Hibeta Amp' 'Teta Power'...
%                        'SMR Power' 'Hibeta Power' 'Alfa Power', ...
%                        'Surrgate', 'EntrPerm', 'SurrEntrPerm'};
%             ylims = {[-1 1] [2 8] [3 5] [1 7] [4 21] [5 9] [1 10] [14 25] [-1 1] [0.8 1] [0.8 1]};

        for i=1:length(chans)
            figure('visible','off');

            pop_erpimage(EEG,1, chans(i),[[]],band_labels{i}, smooth,1,...
                         {RWD_label},[],'', 'yerplabel','\muV', ...
                         'erp',1,'cbar','on', 'erpstd', 'on', 'align', [0], ...
                         'limits',[NaN NaN NaN NaN NaN NaN NaN], 'renorm','yes'); %,...
%                              'plotamps','on');
        end
    end
% -------------------------------------------------------------------------
    function fig = show_pebaselines(EEG, metrics, lang)

        NFB = 1;
        BS1 = 2;
        BS2 = 3;

        corbs = [ metrics.corr{1,NFB}(1,2)
                  metrics.corr{1,BS1}(1,2)
                  metrics.corr{1,BS2}(1,2) ];
        wass = metrics.wass;
    % ----------------------------------------------------------------------------------
        if (nargin < 5 || strcmp(lang, 'ptbr'))
            y_label = 'Entropia de Permutação';
            x_label = 'Tempo (segundo)';
            data = 'Dado';
            surr = '\it{Surrogate}';
            nfbtitle = 'Tarefa (NFB)';
        elseif strcmp(lang, 'en')
            y_label = 'Permutation Entropy';
            x_label = 'Time (ms)';
            data = 'Data';
            surr = 'Surrogate';
            nfbtitle = 'NFB Task';
        end
        
        fig = figure('visible','off');
        fig.Position = [200 200 1080 720];

        cursbplot = subplot(2,2,1);
        plot(EEG.bs1.times./1000, EEG.bs1.data(12,:,:), 'linewidth', 1.2); hold on;
        plot(EEG.bs1.times./1000, EEG.bs1.data(13,:,:), 'linewidth', 1.2); hold off;

        title('Baseline - Pre');
        xlabel(x_label);
        ylabel(y_label);
        lgnd = legend(data, surr);
        set(lgnd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1 1 1 0.5]'));
        infodata = {[ sprintf("Corr: %7.4f", string(corbs(BS1))), ...
                      sprintf("Wass: %7.4f", string(wass(BS1))) ], 1};
        setannon__(cursbplot, infodata); %corbs(BS1), wass(BS1));
    % ----------------------------------------------------------------------------------
        cursbplot = subplot(2,2,2);
        plot(EEG.bs2.times./1000, EEG.bs2.data(12,:,:), 'linewidth', 1.2); hold on;
        plot(EEG.bs2.times./1000, EEG.bs2.data(13,:,:), 'linewidth', 1.2); hold off;

        title('Baseline - Pos');
        xlabel(x_label);
        lgnd = legend(data, surr);
        set(lgnd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1 1 1 0.5]'));
        infodata = {[ sprintf("Corr: %7.4f", string(corbs(BS2))), ...
                      sprintf("Wass: %7.4f", string(wass(BS2))) ], 1};
        setannon__(cursbplot, infodata); %corbs(BS2), wass(BS2));
    % ----------------------------------------------------------------------------------
        cursbplot = subplot(2,2,[3,4]);
        plot(EEG.nfb.times./1000, EEG.nfb.data(12,:,:), 'linewidth', 1.2); hold on;
        plot(EEG.nfb.times./1000, EEG.nfb.data(13,:,:), 'linewidth', 1.2); hold off;

        title(nfbtitle);
        xlabel(x_label);
        ylabel(y_label);
        lgnd = legend(data, surr);
        set(lgnd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1 1 1 0.5]'));
        infodata = {[ sprintf("Corr: %7.4f", string(corbs(NFB))), ...
                      sprintf("Wass: %7.4f", string(wass(NFB))) ], 1};
        setannon__(cursbplot, infodata); %corbs(NFB), wass(NFB));
    end
% -------------------------------------------------------------------------
    function savefigure(fig, rootpath, filename, figlabel, subfolder)

        figDir = join([rootpath, "fig"], "\");
%             if ~isfolder(figDir)
%                 mkdir(char(figDir));
%             end

        sbjName = strsplit(filename, "-");
        sbjName = sbjName{1};
        if nargin < 5
            sbjFigDir = join([figDir, sbjName], "\");
        else
            sbjFigDir = join([figDir, sbjName, subfolder], "\");
        end

        if ~isfolder(sbjFigDir)
            mkdir(char(sbjFigDir));
        end

        fullfilename = filename + "_" + figlabel;
        fullName = join([sbjFigDir, fullfilename], "\");
        figName = join([fullName "png"], ".");

        saveas(fig,figName);

        disp("Figure save at " + figName);
    end
% -------------------------------------------------------------------------
    function boxfig = gen_boxchart(SBJ_DT, bandlabels, lang)
        % BOXPLOT
        bndPWR = SBJ_DT.bandspower;
        lenband = length(bandlabels);
        boxfig = cell(1,length(bandlabels));
        
        if (nargin < 3 || strcmp(lang, 'ptbr'))
            y_label = 'Potência Relativa';
            x_label = 'Sessões';
        elseif strcmp(lang, 'en')
            y_label = 'Relative Band Power';
            x_label = 'Sessions';
        end
        
        for band=1:lenband
            banDatas = cellfun(@(x) x(:,band), bndPWR,'UniformOutput',false);
            bandsvec = vertcat(banDatas{:});

            g = [];
            for sess=1:5
                s = split(SBJ_DT.filespath{sess}, "_");
                g = [g; repmat(s(end),size(banDatas{sess},1),1) ];
            end
            
            boxfig{band} = figure('visible','off');
            boxchart(categorical(g), bandsvec); %,'notch','on');
            title(bandlabels{band});
            xlabel(x_label);
            ylabel(y_label);
        end
    end

% -------------------------------------------------------------------------
function boxfig = gen_coefsboxchart(SBJ_DT, lang)
        % BOXPLOT
%         bndPWR = SBJ_DT.bandspower;
%         lenband = length(bandlabels);
%         boxfig = cell(1,length(labels));

        if (nargin < 2 || strcmp(lang, 'ptbr'))
            y_label = 'Coeficiente Angular';
            x_label = 'Sessões';
            cat = {'Antes Reforço','Depois Reforço'};
        elseif strcmp(lang, 'en')
            y_label = 'Slope Coefficients';
            x_label = 'Sessions';
            cat = {'Non-Reward','Reward'};
        end
        
        slpvec = cell(1,5);
        g = [];
        coefvals = [];
        sessions = [];
        for sess=1:5
            slpcoefs = { SBJ_DT(1).metrics(sess).pslopecoefs;
                         SBJ_DT(1).metrics(sess).slopecoefs  }';
                    
            coefDatas = cellfun(@(x) x(:,1), slpcoefs(1,:),'UniformOutput',false);
            slpvec = vertcat(coefDatas{:});
            coefvals = [coefvals; slpvec];
            
            lenslps = length(slpcoefs);
            sesslabel = split(SBJ_DT(1).filespath{sess}, "_");

            for slp=1:lenslps
                g = [g; repmat(cat(slp),size(coefDatas{slp},1),1) ];
                sessions = [sessions; repmat(sesslabel(2),size(coefDatas{slp},1),1)];
            end
        end
        
        hexcolors = {'#D95319', '#0072BD'};
        
        boxfig = figure('Position', [230,340,960,640]); %('visible','off');
        bc = boxchart(categorical(sessions), coefvals, ...
                      'GroupByColor', g, ...
                      'BoxFaceColorMode', 'manual', ...
                      'MarkerColorMode',  'manual', ...
                      'JitterOutliers','on'); %,'notch','on');
        
        bc(1).BoxFaceColor = hexcolors{1};
        bc(1).MarkerColor  = hexcolors{1};
%         bc(1).SeriesIndex  = 1;
        bc(2).BoxFaceColor = hexcolors{2};
        bc(2).MarkerColor  = hexcolors{2};
%         bc(2).SeriesIndex  = 2;
%             title(labels{slp});
        xlabel(x_label, 'FontSize',16);
        ylabel(y_label, 'FontSize',16);
        legend('show', 'Location', 'northeast', ...
               'FontWeight','bold', ...
               'FontSize',10);
    end

% -------------------------------------------------------------------------
end
end

function setannon__(fig, infodata ) % corbs, wass)

%     dim = [.2 .5 .3 .3];
%     infodata = [ sprintf("Correlation: %7.4f", string(corbs)), ...
%                  sprintf("Wass Dist: %7.4f", string(wass)) ];

    infotxt = infodata{1};
    figposition = infodata{2};
    
    if figposition
        annotation('textbox', 'String', infotxt, ...
                   'Vert','bottom', ...
                   'Horiz', 'right', ...
                   'FaceAlpha', .5, ...
                   'BackgroundColor', 'white', ...
                   'Position', fig.Position, ...
                   'FitBoxToText','on');
    else
        box = infodata{3};
        annotation('textbox', box, ... %[.89, 0.3, 0.1, 0.1], ...
                   'String', infotxt, ...
                   'Vert','bottom', ...
                   'Horiz', 'right', ...
                   'FaceAlpha', .5, ...
                   'BackgroundColor', 'white', ...
                   'FitBoxToText','on', ...
                   'FontSize', 10);
    end
end


function mM = check_minmax(dtmeansur, mM)
    
    curmM = minmax(dtmeansur);
    if curmM(1) < mM(1)
        mM(1) = curmM(1);
    end
    if curmM(2) > mM(2)
        mM(2) = curmM(2);
    end
end

