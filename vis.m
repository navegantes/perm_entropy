classdef vis
    methods (Static)
        function fig = show_hspecperm(EEG)
            
            [te, specH, dwnPE] = fntools.getHSpecPerm__(EEG);
            ccoef = corrcoef(dwnPE(1:end-1), specH);
            
            fig = figure; %('visible','off');
            plot(te,specH, 'LineWidth',2 );
            hold on;
            plot(te,dwnPE(1:end-1), 'LineWidth',2);
            
            legend({'Entropia Espectral','Entropia Permutação'})
            infodata = sprintf("Correlation: %7.4f", string(ccoef(1,2)));
            annotation('textbox', ...
                       'String',infodata, ...
                       'Vert','bottom', ...
                       'Horiz', 'right', ...
                       'FaceAlpha', .5, ...
                       'BackgroundColor', 'white', ...
                       'FitBoxToText','on');
            hold off;
            
        end
        
        function [time2plot, fig] = show_erpspectrum(EEG, fs, t_range, f_range, chans, infodata)
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
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
%                 if ~isempty(figtitle{d})
                title(figtitle{d});
%                 end
                caxis([mM(1,1), mM(2,2)]);
                colorbar;
            end
            hold off;
        end
% -------------------------------------------------------------------------
        function [fig] = show_evtrials(EEGev, ev_range)
            
            Hperm = permute(EEGev.data(12,:,:), [2 3 1]);
            HpermSur = permute(EEGev.data(13,:,:), [2 3 1]);
            
            mM = minmax(minmax([Hperm; HpermSur])');

            [len_time, trials] = size(Hperm);
            time = ev_range(1):1/EEGev.srate:ev_range(2);
            time = time(1:len_time)*1000;
            frex = 1:trials; %logspace(log10(1), log10(trials), trials);

            fig = figure('visible','off');
            
            fig.Position = [200 200 1080 720];
            subplot(1,2,1)
            contourf(time, frex, Hperm', 200, 'linecolor', 'none');
            xlabel('Time (ms)');
            ylabel('Trials');
            title('Permutation Entropy ERP');
            caxis([mM(1,1), mM(2,2)]);
            colorbar;

            % figure;
            subplot(1,2,2)
            contourf(time, frex, HpermSur', 200, 'linecolor', 'none');
            caxis([mM(1,1), mM(2,2)]);
            colorbar;
            xlabel('Time (ms)');
            ylabel('Trials');
            title('PE Surrogate ERP');
            % [time, frex] = fntools.show_evspec(permute(Hperm, [2 1]), EEGev.srate, ev_range, f_range, 'Reward Data');
        end
% -------------------------------------------------------------------------
        function EEG = show_erp_events(EEG, RWD_label)

            smooth = 5;
            chans = [1 4 5 6 7 8 9 10 11 12 13];
            band_labels = {'EEG' 'Theta Amp' 'SMR Amp' 'Hibeta Amp' 'Theta Power'...
                           'SMR Power' 'Hibeta Power' 'Alpha Power', ...
                           'Surrgate', 'EntrPerm', 'SurrEntrPerm'};
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
        function fig = show_pebaselines(EEG, metrics)
            
            NFB = 1;
            BS1 = 2;
            BS2 = 3;
            
%             cornfb = metrics.corr{1,NFB}(1,2);
%             corbs2 = metrics.corr{1,BS2}(1,2);
            corbs = [ metrics.corr{1,NFB}(1,2)
                      metrics.corr{1,BS1}(1,2)
                      metrics.corr{1,BS2}(1,2) ];
            wass = metrics.wass;
% ----------------------------------------------------------------------------------
            fig = figure('visible','off');
            fig.Position = [200 200 1080 720];

            cursbplot = subplot(2,2,1);
            plot(EEG.bs1.times./1000, EEG.bs1.data(12,:,:), 'linewidth', 1.2); hold on;
            plot(EEG.bs1.times./1000, EEG.bs1.data(13,:,:), 'linewidth', 1.2); hold off;
            
            title('Baseline - Pre');
            xlabel('Time (s)');
            ylabel('Permutation Entropy');
            lgnd = legend('Data','Surrogate');
            set(lgnd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1 1 1 0.5]'));
            setannon__(cursbplot, corbs(BS1), wass(BS1));
% ----------------------------------------------------------------------------------
            cursbplot = subplot(2,2,2);
            plot(EEG.bs2.times./1000, EEG.bs2.data(12,:,:), 'linewidth', 1.2); hold on;
            plot(EEG.bs2.times./1000, EEG.bs2.data(13,:,:), 'linewidth', 1.2); hold off;
            
            title('Baseline - Pos');
            xlabel('Time (s)');
            lgnd = legend('Data','Surrogate');
            set(lgnd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1 1 1 0.5]'));
            setannon__(cursbplot, corbs(BS2), wass(BS2));
% ----------------------------------------------------------------------------------
            cursbplot = subplot(2,2,[3,4]);
            plot(EEG.nfb.times./1000, EEG.nfb.data(12,:,:), 'linewidth', 1.2); hold on;
            plot(EEG.nfb.times./1000, EEG.nfb.data(13,:,:), 'linewidth', 1.2); hold off;
            
            title('NFB Task');
            xlabel('Time (s)');
            ylabel('Permutation Entropy');
            lgnd = legend('Data','Surrogate');
            set(lgnd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1 1 1 0.5]'));
            setannon__(cursbplot, corbs(NFB), wass(NFB));
        end
% -------------------------------------------------------------------------
        function savefigure(fig, rootpath, filename, figlabel)
            
            figDir = join([rootpath, "fig"], "\");
%             if ~isfolder(figDir)
%                 mkdir(char(figDir));
%             end
            
            sbjName = strsplit(filename, "-");
            sbjName = sbjName{1};
            sbjFigDir = join([figDir, sbjName], "\");
            
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
        function boxfig = gen_boxchart(SBJ_DT, bandlabels)
            % BOXPLOT
            bndPWR = SBJ_DT.bandspower;
            lenband = length(bandlabels);
            boxfig = cell(1,length(bandlabels));

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
                xlabel('Sessions');
                ylabel('Relative Band Power');
            end
        end

% -------------------------------------------------------------------------
    end
end

function setannon__(fig, corbs, wass)

%     dim = [.2 .5 .3 .3];
    infodata = [ sprintf("Correlation: %7.4f", string(corbs)), ...
                 sprintf("Wass Dist: %7.4f", string(wass)) ];
    annotation('textbox', ... %dim, ...
               'String',infodata, ...
               'Vert','bottom', ...
               'Horiz', 'right', ...
               'FaceAlpha', .5, ...
               'BackgroundColor', 'white', ...
               'Position', fig.Position, ...
               'FitBoxToText','on');   
end
