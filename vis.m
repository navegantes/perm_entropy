classdef vis
    methods (Static)
        function [time, frex, fig] = show_erpspectrum(EEG, fs, t_range, f_range, chans)
            
            data = { permute(EEG.data(chans(1),:,:), [2 3 1]), ...
                     permute(EEG.data(chans(2),:,:), [2 3 1])};
            
            figtitle = { 'Original Data', ...
                       'Surrogates Data'};
       
            len_time= size(data{1},1);
            time = t_range(1):1/fs:t_range(2);
            time = time(1:len_time).*1000; % to ms
%             time = linspace(t_range(1), t_range(2), fs);
            num_frex= round(f_range(2)-f_range(1));
            frex = logspace(log10(f_range(1)), log10(f_range(2)), num_frex);
            
            d_spec = cell(1,2);
            for i=1:length(data)
               d_spec{i} = wav_timeFreq(data{i}, fs, f_range(1), f_range(2)); 
            end
            
            mM = minmax(minmax(cell2mat(d_spec))');
            
            fig = figure('visible','off');
%             set(fig, 'Visible', 'off');
            
            fig.Position = [100 100 1080 720];
            for d=1:length(data)
%                 d_spec= wav_timeFreq(data{d}, fs, f_range(1), f_range(2));
                subplot(1,2,d)
                hold on;
                contourf(time, frex, d_spec{d}, 20, 'linecolor', 'none')
                plot([0 0], [1 frex(end)], '--r', 'LineWidth',2);
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                if ~isempty(figtitle{d})
                    title(figtitle{d});
                end
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
%             RWD_label = 'RWD-250';
%             pop_eegplot( EEG, 1, 1, 1);

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
        function fig = show_pebaselines(EEG)
            
            EEGbs1 = pop_select( EEG,'time',[0 60] );
            EEGbs2 = pop_select( EEG,'time',[1860 1920] );
            
            bs1data = { EEGbs1.data(12,:,:), EEGbs1.data(13,:,:) };
            bs2data = { EEGbs2.data(12,:,:), EEGbs2.data(13,:,:) };
            
            corbs1 = corrcoef(bs1data{1}, bs1data{2});
            corbs2 = corrcoef(bs2data{1}, bs2data{2});
            wsbs1 = ws_distance(bs1data{1}, bs1data{2});
            wsbs2 = ws_distance(bs2data{1}, bs2data{2});

            fig = figure('visible','off');
            
            fig.Position = [200 200 1080 720];

            subplot(1,2,1)
            plot(EEGbs1.times./1000, EEGbs1.data(12,:,:));
            hold on;
            plot(EEGbs1.times./1000, EEGbs1.data(13,:,:));
            title(['Baseline - Pre', ...
                   'Correlation: ' + string(corbs1(1,2)), ...
                   'WS Distance: ' + string(wsbs1)] );
            hold off;
            xlabel('Time (s)');
            ylabel('Permutation Entropy');
            legend('Data','Surrogate');

            subplot(1,2,2)
            plot(EEGbs2.times./1000, EEGbs2.data(12,:,:));
            hold on;
            plot(EEGbs2.times./1000, EEGbs2.data(13,:,:));
            title(['Baseline - Pos', ...
                   'Correlation: ' + string(corbs2(1,2)), ...
                   'WS Distance: ' + string(wsbs2)]);
            hold off;
            xlabel('Time (s)');
            legend('Data','Surrogate');
            
            suptitle('PE Baselines');
        end
% -------------------------------------------------------------------------
        function savefigure(fig, rootpath, filename, figlabel)
            
%             currpath = pwd;
%             JLC-240120_S5
            
            sbjName = strsplit(filename, "-");
            sbjName = sbjName{1};
            figDir = join([rootpath, "fig", sbjName], "\");
            
            fullfilename = filename + "_" + figlabel;
            fullName = join([figDir, fullfilename], "\");
            figName = join([fullName "png"], ".");
            
            saveas(fig,figName);
            
            disp("Figure save at " + figName);
        end
% -------------------------------------------------------------------------
    end
end