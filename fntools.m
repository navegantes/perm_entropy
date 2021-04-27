classdef fntools
    methods (Static)
        function [Hdatas, EEG] = hsurdata(EEG) %, delay, order, windowSize)
            
            [~, npnts] = size(EEG.data);
            Hdatas = struct();
            permH = zeros(1,npnts);
            surpermH = zeros(1,npnts);
            surdata = zeros(1,npnts);
            labels = {'Surrogate', 'PermH', 'SurPermH'};

            delay = 1; % delay 1 between points in ordinal patterns (successive points)
            order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
            windowSize = 5*EEG.srate;
            
            surdata = IAAFTsur(EEG.data(1,:, :), 1);
            permH(1:npnts-windowSize-2) = PE(EEG.data(1,:, :)', delay, order, windowSize);
            surpermH(1:npnts-windowSize-2) = PE(surdata', delay, order, windowSize);
            
            data = {surdata, permH, surpermH};

            for i=1:3
                EEG.data(end+1,:) = data{i};
                EEG.nbchan = EEG.nbchan + 1;
                EEG.chanlocs(end+1).labels = labels{i};
            end
            
            Hdatas.surdata = surdata;
            Hdatas.permH = permH;
            Hdatas.surpermH = surpermH;
%             EEG.data(end+1,:) = Hdatas.surdata;
%             EEG.data(end+1,:) = Hdatas.permH;
%             EEG.data(end+1,:) = Hdatas.surpermH;
        end
            
% -------------------------------------------------------------------------
        function EEG = show_erp_events(EEG, RWD_label)
%             RWD_label = 'RWD-250';
            pop_eegplot( EEG, 1, 1, 1);

            smooth = 5;
            chans = [1 4 5 6 7 8 9 10];
            band_labels = {'EEG' 'Theta Amp' 'SMR Amp' 'Hibeta Amp' 'Theta Power'...
                           'SMR Power' 'Hibeta Power' 'Alpha Power'};
            ylims = {[-1 1] [2 8] [3 5] [1 7] [4 21] [5 9] [1 10] [14 25]};

            for i=1:length(chans)
                figure;
                pop_erpimage(EEG,1, chans(i),[],band_labels{i}, ...
                             smooth,1,{RWD_label},[],'epoch' ,...
                             'yerplabel','\muV','erp','on','cbar','on', ...
                             'limits',[NaN NaN ylims{i} NaN NaN NaN NaN] );
            end
        end
% -------------------------------------------------------------------------
        function EEG = create_events(EEG, chan, latency, RWD_label, ev_duration, ev_range)

            if isempty(RWD_label)
                RWD_label = 'RWD-250';
            end
            
            if isempty(ev_duration)
                ev_duration = '.250<=2';
            end
            
%             if length(ev_range)<1
%                 ev_range = [-0.5 0.5];
%             end
            
            if isempty(latency)
                latency = '60<=1860';
            end
            
            
            % define os eventos  disp([newline 'Def events']);
            EEG = pop_chanevent( EEG, chan, 'edge', 'leading', 'duration', 'on',...
                                 'typename', 'RWD', 'delchan', 'off', 'edgelen', 1);

            % Seleciona eventos duracao maior que 250 ms
            EEG = pop_selectevent( EEG, 'latency',latency,'duration',ev_duration,...
                                   'renametype',RWD_label,'oldtypefield','RWD','deleteevents','on');
            EEG = pop_epoch( EEG, { RWD_label }, ev_range, ...
                             'newname', 'JLC-210120_S2-epochs', 'epochinfo', 'yes'); 
        end
% -------------------------------------------------------------------------
        function H_perm = evpermentropy(EEG)
            addpath('PE');
            
            [~, len, numevnts] = size(EEG.data);
            delay = 1; % delay 1 between points in ordinal patterns (successive points)
            order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
            windowSize = round(EEG.srate/12);
            
            H_perm = zeros(numevnts, len);
            for iev=1:numevnts
                h = PE(EEG.data(1,:, iev)', delay, order, windowSize);
                H_perm(iev, 1:length(h)) = h;
                %Computing Entropy
                %EEG.data(end+1,:) = H_perm; 
            end
        end
% -------------------------------------------------------------------------
        function [time, frex] = show_evspec(data, fs, t_range, f_range)
            
%             p_data = permute(data(1,:,:), [2 3 1]);
            len_time= size(data,1);
            time = t_range(1):1/fs:t_range(2);
            time = time(1:len_time);
%             time = linspace(t_range(1), t_range(2), fs);
            num_frex= round(f_range(2)-f_range(1));
            frex = logspace(log10(f_range(1)), log10(f_range(2)), num_frex);
            
            d_spec= wav_timeFreq(data, fs, f_range(1), f_range(2));
            
            figure;
            contourf(time, frex, d_spec, 20, 'linecolor', 'none')
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            colorbar;
        end
    end
end