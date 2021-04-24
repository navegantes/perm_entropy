classdef fn_tools
    methods (Static)
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
                pop_erpimage(EEG,1, chans(i),[],band_labels{i},smooth,1,{RWD_label},[],'epoch' ,...
                    'yerplabel','\muV','erp','on','cbar','on', 'limits',[NaN NaN ylims{i} NaN NaN NaN NaN] );
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
            
            if length(ev_range)<1
                ev_range = [-0.3 0.5];
            end
            
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
    end
end