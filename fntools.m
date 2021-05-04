classdef fntools
    methods (Static)
        function [Hdatas, EEG] = pesurdata(EEG, filename) %, delay, order, windowSize)
            % JLC-240120_S5
%             surfilename = filename{end};
            sbjname = strsplit(filename, {'-'});
            sbjname = sbjname(1);
%             session = strsplit(filename, '_');
%             session = session(2);
            surfolder = "datas/" + sbjname + "/sur/";

            savesurdata = false;
            istherefile = true;
            
            [~, npnts] = size(EEG.data);
            Hdatas = struct();
            permH = zeros(1,npnts);
            surpermH = zeros(1,npnts);
            surdata = zeros(1,npnts);
            labels = {'Surrogate', 'PermH', 'SurPermH'};

            delay = 1; % delay 1 between points in ordinal patterns (successive points)
            order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
            windowSize = floor(2*EEG.srate); %5*EEG.srate;
            
            if ~isfile(surfolder + filename + ".mat")
                mkdir(char(surfolder));
                istherefile = false;
                savesurdata = true;
            end
            
            if istherefile % || updatesurdata
%                 disp('>> Surrogates already loaded...');
                disp('>> Loading surrogates from file...')
                load(surfolder + filename + ".mat", 'surdata');
            else
                disp('>> Calculating Surrogates...');
                surdata = IAAFTsur(EEG.data(1,:, :), 1);
            end
            
%             surdata = IAAFTsur(EEG.data(1,:, :), 1);
            permH(1:npnts-windowSize-2) = PE(EEG.data(1,:, :)', delay, order, windowSize);
            surpermH(1:npnts-windowSize-2) = PE(surdata', delay, order, windowSize);
            
            data = {surdata, permH, surpermH};
            
            [row, ~] = size(EEG.data);
            if row < 11
                for i=1:3
                    EEG.data(end+1,:) = data{i};
                    EEG.nbchan = EEG.nbchan + 1;
                    EEG.chanlocs(end+1).labels = labels{i};
                end
            end
            
            Hdatas.surdata = surdata;
            Hdatas.permH = permH;
            Hdatas.surpermH = surpermH;
            
            if savesurdata
                disp('>> Saving surDatas...')
                save(surfolder + filename, 'surdata');
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
        function perm_entr = pe_bytrials(EEG, chan)
            addpath('PE');
            
            [~, len, numevnts] = size(EEG.data);
            delay = 1; % delay 1 between points in ordinal patterns (successive points)
            order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
            windowSize = floor(.25*EEG.srate); %floor(EEG.srate/12);
            
            perm_entr = zeros(numevnts, len);
            for iev=1:numevnts
                h = PE(EEG.data(chan,:, iev)', delay, order, windowSize);
                perm_entr(iev, 1:length(h)) = h;
                %Computing Entropy
                %EEG.data(end+1,:) = H_perm; 
            end
        end
% -------------------------------------------------------------------------
        function savefigure(gcf, edf_file, figlabel)
            
            currpath = pwd;
            
            sbjName = strsplit(edf_file, "\");
            sbjName = sbjName{5};
            figDir = join([currpath, "fig", sbjName], "\");
            session = strsplit(edf_file(end-5:end), '.');
            session = session(1);

            fileName = join([sbjName, session{1}, figlabel], "_");
            fullName = join([figDir, fileName], "\");
            figName = join([fullName "png"], ".");
            saveas(gcf,figName);
            
            disp("Figure save at " +figName);
        end
% -------------------------------------------------------------------------
    end
end