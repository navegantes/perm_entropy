classdef fntools
    methods (Static)
%         function [Hdatas, EEG] = pesurdata(EEG, filename, istherefile, savesurdata) %, delay, order, windowSize)
        function [EEG] = pesurdata(EEG, filename, istherefile, savesurdata)
            
            if nargin < 3   % if the number of inputs equals 2
              istherefile = true; % then make the third value, z, equal to my default value, 5.
            end
            if nargin < 4
              savesurdata = false;
            end
            
            sbjname = strsplit(filename, {'-'});
            sbjname = sbjname(1);
            surfolder = "datas/" + sbjname + "/sur/";
            
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
                if ~isfolder(surfolder)
                    mkdir(char(surfolder));
                end
                istherefile = false;
            end
            
            if istherefile
                clear surdata;
                disp(">> Loading surrogates from file..." + filename)
                load(surfolder + filename + ".mat", 'surdata');
            else
                disp(">> Calculating Surrogates..." + filename);
                surdata = IAAFTsur(EEG.data(1,:, :), 1);
                savesurdata = true;
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
            
%             Hdatas.surdata = surdata;
%             Hdatas.permH = permH;
%             Hdatas.surpermH = surpermH;
            
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
        function sbj_dt = gendatastruct(dataPath, sbj_dt, subj_list)
            
            for indx=1:length(subj_list)
                sujname = subj_list(indx);
                nfbpath = dataPath + sujname + '\NFB\';
                dir_info = dir(char(nfbpath));
                sess_dir = {dir_info.name};

                sbj_dt(indx).name = sujname;
                sbj_dt(indx).nfb_path = nfbpath;
                sbj_dt(indx).folders = sess_dir(3:end)';
                sbj_dt(indx).dir_info = dir_info;

                for sess=1:length(sbj_dt(indx).folders)
                    folder = sbj_dt(indx).folders(sess);
                    splt = split(folder, '-');
                    ffilename = sbj_dt(indx).nfb_path + folder{1}+'\'+sujname+...
                                '-'+splt(2)+'_S'+num2str(sess);
                    sbj_dt(indx).filespath{sess} = ffilename;
                end
                sbj_dt(indx).filespath = sbj_dt(indx).filespath';
            end
        end
% -------------------------------------------------------------------------
    end
end