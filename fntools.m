classdef fntools
    methods (Static)
%         function [Hdatas, EEG] = pesurdata(EEG, filename, istherefile, savesurdata) %, delay, order, windowSize)
        function [EEG] = pesurdata(EEG, filename, savesurpath, loadsurdata, savesurdata)
            
            if nargin < 4   % if the number of inputs equals 2
              loadsurdata = true; % then make the third value, z, equal to my default value, 5.
            end
            if nargin < 5
              savesurdata = false;
            end
            
            sbjname = strsplit(filename, {'-'});
            sbjname = sbjname(1);
%             surfolder = "datas/" + sbjname + "/sur/";
%             surfolder = string(savesurpath) + sbjname + "\sur\";
            surfolder = join([savesurpath, sbjname, "sur"], "\");
            
            [~, npnts] = size(EEG.data);
%             Hdatas = struct();
            permH = zeros(1,npnts);
            surpermH = zeros(1,npnts);
            surdata = zeros(1,npnts);
            labels = {'Surrogate', 'PermH', 'SurPermH'};

            delay = 1; % delay 1 between points in ordinal patterns (successive points)
            order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
            windowSize = floor(2*EEG.srate); %5*EEG.srate;
            
            surdatapath = join([surfolder, filename], "\");
            disp(newline + ">> " + surdatapath  + ".mat");
            if ~isfile(surdatapath + ".mat")
                if ~isfolder(surfolder)
                    mkdir(char(surfolder));
                end
                loadsurdata = false;
            end
            
            if loadsurdata
                clear surdata;
                disp(">> Loading surrogates from file..." + filename)
                load(surdatapath + ".mat", 'surdata');
            else
                disp(">> Calculating Surrogates..." + filename);
                surdata = IAAFTsur(EEG.data(1,:), 1);
%                 savesurdata = true;
            end
            
            permH(1:npnts-windowSize-2) = PE(EEG.data(1,:)', delay, order, windowSize);
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
                disp(">> Saving surDatas to..." +newline+ surdatapath);
%                 save(surfolder + filename, 'surdata');
                save(surdatapath, 'surdata');
            end
        end
% -------------------------------------------------------------------------
        function [EEG] = create_events(EEG, chan, latency, RWD_label, ...
                                       ev_duration, ev_range, ...
                                       filename, rejspecevent)

            if isempty(RWD_label)
                RWD_label = 'RWD-250';
            end
            
            if isempty(ev_duration)
                ev_duration = '.250<=2';
            end
            
            if isempty(latency)
%                 latency = '60<=1860';
                latency = '0<=1800';
            end
            
            % define os eventos  disp([newline 'Def events']);
            EEG = pop_chanevent( EEG, chan, ...
                                 'edge', 'leading', ...
                                 'duration', 'on',...
                                 'typename', 'RWD', ...
                                 'delchan', 'off', ...
                                 'edgelen', 1);

            % Seleciona eventos duracao maior que 250 ms
            EEG = pop_selectevent( EEG, 'latency', latency, ...
                                   'duration', ev_duration,...
                                   'renametype', RWD_label, ...
                                   'oldtypefield', 'RWD', ...
                                   'deleteevents', 'on');
            EEG = pop_epoch( EEG, { RWD_label }, ev_range, ...
                             'newname', [filename '-epochs'], ...
                             'epochinfo', 'yes');
            if rejspecevent
                [EEG, rejindx] = pop_rejspec( EEG, 1, ...
                                        'elecrange', [1:1], ...
                                        'threshold', [-30 30], ...
                                        'freqlimits', [1 100], ...
                                        'method', 'fft', ...
                                        'eegplotreject', 1, ...
                                        'eegplotplotallrej', 0);
                EEG.rejindices = rejindx;
            end
        end
% -------------------------------------------------------------------------
        function perm_entr = pe_bytrials(EEG, chan)
%             addpath('PE');
            
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

                sbj_dt(indx).names = sujname;
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
        function bandsPWR = calc_bandpower(EEGev, frange)
            
            if nargin<2 || isempty(frange)
                frange = {[4 7] [12 15] [20 30] [8 12]};
            end

            numBands = length(frange);
            bandsPWR = cell(1, numBands);
            dt = permute(EEGev.data(1,:,:), [2 3 1]);

            for bnd=1:numBands
                bandsPWR{bnd} = bandpower(dt, EEGev.srate, frange{bnd});
%                 disp(["BAND: " + bnd]);
            end
            bandsPWR = vertcat(bandsPWR{:})';
        end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    end
end