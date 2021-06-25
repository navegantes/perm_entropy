classdef fntools
    methods (Static)
        function [te, se] = calcSpecH(EEG) %te, specH, rhpe, ccoef)
    
            dt       = EEG.nfb.data( 1,:,:);
            times    = EEG.nfb.times;
            srate    = EEG.nfb.srate;
            numpnts  = EEG.nfb.pnts;
            timeres  = (numpnts/srate)/10
            
            [p,fp,tp] = pspectrum(dt, times,'spectrogram', ...
                                  'TimeResolution', timeres, ...
                                  'Leakage',0.85, ... %hann window
                                  'OverlapPercent', 50, ...
                                  'FrequencyLimits', [0 100]);
                              
            [se,te] = pentropy(p,fp,tp);
           
        end
% -------------------------------------------------------------------------
        function [EEG] = pesurdata(EEG, filename, savesurpath, loadsurdata, savesurdata)
            
            if nargin < 4   % if the number of inputs equals 2
              loadsurdata = true; % then make the third value, z, equal to my default value, 5.
            end
            if nargin < 5
              savesurdata = false;
            end
            
            sbjname = strsplit(filename, {'-'});
            sbjname = sbjname(1);
            surfolder = join([savesurpath, 'surdatas', sbjname, "sur"], "\");
            
            [~, npnts] = size(EEG.data); %             Hdatas = struct();
            permH      = zeros(1,npnts);
            surpermH   = zeros(1,npnts);
            surdata    = zeros(1,npnts);
%             SpecH      = zeros(1,npnts);
%             SurSpecH   = zeros(1,npnts);

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
                disp(">> Loading surrogates from file..." +newline+ filename);
                load(surdatapath + ".mat", 'surdata');
            else
                disp(">> Calculating Surrogates..." + filename);
                surdata = IAAFTsur(EEG.data(1,:), 1);
            end
            
            permH(1:npnts-windowSize-2)    = PE(EEG.data(1,:)', delay, order, windowSize);
            surpermH(1:npnts-windowSize-2) = PE(surdata', delay, order, windowSize);
            SpecH    = fntools.interpSpecH__(EEG, 1);  % 1 - data chan
            data = {surdata, permH, surpermH, SpecH};
            labels = {'Surrogate', 'PermH', 'SurPermH', 'SpecH'};
            
            [row, ~] = size(EEG.data);
            if row < 11
                for i=1:4
                    EEG.data(end+1,:) = data{i};
                    EEG.nbchan = EEG.nbchan + 1;
                    EEG.chanlocs(end+1).labels = labels{i};
                end
            end
            
%             SpecH    = interpSpecH(EEG, 1);  % 1 - data
%             SurSpecH = interpSpecH(EEG, 11); % 11 - surrogate chan
%             slabel = {'SpecH', 'SurSpecH'};
            
            EEG.data(end+1,:) = fntools.interpSpecH__(EEG, 11); % 11 - surrogate chan
            EEG.nbchan = EEG.nbchan + 1;
            EEG.chanlocs(end+1).labels = 'SurSpecH';
            
            if savesurdata
                disp(">> Saving surDatas to..." +newline+ surdatapath);
                save(surdatapath, 'surdata');
            end
        end
% -------------------------------------------------------------------------
        function  se = interpSpecH__(EEG, chan) %te, specH, rhpe, ccoef)

            dt       = EEG.data( chan,:,:);
            times    = EEG.times;
            srate    = EEG.srate;
            numpnts  = EEG.pnts;

            [p,fp,tp] = pspectrum(dt, times,'spectrogram', ...
                                  'TimeResolution', numpnts/srate, ...
                                  'Leakage',0.85, ... %hann window
                                  'OverlapPercent', 50, ...
                                  'FrequencyLimits', [0 100]);

            [se,te] = pentropy(p,fp,tp);

            se = interp1(te, se, times, 'spline');
        %             te = times;
        end
% -------------------------------------------------------------------------
        function [EEG] = create_events(EEG, rwdchan, latency, RWD_label, ...
                                       ev_duration, ev_range, ...
                                       filename, rejspecevent)

            if isempty(RWD_label)
                RWD_label = 'RWD-250';
            end
            
            if isempty(ev_duration)
                ev_duration = '.250<=2';
            end
            
            if isempty(latency)
                latency = '0<=1800';
            end
            
            % define os eventos  disp([newline 'Def events']);
            EEG = pop_chanevent( EEG, rwdchan, ...
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
                         
            % Rejection Trends
            % OUTEEG = pop_rejtrend( INEEG, typerej, elec_comp, winsize, maxslope, minR, superpose, reject,calldisp);
            disp("..." +newline+ ">> Trend rejection...")
            winsize = floor(0.3 * EEG.srate);
            rejmarked = 1;
            EEG = pop_rejtrend(EEG, 1, [1], winsize, .6, .35, 1,rejmarked, 0);
            
            disp("..." +newline+ ">> Spec rejection...")
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
        function EEG = rejtrends(EEG)
            % Rejection Trends
            % OUTEEG = pop_rejtrend( INEEG, typerej, elec_comp, winsize, maxslope, minR, superpose, reject,calldisp);
            disp("..." +newline+ ">> Trend rejection...")
            winsize = floor(0.25 * EEG.srate);
            rejmarked = 1;
            EEG = pop_rejtrend(EEG, 1, [1], winsize, .6, .65, 1,rejmarked, 0);
        end

        function EEG = rejspec(EEG)
            disp("..." +newline+ ">> Spec rejection...")

            [EEG, rejindx] = pop_rejspec( EEG, 1, ...
                                          'elecrange', [1:1], ...
                                          'threshold', [-30 30], ...
                                          'freqlimits', [1 100], ...
                                          'method', 'fft', ...
                                          'eegplotreject', 1, ...
                                          'eegplotplotallrej', 0);
            EEG.rejindices = rejindx;
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
        function bandsPWR = calc_bandpower(EEGev, chan, frange)
            
            if nargin<2 || isempty(frange)
                frange = {[4 7] [12 15] [20 30] [8 12]};
            end
            
            numBands      = length(frange);
            bandsPWR      = {cell(1, numBands), cell(1, numBands)};
            norwdbandsPWR = cell(1, numBands);
            
            dt    = permute(EEGev.data(chan,:,:), [2 3 1]);
            times = EEGev.times;
            
            trange = [0 250];
            tindx  = get_evtimeindxs(times, trange);
            tdne   = tindx(1);
            tzero  = tindx(2);
            tend   = tindx(3);
            
            totpwr      = bandpower(dt(tzero:tend,:), EEGev.srate, [1 100]);
            norwdtotpwr = bandpower(dt(tdne:tzero,:), EEGev.srate, [1 100]);

%             for rwd=1:2
            for bnd=1:numBands
                pwrBand          = bandpower(dt(tzero:tend,:), EEGev.srate, frange{bnd});
                norwdpwrBand     = bandpower(dt(tdne:tzero,:), EEGev.srate, frange{bnd});
                bandsPWR{1}{bnd} = norwdpwrBand./norwdtotpwr;
                bandsPWR{2}{bnd} = pwrBand./totpwr;
                
%                 disp(["BAND: " + bnd]);
            end
%             end
            bandsPWR{1}      = vertcat(bandsPWR{1}{:})';
            bandsPWR{2}      = vertcat(bandsPWR{2}{:})';
%             norwdbandsPWR = vertcat(norwdbandsPWR{:})';
        end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    end
end

% -------------------------------------------------------------------------
function tindex = get_evtimeindxs(times, trange)
    t_zero  = dsearchn(times', trange(1));
    t_evend = dsearchn(times', trange(2));
    t_evdne = dsearchn(times', -1*trange(2));
    tindex = [t_evdne t_zero t_evend];
end
% -------------------------------------------------------------------------
