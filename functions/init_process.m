
function [EEG] = init_process(filepath, tminmax, zscore_norm)
    
    csv_file = char(join([filepath, 'csv'], '.'));
    edf_file = char(join([filepath, 'edf'], '.'));

    spltpath = strsplit(filepath, {'\'});
    setname = spltpath(8);
    
    tmin = tminmax(1); %10;
    tmax = tminmax(2); %1939;

    % -------------------------------------------------------------------------
    %Read from edf file
    disp("..." +newline+ "Reading ..." + edf_file);
    EEG = pop_biosig(edf_file, 'importevent','off');
    EEG = pop_select( EEG,'time',[tmin tmax] );
    % -------------------------------------------------------------------------
    %Read from csv file
    file = csvread(csv_file, 2, 0);
%     file = file';

    rawCSV = pop_importdata( 'setname', setname, ...
                             'data', file', ...
                             'nbchan', 13, ...
                             'dataformat', 'array', ...
                             'srate', EEG.srate);
    
    EEGcsv = pop_select( rawCSV, 'time',[tmin+1 tmax+1] );
    % -------------------------------------------------------------------------
    % Muda as labels dos canais
    labels = {'EEG','ECG'};
    for i=1:2
        EEG.chanlocs(i).labels = labels{i};
    end
    
    % Normalização
    if zscore_norm
        % EEG.data(1, :) = normalize(EEG.data(1, :));
        EEG.data(1, :) = zscore(EEG.data(1, :));
    end

    % Fitragem passa-faixa.
    EEG = pop_eegfiltnew(EEG, 1,100,900,0,[],0);
    % Fitragem notch.
    EEG = pop_eegfiltnew(EEG, 58,62,424,1,[],0);

    % ------------------------------------------------------------------------
    % Adiciona os canais (11, 2, 5, 8) RWD e Amplitudes theta SMR hibeta, do csv.
    labels = {'RWD' 'Theta:Amp' 'SMR:Amp' 'Hibeta:Amp' 'Theta:Power'...
              'SMR:Power' 'Hibeta:Power' 'Alpha:Power'};
    chanlist = [11 2 5 8];
    for chan=1:length(chanlist)
        EEG.data(end+1,:) = EEGcsv.data(chanlist(chan), :);
        EEG.nbchan = EEG.nbchan + 1;
        EEG.chanlocs(chan+2).labels = labels{chan};
    end
    % ------------------------------------------------------------------------
    % Adiciona potencias instantaneas theta, SMR, hbeta, alpha
    frange = {[4 7] [12 15] [20 30] [8 12]};
    lenfmM = length(frange);
    for f=1:length(frange)
        freq = frange{f};
        EEGfilt = pop_eegfiltnew(EEG, freq(1),freq(2),900,0,[],0);
        hilb = hilbert(EEGfilt.data(1, :));
        EEG.data(end+1,:) = abs(hilb).^2;
        EEG.nbchan = EEG.nbchan + 1;
        EEG.chanlocs(f+lenfmM+2).labels = labels{f+lenfmM};
    end
    
    disp("..." +newline+ "Init process done...");
end

