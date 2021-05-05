
function [EEG] = init_process(filepath, tminmax)
    
    csv_file = char(join([filepath, 'csv'], '.'));
    edf_file = char(join([filepath, 'edf'], '.'));

    spltpath = strsplit(filepath, {'\'});
    setname = spltpath(8);

    file = csvread(csv_file, 2, 0);
    file = file';

    rawCSV = pop_importdata( 'setname', setname, ...
                             'data', file, ...
                             'nbchan', 13, ...
                             'dataformat', 'array', ...
                             'srate', 256);
    tmin = tminmax(1); %10;
    tmax = tminmax(2); %1939;
    EEGcsv = pop_select( rawCSV, 'time',[tmin+1 tmax+1] );

    % -------------------------------------------------------------------------
    %Read from edf file
    EEG = pop_biosig(edf_file, 'importevent','off');
    EEG = pop_select( EEG,'time',[tmin tmax] );

    % Muda as labels dos canais
    labels = {'EEG','ECG'};
    for i=1:2
        EEG.chanlocs(i).labels = labels{i};
    end

    % Fitragem passa-faixa.
    EEG = pop_eegfiltnew(EEG, 1,100,900,0,[],0);

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
    fminmax = {[4 7] [12 15] [20 30] [8 12]};
    lenfmM = length(fminmax);
    for f=1:length(fminmax)
        freq = fminmax{f};
        EEGfilt = pop_eegfiltnew(EEG, freq(1),freq(2),900,0,[],0);
        hilb = hilbert(EEGfilt.data(1, :));
        EEG.data(end+1,:) = abs(hilb).^2;
        EEG.nbchan = EEG.nbchan + 1;
        EEG.chanlocs(f+lenfmM+2).labels = labels{f+lenfmM};
    end
    
    disp('Init process done...');
end

