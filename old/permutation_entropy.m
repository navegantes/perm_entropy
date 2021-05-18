%% -----------------------------------------------------------------
% - entropia de permutação
% - entropia aproximada
% 	- Surrogate
% Teste multi-repos
% Teste2 multi-repos login
    
close all
clear
clc

addpath('PE');
addpath('functions')

dataPath = 'D:\Users\NFB\Pacientes\';
subj_list = ["JLC"]; % ["EYK", "JLC", "JRJ", "SAJ"];
% subj_len = length(subj_list);
hdatas   = cell(5, length(subj_list));
wess_dist = cell(5, length(subj_list));
% surdatas = cell(5, length(subj_list));
hsurdata = cell(5, length(subj_list));
sbj_dt = struct();

savesurdata = false;
istheresfile = true;

ep_len = 120; % (s)
tmin = 10;
tmax = 1940;

for suj=1:length(subj_list)
    p_list = dataPath + subj_list(suj) + '\NFB\';
    dir_info = dir(char(p_list));
    sess_dir = {dir_info.name};
    sbj_dt(suj).name = subj_list(suj);
    sbj_dt(suj).path = p_list;
    sbj_dt(suj).folders = sess_dir(3:end);
    
    if ~isfile(char('datas/'+sbj_dt(suj).name+'/'+sbj_dt(suj).name+'.mat'))
        mkdir(char('datas/' +sbj_dt(suj).name))
        istheresfile = false;
        savesurdata = true;
    else
        load(char('datas/'+sbj_dt(suj).name+'/'+sbj_dt(suj).name+'.mat'))
        istheresfile = true;
    end

    for sess=1:length(sbj_dt(suj).folders)
        folder = sbj_dt(suj).folders(sess);
        splt = split(folder, '-');
        ffilename = p_list + folder{1}+'\'+subj_list(suj)+...
                    '-'+splt(2)+'_S'+num2str(sess)+'.edf';
        sbj_dt(suj).fullpath{sess} = ffilename;
        
        if ~isfile(ffilename)
            fullPathFName(sess) = ''; % sujeito SAJ sem primeira sessao
%             datas{sess, suj} = '';
        else
            fullPathFName(sess) = ffilename;
            disp(fullPathFName(sess));
            EEG = pop_biosig(char(fullPathFName(sess)), 'importevent','off');
%             [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', 'data-name', 'gui','off');
            
            % Fitragem passa-faixa.
            EEG = pop_eegfiltnew(EEG, 1,100,850,0,[],0);
            % tmin = 10; tmax = 1940;
            EEG = pop_select( EEG,'time',[tmin tmax] );
            
            sbj_dt(suj).tasks.bs1 = pop_select( EEG,'time',[0 60] );
            sbj_dt(suj).tasks.nfb = pop_select( EEG,'time',[60 1860] );
            sbj_dt(suj).tasks.bs2 = pop_select( EEG,'time',[1860 1920] );
            
            % Surrogates
            if ~istheresfile % || updatesurdata
                disp('>> Calculating Surrogates...')
                surdatas{sess, suj} = IAAFTsur(EEG.data(1,:, :), 1);
            else
                disp('>> Surrogates already loaded...')
            end
            
            delay = 1; % delay 1 between points in ordinal patterns (successive points)
            order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
            windowSize = 4*EEG.srate;
            %Computing Entropy
            disp('>> Calculating PE...')
%             hdatas{sess, suj}   = PE(EEG.data(1,:, :)', delay, order, windowSize);
            hdatas{sess, suj}   = PE(sbj_dt(suj).tasks.nfb.data(1,:, :)', delay, order, windowSize);
            disp('>> Calculating PE on surrogates...')
            hsurdata{sess, suj} = PE(surdatas{sess, suj}', delay, order, windowSize);
            disp('>> Calculating Wessertein...')
            wess_dist{sess, suj} = ws_distance(hdatas{sess,suj}, hsurdata{sess,suj}, 1);
        end
    end
    
    if savesurdata
        disp('>> Saving surDatas...')
        save(char('datas/'+sbj_dt(suj).name+'/'+sbj_dt(suj).name), 'surdatas');
    end
    
    show_subj_permH(EEG.times, subj_list(1), hdatas, hsurdata);
end

disp('COMPLETE!!')

% show_subj_permH(EEG.times, subj_list(1), hdatas, hsurdata);
            
%             % split datas in epochs of ep_len seconds
%             EEGepch = eeg_regepochs(EEG, 'recurrence', ep_len, 'limits', [0 ep_len],'eventtype', 'trial');

%             epoch_dt = cell(1800/ep_len, 1);
%             disp('Processing entropy');
%             for ep=1:1800/ep_len
% %                 disp(['Processing entropy. Epoch' num2str(ep)]);
%                 epoch_dt{ep, 1} = PE(EEGep.data(1,:, ep)', delay, order, windowSize);
%             end
% %             hdatas{sess, suj} = PE(EEG.data(1,:)', delay, order, windowSize);
%             hdatas{sess, suj} = mean(cell2mat(epoch_dt));




%% -----------------------------------------------------------------
% Plots 
for suj=1:length(subj_list)
%     figure(suj);
    suj_folder = ['fig/' char(subj_list(suj))];
    if ~exist(suj_folder)
        mkdir(suj_folder);
    end
    
    for i=1:5
        figure(i);
%         subplot(5,1,i);
        if subj_list(suj) == "SAJ" && i==1
            avgH = zeros(1, length(hdatas{1,1}));
        else
            avgH = hdatas{i,suj};
            
            plot(EEG.times(1:length(hdatas{i,suj})).*10^-3, hdatas{i,suj}, 'LineWidth',1.5); hold on;
            plot(EEG.times(1:length(hsurdata{i,suj})).*10^-3, hsurdata{i,suj}, 'LineWidth',1.5); % hold off;
            title(['Permutation Entropy ' '(S' num2str(i) ')']);
            legend('Data', 'Surrogate');
            
%             csvfile = strrep(sbj_dt(suj).fullpath{i}, '.edf', '.csv');
%             file = csvread(csvfile, 2, 0);
%             file = file';
%             EEGcsv = pop_importdata( 'setname', 'csvfile', 'data', file, 'nbchan', 13,...
%                                      'dataformat', 'array', 'srate', 256);
%             EEGcsv = pop_select( EEGcsv,'time',[tmin tmax] );
 
%             plot(EEG.times(1:length(avgH)).*10^-3, EEG.data(1,1:length(avgH)), 'LineWidth',1.5);  hold on;
%             plot(EEG.times(1:length(avgH)).*10^-3, EEGcsv.data(12,1:length(avgH))); hold on;
            
%             %Sombra 60 min resting inicio
%             patch([0 0 60 60], [0 1 1 0], [.5 .5 .5],'FaceAlpha',.2);
%             %Sombra 60 min resting final
%             patch([1860 1860 1920 1920], [0 1 1 0], [.5 .5 .5],'FaceAlpha',.2);
%             title(['S' num2str(i)])
            xlabel('Time (s)');
%             ylabel('Permutation Entropy');
            scrollplot('WindowSizeY', .5, 'WindowSize',30); hold off;
        end
        % Saving figure
%         filename = ['fig/' char(subj_list(suj)) '/' char(subj_list(suj)) '-S'... 
%                     num2str(i) '_avg_H-' num2str(ep_len) '.png'];
%         saveas(gcf, filename);
    end
end
            
            % Read from csv file
%             csvfile = strrep(sbj_dt(suj).fullpath{i}, '.edf', '.csv');
%             file = readtable(csvfile, 'HeaderLines',1);
%             rwd = file.REWARD;
%             plot(EEG.times(1:length(avgH))*10^-3, rwd(1:length(avgH))./20.4082 ); hold on;
            
%             patch([0 0 60 60], [0 1 1 0], [.5 .5 .5],'FaceAlpha',.5);
%             hold off;
%% -----------------------------------------------------------------

% disp()
% disp('Calculating Entropies...')
% for suj=1:length(subj_list)
%     for sess=1:length(sbj_dt(suj).folders)
%         datas(sess, suj) = PE(data, delay, order, windowSize);
%     end
% end


%% -----------------------------------------------------------------
% close all
% clc
% clear
% 
% % for fname = sess_dir
% %   fpath = '';
% % end
% 
% file = char('D:\Users\NFB\Pacientes\SAJ\NFB\nfb-221019\SAJ-221019_S2.edf');
% EEG = pop_biosig(file, 'importevent','off');
% 
% data = EEG.data;
% fs = EEG.srate; 
% 
% delay = 1; % delay 1 between points in ordinal patterns (successive points)
% order = 3; % order 3 of ordinal patterns (4-points ordinal patterns)
% % windowSize = .390625*fs; % 100 miliseconds sliding window - fs = 256
% windowSize = 2*fs;
% 
% perm_H = PE(data, delay, order, windowSize);
% % time = EEG.times(1:size(perm_H,2)-1)*10^-3;
% time = 1/fs:1/fs:length(EEG.times)/fs;
% len_time = length(time);
% 
% figure(1);
% plot(time, perm_H(1:length(time)));

%% -----------------------------------------------------------------