function show_subj_permH(times, suj, hdatas, hsurdata)


    suj_folder = ['fig/' char(suj)];
    if ~exist(suj_folder)
        mkdir(suj_folder);
    end
    
    for i=1:5
        figure(i);
%         subplot(5,1,i);
        if suj == "SAJ" && i==1
            avgH = zeros(1, length(hdatas{1,1}));
        else
            avgH = hdatas{i};
            
            plot( times(1:length(hdatas{i, 1})).*10^-3, hdatas{i, 1}, 'LineWidth',1.5); hold on;
            plot( times(1:length(hsurdata{i, 1})).*10^-3, hsurdata{i, 1}, 'LineWidth',1.5); % hold off;
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
            
            % Read from csv file
%             csvfile = strrep(sbj_dt(suj).fullpath{i}, '.edf', '.csv');
%             file = readtable(csvfile, 'HeaderLines',1);
%             rwd = file.REWARD;
%             plot(EEG.times(1:length(avgH))*10^-3, rwd(1:length(avgH))./20.4082 ); hold on;
            
%             patch([0 0 60 60], [0 1 1 0], [.5 .5 .5],'FaceAlpha',.5);
%             hold off;