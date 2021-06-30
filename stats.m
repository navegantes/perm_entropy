classdef stats
methods (Static)

function [p, h, S] = calc_stats(trialsdata, numtrials, tindx, test)
    
    p = zeros(numtrials, 1);
    h = zeros(numtrials, 1);
    S = cell(numtrials, 1);
    
    dnevet = tindx(1); % samp: -250
    tzero  = tindx(2); % samp: 0
    tevend = tindx(3); % samp: 250
    
    for trial=1:numtrials
        x = trialsdata(dnevet:tzero-1,trial);
        y = trialsdata(tzero+1:tevend,trial);
        if strcmp(test, 'ranksum')
            [p(trial), h(trial), S{trial}] = ranksum(x, y, 'tail', 'right');
        end
        if strcmp(test, 'ttest')
            [h(trial), p(trial), ~, S{trial}] = ttest(x, y,'tail','right');
        end
    end
end
% -----------------------------------------------------------------------------------
% numtrials = EEGev(sess).trials;
% pe        = permute(EEGev(sess).data(12,:,:), [2 3 1]);
% pesu      = permute(EEGev(sess).data(13,:,:), [2 3 1]);
% se        = permute(EEGev(sess).data(14,:,:), [2 3 1]);
% sesu      = permute(EEGev(sess).data(15,:,:), [2 3 1]);
% entrData  = {pe, pesu, se, sesu};
function [p, h, S] = cross_stats(xdata, ydata, tindx, test)

    dnevet = tindx(1); % samp: -250
    tzero  = tindx(2); % samp: 0
    tevend = tindx(3); % samp: 250
    
    numtrials = size(xdata, 2);
    
    p = zeros(numtrials, 2); %, zeros(numtrials, 1)};
    h = zeros(numtrials, 2); %, zeros(numtrials, 1)};
    S = {cell(numtrials, 1),cell(numtrials, 1)};  %, cell(numtrials, 1)};
    
    for trial=1:numtrials
        x = {xdata(dnevet:tzero, trial), xdata(tzero:tevend, trial)};
        y = {ydata(dnevet:tzero, trial), ydata(tzero:tevend, trial)};
        
        if strcmp(test, 'ranksum')
            for rwd=1:2
                [p(trial, rwd), h(trial, rwd), S{rwd}{trial}] = ranksum(x{rwd}, y{rwd}, 'tail', 'right');
            end
        elseif strcmp(test, 'ttest')
            for rwd=1:2
                [h(trial, rwd), p(trial, rwd), ~, S{rwd}{trial}] = ttest(x{rwd}, y{rwd},'tail','right');
            end
        else
            error('Error. \nInvalid stats method.');
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
