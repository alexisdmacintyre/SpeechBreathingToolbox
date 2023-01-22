% This function identifies inhalation events (inhalation beginnings and
% endings) from the breath belt signal. Requires MATLAB 2016 or later. 
%
% Usage: [onsets,offsets] =
% breathTimes(vector,Fs,'WinSz',20,'MinDur',100,'MinHeight',0.1,'Plot',1)
%
% Required arguments: vector (Nx1 breath signal); Fs (sample rate)
%
% Optional name pair arguments: 'WinSz' (for moving mean smoothing, default is 20 ms);
% 'MinDur' (minimum duration between onset and offset, default is 100 ms);
% 'MinHeight' (minimum vertical distance between onset and offset, default
% is 0.075 A.U.); 'Plot' (set to 1 to view results, default is 0)
%
% Suggestions: Some breath belt patterns look like inhalations but are
% associated with vocal exhalation, so make sure you corroborate this
% function's output with known speech onsets, etc.
%
% Alexis Deighton MacIntyre
% a.macintyre.17@ucl.ac.uk

function [onsets,offsets] = breathTimes(vector,Fs,varargin)
defaultWin = 20;
defaultDuration = 100;
defaultHeight = 0.075;
defaultPlot = 0;
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validPlot = @(x) (x == 0) || (x == 1);
addRequired(p,'vector');
addRequired(p,'Fs',validScalarPosNum);
addParameter(p,'WinSz',defaultWin,validScalarPosNum);
addParameter(p,'MinDur',defaultDuration,validScalarPosNum);
addParameter(p,'MinHeight',defaultHeight,validScalarPosNum);
addParameter(p,'Plot',defaultPlot,validPlot);
parse(p,vector,Fs,varargin{:});
win = p.Results.WinSz;
minDur = (p.Results.MinDur/1000)*Fs;
minHt = p.Results.MinHeight;
plotResults = p.Results.Plot;
timeChunk = round(Fs/30); % Interval to check linearity of postive-going slope (i.e., inhalations) every ~33 ms
shortChunk = round(0.5*timeChunk); % More precise interval for fine-tuning
minPeakDist = 0.3*Fs; % Separation in ms between peaks
%splitPeakDist = 0.15*Fs; % Threshold at which to join rather than delete breaths
% Smooth
vector_smooth = movmean(vector,Fs*(win/1000));
vector_smooth = rescale(vector_smooth,0,1);
% Initial round-up of peaks
meanSlope = [];
c = 1;
for iii = 1:round((numel(vector_smooth)/timeChunk))-1
    
    t1_loop = timeChunk*iii;
    
    if iii == round((numel(vector_smooth)/timeChunk))-1
        t2_loop = numel(vector_smooth);
    else
        t2_loop = t1_loop+timeChunk;
    end
    vOut = mean(diff(vector_smooth(t1_loop:t2_loop)))*1e+4;
    
    if vOut < -0.25
        vOut(2) = 0;
        vOut(3) = 0;
        c = 1;
    else % Is a continuous, positive-going slope
        vOut(2) = c;
        if c == 1
            vOut(3) = vOut(1);
        elseif c == 2
            vOut(3) = meanSlope(iii-1,1)+vOut(1);
        elseif c > 2
            vOut(3) = meanSlope(iii-1,3)+vOut(1);
        end
        c = c + 1;
    end
    meanSlope = [meanSlope ; vOut];
end
% Join split slopes
for j = 2:size(meanSlope,1)-1
    if meanSlope(j,2) == 0 & ...
            meanSlope(j-1,2) ~= 0 & ...
            find(meanSlope(j:end,2) ~= 0,1,'first') < 4
        
        meanSlope(j,2) = meanSlope(j-1,2) + 1;
        meanSlope(j,3) = meanSlope(j-1,3);
        jj = 1;
        while meanSlope(j+jj,2) ~= 0
            meanSlope(j+jj,2) = meanSlope(j+jj,2) + meanSlope(j,2);
            meanSlope(j+jj,3) = meanSlope(j+jj,3) + meanSlope(j-1,3);
            jj = jj + 1;
            if j+jj > size(meanSlope,1)
                break
            end
        end
    end
end
% Ignore isolated slopes < minimum duration
minSlope = floor(minDur/timeChunk);
for i = 1:size(meanSlope,1)-1
    if (meanSlope(i,2) < minSlope & meanSlope(i,2) > 0) && ...
            meanSlope(i+1,2) == 0
        if meanSlope(i,2) == 1
            meanSlope(i,:) = [0 0 0];
        else
            bg = find(meanSlope(1:i,2)==1,1,'last');
            meanSlope(bg:i,:) = zeros(i-bg+1,size(meanSlope,2));
        end
    end
end
idx = zeros(size(meanSlope,1),1);
idx(meanSlope(:,2)>0) = 1;
idx(1) = 0;
idx(end) = 0;
% Ignore successive slopes of a lesser height
idx_st = strfind(idx',[0 1]);
idx_en = strfind(idx',[1 0]);
idx = idx_en'-idx_st';
ht = vector_smooth(idx_en);
for i = 2:numel(idx)
    if idx(i)*timeChunk < minDur && ...
            ht(i) < ht(i-1)
        idx_en(i) = NaN;
    end
end
idx_en(isnan(idx_en)) = [];
% Determine the local peaks of positive slopes
peaks = [];
for ii = 1:numel(idx_en)
    if ii == 1
        t1 = 1;
        t2 = (timeChunk*idx_en(ii))+timeChunk;
    elseif ii == numel(idx_en)
        t1 = (timeChunk*idx_en(ii))-timeChunk;
        t2 = numel(vector_smooth);
    else
        t1 = (timeChunk*idx_en(ii))-timeChunk;
        if (timeChunk*idx_en(ii))+timeChunk > numel(vector_smooth)
            t2 = numel(vector_smooth);
        else
            t2 = (timeChunk*idx_en(ii))+timeChunk;
        end
    end
    
    [~,m] = max(vector_smooth(t1:t2));
    peaks = [peaks ; t1+m];
    
end
peaks(peaks>=numel(vector_smooth)) = numel(vector_smooth);
peaks(isnan(peaks)) = [];
% Within established slopes, determine more precise inhalation beginning/end points
onsets = [];
offsets = [];
for ii = 1:numel(peaks)
    
    if ii == 1
        t1 = find(vector_smooth(1:peaks(ii))<quantile(vector_smooth(1:peaks(ii)),0.05),1,'last');
        t2 = peaks(ii);
    else
        t1 = find(vector_smooth(peaks(ii-1):peaks(ii))<=quantile(vector_smooth(...
            peaks(ii-1):peaks(ii)),0.025),1,'last');
        t1 = t1 + peaks(ii-1) - 1;
        t2 = peaks(ii);
    end
    
    t1_new = t1;
    t2_new = t2;
    vec = vector_smooth(t1:t2);
    
    if numel(vec) >= minDur
        
        % Analyse slope and trim flat segments preceding true onset/after
        % true peak using more precise window
        meanSlope = [];
        c = 1;
        
        for iii = 1:floor((numel(vec)/shortChunk))
            
            t1_loop = shortChunk*iii;
            
            if iii == floor((numel(vec)/shortChunk))
                t2_loop = numel(vec);
            else
                t2_loop = t1_loop+shortChunk;
            end
            
            vOut = mean(diff(vec(t1_loop:t2_loop)))*1e+4;
            
            if vOut < 0.5
                vOut(2) = 0;
                vOut(3) = 0;
                c = 1;
            else
                vOut(2) = c;
                if c == 1
                    vOut(3) = vOut(1);
                elseif c == 2
                    vOut(3) = meanSlope(iii-1,1)+vOut(1);
                elseif c > 2
                    vOut(3) = meanSlope(iii-1,3)+vOut(1);
                end
                c = c + 1;
            end
            meanSlope = [meanSlope ; vOut];
        end
 
        % Join splits
        isPos = (meanSlope(:,2)>1)';
        isPos(1) = 1;
        isPos(end) = 1;
        a = strfind(isPos,[1 0]);
        b = strfind(isPos,[0 1]);
        c = b-a;
        c = c/size(meanSlope,1) < 0.06; % Join if split is <= 5% of total slope
        a = a(c);
        b = b(c);
        
        for j = 1:numel(a)
            if (a(j) > 0.15*numel(isPos)) & ... % Don't join if closeby edges
                    (numel(isPos)-b(j) > 0.15*numel(isPos))
                
                firstPos = find(meanSlope(1:a(j),2)==1,1,'last');
                lastPos = find(meanSlope(b(j):end,2)<1,1,'first');
                lastPos = b(j) + lastPos - 2;
                
                if isempty(firstPos)
                    firstPos = 1;
                end
                if isempty(lastPos)
                    lastPos = size(meanSlope,1);
                end
                
                meanSlope(firstPos:lastPos,2) = ...
                    1:(lastPos-firstPos)+1;
                meanSlope(firstPos:lastPos,3) = ...
                    cumsum(meanSlope(firstPos:lastPos,1));
            end
        end
        % Narrow down begin/end point windows
        while size(meanSlope,1) > 0
            maxTot = find(meanSlope(:,3)==max(meanSlope(:,3)),1,'last'); % End of greatest pos. slope
            pb = find(meanSlope(1:maxTot,2)==1,1,'last'); % Beginning of biggest slope
            
            [t1_new,t2_new,pb,maxTot,meanSlope,vec] = ...
                updateFunc(t1_new,pb,maxTot,meanSlope,vector_smooth,shortChunk);
            
            % Check absolute grade of slope at edges
            seg = round(0.2*size(meanSlope,1));
            if seg >= 1
                n = find(meanSlope(1:seg,1)<1,1,'last');
                if ~isempty(n)
                    pb = pb + n;
                    [t1_new,t2_new,pb,maxTot,meanSlope,vec] = ...
                        updateFunc(t1_new,pb,maxTot,meanSlope,vector_smooth,shortChunk);
                end
            end
            
            seg = round(0.7*size(meanSlope,1));
            if seg <= size(meanSlope,1) && seg > 0
                n = find(meanSlope(seg:end,1)<1,1,'first');
                if ~isempty(n)
                    if n == 1
                        maxTot = seg - 1;
                    else
                        maxTot = seg + n - 2;
                    end
                    [t1_new,t2_new,pb,maxTot,meanSlope,vec] = ...
                        updateFunc(t1_new,pb,maxTot,meanSlope,vector_smooth,shortChunk);
                end
            end
            
            % Check relative grade of slope at edges
            meanSlope(:,3) = cumsum(meanSlope(:,1));
            rt = meanSlope(:,3)/max(meanSlope(:,3));
            if ~isempty(find(rt<0.06,1,'last'))
                pb = find(rt<0.06,1,'last') + 1; % trim any preceding flat section
                
                [t1_new,t2_new,pb,maxTot,meanSlope,vec] = ...
                    updateFunc(t1_new,pb,maxTot,meanSlope,vector_smooth,shortChunk);
                meanSlope(:,3) = cumsum(meanSlope(:,1));
                rt = meanSlope(:,3)/max(meanSlope(:,3));
            end
            
            if ~isempty(find(rt>=0.95,1,'first'))
                if find(rt>=0.95,1,'first') < numel(rt) % trim any subsequent flat section
                    maxTot = find(rt>=0.95,1,'first');
                    [t1_new,t2_new,pb,maxTot,meanSlope,vec] = ...
                        updateFunc(t1_new,pb,maxTot,meanSlope,vector_smooth,shortChunk);
                end
            end
            
            % Specify inhale begin
            p = quantile(vec,0.01);
            p = find(vec<=p,1,'last');
            
            new_beg = t1_new + p - 1;
            
            % Specify inhale end
            vec = vector_smooth(new_beg:t2_new);
            
            [~,p] = max(vec);
            
            if p == numel(vec)
                m = 0;
                while p == numel(vec) && t2_new+m <= numel(vector_smooth)
                    vec = vector_smooth(new_beg:t2_new+m);
                    [~,p] = max(vec);
                    m = m + 1;
                end
            end
            
            new_en = new_beg + p - 1;
            
            if new_en - new_beg >= minDur & ...
                mean(meanSlope(:,1)) > 1
                onsets = [onsets ; new_beg];
                offsets = [offsets ; new_en];
            end
            break
        end
    end
end
onsets(offsets>=numel(vector_smooth)) = [];
offsets(offsets>=numel(vector_smooth)) = [];
% Join any remaining split breaths
% If two separated peaks are close together, choose the first one (based on
% comparison to acoustic signal, it seems that second/later peaks, even if larger, tend to be
% speech and not inhalation)
for ii = 1:numel(offsets)-1
    if ~(isnan(onsets(ii)) & isnan(offsets(ii)))
        
        t0 = onsets(ii);
        
        if isnan(t0)
            j = ii;
            while isnan(t0)
                t0 = onsets(j);
                j = j - 1;
            end
        end
        
        t1 = offsets(ii);
        t2 = onsets(ii+1);
        t3 = offsets(ii+1);
        
        if t2-t1 < 1.5*minDur
            if t2-t1 >= 1.25*timeChunk
                % Keep split
                ht = [vector_smooth(t1)-vector_smooth(t0) ...
                    vector_smooth(t3)-vector_smooth(t0)];
                
                if ht(1)/ht(2) > 0.5
                    % only keep first slope
                    onsets(ii+1) = NaN;
                    offsets(ii+1) = NaN;
                else
                    % only keep second slope
                    [~,new_beg] = min(vector_smooth(t1:t2));
                    new_beg = new_beg + t1 - 1;
                    onsets(ii+1) = new_beg;
                    
                    onsets(ii) = NaN;
                    offsets(ii) = NaN;
                end
                
            else
                % Join
                onsets(ii+1) = NaN;
                offsets(ii) = NaN;
            end
        end
    end
end
onsets(isnan(onsets)) = [];
offsets(isnan(offsets)) = [];
% Remove sub-threshold breaths
z_ht = vector_smooth(offsets)-vector_smooth(onsets); % By height
onsets(z_ht<minHt) = [];
offsets(z_ht<minHt) = [];
z_ht = vector_smooth(offsets)-vector_smooth(onsets); % By grade
z_dur = offsets-onsets;
z_grade = z_ht./z_dur;
lt = prctile(z_grade,25)-2*iqr(z_grade); %prctile(iqr(z_grade),25)-1.5*iqr(z_grade); 
onsets(z_grade<lt) = []; 
offsets(z_grade<lt) = [];
if plotResults == 1
    figure;
    plot(vector_smooth,'LineWidth',0.75,'Color',[0.3010 0.7450 0.9330]);
    hold on
    plot(onsets,vector_smooth(onsets),'gx','MarkerSize',8)
    plot(offsets,vector_smooth(offsets),'rx','MarkerSize',8)
    yticks(-1:0.5:1)
    
    xlabel('Time (seconds)');
    xticks(1:round(numel(vector)/20):numel(vector)+1);
    t = numel(vector)/Fs;
    xticklabels(0:round(t/19):t);
    axis tight
    grid on
    set(gcf, 'Position', [500, 500, 900, 350])
    legend({'Breath Signal','Initiation','Maxima'},'Location','northeastoutside')
    title('Inhalation Events')
    set(gcf,'color','w')
end
end
function [t1_new,t2_new,pb,maxTot,meanSlope,vec] = ...
    updateFunc(t1_new,pb,maxTot,meanSlope,vector_smooth,shortChunk)
t2_new = t1_new + (maxTot*shortChunk) + shortChunk;
if t2_new > numel(vector_smooth)
    t2_new = numel(vector_smooth);
end
if pb ~= 1
    t1_new = t1_new + (pb*shortChunk) - shortChunk;
end
vec = vector_smooth(t1_new:t2_new);
meanSlope = meanSlope(pb:maxTot,:);
pb = 1;
maxTot = size(meanSlope,1);
end