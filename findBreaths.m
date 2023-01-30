% Automatic breath belt annotations

% Note that this toolbox is under development.

% This function accepts a breath belt vector, its sampling rate, and
% various optional parameters.

% breathBelt = breath belt signal (of dimension nSamp x 1, for now)
% breathFs = breath belt sampling rate

% Optional name-value paired arguments:

% 'WinSz' (for moving mean smoothing, default is 20 ms)
% 'MinDur' (minimum duration between onset and offset, default is 150 ms)
% 'MinHeight' (minimum vertical distance between onset and offset, default
% is 0.075 A.U.)

% Example usage:
% [beg,en] = findBreaths(breathBelt,1000,'MinDur',200);   

function [beg,en] = findBreaths(breathBelt,breathFs,varargin)

defaultWin = 20;
defaultDuration = 150;
defaultHeight = 0.075;

argIn = inputParser;
validParam = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(argIn,'breathBelt');
addRequired(argIn,'breathFs',validParam);
addParameter(argIn,'WinSz',defaultWin,validParam);
addParameter(argIn,'MinDur',defaultDuration,validParam);
addParameter(argIn,'MinHeight',defaultHeight,validParam);

parse(argIn,breathBelt,breathFs,varargin{:});

if breathFs ~= 1000 % Work with breath belt data in milliseconds
    FsOut = breathFs;
    [p,q] = rat(1000/breathFs);
    breathBelt = resample(breathBelt,p,q);
    breathFs = 1000;
else
    FsOut = breathFs;
end

win = argIn.Results.WinSz;
minHt = argIn.Results.MinHeight;
minDur = argIn.Results.MinDur;

timeChunk = round(breathFs/30); % Interval to check linearity of postive-going slope (i.e., inhalations) every ~33 ms
shortChunk = round(0.5*timeChunk); % More precise interval for fine-tuning
%minPeakDist = 0.3*breathFs; % Separation in ms between peaks
%splitPeakDist = 0.15*Fs; % Threshold at which to join rather than delete breaths

% Smooth
vector_smooth = movmean(breathBelt,breathFs*(win/1000));
vector_smooth = rescale(vector_smooth - mean(vector_smooth),0,1);

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
beg = [];
en = [];

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
                beg = [beg ; new_beg];
                en = [en ; new_en];
            end
            break
        end
    end
end

beg(en>=numel(vector_smooth)) = [];
en(en>=numel(vector_smooth)) = [];


% Join any remaining split breaths

% If two separated peaks are close together, choose the first one (based on
% comparison to acoustic signal, it seems that second/later peaks, even if larger, tend to be
% speech and not inhalation)

for ii = 1:numel(en)-1
    if ~(isnan(beg(ii)) & isnan(en(ii)))
        
        t0 = beg(ii);
        
        if isnan(t0)
            j = ii;
            while isnan(t0)
                t0 = beg(j);
                j = j - 1;
            end
        end
        
        t1 = en(ii);
        t2 = beg(ii+1);
        t3 = en(ii+1);
        
        if t2-t1 < 1.5*minDur
            if t2-t1 >= 1.25*timeChunk
                % Keep split
                ht = [vector_smooth(t1)-vector_smooth(t0) ...
                    vector_smooth(t3)-vector_smooth(t0)];
                
                if ht(1)/ht(2) > 0.5
                    % only keep first slope
                    beg(ii+1) = NaN;
                    en(ii+1) = NaN;
                else
                    % only keep second slope
                    [~,new_beg] = min(vector_smooth(t1:t2));
                    new_beg = new_beg + t1 - 1;
                    beg(ii+1) = new_beg;
                    
                    beg(ii) = NaN;
                    en(ii) = NaN;
                end
                
            else
                % Join
                beg(ii+1) = NaN;
                en(ii) = NaN;
            end
        end
    end
end

beg(isnan(beg)) = [];
en(isnan(en)) = [];

% Remove sub-threshold breaths
z_ht = vector_smooth(en)-vector_smooth(beg); % By height
beg(z_ht<minHt) = [];
en(z_ht<minHt) = [];

z_ht = vector_smooth(en)-vector_smooth(beg); % By grade
z_dur = en-beg;
z_grade = z_ht./z_dur;
lt = prctile(z_grade,25)-2*iqr(z_grade); %prctile(iqr(z_grade),25)-1.5*iqr(z_grade); 
beg(z_grade<lt) = []; 
en(z_grade<lt) = [];

% Convert back to original Fs
if FsOut ~= breathFs
    f = (FsOut/breathFs);
    beg = round(f*beg);
    en = round(f*en);    
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