% Acoustic calibration for breath belt annotations

% Note that this toolbox is under development.

% This function requires the output breath beginning and ends generated by
% breathTimes, breath belt and its Fs, plus the corresponding acoustic recording (and Fs), 
% which is used to reject spurious inhalation events and adjust the inhalation
% onsets and ends as needed.

% The name-value paired argument envMethod determines which envelope extraction
% menthod to use. The default (1) method filters the signal into critical
% bands based on the bark scale and is slower but more accurate. The
% alternative (2) method is the Hilbert transform and is faster but
% untested and potentially less accurate.

% beg = vector of breath onsets
% en = vector of breath ends
% breathBelt = breath belt signal (nSamp x 1)
% breathFs = sample rate of breath belt/annotations
% audioData = acoustic data
% audioFs  = acoustic sample rate

% Optional name-value paired arguments:
% 'EnvMethod' (1 or 2, default is 1)
% 'MinDur' (minimum duration between onset and offset, default is 150 ms)
% 'MinHeight' (minimum vertical distance between onset and offset, default
% is 0.075 A.U.)

% Example usage:
% [beg,en] = breathSpeechCompare(beg,en,breathBelt,1000,speechRecording,44100,2,'MinHeight',0.10);                


function [beg,en] = breathSpeechCompare(beg,en,breathBelt,breathFs,audioData,audioFs,varargin)

argIn = inputParser;

validParam = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validEnv = @(x) x == 1 || x == 2;
defaultEnvelope = 1;
defaultDuration = 150;
defaultHeight = 0.075;

addRequired(argIn,'beg');
addRequired(argIn,'en');
addRequired(argIn,'breathBelt');
addRequired(argIn,'breathFs',validParam);
addRequired(argIn,'audioData');
addRequired(argIn,'audioFs');
addParameter(argIn,'EnvMethod',defaultEnvelope,validEnv);
addParameter(argIn,'MinDur',defaultDuration,validParam);
addParameter(argIn,'MinHeight',defaultHeight,validParam);

parse(argIn,beg,en,breathBelt,breathFs,audioData,audioFs,varargin{:});

if breathFs ~= 1000 % Work in milliseconds
    FsOut = breathFs;
    [p,q] = rat(1000/breathFs);
    breathBelt = resample(breathBelt,p,q);
    breathFs = 1000;
    beg = beg*(p/q);
    en = en*(p/q);
else
    FsOut = breathFs;
end

env = argIn.Results.EnvMethod;
minDur = argIn.Results.MinDur;
minHt = argIn.Results.MinHeight;

breathBelt = rescale(breathBelt - mean(breathBelt),0,1);

% Set up parameters
maxBlip = 45*breathFs;
Rp = 1;
Rs = 50;
Fn = audioFs/2;
Wp_sp = [2.1  995.0]/Fn;
Ws_sp = [1.8  997.0]/Fn;
[n,Ws_sp]  = cheb2ord(Wp_sp,Ws_sp,Rp,Rs);
[z,p,k] = cheby2(n,Rs,Ws_sp);
[sos_sp,g_sp] = zp2sos(z,p,k);

CV = @(alpha)-sqrt(2)*erfcinv(2*alpha);
alpha = 0.6;
zs = CV([(1-alpha)/2   1-(1-alpha)/2]);

% Preprocess sound
if size(audioData,2) > size(audioData,1)
    audioData = audioData';
end
if size(audioData,2) == 2
    audioData = audioData(:,1);
end

audioData = audioData - mean(audioData);
audioData = audioData./(max(abs(audioData)));
sigLim = mean(audioData) + zs*std(audioData);
audioData(audioData<sigLim(1)) = sigLim(1);

audioData(audioData>sigLim(2)) = sigLim(2);
audioData = filtfilt(sos_sp, g_sp, audioData);

% Extract envelope

if env == 1
    spOut = env1(audioData,audioFs,breathFs);
elseif env == 2
    spOut = env2(audioData,audioFs,breathFs);    
end

spOut(spOut>prctile(spOut,95)) = prctile(spOut,95);
spOut = rescale(spOut);
thresh = quantile(movmean(spOut,10),50);
thresh = thresh([24 30]);

for ii = 1:numel(beg)
    fl = 0;
    speechVector = spOut(beg(ii):en(ii));
    breathVector = breathBelt(beg(ii):en(ii));
    og_beg = beg(ii);
    og_en = en(ii);
    
    % Simplify speech vector
    spSil = movmean(speechVector,5);
    spSil(spSil <= thresh(1)) = 0;
    spSil(spSil > thresh(1) & spSil <= thresh(2)) = thresh(1);
    spSil(spSil > thresh(2)) = thresh(2);
    
    % Bulk of inhalation should take place during non-speech
    brDer = diff(breathVector);
    winWid = ceil(numel(breathVector)/8);
    if winWid*6 ~= numel(brDer)
        brDer = [brDer ; nan(winWid*8-numel(brDer),1)];
    end
    rBrDer = reshape(brDer, winWid, []);
    mBrDer = nanmean(rBrDer);
    [~,a] = maxk(mBrDer,4);
    spDer = diff(speechVector);
    winWid = ceil(numel(spDer)/8);
    if winWid*6 ~= numel(spDer)
        spDer = [spDer ; nan(winWid*8-numel(spDer),1)];
    end
    rSpDer = reshape(spDer, winWid, []);
    mSpDer = nanmean(rSpDer);
    [~,b] = maxk(mSpDer,4);
    if sum(speechVector>thresh(1))/numel(speechVector) > 0.5 && ...
            sum(ismember(a,b)) > 2
        fl = 1; % flag otherwise
    end
    
    % Decide if 'silence'
    isSilence = [];
    
    for jj = 1:numel(spSil)-1
        
        tn0 = jj-find(spSil(1:jj-1)==0,1,'last');
        if isempty(tn0)
            tn0 = inf;
        end
        tn2 = jj-find(spSil(1:jj-1)==thresh(2),1,'last');
        if isempty(tn2)
            tn2 = inf;
        end
        t0 = find(spSil(jj:end)==0,1,'first');
        if isempty(t0)
            t0 = inf;
        end
        t2 = find(spSil(jj:end)==thresh(2),1,'first');
        if isempty(t2)
            t2 = inf;
        end
        
        if spSil(jj) == 0 % 0, then silence
            isSilence(jj) = 1;
            
        elseif spSil(jj) == thresh(1)
            % if threshold(1) and coming from or leading to
            % threshold(2), then not silence
            
            if t2 < t0
                isSilence(jj) = 0;
            elseif tn2 < tn0
                isSilence(jj) = 0;
            else
                isSilence(jj) = 1;
            end
            
        elseif spSil(jj) == thresh(2) % threshold(2), then not silence
            
            isSilence(jj) = 0;
            
        end
    end
    
    % Permit ingressive sounds (e.g., pops and clicks)
    
    isSilence(1) = 1;
    isSilence(end) = 1;
    a = strfind(isSilence,[1 0]);
    b = strfind(isSilence,[0 1]);
    for jj = 1:numel(a)
        if a(jj) > 0.05*numel(isSilence) && ... % Skip edges
                numel(isSilence)-b(jj) > 0.05*numel(isSilence)
            
            if sum(spSil(a(jj)-5:b(jj)+5) > thresh(1)) <= maxBlip
                isSilence(a(jj):b(jj)) = 1;
            end
            
        end
    end
    
    if sum(isSilence) < 1.25*minDur && ...
            sum(isSilence)/numel(isSilence) < 0.85
        fl = 1; % Flag if extra short breath with non-silence
    end
    
    % Find silence edges
    isSilence(1) = 0;
    isSilence(end) = 0;
    a = strfind(isSilence,[0 1]);
    b = strfind(isSilence,[1 0]);
    
    [~,l] = max(b-a);
    
    new_beg = beg(ii) + a(l);
    new_en = beg(ii) + b(l);
    breathVector = breathBelt(new_beg:new_en);
    %speechVector = spOut(new_beg:new_en);
    
    % Final Flags
    if ~isempty(breathVector)
        
        % More that 75% of slope trimmed
        if (breathBelt(og_en)-breathBelt(og_beg))/(breathBelt(new_en)-breathBelt(new_beg)) < 0.75
            fl = 1;
        end
        
        % More than 66% of duration trimmed
        if (new_en-new_beg)/numel(isSilence) < 0.33
            fl = 1;
        end
        
    else
        fl = 1;
    end
    
    if fl == 0
        
        firstDeriv = movmean(diff(breathVector),5);
        
        if new_beg > og_beg+1
            % New inhale begin
            nM = find(firstDeriv(1:round(0.5*numel(breathVector)))<=...
                quantile(firstDeriv,.05),1,'last');
            
            if isempty(nM)
                [~,nM] = min(breathVector);
            end
            
            new_beg = new_beg + nM - 1;
            
            breathVector = breathBelt(new_beg:new_en);
        end
        
        
        if new_en < og_en-1
            % New inhale end
            id = round(0.8*numel(firstDeriv));
            nM = find(...
                firstDeriv(id:end)<=...
                quantile(firstDeriv(id:end),.2),1,'first');
            if isempty(nM)
                [~,nM] = max(breathVector(id:end));
            end
            
            new_en = new_beg + ...
                id + nM - 1;
        end
        
        
        if new_en-new_beg >= minDur
            
            beg(ii) = new_beg;
            en(ii) = new_en;
        else
            beg(ii) = NaN;
            en(ii) = NaN;
        end
        
    else
        beg(ii) = NaN;
        en(ii) = NaN;
    end
end

beg(isnan(beg)) = [];
en(isnan(en)) = [];

hts = breathBelt(en)-breathBelt(beg); % Remove breaths that are too small/short to be realistic

beg(hts < minHt) = [];
en(hts < minHt) = [];

durs = en-beg;
beg(durs < minDur) = [];
en(durs < minDur) = [];

% Convert back to original Fs
if FsOut ~= breathFs
    f = (FsOut/breathFs);
    beg = round(f*beg);
    en = round(f*en);    
end

end