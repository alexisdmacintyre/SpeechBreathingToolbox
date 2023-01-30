function plotBreaths(beg,en,breathBelt,breathFs,varargin)

% Plot breath belt with annotations

% Note that this toolbox is under development.

% This function requires the output breath beginning and ends generated by
% breathTimes, and the breath belt and its Fs.

% You can also optionally plot the corresponding acoustic recording.

% beg = vector of breath onsets
% en = vector of breath ends
% breathBelt = breath belt signal (nSamp x 1)
% breathFs = sample rate of breath belt/annotations

% Optional additional arguments
% audioData = acoustic data
% audioFs  = acoustic sample rate

% Example usage:
% plotBreaths(beg,en,breathBelt,1000,speechRecording,44100); 


argIn = inputParser;

validParam = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addRequired(argIn,'beg');
addRequired(argIn,'en');
addRequired(argIn,'breathBelt');
addRequired(argIn,'breathFs',validParam);
addOptional(argIn,'audioData',[]);
addOptional(argIn,'audioFs',[]);

parse(argIn,beg,en,breathBelt,breathFs,varargin{:});

audioData = argIn.Results.audioData;
audioFs = argIn.Results.audioFs;

breathBelt = detrend(breathBelt);
breathBelt = rescale(movmean(breathBelt,50),-1,1);

figure('Units', 'normalized', 'Position', [0.1,0.25,0.85, 0.5]);

if ~isempty(audioData)

    if audioFs ~= breathFs
        [p,q] = rat(breathFs/audioFs);
        audioData = resample(audioData,p,q);
    end

    minA = median(audioData)-(3*mad(audioData));
    maxA = median(audioData)+(3*mad(audioData));

    audioData(audioData>maxA) = maxA;
    audioData(audioData<minA) = minA;

    audioData = rescale(audioData,0.5*min(breathBelt),0.5*max(breathBelt));
    
    plot(audioData,'LineWidth',1,'Color',[.88 .925 .9]);
    
    lString = {'Speech','Breath Belt','Breath Onset','End'};
    
else
    lString = {'Breath Belt','Breath Onset','End'};
    
end

hold on
plot(breathBelt,'LineWidth',2,'Color',[0.01 0.01 0.01])
plot(beg,breathBelt(beg),'gx','MarkerSize',12)
plot(en,breathBelt(en),'rx','MarkerSize',12)


yticks(-1:0.25:1)

xlabel('Time (seconds)','FontSize',15);
ylabel('Arbitrary Units','FontSize',15);
xticks(1:round(numel(breathBelt)/20):numel(breathBelt)+1);
t = numel(breathBelt)/breathFs;
xticklabels(0:round(t/19):t);

set(gcf,'color','w');

ylim([-1.1 1.1]);
set(gca,'box','off')
grid on
legend(lString,...
    'FontSize',15,'NumColumns',2,'EdgeColor',[1 1 1],...
    'Location','northeastoutside');

title('Inhalation Onsets and Ends','FontWeight','normal','FontSize',20);

zoom on

end