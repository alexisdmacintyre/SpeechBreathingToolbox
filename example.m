% Example script for speechBreathingToolbox

% Note that this toolbox is under development.

close all
clearvars
clc

localPath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(localPath))

% Load breath belt
[breathBelt,breathFs] = audioread(fullfile(localPath,'exampleBreathBelt.wav'));

% Detect inhalation onsets and ends
[beg,en] = findBreaths(breathBelt,breathFs);

% Load acoustic recording
[audioData,audioFs] = audioread(fullfile(localPath,'exampleAudio.wav'));

% Adjust inhalation onsets and ends based on speech envelope
[beg,en] = breathSpeechCompare(beg,en,breathBelt,breathFs,audioData,audioFs,'EnvMethod',2);

% Plot the results
plotBreaths(beg,en,breathBelt,breathFs,audioData,audioFs);

