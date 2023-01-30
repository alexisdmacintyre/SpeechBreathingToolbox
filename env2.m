function envelopeOut = env2(signalIn,FsIn,FsOut,lowPass)

% This method is simply to compute the analytic signal using the Hilbert
% transform of the broadband speech signal.

% This method is commonly used across speech sciences and is described, for example, in:
% He, L., & Dellwo, V. (2016, September). A Praat-based algorithm to
% extract the amplitude envelope and temporal fine structure using the
% Hilbert transform. ISCA.

if ~exist('lowPass','var')
    lowPass = 10;
end

% Transform
env = hilbert(signalIn);
env = abs(env);

% Low pass filter
[z,p,k] = butter(8,lowPass/(round(FsIn/2)),'low');
%[z,p,k] = butter(8,[1 lowPass]/(FsIn/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

env = filtfilt(sos,g,env);

% Resample
[p,q] = rat(FsOut/FsIn); 
env = resample(env,p,q);

envelopeOut = env;
end