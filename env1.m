function envelopeOut = env1(signalIn,Fs,FsOut)

if exist('FsOut','var')
    if isempty(FsOut)
        FsOut = 1000;
    end
elseif ~exist('FsOut','var')
    FsOut = 1000;
end

% Oganian, Y., & Chang, E. F. (2019). A speech
% envelope landmark for syllable encoding in human superior temporal gyrus.
% Science advances, 5(11).

% This identifies vowel nuclei according to the method of Schotola (1984).
% Written by Eric Edwards, adapted by Yulia Oganian.

signalIn = signalIn(:);
N = length(signalIn);
tN = (0:N-1)'./Fs;
T = 1/Fs;

%Loudness functions will be sampled at X Hz
sr = FsOut;
N1 = fix(N*sr/Fs);
t1 = (0:N1-1)'./sr;

%Critical-band filters are applied in freq-domain
signalIn = fft(signalIn,N);
frqs = (0:N/2)*(Fs/N); %FFT positive freqs
nfrqs = length(frqs);

%Set up for cricial-band filter bank (in freq-domain)
z = 13*atan(.76*frqs/1000) + 3.5*atan(frqs/7500).^2; % Bark (Zwicker & Terhardt 1980)
z = z(2:nfrqs-1); % to set DC and Nyquist to 0

%Set up for simple (RC) smoothing with 1.3-ms time-constant
tau = 0.0013; r = exp(-T/tau);
b1 = 1-r; a1 = poly(r);

F = zeros([N 1],'double'); % will hold critical-band filter shape
czs = 1:22; nczs = length(czs);
Nv = zeros ([N1 nczs],'double'); % will hold the specific loudnesses
for ncz=1:nczs, cz = czs(ncz);
    F(2:nfrqs-1) = 10.^(.7-.75*((z-cz)-.215)-1.75*(0.196+((z-cz)-.215).^2));
    % F = F./sum(F);
    Lev = real(ifft(signalIn.*F,N));
    Lev = Lev.^2; % square-law rectification (Vogel 1975)
    Lev = filter(b1,a1,Lev); % smoothing (1.3-ms time-constant)
    Lev = flipud(filter(b1,a1,flipud(Lev)));
    % the last line makes zero-phase smoothing, comment out to leave causal
    Lev = log(Lev); %logarithmic nonlinearity; this is now "excitation level"
    %The "specific loudness", Nv(t), is made by "delog" of .5*Lev(t)
    Nv(:,ncz) = interp1q(tN,exp(.5*Lev),t1);

end

Nv = max(0,Nv); %in case negative values, set to 0

% Nm is the total modified loudness, Nm(t). This is LPFed by a 3-point
% triangular filter, n=26 times (I do 13 fwd and 13 backwd for zero-phase),
% which results in ~Gaussian smoothing.
gv = ones([nczs 1],'double'); %weights
gv(czs<3) = 0; gv(czs>19) = -1;
Nm = Nv*gv;
b = ones([3 1],'double')./3; n = 13;
for nn = 1:n, Nm = filtfilt(b,1,Nm); end

envelopeOut = Nm;

%REFERENCES:

% Bismarck Gv (1973). Vorschlag fur ein einfaches Verfahren zur
% Klassifikation station?rer Sprachschalle. Acustica 28(3): 186-7.
%
% Zwicker E, Terhardt E, Paulus E (1979). Automatic speech recognition using
% psychoacoustic models. J Acoust Soc Am 65(2): 487-98.
%
% Ruske G, Schotola T (1978). An approach to speech recognition using
% syllabic decision units. ICASSP: IEEE. 3: 722-5.
%
% Ruske G, Schotola T (1981). The efficiency of demisyllable segmentation in
% the recognition of spoken words. ICASSP: IEEE. 6: 971-4
%
% Schotola T (1984). On the use of demisyllables in automatic word
% recognition. Speech Commun 3(1): 63-87.

end
