%This script is used to compute the position, the weight and the time delay for each
%microphone in the array. If the “config” variable does not exist, or it is set to zero, the script
%becomes interactive and prompts the user.

v = 343; % sound speed in air [m/s], 343 m/s at 20degC, 331 m/s at 0degC
% reference equation: D/v = 1/(2*f) = NT*Ts
if ~exist('config','var') || config==0,
    BFtype = input('BeamForming type: 0=broadside(summing), 1=endfire(differential)? ');
    Nmic = input('Number of microphones (1-8, but for endfire: 2-3)? ');
    if BFtype>0, % endfire, D = v*NT*Ts, Ts=1/Fs
        if (Nmic<2) || (Nmic>3), fprintf('Wrong number of mics!\n'); return; end;
        Fs = input('Sampling frequency (typ. 48e3 for PCM, 1e6 for PDM)? '); Ts = 1/Fs;
        f1 = 3e3; N1 = ceil(1/(2*f1*Ts)); % null is at 180deg, null freq range 3-48kHz
        f2 = 48e3; N2 = max(floor(1/(2*f2*Ts)),ceil(0.004/v/Ts)); % mic distance >= 4mm
        for NT = round(linspace(max(N2,1),N1,min(8,N1-N2))),
            fprintf('Delay NT = %d (D = %.1fmm) for null at %.1fkHz\n', ...
                NT,v*NT*Ts*1000,1/(2*NT*Ts)/1e3);
        end;
        NT = input('Delay (N samples, set null freq, typ. 2x freq of interest) ? ');
        D = v*NT*Ts; f = v/(2*D);
    else % broadside
        Ts = 0; NT = 0; % default time delay NT*Ts=0
        if (Nmic<1) || (Nmic>8), fprintf('Wrong number of mics!\n'); return; end;
        f = input('Frequency (Hz, set null freq, typ. 1x freq of interest)? ');
        D = v/(2*f);
    end
    D = input(sprintf('Microphone spacing (%.1fmm for null at %.1fkHz)? ',D*1000,f/1e3));
    D = D/1000; f = v/(2*D); fprintf('Actual null freq: %.1fkHz\n',f/1e3);
    if BFtype>0, f = f/2; end; % endfire: show peak at quarter-wavelength
    config = 1;
else % config exists
    fprintf('BFtype=%d; Nmic=%d; D=%.3f; NT=%d; Ts=1/%.4g;\n',BFtype,Nmic,D,NT,1/Ts);
end;
% y ^
% |
% z .--> x
% set position, weight and time delay of each microphone
if BFtype>0, % endfire (differential),
    BFstr = 'endfire(diff)';
    if Nmic==2, % 1st order
        m = [+D/2, 0; -D/2, 0]; % XY position [m]
        mw = [+1; -1]; % weights (no need to normalize)
        mt = [ 0; NT*Ts]; % time delays [s]
    else %Nmic==3 % 2nd order
        m = [+D, 0; 0, 0; -D, 0]; % XY position [m]
        mw = [+1; -2; +1]; % weights (no need to normalize)
        mt = [ 0; NT*Ts; 2*NT*Ts]; % time delays [s]
    end
else % broadside (summing)
    BFstr = 'broadside(sum)';
    m = zeros(Nmic,2); % XY position [m]
    m(1:Nmic,1) = [0:Nmic-1]*D; % X multiple of D
    m(:,1) = m(:,1) - max(m(:,1))/2; % center around X=0
    mw = ones(Nmic,1); % weights (no need to normalize)
    mt = zeros(Nmic,1); % time delays [s]
end
m = [m(:,1),m(:,2),zeros(size(m,1),1)]; % add 3rd coord z=0
%mw = sum(mw); % no need to normalize weights
maxgain = sum(abs(mw)); % max gain may be used to normalize plots
% other simulation parameters
f = 3e3; % ref freq, 300-3400Hz voice band, 3-4kHz min hearing threshold
f1 = 0.1e3; f2 = 20e3; % lower and upper freq bounds [Hz]
BFstr = sprintf('%s, %dmics, D=%.1fmm, T=%d*%.2fus',BFstr,Nmic,D*1000,NT,Ts*1e6);