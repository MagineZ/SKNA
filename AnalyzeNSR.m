clear ; close all

% load all needed codes.
addpath('C:\Users\magin\OneDrive\Documents\Skin_electric_cardiacsignal_removal\code') ;
addpath('C:\Users\magin\OneDrive\Documents\Skin_electric_cardiacsignal_removal\code\UNSW') ;
addpath(genpath('C:\Users\magin\Downloads\Code Archive3\fecgsyn\fecgsyn-master\subfunctions\extraction-methods'))
addpath(genpath('C:\Users\magin\Downloads\ecg-kit-0.1.6\ecg-kit-0.1.6'))
scrsz = get(0,'ScreenSize') ;
% use SimulateSKNA to generate SKNA
% use high quality ECG segments with nonlocal median to generate high quality clean
% ECG
ggg = dir('./*.abf') ;


% these are indices to evaluate
AR = {};


for Jidx = 1:length(ggg)


    fprintf([ggg(Jidx).name '\n'])

    % all simulated data are saved in these mat files. 
    data = abfload(ggg(Jidx).name);
    SKNAR = data(:,2) - mean(data(:,2)) ;
    Fs = 10000;  % Sampling frequency in Hz

    trend = 0 ;
    SKNAI00 = -data(:,1) - trend ;
    Wn = 1 / (Fs / 2);
    [b, a] = butter(4, Wn, 'high');
    SKNAI0 = filter(b, a, SKNAI00) ;
    

    
    x = SKNAI0;
    fprintf(['length = ' num2str(length(x)/Fs/60), 'min\n'])


    % detect R peaks 
    [locsf0] = RpeakUltraLong(x(1:10:end), 1000) ;
    locsf0 = locsf0*10 ;
    % change +- of x by the sign of R peaks
    if median(x(locsf0))<0
        x= -x; SKNAI0 = -SKNAI0;
    end
    locsf = locsf0 ;

    for jj = 1: length(locsf)
        qidx = max(1, locsf(jj)-500) : min(locsf(jj)+500, length(x)) ;
        [~, b] = max(x(qidx)) ;
        locsf(jj) = max(1, locsf(jj)-500) + b(1) - 1 ;
    end
    %Remove unnormal beats
    locsf = unique(locsf) ;
    ind = ones(1,length(locsf));
    for i = 1:(length(locsf)-1)
        if locsf(i+1)-locsf(i)<100
            if x(locsf(i+1))>x(locsf(i+1))
                ind(i) = 0;
            else
                ind(i+1)=0;
            end
        end
    end
    locsf = locsf(logical(ind));
    % use the traditional BPF (other algorithms need R peaks info)
    % Normalize the passband frequencies by the Nyquist frequency
    Wn = [500 1000] / (Fs / 2) ;
    [b, a] = butter(4, Wn, 'bandpass') ;
    SKNAIbpf = filter(b, a, x) ;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply other algorithms

    Wn = [0.5 150] / ((Fs/20) / 2);
    [b, a] = butter(4, Wn, 'bandpass');
    ECG = resample(filter(b, a, x(1:20:end)), 10000, 500) ;
    ECG = ECG(1:length(x)) ;


    if 0
    t1 = tic ;
    SKNAts_suzanna = FECGSYN_ts_extraction(locsf, x', 'TS-SUZANNA', 0, 20, 2, Fs) ;
    SKNAts_cerutti = FECGSYN_ts_extraction(locsf, x', 'TS-CERUTTI', 0, 20, 2, Fs) ;
    SKNAts_lp = FECGSYN_ts_extraction(locsf, x', 'TS-LP', 0, 20, 2, Fs) ;
    SKNAesn = FECGSYN_adaptfilt_extraction(x', ECG', 'ESN', 0);
    SKNAlms = FECGSYN_adaptfilt_extraction(x', ECG', 'LMS', 0);
    SKNArls = FECGSYN_adaptfilt_extraction(x', ECG', 'RLS', 0);
    SKNAekf = FECGSYN_kf_extraction(locsf, x', 0) ;
    toc(t1) ;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start to apply our OS algorithm to remove ECG

    p1 = ceil((quantile(diff(locsf), .95))*1/2) ;
    p2 = ceil((quantile(diff(locsf), .95))*3/4) ;
    [~, ~, Y, YOSc, YOSc2] = OSdecomp(x, locsf, p1, p2, 0, 1) ;

    % stitch all cycles back to get cardiogenic artifact free cycles.
    [xqECG, CHECK] = stitchALL(YOSc2, locsf, length(x), p1) ;

    SKNAos0 = x - xqECG ;
    trend = smooth(medfilt1(SKNAos0, Fs/10), Fs/8, 'loess') ;

    SKNAos = SKNAos0 - trend ;
    xq = xqECG + trend ;


    keyboard


    EST = [SKNAos(:) SKNAts_suzanna(:) SKNAts_cerutti(:) SKNAts_lp(:) SKNAlms(1:length(SKNAos)),SKNArls(1:length(SKNAos)),SKNAesn(1:length(SKNAos)),SKNAekf(1:length(SKNAos))'];
    
    % Evaluate AR
    [Metrics] = EvaluatePerformance_AR1(x,EST,locsf,Fs);
    AR{Jidx} = Metrics;
    save(string(Jidx) +"_EST",'EST');
end

save('AR','AR');
