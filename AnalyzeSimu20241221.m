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
ggg = dir('./SKNA*.mat') ;


% these are indices to evaluate
NRMSE = [] ;
NRMSE10_500 = [] ;
NRMSE500_1000 = [] ;
NRMSE1000_2000 = [] ;
NRMSE2000_5000 = [] ;
AR = [] ;


for Jidx = 1:length(ggg)


    fprintf([ggg(Jidx).name '\n'])

    % all simulated data are saved in these mat files. Check if they are good
    % visually
    load(['SKNA' num2str(Jidx) '.mat'])
    Fs = 10000;  % Sampling frequency in Hz


    % %%%%%%%%%%%%%%%%%%%%%%
    % generate simulated ECG
    if 0
        [locsf0] = RpeakUltraLong(ECG(1:10:end), 1000) ;
        locsf0 = locsf0*10 ;
        locsf = locsf0 ;
        for jj = 1: length(locsf)
            qidx = max(1, locsf(jj)-500) : min(locsf(jj)+500, length(SKNAI)) ;
            [~, b] = max(SKNAI(qidx)) ;
            locsf(jj) = max(1, locsf(jj)-500) + b(1) - 1 ;
        end

        % determine SQI
        rail_mv = [-5.554198887532222,5.556912223578890];
        ECGdown = ECG(1:20:end) ;

        % get bad segment mask
        [rcmask,hfmask,lpmask,bsmask,finalmask] = ...
            UNSW_ArtifactMask(ECGdown'/1000, rail_mv, 60, Fs/20) ;

        % collect SQI
        [locsf1,RRlist,nRR,mRR,nSections] = UNSW_QRSDetector(ECGdown'/1000, Fs/20, [], 0);
        [locsf2] = RpeakUltraLong(ECGdown, 500) ;
        locsf1 = locsf1(:) ; locsf2 = locsf2(:) ;

        Qn = floor(length(ECGdown)/500/5) ;
        SQI = zeros(Qn, 1) ;
        MASK = [] ;
        for qq = 1: Qn
            SP = max(1, (qq-1)*500*5 - 500*5 + 1) ;
            EP = min(length(ECGdown), (qq-1)*500*5 + 500*5) ;
            idx1 = find(locsf1 >= SP & locsf1 < EP) ;
            idx2 = find(locsf2 >= SP & locsf2 < EP) ;
            [tmp, ~] = bsqi_matlab(locsf1(idx1), locsf2(idx2), 0.1, 500) ;
            SQI(qq) = tmp ;
            if tmp < 0.95 ; MASK = [MASK; idx2] ; end
        end

        MASK = unique(MASK) ;

        % QQ are good beats;
        QQ = setdiff(1:length(locsf2), MASK) ;

        locsf2 = locsf2 * 20 ;
        for jj = 1: length(locsf2)
            qidx = max(1, locsf2(jj)-500) : min(locsf2(jj)+500, length(SKNAI)) ;
            [~, b] = max(ECG(qidx)) ;
            locsf2(jj) = max(1, locsf2(jj)-500) + b(1) - 1 ;
        end


        locsf3 = locsf2(QQ) ;
        Qidx = find(diff(QQ) > 1) ;

        ECG2 = ECG ;
        for jj = length(Qidx): -1: 1
            tmp = locsf3(Qidx(jj)): locsf3(Qidx(jj)+1) ;
            ECG2(tmp) = [] ;
        end

        [locsf0] = RpeakUltraLong(ECG2(1:10:end), 1000) ;
        locsf0 = locsf0*10 ;
        locsf = locsf0 ;
        for jj = 1: length(locsf)
            qidx = max(1, locsf(jj)-500) : min(locsf(jj)+500, length(ECG2)) ;
            [~, b] = max(ECG2(qidx)) ;
            locsf(jj) = max(1, locsf(jj)-500) + b(1) - 1 ;
        end

    end
    % %%%%%%%%%%%%%%%%%%%%%

    % prepare the simulated signal
    %Wn = [0.5 125] / (256 / 2);
    %[b, a] = butter(4, Wn, 'bandpass');
    %Q = filter(b, a, s) ;
    %ECG = resample(Q, Fs, 256) ;

    % the semi-real simulated SKNA and ECG are of different length. Make them the
    % same length before summing them
    n = min(length(ECG), length(SKNA)) ;
    ECG = ECG(1:n) ; ECG = ECG(:) ;
    SKNA = SKNA(1:n) ; SKNA = SKNA(:) ;

    SKNAI = ECG + SKNA ;
    SKNAtrue = SKNA ;

    clear SKNA n

    x = SKNAI ;
    fprintf(['length = ' num2str(length(x)/Fs/60), 'min\n'])



    % use the traditional BPF (other algorithms need R peaks info)
    % Normalize the passband frequencies by the Nyquist frequency
    Wn = [500 1000] / (Fs / 2) ;
    [b, a] = butter(4, Wn, 'bandpass') ;
    SKNAIbpf = filter(b, a, x) ;



    % detect R peaks 
    [locsf0] = RpeakUltraLong(x(1:10:end), 1000) ;
    locsf0 = locsf0*10 ;
    % change +- of x by the sign of R peaks
    if median(x(locsf0))<0
        x= -x; SKNAtrue = -SKNAtrue;
    end
    locsf = locsf0 ;

    for jj = 1: length(locsf)
        qidx = max(1, locsf(jj)-500) : min(locsf(jj)+500, length(x)) ;
        [~, b] = max(x(qidx)) ;
        locsf(jj) = max(1, locsf(jj)-500) + b(1) - 1 ;
    end

    %plot(SKNAI) ; hold on ; plot(locsf, SKNAI(locsf), 'ro') ; pause



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply other algorithms

    Wn = [0.5 150] / ((Fs/20) / 2);
    [b, a] = butter(4, Wn, 'bandpass');
    ECG = resample(filter(b, a, SKNAI(1:20:end)), 10000, 500) ;
    ECG = ECG(1:length(SKNAI)) ;
    
    % Algorithms for comparison
    t1 = tic ;
    SKNAts_suzanna = FECGSYN_ts_extraction(locsf, x', 'TS-SUZANNA', 0, 20, 2, Fs) ;
    SKNAts_cerutti = FECGSYN_ts_extraction(locsf, x', 'TS-CERUTTI', 0, 20, 2, Fs) ;
    SKNAts_lp = FECGSYN_ts_extraction(locsf, x', 'TS-LP', 0, 20, 2, Fs) ;
    SKNAesn = FECGSYN_adaptfilt_extraction(x', ECG', 'ESN', 0);
    SKNAlms = FECGSYN_adaptfilt_extraction(x', ECG', 'LMS', 0);
    SKNArls = FECGSYN_adaptfilt_extraction(x', ECG', 'RLS', 0);
    SKNAekf = FECGSYN_kf_extraction(locsf, x', 0) ;
    toc(t1) ;

    
   
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




    EST = [SKNAos(:) SKNAts_suzanna(:) SKNAts_cerutti(:) SKNAts_lp(:) SKNAlms(1:length(SKNAos)),SKNArls(1:length(SKNAos)),SKNAesn(1:length(SKNAos)),SKNAekf(1:length(SKNAos))'];
    % EvaluatePerformance
    [Metrics] = EvaluatePerformance(SKNAtrue, EST) ;
    % EvaluatePerformance: band pass filters
    [Metrics10_500] = EvaluatePerformance_bp(SKNAtrue, EST,Fs,50,500);
    [Metrics500_1000] = EvaluatePerformance_bp(SKNAtrue, EST,Fs,500,1000);
    [Metrics1000_2000] = EvaluatePerformance_bp(SKNAtrue, EST,Fs,1000,2000);
    [Metrics2000_5000] = EvaluatePerformance_bp(SKNAtrue, EST,Fs,2000,4999);

    % save all NRMSE
    NRMSE = [NRMSE; Metrics.NRMSE] ;
    NRMSE10_500 = [NRMSE10_500;Metrics10_500.NRMSE] ;
    NRMSE500_1000 = [NRMSE500_1000;Metrics500_1000.NRMSE] ;
    NRMSE1000_2000 = [NRMSE1000_2000;Metrics1000_2000.NRMSE] ;
    NRMSE2000_5000 = [NRMSE2000_5000;Metrics2000_5000.NRMSE] ;
    % Evaluate AR
    [MetricsAR] = EvaluatePerformance_AR1(x,EST,locsf,Fs);
    AR{Jidx} = MetricsAR;
end
%Save AR
save('AR_simulated','AR');
