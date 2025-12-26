clear ; close all ; clc
% load all needed codes.
addpath('/Users/hautiengwu/Dropbox/___working tex since 20221228/_WORKING SKNA/code') ;
addpath('/Users/hautiengwu/Dropbox/code/fecgsyn-master/subfunctions/extraction-methods') ;
addpath('/Users/hautiengwu/Dropbox/code/fecgsyn-master/subfunctions/extraction-methods/libs') ;
addpath('/Users/hautiengwu/Dropbox/code/fecgsyn-master/subfunctions/extraction-methods/libs/EKF') ;
addpath('/Users/hautiengwu/Dropbox/code/fecgsyn-master/subfunctions/extraction-methods/libs/ESNTOOL') ;
addpath('/Users/hautiengwu/Dropbox/code/fecgsyn-master/subfunctions/extraction-methods/libs/FECGESN') ;

scrsz = get(0,'ScreenSize') ;



%RESULTDIR = 'decomp20241010' ;

DIR = 'mice control' ;

ggg = dir([DIR '/*.edf']) ;


for Jidx = 1: length(ggg)

    fprintf(['Analyze ' DIR '/' ggg(Jidx).name '\n'])

    [profile, data] = edfread([DIR '/' ggg(Jidx).name]) ;
    data = data' ;

    CHidx = -1 ;
    for lll = 1: length(profile.label)
        if strcmp(profile.label{lll}, 'tungstenraw')
            CHidx = lll ;
        end
    end

    if CHidx == -1
        keyboard
    end

    SKNAI00 = data(:, CHidx) - mean(data(:, CHidx)) ;

    Fs = profile.samples(CHidx) ;  % Sampling frequency in Hz

    if Fs ~= 10000
        keyboard
    end

    % Remove abnormal segments in each case
    if Jidx == 1
        SKNAI00(1:100*Fs) = [];
    elseif Jidx == 2
        SKNAI00(1:455*Fs) = [];
    elseif Jidx == 4
        SKNAI00(1:21*Fs) = [];
    end
    Wn = 1 / (Fs / 2);
    [b, a] = butter(4, Wn, 'high');
    SKNAI0 = filter(b, a, SKNAI00) ;

    x = SKNAI0 ;
    % In mice database we need to remove the 60 HZ artifact from AC
    if Jidx == 4
        Wn = [50 70] / (Fs / 2);
        [b, a] = butter(4, Wn, 'stop');
        x = filter(b, a, SKNAI0);
    else
        % Design a notch filter
        F0 = 60;  % Frequency to remove (60 Hz)
        Q = 35;   % Quality factor (higher value = narrower notch)
        [b, a] = iirnotch(F0/(Fs/2), F0/(Fs*Q));
        x = filter(b, a, x);
    end
    fprintf(['length = ' num2str(length(x)/Fs/60), 'min\n'])


    Wn = [500 1000] / (Fs / 2);
    [b, a] = butter(4, Wn, 'bandpass');
    SKNAIbpf = filter(b, a, x) ;


    % detect R peaks (for mice, use 250 to adjust mice's fast heart rate)
    Fmice = Fs/10;
    locsf0 = RpeakUltraLong(x(1:10:end), Fmice/10) ;
    locsf0 = locsf0*10;
    if median(x(locsf0))<0
        x= -x; SKNAI0 = -SKNAI0;
    end
    locsf = locsf0 ;


    for jj = 1: length(locsf)
        qidx = max(1, locsf(jj)-350) : min(locsf(jj)+350, length(x)) ;
        [~, b] = max(x(qidx)) ;
        locsf(jj) = max(1, locsf(jj)-350) + b(1) - 1 ;
    end

    MANUALcorrect ;

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

    % use the traditional BPF
    % Normalize the passband frequencies by the Nyquist frequency
    Wn = [500 1000] / (Fs / 2);

    % Design the Butterworth bandpass filter
    % n is the filter order, you can start with 4 for a moderate filter
    [b, a] = butter(4, Wn, 'bandpass');
    SKNAI = filter(b, a, SKNAI0) ;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply other algorithms

    Wn = [0.5 150] / ((Fs/20) / 2);
    [b, a] = butter(4, Wn, 'bandpass');
    ECG = resample(filter(b, a, x(1:20:end)), 10000, 500) ;
    ECG = ECG(1:length(x)) ;


    t1 = tic ;
    SKNAts_suzanna = FECGSYN_ts_extraction(locsf, x', 'TS-SUZANNA', 0, 20, 2, Fs) ;
    SKNAts_cerutti = FECGSYN_ts_extraction(locsf, x', 'TS-CERUTTI', 0, 20, 2, Fs) ;
    SKNAts_lp = FECGSYN_ts_extraction(locsf, x', 'TS-LP', 0, 20, 2, Fs) ;
    SKNAesn = FECGSYN_adaptfilt_extraction(x', ECG', 'ESN', 0,Fs);
    SKNAlms = FECGSYN_adaptfilt_extraction(x', ECG', 'LMS', 0,Fs);
    SKNArls = SKNAlms ; %FECGSYN_adaptfilt_extraction(x', ECG', 'RLS', 0,Fs);
    SKNAekf = FECGSYN_kf_extraction(locsf, x', 0,Fs) ;
    toc(t1) ;



    % Step 2: run OS to remove ECG
    p1 = ceil((quantile(diff(locsf), .95))*1/2) ;
    p2 = ceil((quantile(diff(locsf), .95))*3/4) ;

    [~, ~, Y, YOSc, YOSc2] = OSdecomp(x, locsf, p1, p2, 0,1) ;



    % stitch all cycles back to get cardiogenic artifact free cycles.
    [xq0, CHECK] = stitchALL(YOSc2, locsf, length(x), p1) ;


    SKNAIos0 = x - xq0 ;
    trend = smooth(medfilt1(SKNAIos0, round(Fs/10)), round(Fs/8), 'loess') ;


    %if strcmp(ggg(Jidx).name, '07503322.abf')
    SKNAIos = SKNAIos0 - trend ;
    xq = xq0 + trend ;

    tmp = SKNAIos(1: floor(length(SKNAIos)/1000) * 1000) ;
    gg = reshape(abs(tmp), 1000, length(tmp)/1000) ;
    SKNAIosAM = mean(gg, 1);



    %save([DIR '-' ggg(Jidx).name(1:end-4)], 'SKNAIos', 'Fs')



    %EST = [SKNAIos(:) SKNAts_suzanna(:) SKNAts_cerutti(:) SKNAts_lp(:) SKNAlms(1:length(SKNAIos)),SKNArls(1:length(SKNAIos)),SKNAesn(1:length(SKNAIos)),SKNAekf(1:length(SKNAIos))'];


    % Evaluate AR for Mice
    %Fmice = 1500;
    %[Metrics] = EvaluatePerformance_ARMice(x,EST,locsf,Fmice);
    %AR{Jidx} = Metrics;



    %% plot results
    if Jidx == 1

        figure('Position',[1 scrsz(4) scrsz(3)*2/3 scrsz(4)]) ;

        MMM = quantile(xq, .998) * 1.2 ;
        mmm = quantile(xq, .002) * 1.2 ;
        MMMos = quantile(SKNAIos, .99) * 1.4 ;
        mmmos = quantile(SKNAIos, .01) * 1.4 ;
        MMMbpf = quantile(SKNAIbpf, .995) * 1.2 ;
        mmmbpf = quantile(SKNAIbpf, .005) * 1.2 ;

        ax(1) = subplot_tight(9,1,1,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, x, 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('SCNA-i') ; ylabel('(a)')

        ax(2) = subplot_tight(9,1,2,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, xq, 'r', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('ECG') ; ylabel('(b)')

        ax(3) = subplot_tight(9,1,3,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAIos, 'b', 'linewidth', 2) ;
        axis([-inf inf mmmos MMMos]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-S3') ; ylabel('(c)')

        ax(4) = subplot_tight(9,1,4,[0.008 0.062]) ;
        plot([1:length(SKNAIbpf)]/Fs, SKNAIbpf, 'b', 'linewidth', 2);
        axis([-inf inf mmmbpf MMMbpf]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-BPF') ; ylabel('(d)')

        MMM = quantile(SKNAts_lp, .999) * 1.2 ;
        mmm = quantile(SKNAts_lp, .001) * 1.2 ;
        ax(5) = subplot_tight(9,1,5,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAts_lp(:), 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-TSlp') ; ylabel('(e)')

        MMM = quantile(SKNAlms(1:length(SKNAIos)), .999) * 1.2 ;
        mmm = quantile(SKNAlms(1:length(SKNAIos)), .001) * 1.2 ;
        ax(6) = subplot_tight(9,1,6,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAlms(1:length(SKNAIos)), 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-LMS') ; ylabel('(f)')

        MMM = quantile(SKNAesn(1:length(SKNAIos)), .999) * 1.2 ;
        mmm = quantile(SKNAesn(1:length(SKNAIos)), .001) * 1.2 ;
        ax(7) = subplot_tight(9,1,7,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAesn(1:length(SKNAIos)), 'k', 'linewidth', 2) ;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-ESN') ; ylabel('(g)')

        MMM = quantile(SKNAekf(1:length(SKNAIos)), .999) * 1.2 ;
        mmm = quantile(SKNAekf(1:length(SKNAIos)), .001) * 1.2 ;
        ax(8) = subplot_tight(9,1,8,[0.008 0.062]) ;
        plot([1:length(SKNAIbpf)]/Fs, SKNAekf(1:length(SKNAIos))', 'k', 'linewidth', 2);
        axis([-inf inf mmm MMM]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-EKF') ; ylabel('(h)')




        xlabel('Time (s)') ;
        linkaxes(ax, 'x') ;

        xlim([215 219])
        export_fig(['MiceComparison'], '-transparent', '-pdf')

    end

keyboard

end

%save('AR_control','AR');
