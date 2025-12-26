clear ; close all

scrsz = get(0,'ScreenSize') ;


%07503322.abf: first file & ignore it
%03213021.abf: pace maker
%0011.abf: Af subject

ggg = dir('./PMandAF/*.abf') ;

% 2=PM, 1=AF
for Jidx = 2 %1:length(ggg)

    fprintf(['./PMandAF/' ggg(Jidx).name '\n'])

    PACEMAKER = 0 ;

    if strcmp(ggg(Jidx).name, '03213021.abf')
        PACEMAKER = 1 ;
        fprintf('Analyze pacemaker\n')
    end

    data = abfload(['./PMandAF/' ggg(Jidx).name]);
    SKNAR = data(:,2) - mean(data(:,2)) ;

    Fs = 10000;  % Sampling frequency in Hz

    % remove the baseline wandering
    %trend = smooth(-data(:,1), Fs*10, 'loess') ;
    trend = 0 ;
    SKNAI00 = -data(:,1) - trend ;

    Wn = 1 / (Fs / 2);
    [b, a] = butter(4, Wn, 'high');
    SKNAI0 = filter(b, a, SKNAI00) ;



    % for the pacemaker signal, the first part contains strong harmonics.
    % Simply remove this stupid segment
    if PACEMAKER
        SKNAI = SKNAI0(9.04e5+1:end) ;
    else
        SKNAI = SKNAI0 ;
        SKNAI(162*Fs+1: 223*Fs) = [] ;
    end



    % use the traditional BPF (other algorithms need R peaks info)
    % Normalize the passband frequencies by the Nyquist frequency
    Wn = [500 1000] / (Fs / 2);
    [b, a] = butter(4, Wn, 'bandpass');
    SKNAIbpf = filter(b, a, SKNAI) ;



    x = SKNAI ;

    % Remove the 60 Hz AC artifact
    low_cutoff = 50; % Lower bound of the bandstop filter
    high_cutoff = 70; % Upper bound of the bandstop filter
    [b, a] = butter(4, [low_cutoff, high_cutoff]/(Fs/2), 'stop');
    x = filter(b, a, x);
    % use the traditional BPF (other algorithms need R peaks info)
    % Normalize the passband frequencies by the Nyquist frequency
    Wn = [500 1000] / (Fs / 2);
    [b, a] = butter(4, Wn, 'bandpass');
    SKNAI = filter(b, a, SKNAI) ;









    if strcmp(ggg(Jidx).name, '0011.abf')
        for jj = 20:-1:1
            [b, a] = iirnotch(60*jj/(Fs/2), 60*jj/(Fs*30));
            SKNAI = filter(b, a, SKNAI);
            SKNAI = filter(b, a, SKNAI);
        end
    end




    % detect R peaks
    [locsf0] = RpeakUltraLong(x(1:10:end), 1000) ;
    locsf0 = locsf0*10 ;
    if median(x(locsf0))<0
        x= -x; SKNAI0 = -SKNAI0;
    end
    locsf = locsf0 ;

    for jj = 1: length(locsf)
        qidx = max(1, locsf(jj)-500) : min(locsf(jj)+500, length(x)) ;
        [~, b] = max(x(qidx)) ;
        locsf(jj) = max(1, locsf(jj)-500) + b(1) - 1 ;
    end

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



    dlocsf = diff(locsf) ; dlocsf = [dlocsf(1); dlocsf] ;

    

    % stupid manual correction
    if strcmp(ggg(Jidx).name, '03213021.abf')
        qidx = find(locsf>1.1332e6 & locsf<1.13324e6) ;
        locsf(qidx) = [] ;
        qidx = find(locsf>1.14226e6 & locsf<1.14227e6) ;
        locsf(qidx) = [] ;
        qidx = find(locsf>1.14766e6 & locsf<1.14768e6) ;
        locsf(qidx) = [] ;
    elseif strcmp(ggg(Jidx).name, '0011.abf')
        qidx = find(locsf>4.849e5 & locsf<4.851e5) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.44e6 & locsf<1.441e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.4575e6 & locsf<1.458e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.47e6 & locsf<1.484e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.54e6 & locsf<1.541e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.63e6 & locsf<1.631e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.64e6 & locsf<1.6405e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>2.432e6 & locsf<2.433e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.0832e7 & locsf<1.0834e7) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.748e6 & locsf<1.785e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.644e6 & locsf<1.645e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.657e6 & locsf<1.658e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.823e6 & locsf<1.824e6) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
        qidx = find(locsf>1.0221e7 & locsf<1.0222e7) ;
        locsf(qidx) = [] ; dlocsf(qidx) = [] ;
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply other algorithms

    Wn = [0.5 150] / ((Fs/20) / 2);
    [b, a] = butter(4, Wn, 'bandpass');
    ECG = resample(filter(b, a, x(1:20:end)), 10000, 500) ;
    ECG = ECG(1:length(x)) ;

    fprintf('Start to run EKF\n')
    t1 = tic ;
    SKNAekf = FECGSYN_kf_extraction(locsf, x', 0) ;
    toc(t1) ;
    fprintf('Start to run TS-suzanna\n')
    t1 = tic ;
    SKNAts_suzanna = FECGSYN_ts_extraction(locsf, x', 'TS-SUZANNA', 0, 20, 2, Fs) ;
    toc(t1) ;
    fprintf('Start to run TS-cerutti\n')
    t1 = tic ;
    SKNAts_cerutti = FECGSYN_ts_extraction(locsf, x', 'TS-CERUTTI', 0, 20, 2, Fs) ;
    toc(t1) ;
    fprintf('Start to run TS-lp\n')
    t1 = tic ;
    SKNAts_lp = FECGSYN_ts_extraction(locsf, x', 'TS-LP', 0, 20, 2, Fs) ;
    toc(t1) ;
    fprintf('Start to run ESN\n')
    t1 = tic ;
    SKNAesn = FECGSYN_adaptfilt_extraction(x', ECG', 'ESN', 0);
    toc(t1)
    fprintf('Start to run LMS\n')
    t1 = tic ;
    SKNAlms = FECGSYN_adaptfilt_extraction(x', ECG', 'LMS', 0);
    toc(t1)
    fprintf('Start to run RLS\n')
    t1 = tic ;
    SKNArls = SKNAlms ; %FECGSYN_adaptfilt_extraction(x', ECG', 'RLS', 0);
    toc(t1)





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% start to apply our algorithm


    fprintf('Start to run S3\n')
    % Step 1: run OS to remove pacemaker (We will have the knowledge of pacemaker, so
    % we can assume we know if a subject has a pacemaker, and on which mode, etc)
    if PACEMAKER

        xq = zeros(size(x)) ;
        diffx = diff(x) ;

        % handle second peak -1
        locsf1 = locsf ;
        for jj = 1: length(locsf1)
            qidx = max(1, locsf(jj)-1500) : locsf(jj) ;
            [~, b] = max(diffx(qidx)) ;
            locsf1(jj) = max(1, locsf(jj)-1500) + b(1) - 1 ;

            qidx = max(1, locsf1(jj)-5) : locsf1(jj)+5 ;
            [~, b] = max(x(qidx)) ;
            locsf1(jj) = max(1, locsf1(jj)-5) + b(1) - 1 ;
        end

        p1 = 10 ;
        p2 = 10 ;
        [xc0, xq0, ~, ~, ~] = OSdecomp(x, locsf1, p1, p2, 20, 1) ;

        % handle second peak -2
        locsf2 = locsf1 ;
        for jj = 1: length(locsf1)
            qidx = locsf1(jj)+100 : locsf1(jj) + 160 ;
            [~, b] = max(diffx(qidx)) ;
            locsf2(jj) = locsf1(jj)+100 + b(1) - 1 ;

            qidx = max(1, locsf2(jj)-5) : locsf2(jj)+5 ;
            [~, b] = min(x(qidx)) ;
            locsf2(jj) = max(1, locsf2(jj)-5) + b(1) - 1 ;
        end

        [xc1, xq1, ~, ~, ~] = OSdecomp(xc0, locsf2, p1, p2, 20, 1) ;


        % first peaks -1
        locsf3 = locsf1 ;
        for jj = 1: length(locsf1)
            qidx = max(1, locsf1(jj)-2500) : locsf1(jj)-10 ;
            [~, b] = max(diffx(qidx)) ;
            locsf3(jj) = max(1, locsf1(jj)-2500) + b(1) - 1 ;

            qidx = max(1, locsf3(jj)-15) : locsf3(jj)+15 ;
            [~, b] = max(x(qidx)) ;
            locsf3(jj) = max(1, locsf3(jj)-15) + b(1) - 1 ;
        end

        % for the pacemaker case,
        % 4.608e6 switch from DDD to VVI
        % 8.297e6 switch from VVI to DDD
        %locsf3(locsf3>4.608e6 & locsf3<8.297e6) = [] ;
        p1 = 10 ;
        p2 = 20 ;
        [xc2, xq2, ~, ~, ~] = OSdecomp(xc1, locsf3, p1, p2, 20, 1) ;



        % first peaks -2
        locsf4 = locsf3 ;
        for jj = 1: length(locsf3)
            qidx = locsf3(jj)+50 : locsf3(jj)+200 ;
            [~, b] = max(abs(diffx(qidx))) ;
            locsf4(jj) = locsf3(jj)+50 + b(1) - 1 ;

            qidx = max(1, locsf4(jj)-10) : locsf4(jj)+10 ;
            [~, b] = min(x(qidx)) ;
            locsf4(jj) = max(1, locsf4(jj)-10) + b(1) - 1 ;
        end

        % 4.608e6 switch from DDD to VVI
        % 8.297e6 switch from VVI to DDD
        %locsf(locsf>4.608e6 & locsf<8.297e6) = [] ;
        p1 = 20 ;
        p2 = 20 ;
        [xc3, xq3, ~, ~, ~] = OSdecomp(xc2, locsf4, p1, p2, 20, 1) ;
        xqPM = xq0 + xq1 + xq2 + xq3 ;

    else

        xc3 = x ;
        xqPM = zeros(size(x)) ;

    end



    % Step 2: run OS to remove ECG

    dlocsf = diff(locsf) ; dlocsf = [dlocsf(1); dlocsf] ;
    % collect heart cycles
    
    if ~strcmp(ggg(Jidx).name, '0011.abf') 

        % this is for PM case
        p1 = ceil((quantile(diff(locsf), .95))*1/2) ;
        p2 = ceil((quantile(diff(locsf), .95))*3/4) ;
        % test is taking the manifold denoising is better. This part hasn't been fine
        % tuned since it is not the focus. But you can play with it.
        % run a simple eOptShrink
        [~, ~, Y, YOSc, YOSc2] = OSdecomp(xc3, locsf, p1, p2, 0, 1) ;

        % run ROSDOS
        %[~, ~, ~, ~, YlOSc2] = OSdecomp(xc3, locsf, p1, p2, 0, 0) ;
        YlOSc2 = YOSc2 ;
        
    else
    
        % this is AF case
        % for the Af subject, need to take care of the heart rate. Not an ideal
        % solution for sure!!!!!! We may need to improve it, maybe by some rescaling
        % technique.
        % this is the easy version, with simply kmeans to split
        % the data into 2 groups.
        [cidx, ctrs] = kmeans(dlocsf, 2) ;
        Qidx1 = find(cidx==1) ;
        p11 = ceil(quantile(dlocsf(Qidx1), .95) * 1/2) ;
        p12 = ceil(quantile(dlocsf(Qidx1), .95) * 3/4) ;
        [~, ~, Y1, YOSc1, YOSc21] = OSdecomp(xc3, locsf(Qidx1), p11, p12, 500, 1, 20) ;

        % this is ROSDOS (not yet tuned)
        %[~, ~, ~, ~, YlOSc21] = OSdecomp(xc3, locsf(Qidx1), p11, p12, 500, 0, 20) ;
        Qidx2 = find(cidx==2) ;
        p21 = ceil(quantile(dlocsf(Qidx2), .95) * 1/2) ;
        p22 = ceil(quantile(dlocsf(Qidx2), .95) * 3/4) ;
        [~, ~, Y2, YOSc2, YOSc22] = OSdecomp(xc3, locsf(Qidx2), p21, p22, 500, 1, 20) ;
        %[~, ~, ~, ~, YlOSc22] = OSdecomp(xc3, locsf(Qidx2), p21, p22, 500, 0, 20) ;
        p1 = max(p11, p21) ; p2 = max(p12, p22) ;
        Y = zeros(p1+p2+1, length(dlocsf)) ;
        Y(p1-p11+1: p1+p12+1, Qidx1) = Y1 ;
        Y(p1-p21+1: p1+p22+1, Qidx2) = Y2 ;
        YOSc = zeros(p1+p2+1, length(dlocsf)) ;
        YOSc(p1-p11+1: p1+p12+1, Qidx1) = YOSc1 ;
        YOSc(p1-p21+1: p1+p22+1, Qidx2) = YOSc2 ;
        YOSc2 = zeros(p1+p2+1, length(dlocsf)) ;
        YOSc2(p1-p11+1: p1+p12+1, Qidx1) = YOSc21 ;
        YOSc2(p1-p21+1: p1+p22+1, Qidx2) = YOSc22 ;

        %YlOSc2 = zeros(p1+p2+1, length(dlocsf)) ;
        %YlOSc2(p1-p11+1: p1+p12+1, Qidx1) = YlOSc21 ;
        %YlOSc2(p1-p21+1: p1+p22+1, Qidx2) = YlOSc22 ;
        YlOSc2 = YOSc2 ;

    end





    % stitch all cycles back to get cardiogenic artifact free cycles.
    [xqECG, CHECK] = stitchALL(YOSc2, locsf, length(x), p1) ;
    [xqECGl, CHECK] = stitchALL(YlOSc2, locsf, length(x), p1) ;


    xq0 = xqPM + xqECG ;

    SKNAIos0 = x - xq0 ;
    trend = smooth(medfilt1(SKNAIos0, Fs/10), Fs/8, 'loess') ;

    SKNAIos = SKNAIos0 - trend ;
    xq = xq0 + trend ;

    % get ROSDOS result (not yet fine tuned)

    xq0l = xqPM + xqECGl ;

    SKNAIlos0 = x - xq0l ;
    trendl = smooth(medfilt1(SKNAIlos0, Fs/10), Fs/8, 'loess') ;

    SKNAIlos = SKNAIlos0 - trend ;
    xql = xq0l + trendl ;

    %EST = [SKNAIlos(:) SKNAts_suzanna(:) SKNAts_cerutti(:) SKNAts_lp(:)...
    %    SKNAlms(1:length(SKNAIlos)) SKNArls(1:length(SKNAIlos)) SKNAesn(1:length(SKNAIlos)) ...
    %    SKNAekf(1:length(SKNAIlos))'];

    %save(string(Jidx) +"_EST",'EST');


    % Evaluate AR
    if 0
        if PACEMAKER
            %x - xqPM remove pacemaker artificats
            [Metrics_PC] = EvaluatePerformance_AR1(x-xqPM,EST,locsf,Fs);
            save('AR_PC','Metrics_PC');
        else
            [Metrics_AF1] = EvaluatePerformance_AR1(x,EST,locsf(Qidx1),Fs);
            [Metrics_AF2] = EvaluatePerformance_AR1(x,EST,locsf(Qidx2),Fs);
            save('AR_AF','Metrics_AF1','Metrics_AF2');
        end
    end









    if 0
        % calculate SKNA strength
        tmp = SKNAIos(1: floor(length(SKNAIos)/1000) * 1000) ;
        gg = reshape(abs(tmp), 1000, length(tmp)/1000) ;
        SKNAIosAM10 = mean(gg, 1);

        tmp = SKNAIbpf(1: floor(length(SKNAIos)/1000) * 1000) ;
        gg = reshape(abs(tmp), 1000, length(tmp)/1000) ;
        SKNAIbpfAM10 = mean(gg, 1);

        tmp = SKNAIos(1: floor(length(SKNAIos)/Fs) * Fs) ;
        gg = reshape(abs(tmp), Fs, length(tmp)/Fs) ;
        SKNAIosAM1 = resample(mean(gg, 1), 10, 1);

        tmp = SKNAIbpf(1: floor(length(SKNAIos)/Fs) * Fs) ;
        gg = reshape(abs(tmp), Fs, length(tmp)/Fs) ;
        SKNAIbpfAM1 = resample(mean(gg, 1), 10, 1);
    end


    %% plot results
    figure('Position',[1 scrsz(4) scrsz(3)*2/3 scrsz(4)]) ;

    MMM = quantile(xqECG+trend, .998) * 1.2 ;
    mmm = quantile(xqECG+trend, .002) * 1.2 ;
    MMMos = quantile(SKNAIos, .99) * 1.4 ;
    mmmos = quantile(SKNAIos, .01) * 1.4 ;
    

    if PACEMAKER
        ax(1) = subplot_tight(11,1,1,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, x, 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('SKNA-I') ; ylabel('(a)')

        ax(2) = subplot_tight(11,1,2,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, xqPM, 'r', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('PM') ; ylabel('(b)')


        ax(3) = subplot_tight(11,1,3,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, xqECG+trend, 'r', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('ECG') ; ylabel('(c)')

        ax(4) = subplot_tight(11,1,4,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, xq, 'color', [.6 .6 .6], 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('ECG+PM') ; ylabel('(d)')

        ax(5) = subplot_tight(11,1,5,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAIos, 'b', 'linewidth', 2) ;
        axis([-inf inf mmmos MMMos]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-S3') ; ylabel('(e)')

        MMMbpf = quantile(SKNAIbpf, .995) * 1.2 ;
        mmmbpf = quantile(SKNAIbpf, .005) * 1.2 ;
        ax(6) = subplot_tight(11,1,6,[0.008 0.062]) ;
        plot([1:length(SKNAIbpf)]/Fs, SKNAIbpf, 'b', 'linewidth', 2);
        axis([-inf inf mmmbpf MMMbpf]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-BPF') ; ylabel('(f)')

        MMM = quantile(SKNAts_lp, .999) * 1.2 ;
        mmm = quantile(SKNAts_lp, .001) * 1.2 ;
        ax(7) = subplot_tight(11,1,7,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAts_lp(:), 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-TSlp') ; ylabel('(g)')

        MMM = quantile(SKNAlms(1:length(SKNAIlos)), .999) * 1.2 ;
        mmm = quantile(SKNAlms(1:length(SKNAIlos)), .001) * 1.2 ;
        ax(8) = subplot_tight(11,1,8,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAlms(1:length(SKNAIlos)), 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-LMS') ; ylabel('(h)')

        MMM = quantile(SKNAesn(1:length(SKNAIlos)), .999) * 1.2 ;
        mmm = quantile(SKNAesn(1:length(SKNAIlos)), .001) * 1.2 ;
        ax(9) = subplot_tight(11,1,9,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAesn(1:length(SKNAIlos)), 'k', 'linewidth', 2) ;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-ESN') ; ylabel('(i)')

        MMM = quantile(SKNAekf(1:length(SKNAIlos)), .999) * 1.2 ;
        mmm = quantile(SKNAekf(1:length(SKNAIlos)), .001) * 1.2 ;
        ax(10) = subplot_tight(11,1,10,[0.008 0.062]) ;
        plot([1:length(SKNAIbpf)]/Fs, SKNAekf(1:length(SKNAIlos))', 'k', 'linewidth', 2);
        axis([-inf inf mmm MMM]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-EKF') ; ylabel('(j)')

        xlabel('Time (s)') ;
        linkaxes(ax, 'x') ;

        xlim([515 523])
        export_fig(['PMComparisonVVI'], '-transparent', '-pdf')

        xlim([130 138])
        export_fig(['PMComparisonDDD'], '-transparent', '-pdf')


    else

        ax(1) = subplot_tight(9,1,1,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, x, 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('SKNA-I') ; ylabel('(a)')

        ax(2) = subplot_tight(9,1,2,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, xqECG+trend, 'r', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('ECG') ; ylabel('(b)') 


        ax(3) = subplot_tight(9,1,3,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAIos, 'b', 'linewidth', 2) ;
        axis([-inf inf mmmos MMMos]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-S3') ; ylabel('(c)')

        MMMbpf = quantile(SKNAIbpf, .999) * 1.3 ;
        mmmbpf = quantile(SKNAIbpf, .001) * 1.3 ;
        ax(4) = subplot_tight(9,1,4,[0.008 0.062]) ;
        plot([1:length(SKNAIbpf)]/Fs, SKNAIbpf, 'b', 'linewidth', 2);
        axis([-inf inf mmmbpf MMMbpf]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-BPF') ; ylabel('(d)')

        MMM = quantile(SKNAts_lp, .99) * 1.1 ;
        mmm = quantile(SKNAts_lp, .01) * 1.1 ;
        ax(5) = subplot_tight(9,1,5,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAts_lp(:), 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-TSlp') ; ylabel('(e)')

        MMM = quantile(SKNAlms(1:length(SKNAIos)), .99) * 1.1 ;
        mmm = quantile(SKNAlms(1:length(SKNAIos)), .01) * 1.1 ;
        ax(6) = subplot_tight(9,1,6,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAlms(1:length(SKNAIos)), 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-LMS') ; ylabel('(f)')

        MMM = quantile(SKNAesn(1:length(SKNAIos)), .995) * 1.2 ;
        mmm = quantile(SKNAesn(1:length(SKNAIos)), .005) * 1.2 ;
        ax(7) = subplot_tight(9,1,7,[0.008 0.062]) ;
        plot([1:length(x)]/Fs, SKNAesn(1:length(SKNAIos)), 'k', 'linewidth', 2) ;
        axis([-inf inf mmm MMM]) ; set(gca,'xtick',[]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-ESN') ; ylabel('(g)')

        MMM = quantile(SKNAekf(1:length(SKNAIos)), .995) * 1.2 ;
        mmm = quantile(SKNAekf(1:length(SKNAIos)), .005) * 1.2 ;
        ax(8) = subplot_tight(9,1,8,[0.008 0.062]) ;
        plot([1:length(SKNAIbpf)]/Fs, SKNAekf(1:length(SKNAIos))', 'k', 'linewidth', 2);
        axis([-inf inf mmm MMM]) ;
        set(gca, 'fontsize', 18) ; legend('rSNA-EKF') ; ylabel('(h)')

        xlabel('Time (s)') ;
        linkaxes(ax, 'x') ;

        xlim([1016 1040])
        export_fig(['AFComparison'], '-transparent', '-pdf')

    end

    

    if 0
        for jjj = 59 %1: floor(length(SKNAIbpf)/Fs/20)
            xlim([(jjj-1)*20 jjj*20]) ; pause(0.1)
            export_fig(['DEC-' num2str(jjj) '-' ggg(Jidx).name(1:end-4)], '-transparent', '-pdf')
        end
    end



    


    %% plot Figure 1 & Figure 2

    if 0
        figure('Position',[1 scrsz(4) scrsz(3)*3/4 scrsz(4)/2]) ;


        ax(1) = subplot(311) ;
        plot([1:length(x)]/Fs, x, 'k', 'linewidth', 2) ; hold on;
        axis([-inf inf mmm MMM])
        set(gca, 'fontsize', 18) ; legend('SKNAI')


        ax(2) = subplot(312) ;
        plot([1:length(SKNAIbpf)]/Fs, SKNAIbpf, 'b', 'linewidth', 2);
        %hold on ;
        %plot([1:length(SKNAIbpfAM10)]/10, SKNAIbpfAM10, 'r', 'linewidth', 2);
        %plot([1:length(SKNAIbpfAM1)]/10, SKNAIbpfAM1, 'm', 'linewidth', 2);
        axis([-inf inf mmmbpf MMMbpf])
        set(gca, 'fontsize', 18) ; legend('SKNAI-BPF')

        ax(3) = subplot(313) ;
        plot([1:length(x)]/Fs, SKNAIos, 'b', 'linewidth', 2) ;
        %hold on ;
        %plot([1:length(SKNAIosAM10)]/10, SKNAIosAM10, 'r', 'linewidth', 2);
        %plot([1:length(SKNAIosAM10)]/10, SKNAIosAM10, 'm', 'linewidth', 2);
        axis([-inf inf mmmos MMMos])
        set(gca, 'fontsize', 18) ; legend('SKNA')
        xlabel('Time (s)') ;
        linkaxes(ax, 'x') ;

        if PACEMAKER
            xlim([250 255]) ;
        else

        end


        export_fig(['DEC-QQ-' ggg(Jidx).name(1:end-4)], '-transparent', '-pdf')



        %% plot Figure 2
        if 0 && PACEMAKER
            figure('Position',[1 scrsz(4) scrsz(3)*3/4 scrsz(4)]) ;
            ax(1) = subplot(411) ;
            plot([1:length(x)]/Fs-632, xqPM, 'r', 'linewidth', 2) ; hold on;
            axis([-inf inf -300 500]) ; set(gca, 'xtick', [])
            set(gca, 'fontsize', 18) ; legend('Pacing artifact')

            ax(2) = subplot(412) ;
            plot([1:length(x)]/Fs-632, xqECG+trend, 'r', 'linewidth', 2) ; hold on;
            axis([-inf inf -300 500]) ; set(gca, 'xtick', [])
            set(gca, 'fontsize', 18) ; legend('ECG')

            ax(3) = subplot(413) ;
            plot([1:length(x)]/Fs-632, SKNAIos, 'b', 'linewidth', 2) ;
            axis([-inf inf -80 80]) ; set(gca, 'xtick', [])
            set(gca, 'fontsize', 18) ; legend('SKNA')

            ax(4) = subplot(414) ;
            plot([1:length(x)]/Fs-632, x, 'k', 'linewidth', 2) ; hold on;
            axis([-inf inf -300 500]) ;
            set(gca, 'fontsize', 18) ; legend('SKNAI')
            xlabel('Time (s)') ;

            linkaxes(ax, 'x') ;
            xlim([0 6]) ;

            export_fig(['Figure1'], '-transparent', '-pdf')
        end
    end



end
