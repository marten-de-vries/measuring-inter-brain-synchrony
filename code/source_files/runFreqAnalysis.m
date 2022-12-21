function freq = runFreqAnalysis(cfg, data)
    % mostly used to let the user drop EOG electrodes, but could
    % theoretically be used more widely.
    data_selected = ft_selectdata(cfg, data);
    freqcfg = [];
    freqcfg.channel = 'all';
    freqcfg.output = 'fourier';
    freqcfg.method = 'mtmconvol'; % implements multitaper time-frequency transformation; multiplication in frequency domain

    freqcfg.keeptrials = 'yes';
    freqcfg.pad = 'nextpow2';  % efficient fourier transform calculation
    % set frequencies of interest
    switch cfg.band
        case 'alpha'
            freqcfg.foi = [9:14];
        case 'theta'
            freqcfg.foi = [4:7];
        otherwise
            error('Unknown frequency band');
    end
    switch cfg.method
        case 'hann'
            freqcfg.taper = 'hanning';
        case 'multitapers'
            freqcfg.taper = 'dpss';
            freqcfg.tapsmofrq = freqcfg.foi * 0.4; % amount of smoothing
        otherwise
            error('Unknown fourier analysis tapering method');
    end
    freqcfg.t_ftimwin = 4./freqcfg.foi; % 4 cycles per time-window
    freqcfg.toi = cfg.toi;
    freq = ft_freqanalysis(freqcfg, data_selected);
end
