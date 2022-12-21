%% Setup
% Reset variable workspace, close all windows, and clear command window
clear all;
close all;
clc;


% set working directory to current location of this script
fullres = matlab.desktop.editor.getActive;
cd(fileparts(fullres.Filename));
addpath(genpath(pwd)) % add subdirectories to include source files

% load fieldtrip
addpath 'Y:\staff\LowCost\fse\bi\cogmod\vanvugt\Lionel\tacit coordination\analyses\EEG\fieldtrip-20211102'
addpath 'Y:\staff\LowCost\fse\bi\cogmod\vanvugt\Lionel\tacit coordination\analyses\EEG\Marten\circstat-matlab-master'
ft_defaults

% DATA PREPROCESSING:
fileroots = {'2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', ...
    '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', ...
    '27', '28', '29', '30', '31', '33', '35', '36', '37', '38', '39', ...
    '40', '41', '42'}; % session numbers to analyze
% Not included due to recording problems etc.: 1, 26, 32, 34

paths.rawDataDir = 'Y:\staff\LowCost\fse\bi\cogmod\vanvugt\Lionel\tacit coordination\data\main\EEG\raw_data';
paths.cleanDataDir = 'Y:\staff\LowCost\fse\bi\cogmod\vanvugt\Lionel\tacit coordination\data\main\EEG\processed_data\temp';
paths.resultsDir = 'Y:\staff\LowCost\fse\bi\cogmod\vanvugt\Lionel\tacit coordination\analyses\EEG\Marten\results';
paths.behavioralDataDir = 'Y:\staff\LowCost\fse\bi\cogmod\vanvugt\Lionel\tacit coordination\data\main\behavioral';
data_in = loadData(fileroots, paths);

% condition triggers:
% low color: 101
% low shape: 102
% high color: 103
% high shape: 104

%% Load behavioural data

behavData = readtable([paths.behavioralDataDir '\behavior_data.csv']);
% reconstruct original (and EEG-matching) trial numbers
behavData.rawTrial = (behavData.block - 1) * 90 + behavData.trial;
% uniquely identify each stimulus in a single column
behavData.item = strcat(behavData.stim_type, num2str(behavData.image_set));
for sessionNum = fileroots
    session = ['sess' sessionNum{1}];

    % get the data for the current session only
    sessData = behavData(behavData.session == str2double(sessionNum{1}), :);
    % put the table in the same order as 'trialinfo', dropping rows for
    % which no EEG data exists
    info = data_in.(session).subj1.trialinfo;
    [t, i] = ismember(info(:, 1), sessData.rawTrial);
    % attach the resulting table to the EEG data for later use
    data_in.(session).info = sessData(i(t), :);
end

%% Frequency analysis of the actual data

% Run time-frequency analysis for each session & subject & electrode
cfg = [];
cfg.band = 'alpha'; % alpha or theta
cfg.method = 'hann'; % hann or multitapers
% no NaNs with this range. And it nicely matches the 1s wait before the user can answer.
cfg.toi = 0:0.01:0.99; % 0.01 (100 samples) or 0.02 (50 samples)
% channels of special interest:
%   '1', '2', '2','3','4','27','28','29','30','31' ... % frontal
% , '21' ... % right temporo-parietal
% , '9','22','13' ... % centro-parietal
% , '14' ... % left parieto-occipital
% , '18','17' ... % right occipital
cfg.channel =  {'*A1', '*A2', '*A3', '*A4', '*A5', '*A6', '*A7', '*A8', ...
                '*A9', '*A10', '*A11', '*A12', '*A13', '*A14', '*A15', ...
                '*A16', '*A17', '*A18', '*A19', '*A20', '*A21', '*A22', ...
                '*A23', '*A24', '*A25', '*A26', '*A27', '*A28', '*A29', ...
                '*A30', '*A31', '*A32'};;
cfg.trialshuffle = 'time'; % time or spectrum
cfg.calc_robust = true;

for sessionNum = fileroots
    session = ['sess' sessionNum{1}]
    % the code below is useful if you want to do a full res frequency
    % analysis. Relatively slow.
    %
    % fullres = data_in.(session).subj1.time{1};
    % cfg.toi = fullres(fullres >= 0 & fullres < 1);  % match original sampling rate

    freq.(session).subj1 = runFreqAnalysis(cfg, data_in.(session).subj1);
    freq.(session).subj2 = runFreqAnalysis(cfg, data_in.(session).subj2);
    % make the behavioural data (the item shown, mostly) accessable from
    % the frequency spectrum object too.
    freq.(session).info = data_in.(session).info;
end


%% 'single trial ERPs'
savefilename = prepSave(paths.resultsDir, ['singletrialerp.dat']);
for sessionNum = fileroots
    session = ['sess' sessionNum{1}];
    result = [];
    result.p3subj1 = singletrialerp(cfg.channel, data_in.(session).subj1);
    result.p3subj2 = singletrialerp(cfg.channel, data_in.(session).subj2);

    sessResult = createSessionTable(result, data_in.(session).subj1.trialinfo, cfg.channel);
    sessResult.session = repelem(sessionNum, size(sessResult, 1))';
    sessResult.sample = repelem(NaN, size(sessResult, 1))';
    writetable(sessResult, savefilename, 'Delimiter','\t', 'WriteMode', 'Append', 'WriteVariableNames', false);
end

%% Calculate inter-brain synchrony measures for the actual data

% i.e. for every trial within a session between homologous electrodes of
% both participants.

savefilename = prepSave(paths.resultsDir, ['synch_' cfg.band '.prototype.dat']);
for sessionNum = fileroots
    session = ['sess' sessionNum{1}];

    % Calculate the measures using the full spectrum for both participants.
    % A spectrum is a 4d array:
    % ~180 trials x 32 channels x 4 frequencies x 100 time points
    measures = calculateMeasures(freq.(session).subj1.fourierspctrm, freq.(session).subj2.fourierspctrm, cfg.calc_robust);
    % put the calculated measures into a 'long' table with one row per
    % measure value.
    sessResult = createSessionTable(measures, freq.(session).subj1.trialinfo, cfg.channel);
    % add the session number to this table
    sessResult.session = repelem(sessionNum, size(sessResult, 1))';
    % and set sample (used to determine which iteration we are in when
    % creating null distributions) to NaN as it's undefined for actual
    % data.
    sessResult.sample = repelem(NaN, size(sessResult, 1))';
    % write the calculated table directly to disk (appending)
    writetable(sessResult, savefilename, 'Delimiter','\t', 'WriteMode', 'Append', 'WriteVariableNames', false);
end


%% 'shuffled trials' null distribution calculation
savefilename = prepSave(paths.resultsDir, ['synch_' cfg.band '_' cfg.trialshuffle 'shuffled.prototype.dat']);

for sessionNum = fileroots
    session = ['sess' sessionNum{1}];
    for sample = 1:200
        {session sample}  % useful to see where we are in the process.

        % enable of of the following:

        switch cfg.trialshuffle
            case 'time'
                % shuffle the original time series data
                subj2 = data_in.(session).subj2;
                subj2.trial = cellfun(@shuffleTrial, subj2.trial, 'UniformOutput', false);
                freq2spctrm = runFreqAnalysis(cfg, subj2).fourierspctrm;
            case 'spectrum'
                % shuffle the already calculated fourier spectrum. Faster, but not
                % quite identical.
                indices = randperm(size(freq.(session).subj2.fourierspctrm, 4));
                freq2spctrm = freq.(session).subj2.fourierspctrm(:, :, :, indices);
            otherwise
                error('Unknown shuffle method');
        end

        % same process as for actually observed data.
        measures = calculateMeasures(freq.(session).subj1.fourierspctrm, freq2spctrm, false);
        sampleResult = createSessionTable(measures, freq.(session).subj1.trialinfo, cfg.channel);
        sampleResult.session = repelem(sessionNum, size(sampleResult, 1))';
        sampleResult.sample = repelem(sample, size(sampleResult, 1))';
        writetable(sampleResult, savefilename, 'Delimiter','\t', 'WriteMode', 'Append', 'WriteVariableNames', false);
    end
end

%% 'random dyads' null distribution calculation
% Here, it's actually feasible to calculate all possible dyad permutations.
% (Mostly because the frequency spectra can be re-used without
% recalculating them.)

savefilename = prepSave(paths.resultsDir, ['synch_' cfg.band '_dyadshuffled.prototype.dat']);
for sessionNumA = fileroots
    sessionA = ['sess' sessionNumA{1}];
    a = freq.(sessionA);

    % filter out the current session. Either that value has already been
    % computed when analysing the real data, or it is trivial (a person's
    % synchrony with themselves)
    for sessionNumB = fileroots(~strcmp(fileroots, sessionNumA))
        sessionB = ['sess' sessionNumB{1}];
        b = freq.(sessionB);
        % indexes_a will filter out rows with items for which no data is
        % available in sessionB.
        % indexes_b will filter out rows with items for which no data is
        % available in sessionA, and will additionally change the order
        % such that all items line up across sessions.
        [indexes_a i] = ismember(a.info.item, b.info.item);
        indexes_b = i(indexes_a);

        [sessionA sessionB]  % to keep track of progress
 
        % there are four possible 'fake dyad' combinations.
        for subjA = {'subj1', 'subj2'}
            for subjB = {'subj1', 'subj2'}
                % actually get the spectra lined up
                aSpectrum = a.(subjA{1}).fourierspctrm(indexes_a, :, :, :);
                bSpectrum = b.(subjB{1}).fourierspctrm(indexes_b, :, :, :);
                % calculate the measures and store the result. Same
                % approach as for the real data above.
                measures = calculateMeasures(aSpectrum, bSpectrum, false);
                fakeResult = createSessionTable(measures, a.subj1.trialinfo(indexes_a, :), cfg.channel);

                id = [sessionA '-' subjA{1} '_' sessionB '-' subjB{1}];
                fakeResult.session = repelem(NaN, size(fakeResult, 1))';
                fakeResult.sample = repelem({id}, size(fakeResult, 1))';
                writetable(fakeResult, savefilename, 'Delimiter','\t', 'WriteMode', 'Append', 'WriteVariableNames', false);
            end
        end
    end
end

%% Export raw frequency analysis data for plotting
% Inspect result
if strcmp(cfg.band, 'alpha')
    subj1 = data_in.sess2.subj1;

    % raw data
    writematrix(subj1.trial{1}(13, :)', prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-raw.csv'));
    m = squeeze(freq.sess2.subj1.fourierspctrm(1, 13, :, :));
    writematrix(abs(m), prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-amplitude.csv'));
    writematrix(angle(m), prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-phase.csv'));

    % cwt on data
    cwtpath = prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-cwt.csv');
    cwt = cwtft({subj1.trial{1}(13, :), 1 / 512}, 'wavelet', 'mexh');
    for scale = 1:length(cwt.scales)
        result = [repelem(cwt.scales(scale), length(subj1.time{1}))' subj1.time{1}' cwt.cfs(scale, :)'];
        writematrix(result, cwtpath, 'Delimiter', '\t', 'WriteMode', 'append');
    end

    % shuffled data
    % TODO: also implement for multitapers & different window sizes?
    subj1.trial = cellfun(@shuffleTrial, subj1.trial, 'UniformOutput', false);
    m2 = squeeze(runFreqAnalysis(cfg, subj1).fourierspctrm(1, 13, :, :));
    indices = randperm(size(freq.sess2.subj1.fourierspctrm, 4));
    m3 = squeeze(freq.sess2.subj2.fourierspctrm(1, 13, :, indices));

    writematrix(abs(m2), prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-timeshuffled-amplitude.csv'));
    writematrix(angle(m2), prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-timeshuffled-phase.csv'));

    writematrix(abs(m3), prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-spectrumshuffled-amplitude.csv'));
    writematrix(angle(m3), prepSave(paths.resultsDir, 'sess2-subj1-alpha-trial1-spectrumshuffled-phase.csv'));

    % useful for the power analysis
    subj1.trial = cellfun(@shuffleTrial, subj1.trial, 'UniformOutput', false);
    spctrm = runFreqAnalysis(cfg, subj1).fourierspctrm;
    % image(reshape(abs(spctrm), [], 100), 'CDataMapping','scaled'); colorbar;
    ampl = squeeze(mean(abs(spctrm), [1, 2, 3]));
    % approximates 'ones()':
    plot(ampl); ylim([0 max(ampl)]);
    mean(ampl)

    phasewrapped = reshape(angle(spctrm), [], 100)';
    % don't limit the phase to [-pi, pi] but keep increasing
    phase = phasewrapped + cumsum([zeros(1, size(phasewrapped, 2)); diff(phasewrapped)] < 0) * 2 * pi;
    image(phase', 'CDataMapping','scaled'); colorbar;

    % on average ~0.75
    slope = mean(phase(100, :) - phase(1, :)) / 100

    % export phase data for both participants in first session, trial & Pz
    % electrode
    a = squeeze(angle(freq.sess2.subj1.fourierspctrm(1, 13, 3, :)));
    b = squeeze(angle(freq.sess2.subj2.fourierspctrm(1, 13, 3, :)));
    writematrix([a b], prepSave(paths.resultsDir, 'sess2-alpha-trial1-phases.csv'));
end

%% Validation of our imaginary coherence & plv impls against fieldtrip

% figure
% plot(calculateMeasures(freq.sess2.subj1.fourierspctrm, freq.sess2.subj2.fourierspctrm).plv(:, 1));
% 
% newfreq = freq.sess2.subj1;
% newfreq.label = [freq.sess2.subj1.label; freq.sess2.subj2.label];
% newfreq.fourierspctrm = cat(2, freq.sess2.subj1.fourierspctrm, freq.sess2.subj2.fourierspctrm);
% 
% cfg = [];
% cfg.method = 'coh'; % or plv
% cfg.channelcmb = {'1-A1'  '2-A1'};
% cfg.complex = 'imag'; % or 'abs' when method 'plv'
% result = ft_connectivityanalysis(cfg, newfreq);
% 
% figure
% if isfield(result, 'cohspctrm')
%     plot(squeeze(mean(result.cohspctrm, 2)))
% else
%     plot(squeeze(mean(result.plvspctrm, 2)))
% end
