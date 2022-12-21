function out = calculateMeasures(spect1, spect2, calc_robust)
    % enable the following two lines if you want to calculate measures over
    % trials instead of time.
    %
    % spect1 = permute(spect1, [4 2 3 1]);
    % spect2 = permute(spect2, [4 2 3 1]);

    dim = 4;
    freqdim = 3;
    % fourierspctrm: trials, channels, freq, time
    phi = angle(spect1); % extract phase
    psi = angle(spect2);

    % Enable the code below to downsample by half using averaging (extra)
    % phi = squeeze(mean(reshape(phi, size(phi, 1), size(phi, 2), size(phi, 3), 2, []), 4));
    % psi = squeeze(mean(reshape(psi, size(psi, 1), size(psi, 2), size(psi, 3), 2, []), 4));

    % Enable the code below to downsample by half by throwing out half the
    % data (extra)
    % phi = phi(:, :, :, (0:size(phi, 4) / 2 - 1) * 2 + 1);
    % psi = psi(:, :, :, (0:size(psi, 4) / 2 - 1) * 2 + 1);

    plv_by_freq = abs(mean(exp(1i*(phi - psi)), dim));
    out.plv = mean(plv_by_freq, freqdim); % average over frequency

    % Formula taken from Burgess (2013).
    %
    % Validated against a loop calling circ_corrcc from the MATLAB circular
    % statistics toolbox (Berens, 2009). Similar implementation, but is
    % vectorized so not a direct copy.
    phi_bar = angle(sum(exp(1i * phi), dim));
    psi_bar = angle(sum(exp(1i * psi), dim));
    sin_phi_diff = sin(phi - phi_bar);
    sin_psi_diff = sin(psi - psi_bar);

    num = sum(sin_phi_diff .* sin_psi_diff, dim);
    denom = sqrt(sum(sin_phi_diff.^2, dim) .* sum(sin_psi_diff.^2, dim));
    ccorr_by_freq = num ./ denom; % note: denom will be zero if a signal is all zero, resulting in NaN.
    out.ccorr = mean(ccorr_by_freq, freqdim); % average over frequency

    % Inspired by:
    % https://stackoverflow.com/questions/41923700/imaginary-part-of-coherence-matlab
    % approximate expected values using mean
    conjSpect2 = conj(spect2);
    Sij = mean(spect1 .* conjSpect2, dim);
    Sii = mean(spect1 .* conj(spect1), dim);
    Sjj = mean(spect2 .* conjSpect2, dim);
    % calculate coherency
    coherency = Sij ./ sqrt(Sii .* Sjj);
    % and take the imaginary part as our measure
    imagcoh_by_freq = imag(coherency);
    out.imagcoh = mean(imagcoh_by_freq, freqdim); % average out frequency

    % concept: Pewsey et al. 2013
    dist = @(a, b) pi - abs(pi - abs(a - b));
    if (calc_robust)
        % validated against circular correlation by showing that the robust
        % measure still correlates with it.
        p = 0.05; % remove 5% of 'outliers'
        n = size(phi, 4);
        target_length = n * (1 - p);
        while n > target_length
            distsum = zeros(size(phi));
            for trial = 1:size(phi, 1)
                for channel = 1:size(phi, 2)
                    for freq = 1:size(phi, 3)
                        for point = 1:size(phi, 4)
                            s = sum(dist(phi(trial, channel, freq, point), phi(trial, channel, freq, :))) + ...
                                sum(dist(psi(trial, channel, freq, point), psi(trial, channel, freq, :)));
                            distsum(trial, channel, freq, point) = s;
                        end
                    end
                end
            end
            phinew = zeros(size(phi, 1), size(phi, 2), size(phi, 3), size(phi, 4) - 1);
            psinew = phinew;  % copy
	        [~, trim_i] = max(distsum, [], 4);
            for trial = 1:size(phi, 1)
                for channel = 1:size(phi, 2)
                    for freq = 1:size(phi, 3)
                        keep = setdiff(1:size(phi, 4), trim_i(trial, channel, freq));
                        phinew(trial, channel, freq, :) = phi(trial, channel, freq, keep);
                        psinew(trial, channel, freq, :) = psi(trial, channel, freq, keep);
                    end
                end
            end

            phi = phinew;
            psi = psinew;
            n = size(phi, 4)
        end
        % robust circular correlation
        phi_bar = angle(sum(exp(1i * phi), dim));
        psi_bar = angle(sum(exp(1i * psi), dim));
        sin_phi_diff = sin(phi - phi_bar);
        sin_psi_diff = sin(psi - psi_bar);
    
        num = sum(sin_phi_diff .* sin_psi_diff, dim);
        denom = sqrt(sum(sin_phi_diff.^2, dim) .* sum(sin_psi_diff.^2, dim));
        rob_ccorr_by_freq = num ./ denom; % note: denom will be zero if a signal is all zero, resulting in NaN.
        out.robccorr = mean(rob_ccorr_by_freq, freqdim); % average over frequency
    end
end
