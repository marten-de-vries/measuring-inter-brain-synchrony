function [out] = shuffleTrial(in)
    % shuffles a single trial (a matrix of electrode x time) over time
    out = in(:, randperm(size(in, 2)));
end
