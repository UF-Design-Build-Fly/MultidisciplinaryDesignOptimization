function plane = MissionGMScore(plane, bestGM_s, yourGM_s)
% Ground Mission normalization for 2026
% GM = best_time / our_time  (higher is better, <= 1 if you're slower than best)
% Inputs:
%   plane      : Airplane object (struct with .performance)
%   bestGM_s   : Best GM time observed among teams (seconds)
%   yourGM_s   : Your GM time (seconds)

    % Fallbacks if arguments are missing/invalid (pull from plane if present)
    if nargin < 2 || ~isfinite(bestGM_s) || bestGM_s <= 0
        if isfield(plane, 'performance') && isfield(plane.performance, 'bestGM_s') ...
                && isfinite(plane.performance.bestGM_s) && plane.performance.bestGM_s > 0
            bestGM_s = plane.performance.bestGM_s;
        else
            bestGM_s = 30; % safe default
        end
    end
    if nargin < 3 || ~isfinite(yourGM_s) || yourGM_s <= 0
        if isfield(plane, 'performance') && isfield(plane.performance, 'yourGM_s') ...
                && isfinite(plane.performance.yourGM_s) && plane.performance.yourGM_s > 0
            yourGM_s = plane.performance.yourGM_s;
        else
            yourGM_s = 45; % safe default
        end
    end

    gm = bestGM_s / max(yourGM_s, eps);    % 0..1+ (cap for safety)
    gm = max(0, min(gm, 2));               % clamp to a reasonable range

    plane.performance.GM      = gm;
    plane.performance.scoreGM = gm;        % legacy alias, if other scripts read it
end
