function plane = ComputeCompetitionScore2026(plane, proposalScore, designScore, participationLevel, M1_success)
% Aggregates mission scores + reports/participation into one number.

    % Ensure M1 exists (1 = completed, 0 = not)
    if nargin >= 5 && ~isempty(M1_success)
        plane.performance.M1 = double(logical(M1_success));
    elseif ~isfield(plane.performance,'M1') || ~isfinite(plane.performance.M1)
        plane.performance.M1 = 1;
    end

    % Defaults for safety
    if ~isfield(plane.performance,'M2') || ~isfinite(plane.performance.M2), plane.performance.M2 = 0; end
    if ~isfield(plane.performance,'M3') || ~isfinite(plane.performance.M3), plane.performance.M3 = 0; end
    if ~isfield(plane.performance,'GM') || ~isfinite(plane.performance.GM), plane.performance.GM = 0; end

    % Mission subtotal
    plane.performance.TotalMissionScore = ...
        max(0,plane.performance.M1) + max(0,plane.performance.M2) + ...
        max(0,plane.performance.M3) + max(0,plane.performance.GM);

    % Reports normalization (0..1)
    if ~isfinite(proposalScore), proposalScore = 0; end
    if ~isfinite(designScore),   designScore   = 0; end
    plane.performance.TotalReportScore = (proposalScore + designScore) / 200;

    % Participation bump (1..3 -> 0, 0.5, 1.0)
    if ~isfinite(participationLevel), participationLevel = 1; end
    participationLevel = max(1, min(3, round(participationLevel)));
    partBump = (participationLevel - 1) * 0.5;

    % Final competition score (keeps ~1.0e3 magnitude like last year)
    plane.performance.CompetitionScore = ...
          1000 * plane.performance.TotalMissionScore ...
        +   20 * plane.performance.TotalReportScore ...
        +   10 * partBump;
end
