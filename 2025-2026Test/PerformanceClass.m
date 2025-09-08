classdef PerformanceClass
  properties
    % Mass mapping (2026 payloads)
    massPerPassenger = 0.0185     % kg per rubber duck
    massPerCargo     = 0.17010    % kg per hockey puck
    x1_addon_kg      = 0.249      % keep or set 0 if unused

    % Payload selections
    passengers = 9; cargo = 2
    declaredPmax = 9; declaredCmax = 2

    % Banner geometry
    bannerLength_in = 20; bannerHeight_in = 4

    % Weights
    m2Weight = -1
    totalEmptyWeight = -1; totalWeight2 = -1; totalWeight3 = -1

    % Aeroperf
    velocity2 = -1; velocity3 = -1
    takeoffDist2 = -1; takeoffDist3 = -1
    landingSpeed2 = -1; landingSpeed3 = -1
    dynamicThrust = -1
    time2 = -1; time3 = -1
    numLaps2 = -1; numLaps3 = -1; laps2 = -1; laps3 = -1

    % Drag snapshot
    inducedDrag = -1; parasiticDrag = -1; totalDrag = -1; skinDrag = -1
    wingPara=-1; hStabPara=-1; vStabPara=-1; fusePara=-1; gearPara=-1
    drag2=-1; drag3=-1

    % Scoring (2026)
    peerMaxNetIncome = 500; peerMaxBanner = 80
    bestGM_s = 30; yourGM_s = 45
    proposalScore = 90; designReportScore = 95; participationLevel = 3
    M1 = 1; M2 = 0; M3 = 0; GM = 0
    score2 = -1; score3 = -1; scoreGM = -1
    TotalMissionScore = -1; TotalReportScore = -1; CompetitionScore = -1
    scoreTotal = -1 % legacy alias
  end
end
