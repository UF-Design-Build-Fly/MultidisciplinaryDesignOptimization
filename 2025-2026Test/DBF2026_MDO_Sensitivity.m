% DBF2026_MDO_SENSITIVITY.m — Multiairfoil MDO + sensitivity (Wichita 2026)
%%  Setup 
clc; clear; close all; warning('off','all');
addpath(genpath(pwd));         % ensure project files are on path

% ===== Strategy knobs =====
DO_REFINE         = true;      % run a second pass around winners
KEEP_TOP_COARSE   = 200;       % carry this many into refine pass

% ===== Site (Wichita) =====
site.name    = 'Wichita, Kansas'; %#ok<NASGU>
site.alt_m   = 406;                               % ~1330 ft
site.temp_K  = 298;                               % 25 °C
site.wind_ms = 4.0;                               % nominal wind
site.rho     = AtmosISA(site.alt_m, site.temp_K); % density at ambient T

% ===== Rules/knobs =====
time.window_s       = 300;    % 5 min
mission.turnMult    = 0.8;    % turn-speed multiplier

rules.span_ft_min   = 3;
rules.span_ft_max   = 5;
rules.banner_min_in = 10;
rules.banner_AR_max = 5;      % L/H <= 5 -> H >= L/5 and H >= 2"
rules.battery_Wh_max= 100;
rules.p_ge_3c       = true;   % passengers >= 3*cargo

% Scoring normalization (update at COMP)
norm.peerMaxNetIncome = 500;
norm.peerMaxBanner    = 80;
norm.bestGM_s         = 30;    % best ground mission time among all teams
norm.yourGM_s         = 45;    % your current GM time guess

% Reports / participation
report.proposal  = 90;
report.design    = 95;
participation    = 3;          % 1 attend, 2 tech pass, 3 attempted flight
M1_success       = 1;

%  Payload mapping (2026) 
payload.massPerPassenger = 0.0185;   % kg per rubber duck
payload.massPerCargo     = 0.17010;  % kg per hockey puck
payload.x1_addon_kg      = 0.249;    % legacy X1 mass (set 0 if unused)

%  Motor table 
load('MotorSpreadsheet2025.mat');     % -> MotorSpreadsheet (table)
% Hygiene filter & sort
MotorSpreadsheet = MotorSpreadsheet( ...
      MotorSpreadsheet.PowerConsumptionW > 0 ...
   &  MotorSpreadsheet.ThrustKg          > 0 ...
   &  MotorSpreadsheet.PitchSpeedms      > 10 , :);
MotorSpreadsheet = sortrows(MotorSpreadsheet, 'SortScore', 'descend');

%%  PASS 1 (coarse) 
fprintf('\n=== PASS 1 (coarse) ===\n');

% Coarse grids (broad coverage)
aspectRatios      = 4:1:8;           % 4,5,6,7,8
spans_ft          = 3:0.5:5;         % 3.0,3.5,4.0,4.5,5.0
spans_m           = spans_ft * 0.3048;

pOptions          = 6:12;            % ducks (coarse)
cOptions          = 0:2;             % pucks
bannerLengths_in  = 12:4:28;         % banner length (H enforced later)

numPowerSystems   = min(height(MotorSpreadsheet), 24);  % speed

% Wing lookup (all airfoils supported by GenWingData)
fprintf('Generating wing lookup table (coarse)...\n');
try
    [wingsLUT, airfoilTable] = GenWingData(aspectRatios, spans_m);
catch
    wingsLUT = GenWingData(aspectRatios, spans_m);
    airfoilTable = [];
end
AFn = size(wingsLUT,1);
airfoilIndices = 1:AFn;
fprintf('Airfoils evaluated: %d\n', AFn);

% Preallocate (conservative)
ARn = numel(aspectRatios); Sn = numel(spans_m);
PnP = numel(pOptions); Cn = numel(cOptions); Ln = numel(bannerLengths_in);
upperBound = ARn*Sn*AFn*numPowerSystems*PnP*Cn*Ln;

planes(upperBound,1) = struct(AirplaneClass);
scores(upperBound,1) = -inf;
k = 0;

wb = waitbar(0,'PASS 1 (coarse)...');

t0 = tic;
for iAR = 1:ARn
  for iS = 1:Sn
    WS_ft = spans_m(iS)*3.28084;
    if WS_ft < rules.span_ft_min || WS_ft > rules.span_ft_max, continue; end

    for iAF = airfoilIndices
      for iPS = 1:numPowerSystems

        % progress
        baseFrac = ((iAR-1)*Sn*AFn*numPowerSystems + (iS-1)*AFn*numPowerSystems + (iAF-1)*numPowerSystems + iPS) ...
                 /(ARn*Sn*AFn*numPowerSystems);
        if mod(iPS,3)==1
            if isempty(airfoilTable)
                afLabel = sprintf('AF %d/%d', iAF, AFn);
            else
                afLabel = airfoilTable(iAF).name;
            end
            waitbar(baseFrac, wb, sprintf('Wing/Power %.0f%% — AR=%.1f Span=%.2fm %s', 100*baseFrac, aspectRatios(iAR), spans_m(iS), afLabel));
        end

        % ---------- Assemble wing/power/static pieces ----------
        template = struct(AirplaneClass);

        template.wing.span        = spans_m(iS);
        template.wing.aspectRatio = aspectRatios(iAR);
        template.wing             = WingClass.SetWingData(template.wing, wingsLUT, iAF, iAR, iS);

        template.powerSystem = PowerClass.SetPowerSystemData(template.powerSystem, MotorSpreadsheet, iPS);
        if ~isfinite(template.powerSystem.batteryCapacity) || template.powerSystem.batteryCapacity <= 0
            template.powerSystem.batteryCapacity = rules.battery_Wh_max;
        end
        template.powerSystem.batteryCapacity = min(rules.battery_Wh_max, template.powerSystem.batteryCapacity);

        template.fuselage  = FuselageClass.CalcFuselageData(template);
        template.fuselage  = FuselageClass.GenLandingGear(template);
        template.empennage = EmpennageClass.GenEmpennage(template, 4, 2);

        % Cheap constants for early pruning
        rho = site.rho; g = 9.80665;
        S = template.wing.planformArea;
        CLmax = max(template.wing.clFlap, template.wing.clm);
        pitchSpeed = max(1, template.powerSystem.propSpeed);   % m/s
        Tstatic_N  = 9.80665 * template.powerSystem.thrust;

        % ---------- Payload/banner sweeps ----------
        for p = pOptions
          for c = cOptions
            if rules.p_ge_3c && p < 3*c, continue; end

            mass_ducks = p * payload.massPerPassenger;
            mass_pucks = c * payload.massPerCargo;

            for L_in = bannerLengths_in
              pl = template;

              % Payload & banner
              pl.performance.passengers   = p;
              pl.performance.cargo        = c;
              pl.performance.declaredPmax = p;
              pl.performance.declaredCmax = c;

              pl.performance.bannerLength_in = L_in;
              pl.performance.bannerHeight_in = max(2, L_in/5); % H >= max(2, L/5)

              pl.performance.massPerPassenger = payload.massPerPassenger;
              pl.performance.massPerCargo     = payload.massPerCargo;
              pl.performance.x1_addon_kg      = payload.x1_addon_kg;

              % Scoring context
              pl.performance.peerMaxNetIncome   = norm.peerMaxNetIncome;
              pl.performance.peerMaxBanner      = norm.peerMaxBanner;
              pl.performance.bestGM_s           = norm.bestGM_s;
              pl.performance.yourGM_s           = norm.yourGM_s;
              pl.performance.proposalScore      = report.proposal;
              pl.performance.designReportScore  = report.design;
              pl.performance.participationLevel = participation;

              % Totals
              pl = FindTotalWeight(pl);

              % ===== Early pruning =====
              W2_kg = pl.performance.totalEmptyWeight + mass_ducks + mass_pucks + payload.x1_addon_kg;
              W2_N  = W2_kg * g;

              Vmin_lvl = sqrt( 2*W2_N / (rho*S*CLmax) );
              if Vmin_lvl > 0.95 * pitchSpeed, continue; end

              if (Tstatic_N / W2_N) < 0.24, continue; end

              % ===== Full performance =====
              pl = GenVelocityTest(pl, 2, site.rho, site.temp_K);
              if (pl.performance.velocity2/1.5 < pl.performance.landingSpeed2) || (pl.performance.velocity2==-1)
                  continue;
              end

              pl = GenVelocityTest(pl, 3, site.rho, site.temp_K);
              if (pl.performance.velocity3/1.5 < pl.performance.landingSpeed3) || (pl.performance.velocity3==-1)
                  continue;
              end

              if (pl.powerSystem.batteryCapacity > rules.battery_Wh_max)
                  pl.powerSystem.batteryCapacity = rules.battery_Wh_max;
              end
              if rules.p_ge_3c && (pl.performance.declaredPmax < 3*pl.performance.declaredCmax)
                  continue;
              end

              try
                  pl = Mission2Score(pl, site.wind_ms, mission.turnMult, site.rho, time.window_s, norm.peerMaxNetIncome);
              catch
                  pl = Mission2Score(pl, site.wind_ms, mission.turnMult, site.rho);
              end
              if pl.performance.time2 > time.window_s, continue; end

              try
                  pl = Mission3Score(pl, site.wind_ms, mission.turnMult, site.rho, time.window_s, norm.peerMaxBanner);
              catch
                  pl = Mission3Score(pl, site.wind_ms, mission.turnMult, site.rho);
              end

              try
                  pl = MissionGMScore(pl, norm.bestGM_s, norm.yourGM_s);
              catch
                  pl = MissionGMScore(pl);
              end

              try
                  pl = ComputeCompetitionScore2026(pl, report.proposal, report.design, participation, M1_success);
              catch
                  pl.performance.CompetitionScore = pl.performance.score2Normalized + pl.performance.score3Normalized + pl.performance.scoreGMNormalized;
              end

              k = k + 1;
              if k > numel(scores)
                  % grow in chunks
                  scores(end+10000,1) = -inf; %#ok<AGROW>
                  planes(end+10000,1) = struct(AirplaneClass); %#ok<AGROW>
              end
              planes(k) = pl;
              scores(k) = pl.performance.CompetitionScore;

            end % L_in
          end % c
        end % p

      end % motor
    end % airfoil
  end % span
end % AR

elapsed1 = toc(t0);
try, close(wb); end
fprintf('PASS 1 done in %.1f s. Kept %d designs.\n', elapsed1, k);

% Trim to count
planes = planes(1:k); scores = scores(1:k);
if isempty(scores) || ~any(isfinite(scores))
    error('No feasible designs in PASS 1 — relax pruning or grids.');
end

% Rank & keep top for refine
[topScores1, order1] = sort(scores, 'descend'); %#ok<ASGLU>
K1 = min(KEEP_TOP_COARSE, numel(order1));
topCoarse = planes(order1(1:K1));

%%  PASS 2 (refine) 
if DO_REFINE
fprintf('\n=== PASS 2 (refine) around top %d ===\n', K1);

% Neighborhoods around coarse winners (simple & robust)
% Spans (ft)
spans_ft_coarse = zeros(1, K1);
for i=1:K1, spans_ft_coarse(i) = topCoarse(i).wing.span*3.28084; end
baseSpan = round(spans_ft_coarse*4)/4; % to nearest 0.25 ft
span_set_ft = unique([baseSpan-0.25, baseSpan, baseSpan+0.25]);
span_set_ft = span_set_ft(span_set_ft >= rules.span_ft_min & span_set_ft <= rules.span_ft_max);
span_set_m  = span_set_ft * 0.3048;

% Aspect ratios
AR_seed = zeros(1, K1);
for i=1:K1, AR_seed(i) = topCoarse(i).wing.aspectRatio; end
baseAR = round(AR_seed*2)/2;          % to nearest 0.5
AR_set = unique([baseAR-0.5, baseAR, baseAR+0.5]);
AR_set = AR_set(AR_set >= 4 & AR_set <= 8);

% Payload/banner
p_seed = zeros(1, K1); c_seed = zeros(1, K1); L_seed = zeros(1, K1);
for i=1:K1
    p_seed(i) = topCoarse(i).performance.passengers;
    c_seed(i) = topCoarse(i).performance.cargo;
    L_seed(i) = topCoarse(i).performance.bannerLength_in;
end
p_set    = unique(min(14, max(3, [p_seed-2, p_seed, p_seed+2])));
c_set    = unique(min(3,  max(0, [c_seed-1, c_seed, c_seed+1])));
L_set_in = unique(min(30, max(10, [L_seed-4, L_seed, L_seed+4])));

% Motors used most frequently in coarse winners
motorNames = strings(1, K1);
for i=1:K1, motorNames(i) = string(topCoarse(i).powerSystem.motorName); end
[unames,~,idxc] = unique(motorNames);
counts = accumarray(idxc,1);
[~,ordMot] = sort(counts,'descend');
primeMotorNames = unames(ordMot(1:min(numel(unames), 12)));
candIdx = find(ismember(string(MotorSpreadsheet.Motor), primeMotorNames));
if isempty(candIdx), candIdx = 1:min(height(MotorSpreadsheet), 20); end

% Wing LUT for refine neighborhood
try
    [wingsLUT2, airfoilTable2] = GenWingData(AR_set, span_set_m); %#ok<NASGU>
catch
    wingsLUT2 = GenWingData(AR_set, span_set_m);
end
AFn2 = size(wingsLUT2,1);
airfoilIndices2 = 1:AFn2;

% Preallocate (approx)
ARn2 = numel(AR_set); Sn2 = numel(span_set_m);
PnP2 = numel(p_set);  Cn2 = numel(c_set); Ln2 = numel(L_set_in);
upper2 = max(1, ARn2*Sn2*AFn2*numel(candIdx)*PnP2*Cn2*Ln2);

planes2(upper2,1) = struct(AirplaneClass);
scores2(upper2,1) = -inf;
k2 = 0;

wb2 = waitbar(0,'PASS 2 (refine)...');
t1 = tic;

for iAR = 1:ARn2
  for iS = 1:Sn2
    WS_ft = span_set_m(iS)*3.28084;
    if WS_ft < rules.span_ft_min || WS_ft > rules.span_ft_max, continue; end

    for iAF = airfoilIndices2
      for ii = 1:numel(candIdx)
        iPS = candIdx(ii);

        baseFrac = ((iAR-1)*Sn2*AFn2*numel(candIdx) + (iS-1)*AFn2*numel(candIdx) + (iAF-1)*numel(candIdx) + ii) ...
                 /(ARn2*Sn2*AFn2*numel(candIdx));
        if mod(ii,3)==1
            waitbar(baseFrac, wb2, sprintf('Refine %.0f%% — AR=%.1f Span=%.2fm', 100*baseFrac, AR_set(iAR), span_set_m(iS)));
        end

        % Template
        template = struct(AirplaneClass);
        template.wing.span        = span_set_m(iS);
        template.wing.aspectRatio = AR_set(iAR);
        template.wing             = WingClass.SetWingData(template.wing, wingsLUT2, iAF, iAR, iS);

        template.powerSystem = PowerClass.SetPowerSystemData(template.powerSystem, MotorSpreadsheet, iPS);
        if ~isfinite(template.powerSystem.batteryCapacity) || template.powerSystem.batteryCapacity <= 0
            template.powerSystem.batteryCapacity = rules.battery_Wh_max;
        end
        template.powerSystem.batteryCapacity = min(rules.battery_Wh_max, template.powerSystem.batteryCapacity);

        template.fuselage  = FuselageClass.CalcFuselageData(template);
        template.fuselage  = FuselageClass.GenLandingGear(template);
        template.empennage = EmpennageClass.GenEmpennage(template, 4, 2);

        % Cheap constants for pruning
        rho = site.rho; g = 9.80665;
        S = template.wing.planformArea;
        CLmax = max(template.wing.clFlap, template.wing.clm);
        pitchSpeed = max(1, template.powerSystem.propSpeed);
        Tstatic_N  = 9.80665 * template.powerSystem.thrust;

        for p = p_set
          for c = c_set
            if rules.p_ge_3c && p < 3*c, continue; end

            mass_ducks = p * payload.massPerPassenger;
            mass_pucks = c * payload.massPerCargo;

            for L_in = L_set_in
              pl = template;

              pl.performance.passengers   = p;
              pl.performance.cargo        = c;
              pl.performance.declaredPmax = p;
              pl.performance.declaredCmax = c;

              pl.performance.bannerLength_in = L_in;
              pl.performance.bannerHeight_in = max(2, L_in/5);

              pl.performance.massPerPassenger = payload.massPerPassenger;
              pl.performance.massPerCargo     = payload.massPerCargo;
              pl.performance.x1_addon_kg      = payload.x1_addon_kg;

              pl.performance.peerMaxNetIncome   = norm.peerMaxNetIncome;
              pl.performance.peerMaxBanner      = norm.peerMaxBanner;
              pl.performance.bestGM_s           = norm.bestGM_s;
              pl.performance.yourGM_s           = norm.yourGM_s;
              pl.performance.proposalScore      = report.proposal;
              pl.performance.designReportScore  = report.design;
              pl.performance.participationLevel = participation;

              pl = FindTotalWeight(pl);

              % Pruning
              W2_kg = pl.performance.totalEmptyWeight + mass_ducks + mass_pucks + payload.x1_addon_kg;
              W2_N  = W2_kg * g;
              Vmin_lvl = sqrt( 2*W2_N / (rho*S*CLmax) );
              if Vmin_lvl > 0.95 * pitchSpeed, continue; end
              if (Tstatic_N / W2_N) < 0.24, continue; end

              % Performance
              pl = GenVelocityTest(pl, 2, site.rho, site.temp_K);
              if (pl.performance.velocity2/1.5 < pl.performance.landingSpeed2) || (pl.performance.velocity2==-1)
                  continue;
              end
              pl = GenVelocityTest(pl, 3, site.rho, site.temp_K);
              if (pl.performance.velocity3/1.5 < pl.performance.landingSpeed3) || (pl.performance.velocity3==-1)
                  continue;
              end

              if (pl.powerSystem.batteryCapacity > rules.battery_Wh_max)
                  pl.powerSystem.batteryCapacity = rules.battery_Wh_max;
              end
              if rules.p_ge_3c && (pl.performance.declaredPmax < 3*pl.performance.declaredCmax)
                  continue;
              end

              try
                  pl = Mission2Score(pl, site.wind_ms, mission.turnMult, site.rho, time.window_s, norm.peerMaxNetIncome);
              catch
                  pl = Mission2Score(pl, site.wind_ms, mission.turnMult, site.rho);
              end
              if pl.performance.time2 > time.window_s, continue; end

              try
                  pl = Mission3Score(pl, site.wind_ms, mission.turnMult, site.rho, time.window_s, norm.peerMaxBanner);
              catch
                  pl = Mission3Score(pl, site.wind_ms, mission.turnMult, site.rho);
              end
              try
                  pl = MissionGMScore(pl, norm.bestGM_s, norm.yourGM_s);
              catch
                  pl = MissionGMScore(pl);
              end
              try
                  pl = ComputeCompetitionScore2026(pl, report.proposal, report.design, participation, M1_success);
              catch
                  pl.performance.CompetitionScore = pl.performance.score2Normalized + pl.performance.score3Normalized + pl.performance.scoreGMNormalized;
              end

              k2 = k2 + 1;
              if k2 > numel(scores2)
                  scores2(end+5000,1) = -inf; %#ok<AGROW>
                  planes2(end+5000,1) = struct(AirplaneClass); %#ok<AGROW>
              end
              planes2(k2) = pl;
              scores2(k2) = pl.performance.CompetitionScore;

            end % L
          end % c
        end % p

      end % motor
    end % airfoil
  end % span
end % AR

elapsed2 = toc(t1);
try, close(wb2); end
fprintf('PASS 2 done in %.1f s. Kept %d designs.\n', elapsed2, k2);

planes2 = planes2(1:k2); scores2 = scores2(1:k2);
planes = [planes(:); planes2(:)];
scores = [scores(:); scores2(:)];
end % DO_REFINE

%%  Save / Rank / Plot 
if isempty(scores) || ~any(isfinite(scores))
    error('No feasible designs — adjust grids or pruning.');
end

[topScores, order] = sort(scores, 'descend');
K = min(100, numel(order));
topPlanes = planes(order(1:K)); %#ok<NASGU>
topScores = topScores(1:K);     %#ok<NASGU>
save('topPlanes_2026.mat','topPlanes','topScores');

fprintf('\nTop 5 scores:\n'); disp(topScores(1:min(5,K)));

% Optional CSV
if exist('SaveTopPlanesCSV.m','file')==2
    try
        SaveTopPlanesCSV(topPlanes, 'topPlanes_2026.csv');
    catch ME
        warning(ME.identifier, 'CSV export skipped: %s', ME.message);
    end
end

% Sensitivity plots for the best design
try
    best = topPlanes(1);
    fprintf('Generating sensitivity plots...\n');
    MakeSensitivityPlots(best, site, mission, time, norm, rules);
catch ME
    warning(ME.identifier, 'Sensitivity plotting skipped: %s', ME.message);
end

fprintf('\nAll done.\n');
