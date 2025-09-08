function MakeSensitivityPlots(best, site, mission, time, norm, rules)
% MakeSensitivityPlots — M2/M3 sensitivity figures for a saved plane.

PI  = 3.141592653589793;      % robust against shadowed built-ins
EPS = 2.220446049250313e-16;

% defaults
if ~isfield(site,'rho'),        site.rho = 1.225; end
if ~isfield(site,'wind_ms'),    site.wind_ms = 0; end
if ~isfield(mission,'turnMult'), mission.turnMult = 0.8; end
if ~isfield(time,'window_s'),   time.window_s = 300; end
if nargin < 6 || isempty(rules)
    rules.span_ft_min = 3; rules.span_ft_max = 5; rules.banner_min_in = 10;
end
if ~isfield(norm,'peerMaxNetIncome'), norm.peerMaxNetIncome = 500; end
if ~isfield(norm,'peerMaxBanner'),    norm.peerMaxBanner    = 80;  end

rho = site.rho;
V2  = best.performance.velocity2;
V3  = best.performance.velocity3;

% Lap-time model
S     = best.wing.planformArea;
CLmax = max(best.wing.clFlap, best.wing.clm);
maxLift = @(V) 0.5*rho*(V*mission.turnMult).^2 * S * CLmax;
turnAcc = @(V,WN) max(maxLift(V)./WN, 1e-6);
turnRad = @(V,WN) V.^2 ./ turnAcc(V,WN);
turnTime= @(V,WN) 4*PI*turnRad(V,WN)./(V*mission.turnMult);
down    = @(V) 304.8./(V + site.wind_ms);
up      = @(V) 304.8./max(V - site.wind_ms, 0.1);
lapTime = @(V,WN) turnTime(V,WN) + down(V) + up(V);

g  = 9.80665;
WN2 = best.performance.totalWeight2 * g;
WN3 = best.performance.totalWeight3 * g;

% ---------------- M2 sensitivity ----------------
p   = best.performance.passengers; 
c   = best.performance.cargo;
Wh  = min(100, max(0, best.powerSystem.batteryCapacity));
EF  = Wh/100;

m2FromLaps = @(L,p,c,EF) 1 + ((p*(6+2*L)+c*(10+8*L)) - L*(10+0.5*p+2*c)*EF) / max(norm.peerMaxNetIncome,EPS);

N    = 800;

% 1) speed
Vvec  = linspace(0.7*V2, 1.3*V2, N);
Lvec  = floor(time.window_s ./ lapTime(Vvec, WN2));
M2_sp = arrayfun(@(L) m2FromLaps(L,p,c,EF), Lvec);

% 2) empty mass (changes WN2 only)
m0    = best.performance.totalEmptyWeight;
mVec  = linspace(0.7*m0, 1.3*m0, N);
WN2v  = (mVec + (best.performance.totalWeight2 - m0)) * g;
LE    = floor(time.window_s ./ lapTime(V2, WN2v));
M2_em = arrayfun(@(L) m2FromLaps(L,p,c,EF), LE);

% 3) battery Wh (EF)
WhVec = linspace(max(0,Wh*0.7), 100, N);
EFvec = WhVec/100;
baseL2 = floor(time.window_s / lapTime(V2, WN2));
M2_bw  = 1 + ((p*(6+2*baseL2)+c*(10+8*baseL2)) - baseL2*(10+0.5*p+2*c).*EFvec) / max(norm.peerMaxNetIncome,EPS);

% 4) passengers (discrete)
pVec  = unique(max(3*c, p-3):(p+3));
M2_pa = zeros(size(pVec));
for i=1:numel(pVec)
    pi = pVec(i);
    M2_pa(i) = 1 + ((pi*(6+2*baseL2)+c*(10+8*baseL2)) - baseL2*(10+0.5*pi+2*c)*EF)/max(norm.peerMaxNetIncome,EPS);
end

M2_base = m2FromLaps(baseL2,p,c,EF);
pct = @(x,x0) 100*(x - x0) ./ max(abs(x0), 1e-9);

figure; hold on; grid on;
plot(100*(Vvec - V2)/V2,                 pct(M2_sp,M2_base), 'DisplayName','Plane Speed');            %;
plot(100*(mVec - m0)/m0,                 pct(M2_em,M2_base), 'DisplayName','Empty Mass');             %;
plot(100*(WhVec - max(Wh,1))/max(Wh,1),  pct(M2_bw,M2_base), 'DisplayName','Battery Wh');             %;
plot(100*(pVec - max(p,1))/max(p,1),     pct(M2_pa,M2_base), 'o-','DisplayName','Passengers');        %;
xlabel('% Change in Parameter'); ylabel('% Change in M2 Score');
title('M2 Score Sensitivity (2026 rules)'); legend('Location','best'); xlim([-30 30]);

% ---------------- M3 sensitivity ----------------
WS_in      = best.wing.span / 0.0254;
WS_in_meas = round(WS_in);
WS_ft_rec  = WS_in_meas/12;
RAC        = max(0.90, 0.05*WS_ft_rec + 0.75);

L_in   = best.performance.bannerLength_in;
m3From = @(Laps, RACx, Lin) 2 + (Laps .* (Lin ./ RACx)) / max(norm.peerMaxBanner,EPS);

% 1) speed
Vvec3 = linspace(0.7*V3, 1.3*V3, N);
Laps3 = floor(time.window_s ./ lapTime(Vvec3, WN3));
M3_sp = m3From(Laps3, RAC, L_in);

% 2) wingspan -> RAC (rounded inch)
WSftVec = min(rules.span_ft_max, max(rules.span_ft_min, linspace(0.7*WS_ft_rec, 1.3*WS_ft_rec, N)));
WSinVec = round(WSftVec*12);
RACvec  = max(0.90, 0.05*(WSinVec/12) + 0.75);
baseL3  = floor(time.window_s / lapTime(V3, WN3));
M3_ws   = m3From(baseL3*ones(size(RACvec)), RACvec, L_in);

% 3) banner length
LvecBL  = linspace(max(rules.banner_min_in, 0.7*L_in), 1.3*L_in, N);
M3_bl   = m3From(baseL3*ones(size(LvecBL)), RAC, LvecBL);

M3_base = m3From(baseL3, RAC, L_in);

figure; hold on; grid on;
plot(100*(Vvec3 - V3)/V3,                         pct(M3_sp,M3_base), 'DisplayName','Plane Speed');           %;
plot(100*(WSftVec - max(WS_ft_rec,1))/max(WS_ft_rec,1), pct(M3_ws,M3_base), 'DisplayName','Wingspan (via RAC)'); %;
plot(100*(LvecBL - max(L_in,1))/max(L_in,1),      pct(M3_bl,M3_base), 'DisplayName','Banner Length');         %;
xlabel('% Change in Parameter'); ylabel('% Change in M3 Score');
title('M3 Score Sensitivity (2026 rules)'); legend('Location','best'); xlim([-30 30]);
end
