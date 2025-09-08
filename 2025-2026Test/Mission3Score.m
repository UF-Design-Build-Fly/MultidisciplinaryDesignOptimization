function plane = Mission3Score(plane, windSpeed, turnMult, rho, window_s, peerMaxBanner)
    V = plane.performance.velocity3;
    if V<=0 || ~isfinite(V), plane.performance.M3=0; plane.performance.score3=0; return; end

    % Lap time with turning at turnMult
    maxLift = WingClass.FindMaxLift(plane.wing, V*turnMult, rho);
    turnAcc = maxLift / max(1e-6, plane.performance.totalWeight3*9.80665);
    turnRad = V^2 / max(1e-6, turnAcc);
    turnTime = 2*2*pi*turnRad / (V*turnMult);
    lapTime = turnTime + 304.8/(V+windSpeed) + 304.8/max(V-windSpeed,0.1);

    laps = floor(window_s / lapTime);
    plane.performance.time3 = laps*lapTime; plane.performance.laps3=laps; plane.performance.numLaps3=laps;
    if laps<=0, plane.performance.M3=0; plane.performance.score3=0; return; end

    % RAC from rounded span (nearest inch), RAC >= 0.90
    WS_in = plane.wing.span / 0.0254;
    WS_in_meas = round(WS_in);
    WS_ft_rec  = WS_in_meas/12;
    RAC = max(0.90, 0.05*WS_ft_rec + 0.75);

    % Enforce banner constraints (H >= 2", L/H <= 5)
    L_in = plane.performance.bannerLength_in;
    H_in = plane.performance.bannerHeight_in;
    H_in = max(H_in, max(2, L_in/5));
    plane.performance.bannerHeight_in = H_in; % keep consistent

    metric = laps * (L_in / RAC);
    plane.performance.score3 = metric;
    plane.performance.M3 = 2 + metric/max(peerMaxBanner,eps);
end
