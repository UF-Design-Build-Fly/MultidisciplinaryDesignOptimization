function plane = Mission2Score(plane, windSpeed, turnMult, rho, window_s, peerMaxNI)
    V = plane.performance.velocity2;
    if V<=0 || ~isfinite(V), plane.performance.M2=0; plane.performance.score2=0; return; end

    maxLift = WingClass.FindMaxLift(plane.wing, V*turnMult, rho);
    turnAcc = maxLift / max(1e-6, plane.performance.totalWeight2*9.80665);
    turnRad = V^2 / max(1e-6, turnAcc);
    turnTime = 2*2*pi*turnRad / (V*turnMult);
    lapTime = turnTime + 304.8/(V+windSpeed) + 304.8/max(V-windSpeed,0.1);

    laps = floor(window_s / lapTime);
    plane.performance.time2 = laps*lapTime; plane.performance.laps2=laps; plane.performance.numLaps2=laps;
    if laps<=0, plane.performance.M2=0; plane.performance.score2=0; return; end

    p=plane.performance.passengers; c=plane.performance.cargo;
    Wh = max(0, min(100, plane.powerSystem.batteryCapacity));
    EF = Wh/100;
    income = p*(6+2*laps) + c*(10+8*laps);
    cost   = laps*(10 + 0.5*p + 2*c) * EF;
    netIncome = income - cost;

    plane.performance.score2 = netIncome;
    plane.performance.M2 = 1 + netIncome/max(peerMaxNI,eps);
end
