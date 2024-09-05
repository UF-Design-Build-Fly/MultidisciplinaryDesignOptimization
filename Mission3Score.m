function plane = Mission3Score(plane)

    %avgTurnG, turnVelMultiplier, windVel are defined in main
    global avgTurnG;
    global turnVelMultiplier;
    global windVel;

    turnAcceleration = 32*avgTurnG;
    turnRadius = (plane.performance.velocity3*turnVelMultiplier)^2/turnAcceleration;
    %2 360 turns per lap
    turnDist = 2*2*pi*turnRadius;
    turnTime = turnDist/(plane.performance.velocity3*turnVelMultiplier);
    
    downWindDist = 500*2;
    downWindTime = downWindDist/(plane.performance.velocity3+windVel);
    upWindDist = 500*2;
    upWindTime = upWindDist/(plane.performance.velocity3-windVel);

    lapTime = turnTime + downWindTime + upWindTime;
    
    plane.performance.time3 = lapTime;
    plane.performance.numLaps3 = min(plane.powerSystem.time, 300)/lapTime;
    plane.performance.score3 = plane.performance.numLaps3 + 2.5;

end