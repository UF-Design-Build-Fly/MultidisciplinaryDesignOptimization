function plane = Mission2Score(plane)

    %avgTurnG, turnVelMultiplier, windVel are defined in main

    turnAcceleration = 32*avgTurnG;
    turnRadius = (plane.performance.velocity2*turnVelMultiplier)^2/turnAcceleration;
    %2 360 turns per lap
    turnDist = 2*2*pi*turnRadius;
    turnTime = turnDist/(plane.performance.velocity2*turnVelMultiplier);
    
    downWindDist = 500*2;
    downWindTime = downWindDist/(plane.performance.velocity2+windVel);
    upWindDist = 500*2;
    upWindTime = upWindDist/(plane.performance.velocity2-windVel);

    lapTime = turnTime + downWindTime + upWindTime;

    plane.performance.time2 = 3*lapTime;
    plane.performance.score2 = plane.performance.m2Weight/lapTime;

end