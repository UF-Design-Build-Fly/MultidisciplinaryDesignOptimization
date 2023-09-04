function plane = Mission2Score(plane)

    turnRadius = plane.performance.velocity2^2/plane.performance.turnAcceleration;
    %2 360 turns per lap
    lapDist = (500*4)+(2*2*pi*turnRadius);

    lapTime = lapDist/plane.performance.velocity2;    
    plane.performance.time2 = 3*lapTime;
    plane.performance.score2 = plane.performance.m2Weight/lapTime;

end