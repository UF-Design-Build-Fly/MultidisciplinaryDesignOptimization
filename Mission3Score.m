function plane = Mission3Score(plane)

    turnAcceleration = 32*1.5;
    turnRadius = plane.performance.velocity2^2/turnAcceleration;
    %2 360 turns per lap
    lapDist = (500*4)+(2*2*pi*turnRadius);

    lapTime = lapDist/plane.performance.velocity3;    

    plane.performance.time3 = lapTime;
    plane.performance.numLaps3 = min(plane.powerSystem.time, 300)/lapTime;
    plane.performance.score3 = plane.performance.numLaps3 * plane.performance.numPassengers * plane.powerSystem.efficiency / plane.performance.drag3;

end