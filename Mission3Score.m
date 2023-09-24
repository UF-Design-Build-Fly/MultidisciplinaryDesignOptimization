function plane = Mission3Score(plane)

    turnRadius = plane.performance.velocity3^2/plane.performance.turnAcceleration;
    %2 360 turns per lap
    lapDist = (500*4)+(2*2*pi*turnRadius);

    lapTime = lapDist/plane.performance.velocity3;    

    plane.performance.time3 = lapTime;
    plane.performance.numLaps3 = floor(min(plane.powerSystem.time, 300)/lapTime);
    plane.performance.score3 = plane.performance.numLaps3 * plane.performance.numPassengers/plane.powerSystem.batteryCapacity;

end