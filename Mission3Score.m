function plane = Mission3Score(plane, windSpeed, turnSpeedMultiplier)

	maxLift = WingClass.FindMaxLift(plane.wing, plane.performance.velocity3*turnSpeedMultiplier);
	turnAcceleration = 32*maxLift/plane.performance.totalWeight2;
	turnRadius = plane.performance.velocity2^2/turnAcceleration;
	
	%2 360 turns per lap
	turnTime = 2*2*pi*turnRadius/(plane.performance.velocity2*turnSpeedMultiplier);
	
	downWindTime = 1000/(plane.performance.velocity2+windSpeed);
	upWindTime = 1000/(plane.performance.velocity2-windSpeed);

	lapTime = turnTime + downWindTime + upWindTime;
	
	plane.performance.time3 = lapTime;
	plane.performance.numLaps3 = min(plane.powerSystem.time, 300)/lapTime;
	plane.performance.score3 = plane.performance.numLaps3 + 2.5;

end