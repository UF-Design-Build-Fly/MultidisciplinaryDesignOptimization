function plane = Mission2Score(plane, windSpeed, turnSpeedMultiplier)

	maxLift = WingClass.FindMaxLift(plane.wing, plane.performance.velocity2*turnSpeedMultiplier);
	turnAcceleration = 32*maxLift/plane.performance.totalWeight2;
	turnRadius = plane.performance.velocity2^2/turnAcceleration;
	
	%2 360 turns per lap
	turnTime = 2*2*pi*turnRadius/(plane.performance.velocity2*turnSpeedMultiplier);
	
	downWindTime = 1000/(plane.performance.velocity2+windSpeed);
	upWindTime = 1000/(plane.performance.velocity2-windSpeed);

	lapTime = turnTime + downWindTime + upWindTime;

	plane.performance.time2 = 3*lapTime;
	plane.performance.score2 = plane.performance.m2Weight/lapTime;

end