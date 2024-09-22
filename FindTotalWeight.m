function plane = FindTotalWeight(plane)
	
	%DEBUG - Previous runs of analysis seemed to 
	%underestimate fuselage weight so just add a constant multiple here. 
	%Look into fixing root cause.
	
	plane.performance.totalEmptyWeight = plane.empennage.HSweight + ...
		plane.empennage.VSweight + plane.fuselage.weight + plane.fuselage.gearWeight + ...
		plane.fuselage.wheelWeight + plane.powerSystem.weight + plane.wing.weight;
	
	plane.performance.totalWeight2 = plane.performance.totalEmptyWeight + ...
		plane.performance.m2Weight + 0.249; % 0.55lb for X1 Test Vehicle
	
	plane.performance.totalWeight3 = plane.performance.totalEmptyWeight + 0.249; % 0.55lb for X1 Test Vehicle
end