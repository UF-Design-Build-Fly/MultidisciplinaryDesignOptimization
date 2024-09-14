function plane = FindTotalWeight(plane)
	%weight in pounds (verify this with other functions)
	
	%DEBUG - Previous runs of analysis seemed to 
	%underestimate fuselage weight so just add a constant multiple here. 
	%Look into fixing root cause.
	
	plane.performance.totalEmptyWeight = plane.empennage.HSweight + ...
		plane.empennage.VSweight + plane.fuselage.weight + plane.fuselage.gearWeight + ...
		plane.fuselage.wheelWeight + plane.powerSystem.weight + plane.wing.weight;
	
	plane.performance.totalWeight2 = plane.performance.totalEmptyWeight + ...
		plane.performance.m2Weight + 0.5; % 0.5lb for X1 Test Vehicle
	
	plane.performance.totalWeight3 = plane.performance.totalEmptyWeight;
end