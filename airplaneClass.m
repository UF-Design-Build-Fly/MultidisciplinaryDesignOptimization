classdef airplaneClass
	
	properties
		wing = struct(wingClass);
		empennage = struct(empennageClass);
		powerSystem = struct(powerClass);
		fuselageClass = struct(fuselageClass); 
		performance = struct(performanceClass);
		sanityFlag = true;
	end
	
	methods
		%function obj = airplaneClass(AR, powerSystem, aifoil, sensorAR, sensorWeight, num_cargo) %TODO: define constructor in a way that makes sense
			%constructor should set all properties and call other constructors to make blank slate - no holdovers from the last object
		%end
	end	

end