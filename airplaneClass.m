classdef airplaneClass
	
	properties
		wing = struct(wingClass);
		empennage = struct(empennageClass);
		power = struct(powerClass);
		fuselage = struct(fuselageClass); 
		performance = struct(performanceClass);
        failureReason = "Not checked";
        sanityFlag = 1;
        volSanityFlag = 1;
        takeoffFail = 0;
        momentFail = 0;
        convergeFail = 0;
        spaceFail = 0;
        epFail = 0;
	end
	
	methods
		%function obj = airplaneClass(AR, powerSystem, aifoil, sensorAR, sensorWeight, num_cargo) %TODO: define constructor in a way that makes sense
			%constructor should set all properties and call other constructors to make blank slate - no holdovers from the last object
		%end
	end	

end