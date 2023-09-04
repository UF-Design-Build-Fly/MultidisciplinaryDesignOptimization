classdef AirplaneClass
	
	properties
		wing = struct(WingClass);
		empennage = struct(EmpennageClass);
		powerSystem = struct(PowerClass);
		fuselage = struct(FuselageClass); 
		performance = struct(PerformanceClass);
        failureReason = "Not checked";
        sanityFlag = 1;
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