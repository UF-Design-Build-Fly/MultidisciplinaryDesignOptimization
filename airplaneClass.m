classdef airplane
	
	properties
		%I only called these "class" because we already had files with the same name
		wingClass 
		empennageClass
		powerSystemClass
		fuselageClass
		performanceClass
	end
	
	methods
		function obj = airplane(AR, powerSystem, aifoil, sensorAR, sensorWeight, num_cargo) %TODO: define constructor in a way that makes sense
			%constructor should set all properties and call other constructors to make blank slate - no holdovers from the last object!
		end
	end	

end