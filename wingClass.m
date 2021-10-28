classdef wingClass

	properties
		clw = -1; %coefficient of lift for the wing
		clm = -1; %cl max
		cd = -1; %coefficient of drag for the wing
		weight = -1;
		name = -1;
		surfaceArea = -1;
		chord = -1;
        span = 8; %max out length for competition rules this year.
		aspectRatio = -1; %ratio between length and width of wing.
		planformArea = -1;
		clFlap = -1;
		%wings %matrix returned by wing.m. %DEBUG - update this to break out to list of properties that is easier to read
			  %maybe break wing.m into two functions. one does all the calcs then have a translator function to return properties
			  %of the wing given just an index as input
		
	end

	methods
    
	end

end