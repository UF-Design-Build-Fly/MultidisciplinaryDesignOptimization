classdef wingClass

	properties
		cl = -1; %coefficient of lift
		clm = -1; %cl max
		cd0 = -1; %zero velocity coefficient of drag
		weight = -1;
		name = -1;
		surfaceArea = -1;
		chord = -1;
		aspectRatio = -1; %ratio between length and width of wing.
		%wings %matrix returned by wing.m. %DEBUG - update this to break out to list of properties that is easier to read
			  %maybe break wing.m into two functions. one does all the calcs then have a translator function to return properties
			  %of the wing given just an index as input
		
	end

	methods
    
	end

end