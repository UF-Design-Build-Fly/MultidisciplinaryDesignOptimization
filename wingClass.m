classdef wing

	properties
		cl; %coefficient of lift
		clm; %cl max
		cd0; %zero velocity coefficient of drag
		weight;
		name;
		surfaceArea;
		chord;
		aspectRatio; %ratio between length and width of wing.
		%wings %matrix returned by wing.m. %DEBUG - update this to break out to list of properties that is easier to read
			  %maybe break wing.m into two functions. one does all the calcs then have a translator function to return properties
			  %of the wing given just an index as input
		
	end

	methods
    
	end

end