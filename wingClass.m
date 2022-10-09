classdef wingClass

	properties
        span = -1;
		clw = -1; %coefficient of lift for the wing
		clm = -1; %cl max
		cd = -1; %coefficient of drag for the wing
		weight = -1;
		name = -1;
		surfaceArea = -1;
		chord = -1;
		aspectRatio = -1; %ratio between length and width of wing.
		planformArea = -1;
		clFlap = -1;
		%wings %matrix returned by wing.m. %See wingData.m for function to fill all wing parameters. See winger.m and wing_syntax.m to access wing data from main code.
		
	end

	methods
    
	end

end