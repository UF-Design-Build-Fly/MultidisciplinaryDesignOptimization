classdef WingClass

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
		thickness = -1;
		%wings %matrix returned by wing.m. %See wingData.m for function to fill all wing parameters. See winger.m and wing_syntax.m to access wing data from main code.
		
	end

	methods (Static)

		function wing = SetWingData(wing, wings, airfoilIndex, aspectRatioIndex, spanIndex)

			%Set values from wings matrix into plane. See wingClass.m for property descriptions
			wing.clw = wings(airfoilIndex, 1, aspectRatioIndex, spanIndex);
			wing.clm = wings(airfoilIndex, 2, aspectRatioIndex, spanIndex);	 %cl max
			wing.cd = wings(airfoilIndex, 3, aspectRatioIndex, spanIndex);	 %cd i zero velocity coefficient of drag
			wing.clFlap = wings(airfoilIndex, 4, aspectRatioIndex, spanIndex);  %weight?
			wing.weight = wings(airfoilIndex, 5, aspectRatioIndex, spanIndex);
			wing.chord = wings(airfoilIndex, 6, aspectRatioIndex, spanIndex);
			wing.planformArea = wings(airfoilIndex, 7, aspectRatioIndex, spanIndex);
			wing.surfaceArea = wings(airfoilIndex, 8, aspectRatioIndex, spanIndex);
			wing.name = wings(airfoilIndex, 9, aspectRatioIndex, spanIndex);
			wing.thickness=wings(airfoilIndex, 10, aspectRatioIndex, spanIndex);

		end

		function maxLift = FindMaxLift(wingObj, airSpeed, rho)

			maxLift = 0.5 * rho * wingObj.clm * wingObj.planformArea * airSpeed^2;

		end

	end

end