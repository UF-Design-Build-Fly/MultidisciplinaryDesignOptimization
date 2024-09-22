classdef FuselageClass

	properties
	
		frontalSurfaceArea = -1;
		height = -1;
		width = -1;
		length = -1;
		totalSA = -1;
		
		weight = -1; % weight of carbon fiber to make fuselage + mechanisms
		
		gearWeight = -1;
		gearSA = -1;
		gearFrontalSA = -1;
	
		wheelWidth = -1;
		wheelRadius = -1;
		wheelWeight = -1;
		wheelSA = -1;
		wheelFrontalSA = -1;
		
	end
	
	methods (Static)
	
		function fuselage = CalcFuselageData(plane) 
			fuselage = plane.fuselage;
			% Inputs: Aspect Ratio, Wing Area
			% Must be run after wing data is set
			
			in2m = 0.0254; % in to meters
			
			% Nose section up to bulkhead
			noseLength = (7+4)*in2m; %Space systems needs for their stuff (battery + motor)
			
			fuselageHeight = 4*in2m;
			fuselageWidth = 4*in2m;
			
			% Rectangularish middle section
			planeBayLength = 6*in2m;

			%totalFuselageLength = plane.wing.span*0.75;	%Using 75% rule for fuselage length
			totalFuselageLength = plane.wing.span*0.9;
			
			% Length behind passenger section
			tailLength = totalFuselageLength - noseLength - planeBayLength;
			if (tailLength < 0)
				tailLength = 0;
				totalFuselageLength = noseLength + planeBayLength;
			end
			
			% Calculations
			% This will be the surface area of the main portion of the fuselage, aka the middle part
			
			% Assume: cylinder shape A = 1/2 (a+b)*h
			% Vertical Narrowing is offset
			% Calculate the slant height
			l = sqrt(noseLength^2 + (1.5*in2m - 2*in2m)^2);
			% Calculate the lateral surface area
			surfaceAreaNose =  pi * (1.5*in2m + 2*in2m) * l;
			
			% Assume wedge shape		
			% Calculate the areas of the faces
			top = planeBayLength * fuselageWidth;				% Area of the rectangular base
			sideTriangle = 2 * 0.5 * planeBayLength * fuselageHeight;		% Area of the side triangles
			diagFace = sqrt(planeBayLength^2 + fuselageHeight^2) * fuselageWidth;					% Area of the diagonal face
			
			surfaceAreaPlaneBay = top + sideTriangle + diagFace;

			% Assume: boom shape for
			boomDiameter = 1*in2m;
			surfaceAreaTail = pi*boomDiameter*tailLength;

			% Combine surface areas
			totalSurfaceArea = surfaceAreaNose+surfaceAreaPlaneBay+surfaceAreaTail;
			
			% The volume in m^3, 1/16in is the thickness of the carbon fiber so this will help us find the weight of the fuselage
			carbonFiberVolume = totalSurfaceArea*(1/16*in2m);
			
			rhoCarbon = 1281;	% kg/m^3 1410 found in research, but previous plane was around 1281
			fuselageWeight = carbonFiberVolume*rhoCarbon; % Weight of fuselage in kg
			
			% Set the objects values
			fuselage.frontalSurfaceArea = fuselageWidth*fuselageHeight;
			fuselage.height = fuselageHeight;
			fuselage.width = fuselageWidth;
			fuselage.length = totalFuselageLength;
			fuselage.totalSA = totalSurfaceArea;
			fuselage.weight = fuselageWeight;

		end

		function fuselage = GenLandingGear(plane)
			fuselage = plane.fuselage;
			% This function computes different parameters for the landing gear
			% Inputs: PropDiam = diameter of propellor (inches)
			%		fuselageHeight = height of the fuselage (inches)
			% Outputs: A vector G with the following values:
			%   G = [gearWeight gearParaDrag gearSA wheelSA]
			% Height of landing gear is the difference between prop radius and fuselage height
			%
			%								top
			%		 ^			  --------------------
			%		 |			/ |					   \
			%		 |		   /  |						\
			%		 |		c /   b						 \
			%	height		 /	|						  \ 
			%		 |		/__a__|						   \
			%		 |	   |							   |	|
			%		 |	   |							   |	base
			%		 |	   |							   |	|
			%
			%
			%				  <------------width------------->
			
			%----------------------------Defined Constants----------------------------%
			in2m = 0.0254; % in to meters
			
			numWheels = 2;
			wheelWidth = 0.9063*in2m;
			wheelRadius = 1.5*in2m;
			gearwidth = 1.5*in2m;
			gearThickness = 1/8*in2m;
			base = 1.5*in2m;				%Height of vertical portion
			%rhoAluminum = 0.001122368;	%Density for landing gear material (Al) (slug/in^3)
			rhoCarbon = 1281;	% kg/m^3 1410 found in research, but previous plane was around 1281

			%----------------------Calculation of Gear Weight, SA---------------------%
			width = 0.25*plane.wing.span;
			height = plane.powerSystem.propDiameter*in2m/2 - 2*in2m + 3*in2m; %2" from fuselageHeight + 3" clearence 
			fuselageWidth = plane.fuselage.width;

			a = 0.5*(width-fuselageWidth);
			b = height-base;
			c = sqrt(a^2 + b^2);

			L = fuselageWidth + 2*(c + base);

			fuselage.gearSA = L * gearwidth;
			fuselage.gearFrontalSA = L * gearThickness;
			fuselage.gearWeight = rhoCarbon*(fuselage.gearSA * gearThickness);

			fuselage.wheelWeight = numWheels * 0.0241; % (kg) = 0.85ounces 
			fuselage.wheelSA = numWheels * (2*pi*(wheelRadius)^2 + 2*pi*wheelRadius*wheelWidth);
			fuselage.wheelFrontalSA = numWheels * (2*wheelRadius*wheelWidth);

			fuselage.wheelWidth = wheelWidth;
			fuselage.wheelRadius = wheelRadius;
			
		end

	end
end