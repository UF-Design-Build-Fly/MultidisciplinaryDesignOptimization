classdef FuselageClass

	properties
	
		frontalSurfaceArea = -1;
		height = -1;
		width = -1;
		length = -1; % (ft)
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
			
			% Nose section up to bulkhead
			noseLength = (7+4)/12; %Space systems needs for their stuff (battery + motor)
			
			fuselageHeight = 4/12;
			fuselageWidth = 4/12;
			
			% Rectangularish middle section
			planeBayLength = 6/12;

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
			l = sqrt(noseLength^2 + (1.5 - 2)^2);
			% Calculate the lateral surface area
			surfaceAreaNose =  pi * (1.5/12 + 2/12) * l;
			
			% Assume wedge shape		
			% Calculate the areas of the faces
			top = planeBayLength * fuselageWidth;				% Area of the rectangular base
			sideTriangle = 2 * 0.5 * planeBayLength * fuselageHeight;		% Area of the side triangles
			diagFace = sqrt(planeBayLength^2 + fuselageHeight^2) * fuselageWidth;					% Area of the diagonal face
			
			surfaceAreaPlaneBay = top + sideTriangle + diagFace;

			% Assume: boom shape for
			boomDiameter = 1/12;
			surfaceAreaTail = pi*boomDiameter*tailLength;

			% Combine surface areas
			totalSurfaceArea = surfaceAreaNose+surfaceAreaPlaneBay+surfaceAreaTail;
			
			% The volume in ft^3, 1/16 is the thickness of the carbon fiber in inches so this will help us find the weight of the fuselage
			carbonFiberVolume = totalSurfaceArea*(1/16/12);
			
			rhoCarbon = 80;	% lb/ft^3 88 found in research, but previous plane was around 80
			fuselageWeight = carbonFiberVolume*rhoCarbon; % Weight of fuselage in lb
			
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
			numWheels = 2;
			wheelWidth = 0.9063/12;		%(division by 12 to convert to ft)
			wheelRadius = 1.5/12;		%^^^
			gearwidth = 1.5/12;			%(ft)
			gearThickness = 1/8/12;		%(ft)
			base = 1.5/12;				%Height of vertical portion (ft)
			%rhoAluminum = 0.001122368;	%Density for landing gear material (Al) (slug/in^3)
			rhoCarbon = 120.486;		%lb/ft^3

			%----------------------Calculation of Gear Weight, SA---------------------%
			width = 0.25*plane.wing.span;
			height = plane.powerSystem.propDiameter/12/2 - 2/12 + 3/12; %2" from fuselageHeight + 3" clearence 
			fuselageWidth = plane.fuselage.width;

			a = 0.5*(width-fuselageWidth);
			b = height-base;
			c = sqrt(a^2 + b^2);

			L = fuselageWidth + 2*(c + base);

			fuselage.gearSA = L * gearwidth;
			fuselage.gearFrontalSA = L * gearThickness;
			fuselage.gearWeight = rhoCarbon*(fuselage.gearSA * gearThickness);

			fuselage.wheelWeight = numWheels * 0.85/16; %0.85ounces (lb)
			fuselage.wheelSA = numWheels * (2*pi*(wheelRadius)^2 + 2*pi*wheelRadius*wheelWidth);
			fuselage.wheelFrontalSA = numWheels * (2*wheelRadius*wheelWidth);

			fuselage.wheelWidth = wheelWidth;
			fuselage.wheelRadius = wheelRadius;
			
		end

	end
end