classdef FuselageClass

    properties
    
        frontalSurfaceArea = -1;
        height = -1;
        width = -1;
        length = -1; %feet
        totalSA = -1;
        
        weight = -1; %weight of carbon fiber to make fuselage + mechanisms
        
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
            %Inputs: Aspect Ratio, Wing Area
            %Must be run after wing data is set
            
            %Nose section up to bulkhead
            noseLength = (6+1+2)/12; %Space systems needs for their stuff (battery + motor + fuse/speedcontroller)
            
            fuselageHeight = 4;
            fuselageWidth = 4;
            
            %Rectangularish middle section

            %totalFuselageLength = plane.wing.span*0.75;    %Using 75% rule for fuselage length
            totalFuselageLength = noseLength + 18/12;    %Using 75% rule for fuselage length
            
            %Length behind passenger section
            tailLength = totalFuselageLength - noseLength;
            if (tailLength < 1)
                tailLength = 1;
                totalFuselageLength = noseLength + passengerSection + tailLength;
            end
            
            %Calculations
            %This will be the surface area of the main portion of the fuselage, aka the middle part
            
            %Assume: trapezoidal shape for narrowing A = 1/2 (a+b)*h
            %Vertical Narrowing is offset
            motorMountSize = 3/12; %Narrow to 3x3 for motor mount
            topLength = sqrt((fuselageHeight - motorMountSize)^2 + (noseLength)^2);
            sideLength = sqrt(((fuselageWidth - motorMountSize)/2)^2 + (noseLength)^2);
            surfaceAreaNose = 0.5*(fuselageWidth + motorMountSize)*(topLength + noseLength) + ...
                              motorMountSize^2 + ...
                              (fuselageHeight + motorMountSize)*sideLength; %2 sides cancel the 1/2
            
            %Assume: trapezoidal shape for narrowing A = 1/2 (a+b)*h
            %Vertical Narrowing is offset
            endHeight = 3.5/12; %Narrow to 3.5x1.5in
            endWidth = 1.5/12;
            topLength = sqrt((fuselageHeight - endHeight)^2 + (tailLength)^2);
            sideLength = sqrt(((fuselageWidth - endWidth)/2)^2 + (tailLength)^2);
            surfaceAreaTailcone = 0.5*(fuselageWidth + endWidth)*(topLength + tailLength) + ...
                                  endHeight*endWidth + ...
                                  (fuselageHeight + endHeight)*sideLength; %2 sides cancel the 1/2

            %Combine surface areas
            totalSurfaceArea = surfaceAreaNose+surfaceAreaTailcone;
            
            %The volume in ft^3, 1/16 is the thickness of the carbon fiber in inches so this will help us find the weight of the fuselage
            carbonFiberVolume = totalSurfaceArea*(1/16/12);
            
            rhoCarbon = 80;	%lb/ft^3 88 found in research, but previous plane was around 80
            fuselageWeight = carbonFiberVolume*rhoCarbon; %Weight of fuselage in lb
            
            %Set the objects values
            fuselage.frontalSurfaceArea = surfaceAreaNose;
            fuselage.height = fuselageHeight;
            fuselage.width = fuselageWidth;
            fuselage.length = totalFuselageLength;
            fuselage.totalSA = totalSurfaceArea;
            fuselage.weight = fuselageWeight;

        end

        function fuselage = GenLandingGear(plane)
            fuselage = plane.fuselage;
            %This function computes different parameters for the landing gear
            %Inputs: PropDiam = diameter of propellor (inches)
            %        fuselageHeight = height of the fuselage (inches)
            %Outputs: A vector G with the following values:
            %   G = [gearWeight gearParaDrag gearSA wheelSA]
            %Height of landing gear is the difference between prop radius and fuselage height
            %
            %                                top
            %         ^              -------------------
            %         |            / |                   \
            %         |           /  |                    \
            %         |        c /   b                     \
            %    height         /    |                      \ 
            %         |        /__a__|                       \
            %         |       |                               |    |
            %         |       |                               |    base
            %         |       |                               |    |
            %
            %
            %                  <------------width------------->
            
            %----------------------------Defined Constants----------------------------%
            numWheels = 2;
            wheelWidth = 0.5/12;    %(division by 12 to convert to ft)
            wheelRadius = 1.25/12;   %^^^
            gearwidth = 2/12;       %(ft)
            gearThickness = 1/8/12;    %(ft)
            base = 2/12;            %Height of vertical portion (ft)
            %rhoAluminum = 0.001122368;   %Density for landing gear material (Al) (slug/in^3)
            rhoCarbon = 120.486;	%lb/ft^3

            %----------------------Calculation of Gear Weight, SA---------------------%
            width = 0.25*plane.wing.span;
            height = plane.powerSystem.propDiameter/24 + 1/6; %1.5" from fuselageHeight + 2.5" clearence 
            fuselageWidth = plane.fuselage.width;

            a = 0.5*(width-fuselageWidth);
            b = height-base;
            c = sqrt(a^2 + b^2);

            L = fuselageWidth + 2*(c + base);

            fuselage.gearSA = L * gearwidth;
            fuselage.gearFrontalSA = L * gearThickness;
            fuselage.gearWeight = rhoCarbon*(fuselage.gearSA * gearThickness);

            fuselage.wheelWeight = numWheels * 3.53/16; %3.53ounces (lb)
            fuselage.wheelSA = numWheels * (2*pi*(wheelRadius)^2 + 2*pi*wheelRadius*wheelWidth);
            fuselage.wheelFrontalSA = numWheels * (2*wheelRadius*wheelWidth);

            fuselage.wheelWidth = wheelWidth;
            fuselage.wheelRadius = wheelRadius;
            
        end

    end
end