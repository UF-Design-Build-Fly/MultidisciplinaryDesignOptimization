classdef FuselageClass

    properties
    
    frontalSurfaceArea = -1;
    height = -1;
    width = -1;
    length = -1; %feet
    totalSA = -1;
    
    weight = -1; %weight of carbon fiber to make fuselage + mechanisms
    
    gearWeight = -1;
    gearParaDrag = -1;
    gearSA = -1;

    wheelSA = -1;
    
    end
    
    methods
    
        function fuselage = CalcFuselageData(plane) 
            %Inputs: Aspect Ratio, Wing Area
            %Must be run after wing data is set
            
            %Nose section up to bulkhead
            noseLength = (6+1+2)/12; %Space systems needs for their stuff (battery + motor + fuse/speedcontroller)
            
            %Dimensions of m2 package (ft)
            m2PackageHeight = 3.5/12;
            m2PackageWidth = 3/12;
            m2PackageLength = 10/12; %3in for med cabinet,  5.5in for patient, 0.5in spacing
            
            spaceAroundPackage = 0.5/12;   %Desired room around package (ft)
            
            fuselageHeight = m2PackageHeight+spaceAroundPackage*2;    %Space for package plus spacing
            fuselageWidth = m2PackageWidth+spaceAroundPackage*2+1/12;      %Space for package plus spacing + 1 in
            
            %Mission 3 min size for passengers
            %3 passengers per row, 1in passengers 0.5in spaceing around them
            m3PassengerLength = (plane.performance.numPassengers/3)*1.5/12 + 0.5/12;

            %Rectangularish middle section
            passengerSection = max(m2PackageLength, m3PassengerLength);

            totalFuselageLength = plane.wing.span*0.75;    %Using 75% rule for fuselage length
            
            %Length behind passenger section
            tailLength = totalFuselageLength - noseLength - passengerSection;
            if (tailLength < 1)
                tailLength = 1;
                totalFuselageLength = noseLength + passengerSection + tailLength;
            end
            
            %Calculations
            %This will be the surface area of the main portion of the fuselage, aka the middle part
            surfaceAreaMain = 2*passengerSection*(fuselageHeight + fuselageWidth);
            
            %Assume: trapezoidal shape for narrowing A = 1/2 (a+b)*h
            %Vertical Narrowing is offset
            motorMountSize = 3/12; %Narrow to 3x3 for motor mount
            topLength = sqrt((fuselageHeight - motorMountSize)^2 + (noseLength)^2);
            sideLength = sqrt(((fuselageWidth - motorMountSize)/2)^2 + (noseLength)^2);
            surfaceAreaNose = 0.5*(fuselageWidth + motorMountSize)*(topLength + noseLength) + ...
                              motorMountHeight^2 + ...
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
            totalSurfaceArea = surfaceAreaMain+surfaceAreaNose+surfaceAreaTailcone;
            
            %The volume in ft^3, 1/16 is the thickness of the carbon fiber in inches so this will help us find the weight of the fuselage
            carbonFiberVolume = totalSurfaceArea*(1/16/12);
            
            rhoCarbon = 120.486;	%lb/ft^3
            fuselageWeight = carbonFiberVolume*rhoCarbon; %Weight of fuselage in lb
            
            %Set the objects values
            fuselage.frontalSurfaceArea = surfaceAreaNose;
            fuselage.height = fuselageHeight;
            fuselage.width = fuselageWidth;
            fuselage.length = totalFuselageLength;
            fuselage.totalSA = totalSurfaceArea;
            fuselage.weight = fuselageWeight;

        end

    end

end