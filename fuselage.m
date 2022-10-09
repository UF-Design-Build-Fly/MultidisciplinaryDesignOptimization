function [plane] = fuselage(plane) 
%Inputs: Aspect Ratio, Wing Area

%Dimensions of electronic package (in)
electronicPackageHeight = 3;
electronicPackageWidth = 3;
electronicPackageLength = 6;

avionicsLength=8; %space systems needs for their stuff

fuselageHeight=0;	%fuse height (in)
fuselageWidth=0;	%fuse width (in)
fuselageLength=0;	%fuse length (in)

rhoCarbon=120.486;	%lb/ft^3
spaceAroundPackage=1;   %Desired room around package (in)

fuselageHeight=electronicPackageHeight+spaceAroundPackage*2;    %Space for package plus 2 inches
fuselageWidth=electronicPackageWidth+spaceAroundPackage*2;      %Space for package plus 2 inches

%calculates wigspan using AR and wing area
%aspectRatio = plane.wing.aspectRatio;
%planformArea = plane.wing.planformArea;
wingspan = plane.wing.span;%sqrt(aspectRatio*planformArea);

totalFuselageLength = wingspan*0.75*12;    %Using 75% rule for fuselage length

minMiddleFuselageLength=(electronicPackageLength+avionicsLength)*1.1;   %Room for package and avionics plus 10% (Doesn't include tail or nose)
noseLength=3;	%assumption for the front length of the fuselage.
tailconeLength=3;	%assumption for the back length of the fueslage

middleFuselageLength = totalFuselageLength-noseLength-tailconeLength;   %Middle length is total-tail-nose

%Checks if the 75% rule makes a big enough fuselage
if middleFuselageLength < minMiddleFuselageLength
    middleFuselageLength = minMiddleFuselageLength;
end

fuselageLength=middleFuselageLength+noseLength+tailconeLength; %final length of fuselage

surfaceAreaMain=(2*middleFuselageLength*fuselageHeight)+(2*middleFuselageLength*fuselageWidth); %this will be the surface area of the main portion of the fuselage, aka the middle part
trapezoidDiagonalLengthNose=sqrt((0.5*fuselageHeight)^2+(noseLength)^2); %we will assume a trapezoidal shape for the front of the fuselage if you view it from the side  
surfaceAreaNose=(trapezoidDiagonalLengthNose*fuselageWidth)+(6*fuselageWidth)+((fuselageHeight/2)*fuselageWidth)+2*((0.75*fuselageHeight*6)); %this formula is used to the find the surface area of the front section. So it is bottom area + top area + front area(closed section which is the nose) + side area, which are trapezoids
triangularDiagonalLengthTailcone=sqrt((fuselageHeight)^2+(tailconeLength)^2); %we will assume a triangular shape for the back of the fuselage if you view it from the side. 
surfaceAreaTailcone=(fuselageWidth*tailconeLength)+(fuselageWidth*triangularDiagonalLengthTailcone)+2*((fuselageHeight*tailconeLength)/2); %this formula is used to find the surface area of the back portion of the fueslage. Same method as front except now we have triangular surface area instead of trapezoidal.
totalSurfaceArea=surfaceAreaMain+surfaceAreaNose+surfaceAreaTailcone; %this will be the temporary surface area of the fuselage total in inches

%Switch to square feet
totalSurfaceArea=totalSurfaceArea/144;%Surface Area fuse Total,ft^2
surfaceAreaNose=(fuselageWidth*fuselageHeight)/144;%Frontal Surface Area,ft^2
carbonFiberVolume=totalSurfaceArea*((1/16)/144);%The volume in ft^3, 1/16 is the thickness of the carbon fiber in inches so this will help us find the weight of the fuselage

fuselageWeight=carbonFiberVolume*rhoCarbon; %Weight of fuselage in lb
fuselageHeight=fuselageHeight/12; %convert to feet
fuselageWidth=fuselageWidth/12; %convert to feet
fuselageLength=fuselageLength/12; %convert to feet

%convert these to object orientation
plane.fuselage.frontalSurfaceArea=surfaceAreaNose;
plane.fuselage.height=fuselageHeight;
plane.fuselage.width=fuselageWidth;
plane.fuselage.length=fuselageLength;
plane.fuselage.totalSA=totalSurfaceArea;
plane.fuselage.weight=fuselageWeight;
end