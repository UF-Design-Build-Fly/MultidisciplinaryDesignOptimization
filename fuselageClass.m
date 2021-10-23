classdef fuselageClass

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

    wheelWidth = 0.5;    %width of wheels (in)
    wheelRadius = 1.5; %wheel radius
    wheelSA = 2*pi*(radius_wheel)^2+ pi*2*radius_wheel*width_wheel)/144;

    fuseDrag1 = -1;
    fuseDrag2 = -1;
    fuseDrag3 = -1;

    numVials = -1;
    numSyringes = -1;

  end

 

end