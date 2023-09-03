function [plane] = findTotalWeight(plane, Electronic_Package_Weight, Antenna_Length)
%weight in pounds (verify this with other functions)

%DEBUG - Previous runs of analysis seemed to 
%underestimate fuselage weight so just add a constant multiple here. 
%Look into fixing root cause.
    plane.fuselage.weight = plane.fuselage.weight*3; 
    plane.empennage.HSweight = plane.empennage.HSweight*1.3;
    plane.fuselage.gearWeight = plane.fuselage.gearWeight*4;
    plane.wing.weight = plane.wing.weight*1.15;
    
    plane.performance.totalEmptyWeight = plane.empennage.HSweight + ...
        plane.empennage.VSweight + plane.fuselage.weight + plane.fuselage.gearWeight + ...
        plane.power.weight + plane.wing.weight;
    
    plane.performance.totalWeight2 = plane.performance.totalEmptyWeight + ...
        Electronic_Package_Weight; %0.042 pounds per syringe
    
    plane.performance.totalWeight3 = plane.performance.totalEmptyWeight + ...
        (.16/6)*Antenna_Length; %.16 lbs/ft of pipe (*2 for counterweight)
end