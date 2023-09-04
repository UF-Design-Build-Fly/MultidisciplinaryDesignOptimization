function plane = FindTotalWeight(plane)
    %weight in pounds (verify this with other functions)
    
    %DEBUG - Previous runs of analysis seemed to 
    %underestimate fuselage weight so just add a constant multiple here. 
    %Look into fixing root cause.
    
    plane.performance.totalEmptyWeight = plane.empennage.HSweight + ...
        plane.empennage.VSweight + plane.fuselage.weight + plane.fuselage.gearWeight + ...
        plane.fuselage.wheelWeight + plane.powerSystem.weight + plane.wing.weight;
    
    plane.performance.totalWeight2 = plane.performance.totalEmptyWeight + ...
        plane.performance.m2weight + 6*3/16; %(lbs) 3oz*4crew+stretcher GET ACTUAL VALUES
    
    plane.performance.totalWeight3 = plane.performance.totalEmptyWeight + ...
        (3/16)*(2+plane.performance.numPassengers); %(lbs) 3oz per for crew GET ACTUAL VALUES
end