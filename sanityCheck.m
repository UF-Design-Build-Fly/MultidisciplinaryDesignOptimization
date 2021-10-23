function [check] = sanityCheck(plane)
%returns true if the aircraft makes sense and meets requirements
    %lap distance is ~2160 feet
    %mission 1 velocity must be sufficient to complete 3 laps in 5 minutes
    %i.e, minimum velocity of 21.6 feet per second
    %use the mission two velocity for this check to save calls to velocity solver
    %
    
    if (plane.performance.takeoffDist2 > 25) || (plane.performance.takeoffDist1 > 25)
        sanityCheck = false;
    end
    
    if (plane.fuselage.length > 7.9)
        sanityCheck = false;
    end
end