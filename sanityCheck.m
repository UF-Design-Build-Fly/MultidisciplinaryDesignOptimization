function [plane] = sanityCheck(plane)
%returns true if the aircraft makes sense and meets requirements
    %key point for this function is that it starts out true but functions
    %can only set it false. No function should set it back to true - if one
    %function thinks the plane shouldn't work then the flag should stay
    %false.
    if ((plane.performance.takeoffDist1 > 25) || (plane.performance.takeoffDist2 > 25) || (plane.performance.takeoffDist3 > 25))
        plane.sanityFlag = false;
    end
    if (plane.fuselage.length > 7.9)
         plane.sanityFlag = false;
    end
    if plane.performance.velocity2 == -1 || plane.performance.velocity3  == -1
        disp("Error Mission 2 or 3 velocity not set"); %DEBUG - change to disp
        plane.sanityFlag = false;
    end
    if(plane.performance.time2 < plane.power.time) %must complete 3 laps in 5 minutes
        plane.sanityFlag = false; %for mission 3 not having enough power time will not cause a mission failure, but it is a hard requirement here
    end
    
    %TODO: this should eventually check for -1 values in properties to make
    %sure everything was initialized correctly
    %for now the scoring functions take care of sanity checking the
    %velocity
end