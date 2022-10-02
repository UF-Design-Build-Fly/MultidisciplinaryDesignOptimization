function [plane] = volSanityCheck(plane, Electronic_Package_Weight)
%returns true if the aircraft makes sense and meets requirements
    %key point for this function is that it starts out true but functions
    %can only set it false. No function should set it back to true - if one
    %function thinks the plane shouldn't work then the flag should stay
    %false.
    if Electronic_Package_Weight <= plane.performance.totalweight2*.3
        plane.volSanityFlag = false;
    end
    if ((plane.fuselage.length*12) + 3.7) + (plane.fuselage.width + (plane.wing.thickness*7*12)) + (plane.empennage.VSchord*12) >= 62
        plane.volSanityFlag = false;
    end