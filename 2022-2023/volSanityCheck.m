function [plane] = volSanityCheck(plane, Electronic_Package_Weight)
%returns true if the aircraft makes sense and meets requirements
    %key point for this function is that it starts out true but functions
    %can only set it false. No function should set it back to true - if one
    %function thinks the plane shouldn't work then the flag should stay
    %false.
    plane.failureReason= "volOK";
    if Electronic_Package_Weight <= plane.performance.totalWeight2*.3
        plane.volSanityFlag = false;
        plane.failureReason = "EP not heavy enough";
        plane.epFail = 1;
    end
    if (plane.fuselage.length*12) + 3.7 + (plane.fuselage.width*12) + (plane.wing.thickness*7*12) + (plane.empennage.VSchord*12) >= 62
            plane.volSanityFlag = false;
            plane.failureReason = "Too big";
            plane.spaceFail = 1;
    end
%     if plane.performance.antennaLength > plane.fuselage.length
%         if ((plane.performance.antennaLength) + 3.7) + (plane.fuselage.width + (plane.wing.thickness*7*12)) + (plane.empennage.VSchord*12) >= 62
%             plane.volSanityFlag = false;
%             plane.failureReason = "Too big";
%             plane.spaceFail = 1;
%         end    
%     elseif ((plane.fuselage.length*12) + 3.7) + (plane.fuselage.width + (plane.wing.thickness*7*12)) + (plane.empennage.VSchord*12) >= 62
%             plane.volSanityFlag = false;
%             plane.failureReason = "Too big";
%             plane.spaceFail = 1;
%     end