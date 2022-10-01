function [plane] = sanityCheck(plane)
%returns true if the aircraft makes sense and meets requirements
    %key point for this function is that it starts out true but functions
    %can only set it false. No function should set it back to true - if one
    %function thinks the plane shouldn't work then the flag should stay
    %false.
    plane.failureReason = "OK"; %plane will stay ok if it passes all the checks below
    if ((plane.performance.takeoffDist1 > 25) || (plane.performance.takeoffDist2 > 25) || (plane.performance.takeoffDist3 > 25))
        plane.sanityFlag = false;
        plane.failureReason = "TakeoffDist";
        plane.takeoffFail = 1;
    end
    %Fuselage length will never be too long this year, will be based on a flat multiple of wingspan
%     if (plane.fuselage.length > 7.8)
%          plane.sanityFlag = false;
%          plane.failureReason = "fuselage";
%          plane.fuselageFail = 1;
%     end
    if ((plane.performance.velocity2 == -1) || (plane.performance.velocity3  == -1))
        %disp("Error Mission 2 or 3 velocity solver did not converge");
        plane.failureReason = "noConverge";
        plane.sanityFlag = false;
        plane.convergeFail = 1;
    end
%I don't think this code is applicable this year, might need to instead remove planes which can't complete mission 3 in time 
% %  FIX THIS CODE LATER   (Update 8/2022 - I don't remember what the error was here -- make sure to double check as it may or may not work)
%     if( (plane.performance.time2/60) > plane.power.time) %Must have enough power to complete mission 2. 
%         plane.sanityFlag = false; %for mission 3 not having enough power time will not cause a mission failure, but it is a hard requirement here
%         plane.failureReason = "Power";
%         plane.powerFail = 1;
%     end
    
    %TODO: this should eventually check for -1 values in properties to make
    %sure everything was initialized correctly
    %for now the scoring functions take care of sanity checking the
    %velocity
end