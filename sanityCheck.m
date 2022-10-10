function [plane] = sanityCheck(plane, span, Antenna_Length)
%returns true if the aircraft makes sense and meets requirements
    %key point for this function is that it starts out true but functions
    %can only set it false. No function should set it back to true - if one
    %function thinks the plane shouldn't work then the flag should stay
    %false.
    plane.failureReason = "OK"; %plane will stay ok if it passes all the checks below
    if ((plane.performance.takeoffDist1 > 60) || (plane.performance.takeoffDist2 > 60) || (plane.performance.takeoffDist3 > 60))
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
    %This code will check if the wing is over-torqued, numbers are based on mastercarr spar, assumed that both spars have similar MIz
    %Units are in lb, ft, s, Pa
    %rho=0.763; antennaWidth=0.5/12; V=plane.performance.velocity3; Cdantenna=1.2; g=32.2;
    sparOuterDiam=.393/12; sparInnerDiam=.313/12; Ecf=70000000000;
    %dragAntenna = .5*(Antenna_Length/12)*antennaWidth*rho*(V^2)*Cdantenna*(1/g); %total drag on antenna (lb.)
    momentAntenna = .5*(Antenna_Length/12)*plane.performance.antDrag; %Wingtip moment created by antenna
    MI = pi()*((sparOuterDiam)^4-(sparInnerDiam)^4)*(1/64); %Mass Moment of Inertia of the spar
    reactionSpar = momentAntenna/(2*(plane.wing.chord*(.75/2))); %reaction force of spar due to moment
    dz = (span*reactionSpar)/(6*Ecf*MI); %deflection of the spar
    twistAngle = atan(dz/(plane.wing.chord*(.72/2))); %twist angle of wing due to spar deflection
    plane.performance.antTwist = twistAngle;
    if twistAngle >= 3*(pi()/180)
        plane.sanityFlag = false;
        plane.failureReason = "antennaMoment";
        plane.momentFail = 1;
    end

    
    %TODO: this should eventually check for -1 values in properties to make
    %sure everything was initialized correctly
    %for now the scoring functions take care of sanity checking the
    %velocity
end