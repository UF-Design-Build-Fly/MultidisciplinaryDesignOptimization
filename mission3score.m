function [plane] = mission3score(plane)
%score = 2 + (num_deployments/max_deployments)
%again, for this function, leave normalizing by max score to post procssing
%so that we don't have to go back through and update every aircraft
%retroactively when a plane achieves a new high score
    %calculate lap time from lap profile given in profile.png included in
    %repository
    landingRollTime = plane.performance.totalWeight3 * 1.25; %estimate based on video review of past aircraft
    landingRollTime = min(max(landingRollTime,5),15); %for sanity constain roll time between 5 and 15 seconds. 
    unloadTime = 15; %estimate - update later in design process when time is known
    taxiTime = 5; %average time to get to drop off area if landing wasn't 
                  %perfect. This is really just a margin safety factor to avoid
                  %the analysis making an airplane that is too big to be
                  %practical
    plane.performance.time3_perLap = (((200*pi)+1000)/plane.performance.velocity3) + ...
        (500/(0.75*plane.performance.velocity3)) + (500/plane.performance.landingSpeed3) + ...
        landingRollTime + unloadTime + taxiTime; %measured in seconds
    
    %10 minute window for mission. If power system can't support 10 minutes fly as long as we can
    %careful - flight time here is in minutes but everything else is in seconds!
    %DEBUG - change this in the future. Requires change to motor data table
    if(plane.power.time < 10) 
        flightTime = plane.power.time;
    else
        flightTime = 10; 
    end
    
    %Number of laps = number of laps supported by mission time limit, power
    %system time limit, or number of vials, whichever is less.
    %Only care about "scoreable" laps where we can deploy vials
    plane.performance.Nlaps3 = (flightTime*60)/plane.performance.time3_perLap;
    
    if(plane.fuselage.numVials < floor(plane.performance.Nlaps3)) %if we're not carrying enough vials to use all flight time
        plane.performance.Nlaps3 = plane.fuselage.numVials;
    end
    
    plane.performance.time3 = plane.performance.time3_perLap*plane.performance.Nlaps3;
    plane.performance.score3 = plane.performance.Nlaps3;%DEBUG - leaving out +2 for now just to get a good idea of the ratios being produced 
end