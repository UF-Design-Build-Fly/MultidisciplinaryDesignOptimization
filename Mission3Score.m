function [plane] = Mission3Score(plane, Antenna_Length)

    %M3 score = 2 + [N_(antenna length / mission time) / Max_(antenna length / mission time)]
    %If mission (3 laps) takes longer than 5 minutes, sanity check will throw
    %out plane
    lapdist = (500*4)+(500*pi)+(250*pi); %estimates the overall lap distance in ft
    time3 = (3*(lapdist/plane.performance.velocity3))/60; %Calculates the time taken to fly 3 laps 
    plane.performance.score3 = Antenna_Length/time3; %Calculates mission 3 score


    if plane.powerSystem.time <= 10
        plane.performance.laps2 = (plane.powerSystem.time*60)/(lapDist/plane.performance.velocity2); %Max # of laps at Vmax w/battery. May need unit conversion (assumes plane.performance.velocity2 is in ft/s
    else
        plane.performance.laps2 = 600/(lapDist/plane.performance.velocity2); %Max # of laps at Vmax in 10 minutes
    end

end