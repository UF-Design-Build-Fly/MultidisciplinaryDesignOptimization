function [plane] = mission2score(plane, Electronic_Package_Weight)
    %Ian-9/18/2022-Completed for now, need to debug, need to use
    %profile.png to find more accurate times
    %M2 score = 1 + [N_(payload weight * # laps flown) / Max_(payload weight * # laps flown]. Need to normalize
    %must be done within 10 minutes

    %below is old code for reference
    %plane.performance.time2 = (plane.performance.lapDist*3)/plane.performance.velocity2;
    %plane.performance.score2 = (plane.fuselage.numSyringes/plane.performance.time2); %leave the task of normalizing scores to post-processing   
    %DEBUG - leaving out +1 for now as well just to get a good idea of the ratios being produced


    lapdist = (500*4)+(500*pi)+(250*pi); %estimates the overall lap distance in ft
    if plane.power.time <= 10
        plane.performance.laps2 = (plane.power.time*60)/(lapdist/plane.performance.velocity2); %Max # of laps at Vmax w/battery. May need unit conversion (assumes plane.performance.velocity2 is in ft/s
    else
        plane.performance.laps2 = 600/(lapdist/plane.performance.velocity2); %Max # of laps at Vmax in 10 minutes
    end
    plane.performance.score2 = Electronic_Package_Weight*plane.performance.laps2;
end