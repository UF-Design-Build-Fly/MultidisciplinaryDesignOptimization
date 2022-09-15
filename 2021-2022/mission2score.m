function [plane] = mission2score(plane)
    %M2 score = 1 + (syringes/time)/(max syringes/time). Need to normalize
    %to eliminate the unit time is measured in and to create a score that's
    %realistic to what we'll see in competition.
    %timed for 3 laps
    %must be done within 5 minutes
    plane.performance.time2 = (plane.performance.lapDist*3)/plane.performance.velocity2;
    plane.performance.score2 = (plane.fuselage.numSyringes/plane.performance.time2); %leave the task of normalizing scores to post-processing
                                                                                     %DEBUG - leaving out +1 for now as well just to get a good idea of the ratios being produced
end