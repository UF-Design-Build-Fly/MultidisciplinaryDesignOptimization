## -*- texinfo -*-
## @deftp {Class} performanceClass
##
## @end deftp

classdef performanceClass

  properties
%some of these properties don't make sense to me, they're just added because they were in last year's list
    
%Mission 1: Deployment flight. No payload, complete 3 laps in 5 minutes. Get 1 whole point for not crashing
%Mission 2: Carry syringes. At least 10. Max window of 5 minutes. Score is syringes/flight time
%Mission 3: Carry vial packages. At most 1/10th of syringes. Max window 10 mins. Deploy a package on each lap. Score is num succesful deployments.
    
    emptyWeight;
    weight2;
    weight3;

    velocity1; 
    velocity2;
    velocity3;

    takeoffDist1;
    takeoffDist2;
    takeoffDist3;

    time2;
    time3_perLap;

    Nlaps;

    
    score1;
    score2;
    score3;
  end

  methods
  end

end