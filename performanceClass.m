classdef performanceClass

  properties
%some of these properties don't make sense to me, they're just added because they were in last year's list
    
%Mission 1: Deployment flight. No payload, complete 3 laps in 5 minutes. Get 1 whole point for not crashing
%Mission 2: Carry syringes. At least 10. Max window of 5 minutes. Score is syringes/flight time
%Mission 3: Carry vial packages. At most 1/10th of syringes. Max window 10 mins. Deploy a package on each lap. Score is num succesful deployments.
    
    totalEmptyWeight = -1;
    totalWeight2 = -1;
    totalWeight3 = -1;

    velocity1 = -1; 
    velocity2 = -1;
    velocity3 = -1;

    takeoffDist1 = -1;
    takeoffDist2 = -1;
    takeoffDist3 = -1;

    landingSpeed1 = -1;
    landingSpeed2 = -1;
    landingSpeed3 = -1;
    
    time1 = -1;
    time2 = -1;
    time3_perLap = -1;

    Nlaps1 = -1;
    Nlaps2 = -1;
    Nlaps3 = -1;

    score1 = -1;
    score2 = -1;
    score3 = -1;

    drag1 = -1; %drag at cruise velocity as calculated by gen velocity solver
    drag2 = -1;
    drag3 = -1;

  end

  methods
  end

end