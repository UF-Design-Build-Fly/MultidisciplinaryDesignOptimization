classdef performanceClass

  properties
%some of these properties don't make sense to me, they're just added because they were in last year's list
    
%Mission 1: Deployment flight. No payload, complete 3 laps in 5 minutes. Get 1 whole point for not crashing
%Mission 2: Carry syringes. At least 10. Max window of 5 minutes. Score is syringes/flight time
%Mission 3: Carry vial packages. At most 1/10th of syringes. Max window 10 mins. Deploy a package on each lap. Score is num succesful deployments.
    
    epWeight = -1;
    antennaLength = -1;

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

    dynamicThrust = -1;
    
%     time1 = -1;
    time2 = -1;
    time3 = -1;
%     time3_perLap = -1;

%     Nlaps1 = -1;
%     Nlaps2 = 3;
%     Nlaps3 = -1;

    score1 = -1;
    score2 = -1;
    score3 = -1;

    totalDrag = -1; %DEBUG VALUE
    inducedDrag = -1;
    parasiticDrag = -1;
    skinDrag = -1;
    wingPara = -1;
    hStabPara = -1;
    vStabPara = -1;
    fusePara = -1;
    gearPara = -1;
    antDrag = -1;
    antTwist = -1;
    
%     drag1 = -1; %drag at cruise velocity as calculated by gen velocity solver  
    drag2 = -1;
    drag3 = -1;
    
%     lapDist = 500+200*pi+1000+500; %distances for each velocity profile.
                                   %see profile.png in analysis folder

  end

  methods
  end

end