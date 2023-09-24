function plane = MissionGMScore(plane)
    
    %Only time taken to load (assume same time taken to unload)
    %2 Crew, 2 EMTS, 1 patient, 1 medical cabinet
    timeCrew = 40;
    timePerWingAssembly = 30;
    timePerPassenger = 5;
    
    if (plane.wing.span > 7.5)
        timeCrew = timeCrew + 2*timePerWingAssembly;
    elseif (plane.wing.span > 2.5)
        timeCrew = timeCrew + timePerWingAssembly;
    end
    
    plane.performance.scoreGM = timeCrew + timePerPassenger*plane.performance.numPassengers;

end