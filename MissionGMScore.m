function plane = MissionGMScore(plane)
    
    %Only time taken to load (assume same time taken to unload)
    timeDrop = 60;
    timePerTank = 10;
    
    numTanks = plane.performance.m2Weight/2;
    plane.performance.scoreGM = timeDrop + timePerTank*numTanks;

end