# Multi-disciplinary Design Optimization

Project to determine optimal aircraft configuration for the 2023-2024 AIAA Design, Build, Fly competition.

Performs an exhaustive search over a given list of aircraft parameters. Once the parameters for each aircraft are chosen a velocity solver is used to estimate cruise velocity for that aircraft. Relative performance can be determined from this velocity and number/size of scoring elements chosen. After all combinations have been simulated the top 200 aircraft are then plugged into the scoring formula defined in the rules to estimate the actual score achievable in competition.

Previous programs are included for reference but are no longer functional. main.m is the main project file that outlines the necessary computations. It is actively maintained and will be updated as empirical data is gathered to improve the estimations.

Notes:

- Naming Convention Stuff
  - functions should be CapitalCapital
  - variables should be lowercaseCapital
  - functions should probably contain a verb
  - variables try not to abbreviate but abbreviate
    - number -> num
    - generate -> gen
    - calculate -> calc

TODO:

- move landinggear into seperate class
- remove reassignment of numPowerSystems
- change max saved planes
- get weight of wood things
- debug total weight. last year things needed a weight
    - add wheel weight
- update times for assemblying plane and putting in passengers


Ask HARIS
- Get updated rho value for Whitchita
- Is the frontalSurfaceArea suppose to include all the tapering surface area of the nose?
- Haris/Austin get turn acceleration estimate put into performanceclass change default to -1