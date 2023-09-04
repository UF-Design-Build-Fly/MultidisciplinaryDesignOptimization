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
- Update Motorspreadsheet (Rerun script, recreate spreadsheet, convert to .mat file)
  - update estimated weight to be in pounds time to be in seconds
- remove reassignment of numPowerSystems
- add wheel weight
- get weight of wood things
- debug total weight. last year things needed a weight
- change max saved planes

- Have Haris look at GenWingData function, Add lift from fuselage?, GenEmpennage function
- Have Dylan look at fuselage and landing gear math - switch to carbon fiber
- Haris/Austin get turn acceleration estimate put into performanceclass change default to -1

- spanFailCount may have to be changed to something else or removed
- whats happening with surface area(wheel surface area to be exact)
- Examine GenVelocityTest, TakeoffChecker
