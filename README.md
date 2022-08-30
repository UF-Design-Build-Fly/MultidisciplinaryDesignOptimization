# Sensitivity-Analysis
Project to determine optimal aircraft configuration for the 2021-2022 AIAA Design, Build, Fly competition.

Performs an exhaustive search over a given list of aircraft parameters (defined in xxxClas.m for each major portion of the aircraft). Once the parameters for each aircraft are chosen a velocity solver is used to estimate cruise velocity for that aircraft. Relative performance can be determined from this velocity and number/size of scoring elements chosen. After all combinations have been simulated the top 200 aircraft are then plugged into the scoring formula defined in the rules to estimate the actual score achievable in competition.

Sensitivity_Analysis_2020.m is included for reference but is no longer functional folder structure. Sensitivity_Analysis_2021.m is the main project file that outlines the necessary computations. It is actively maintained and will be updated as empirical data is gathered to improve the estimations.