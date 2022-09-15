# Sensitivity-Analysis
Project to determine optimal aircraft configuration for the 2022-2023 AIAA Design, Build, Fly competition.

Performs an exhaustive search over a given list of aircraft parameters (defined in xxxClass.m for each major portion of the aircraft). Once the parameters for each aircraft are chosen a velocity solver is used to estimate cruise velocity for that aircraft. Relative performance can be determined from this velocity and number/size of scoring elements chosen. After all combinations have been simulated the top 200 aircraft are then plugged into the scoring formula defined in the rules to estimate the actual score achievable in competition.

Sensitivity_Analysis_2021.m is included under the 2021-2022 folder for reference. Sensitivity_Analysis_2021.m is included as well but it is no longer in a functional folder structure. 

Sensitivity_Analysis_xxxx.m is the main project file that outlines the necessary computations. It is actively maintained and will be updated as empirical data is gathered to improve the estimations.
