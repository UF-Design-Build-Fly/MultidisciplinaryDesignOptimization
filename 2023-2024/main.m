%% ============================= Cleanup Environment ========================================
clc;					% Clears the Command Window
clear;					% Clears all variables
close all;				% Closes are figure windows
warning('off','all');	% Using structs the way we do here generates warnings that slow matlab down. Comment this out (and restart matlab) when debugging.

%% ======================== Constants and Search Parameters ==================================
% Aero Equations
%rho = 0.002391;				% Air Density in Whichita, Ks (slug/ft^3)
rho = 0.00235308;				% Air Density in Tuscon, Az (slug/ft^3)
temp = 303;						% Temperature in Tuscon, Az (Kelvin)
windSpeed = 13.2;				% Average Wind Speed in Tucson, Az (ft/s)

% Search Parameters
aspectRatios = 4:1:8;
m2PackageWeight = 4.5:1:11.5;	% (lbs)
wingSpans = 3:1:6;				% (ft)

% Tail Configuration
vertStabAspectRatio = 2;		%From aero calculations done beforehand
horizStabAspectRatio = 4.5;		% ^^^

%% ============================== Import Motor & Wing Data ===================================
load("MotorSpreadsheet2024.mat");
MotorSpreadsheet = sortrows(MotorSpreadsheet, 'PitchSpeedfts', 'descend');
numPowerSystems = height(MotorSpreadsheet);
numPowerSystems = 10;			% DEBUGGING: Search subset to decrease runtime

wingLookupTable = GenWingData(aspectRatios, wingSpans);		% Creates airfoil data lookup table
numAirfoils = 1;				% Airfoils define in GenWingData()

%% ================================== Unused Features ========================================
% Dynamic Thrust Neural Network
%dynamicThrustNet = importKerasNetwork("Fitting_3_model.h5");
%vMean = 49.708;
%vStd = 45.637;
%dMean = 10.523;
%dStd = 4.487;
%pMean = 7.308001;
%pStd = 3.348207;
%rpmMean = 12508.699;
%rpmStd = 9050.014;
%dynamicThrustStats = [vMean, vStd, dMean, dStd, pMean, pStd, rpmMean, rpmStd];

%Setup Fail Table - 2022-2023
%spanFailCount = zeros([length(wingSpans) 6]); %This creates a matrix to check which failure conditions are most prevailent at each span value
%failCountHeader = {'span', 'space', 'ep weight', 'takeoff', 'moment', 'converge'};

%% ==================================== Search Loop ==========================================
%Progress Bar Setup
progressBar = waitbar(0, "Starting, Calculating, Failing");

totalPlanesSearched = length(aspectRatios)*length(wingSpans)*numAirfoils*numPowerSystems*length(m2PackageWeight);
%numSavedPlanes = 10000;        % About 98% of aircraft will fail and be overwritten so maxSavedPlanes does not have to equal max iterations
planes(1:totalPlanesSearched) = struct(AirplaneClass);% Allocate plane storage array for speed

tic %Starts stopwatch
index = 1;
superIteration = 1;
for aspectRatioIndex = 1:length(aspectRatios)
    for spanIndex = 1:length(wingSpans)
        for airfoilIndex = 1:1:numAirfoils
            for powerSystemIndex = 1:numPowerSystems

                superIteration = superIteration + 1;
                percentComplete = superIteration*length(m2PackageWeight)/totalPlanesSearched;
                timeRemaining = toc * (1/percentComplete - 1);
                waitbar(percentComplete, progressBar,  round(timeRemaining)+ "s Remaing - " + 100*round(percentComplete, 3) + "%");

                for m2PackageWeightIndex = 1:length(m2PackageWeight)
                    %quit = 0;
                    %for m3PassengersIndex = 1:length(m3NumPassengers)
                        
                    %Start with a clean slate(overwrite failure flags) in case the previous plane failed
                    planes(index) = struct(AirplaneClass);

                    %Set starting values for each plane
                    planes(index).wing.span = wingSpans(spanIndex);
                    planes(index).wing.aspectRatio = aspectRatios(aspectRatioIndex);
                    planes(index).performance.m2Weight = m2PackageWeight(m2PackageWeightIndex);

                    %Set values from wings matrix into plane
                    planes(index).wing = WingClass.SetWingData(planes(index).wing, wingLookupTable, airfoilIndex, aspectRatioIndex, spanIndex);

                    %Set values from power system table into plane
                    planes(index).powerSystem = PowerClass.SetPowerSystemData(planes(index).powerSystem, MotorSpreadsheet, powerSystemIndex);

                    %Set payload and fuselage configuration
                    planes(index).fuselage = FuselageClass.CalcFuselageData(planes(index));
                    planes(index).fuselage = FuselageClass.GenLandingGear(planes(index));
                    planes(index).empennage = EmpennageClass.GenEmpennage(planes(index), horizStabAspectRatio, vertStabAspectRatio);

                    planes(index) = FindTotalWeight(planes(index));


                    %Simulate mission takeoffs
                    %planes(index) = TakeoffChecker(planes(index), 2, rho);
                    %if (planes(index).performance.takeoffDist2 > 20)
                    %    quit = 1; %M2 Weight Takeoff/Flight Failed
                    %    break;
                    %end
                    %planes(index) = TakeoffChecker(planes(index), 3, rho);
                    %if (planes(index).performance.takeoffDist3 > 20)
                    %    break;
                    %end

                    %Simulate mission velocities
                    planes(index) = GenVelocityTest(planes(index), 2, rho, temp);%, dynamicThrustNet, dynamicThrustStats); %2 signifies mission 2 configuration
                    if (planes(index).performance.velocity2 < planes(index).performance.landingSpeed2 || planes(index).performance.time2 > 300)
                        %quit = 1; %M2 Weight Takeoff/Flight Failed
                        break;
                    end
                    planes(index) = GenVelocityTest(planes(index), 3, rho, temp);%, dynamicThrustNet, dynamicThrustStats); %3 signifies mission 3 configuration
                    if (planes(index).performance.velocity3 < planes(index).performance.landingSpeed3)
                        break;
                    end

                    %Calculate scores
                    planes(index) = Mission2Score(planes(index), windSpeed);
                    planes(index) = Mission3Score(planes(index), windSpeed);
                    planes(index) = MissionGMScore(planes(index));
                    
                    %Make sure all values are calculated
                    %if (SanityCheck(planes(index)))
                    index = index + 1;
                    %end

                    %end %End m3Passengers Loop
                    %if (quit)
                    %    break; %M2 Weight Takeoff/Flight Failed
                    %end
                end %End m2PackageWeight Loop
            end %End powerSystem Loop
        end %End airfoil Loop
    end %End wingSpan Loop
end %End aspectRatio Loop

toc;

%friendlyFailCount = [failCountHeader; num2cell(spanFailCount)]; %Add headers to fail count table

scoresM2 = zeros(1, index-1); %Initilize arrays with 0s
scoresM3 = scoresM2; %Initilize arrays with 0s
scoresGM = scoresM2; %Initilize arrays with 0s

for (i = 1:index-1) %Load data into arrays that are easier to work with
    scoresM2(i) = planes(i).performance.score2;
    scoresM3(i) = planes(i).performance.score3;
    scoresGM(i) = planes(i).performance.scoreGM;
end

for (i = 1:index-1) %Normalizes scores in plane performance object
    planes(i).performance.score2Normalized = planes(i).performance.score2/max(scoresM2);
    planes(i).performance.score3Normalized = planes(i).performance.score3/max(scoresM3);
    planes(i).performance.scoreGMNormalized = min(scoresGM)/planes(i).performance.scoreGM;
    planes(i).performance.scoreTotal = 1 + planes(i).performance.score2Normalized + planes(i).performance.score3Normalized + planes(i).performance.scoreGMNormalized;
end

scoresM2 = scoresM2/max(scoresM2); %Normalize scores against best performers
scoresM3 = scoresM3/max(scoresM3);
scoresGM = min(scoresGM)./scoresGM;

score = 1 + scoresM2 + scoresM3 + scoresGM;

[winners, indices] = maxk(score, 100); %Find the top planes
topPlanes = planes(indices);
%save("winners.mat", "winners"); %Save top results
%save("winners_scores.mat", "scores");
save("topPlanes.mat", "topPlanes");

%clear; %Clear variables to free RAM. RAM usage is the limiting factor in enabling the analysis to run and avoid crashing.