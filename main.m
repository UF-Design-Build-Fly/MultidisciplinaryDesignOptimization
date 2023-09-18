clc; %Clears the Command Window
clear; %Clears all variables
close all; %Closes are figure windows

tic %Starts stopwatch timer

%Using structs the way we do here generates warnings that slow matlab down. Comment this out (and restart matlab) when debugging.
warning('off','all');


%Define constants for aero equations
%rho = 0.00235308; %Density air at Tuscon, Az (slug/ft^3)
rho = 0.002391; %Density air at Whichita, Ks with average climate data from April 2021
temp = 293; %Temperature in kelvin at competition site


%Dynamic thrust Neural Network setup
dynamicThrustNet = importKerasNetwork("Fitting_3_model.h5");
vMean = 49.708;
vStd = 45.637;
dMean = 10.523;
dStd = 4.487;
pMean = 7.308001;
pStd = 3.348207;
rpmMean = 12508.699;
rpmStd = 9050.014;
dynamicThrustStats = [vMean, vStd, dMean, dStd, pMean, pStd, rpmMean, rpmStd];


%Define plane properties to search
aspectRatios = [6 7 8 9 10];                    % more spans??? 4, 5
m2PackageWeight = (4:1:5); %(lbs)                % 4:1:8 for first run 2.5x
m3NumPassengers = 20:5:25;                      % 15:3:30 for first run 2.5x
wingSpans = 4:1:5;                              % 2.5:1.25:5 for first run 1.33x
load("MotorSpreadsheet.mat");
%numPowerSystems = height(MotorSpreadsheet);
numPowerSystems = 20; %DEBUGGING: Only search first 20 to decrease runtime while redesigning
numAirfoils = 8; %Airfoils define in GenWingData()
numSavedPlanes = 100; %About 98% of aircraft will fail and be overwritten so maxSavedPlanes does not have to equal max iterations

vertStabAspectRatio = 2; %From aero calculations done beforehand
horizStabAspectRatio = 4.5; %^^^


wingLookupTable = GenWingData(aspectRatios, wingSpans); %Creates airfoil data lookup table

%Setup Fail Table - 2022-2023
%spanFailCount = zeros([length(wingSpans) 6]); %This creates a matrix to check which failure conditions are most prevailent at each span value
%failCountHeader = {'span', 'space', 'ep weight', 'takeoff', 'moment', 'converge'};

%Progress Bar Setup
totalPlanesSearched = length(aspectRatios)*length(wingSpans)*numAirfoils*numPowerSystems*length(m2PackageWeight)*length(m3NumPassengers);
progressBar = waitbar(0, "Searching, Calculating, Failing");


planes(1:numSavedPlanes) = struct(AirplaneClass); %Create a matrix to hold all the computed planes
index = 1;
iteration = 1;
for aspectRatioIndex = 1:length(aspectRatios)
    for spanIndex = 1:length(wingSpans)

        %spanFailCount(spanIndex, 1) = wingSpans(spanIndex);

        %disp("Aspect Ratio: " + aspectRatioIndex + "/" + length(aspectRatios) + " - Span: " + spanIndex + "/" + length(wingSpans));
        %toc;

        for airfoilIndex = 1:1:numAirfoils
            for powerSystemIndex = 1:numPowerSystems
                for m2PackageWeightIndex = 1:length(m2PackageWeight)
                    for m3PassengersIndex = 1:length(m3NumPassengers)
    
                        iteration = iteration+1;
                        waitbar(iteration/totalPlanesSearched, progressBar);
                

                        %Start with a clean slate(overwrite failure flags) in case the previous plane failed
                        planes(index) = struct(AirplaneClass);

                        %Set starting values for each plane
                        planes(index).wing.span = wingSpans(spanIndex);
                        planes(index).wing.aspectRatio = aspectRatios(aspectRatioIndex);
                        planes(index).performance.m2Weight = m2PackageWeight(m2PackageWeightIndex);
                        planes(index).performance.numPassengers = m3NumPassengers(m3PassengersIndex);

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
                        planes(index) = TakeoffChecker(planes(index), 2, rho);
                        if (planes(index).performance.takeoffDist2 >= 20)
                            continue;
                        end
                        planes(index) = TakeoffChecker(planes(index), 3, rho);
                        if (planes(index).performance.takeoffDist3 >= 20)
                            continue;
                        end

                        %Simulate mission velocities
                        planes(index) = GenVelocityTest(planes(index), 2, rho, temp, dynamicThrustNet, dynamicThrustStats); %2 signifies mission 2 configuration
                        if (planes(index).performance.velocity2 == -1)
                            continue;
                        end
                        planes(index) = GenVelocityTest(planes(index), 3, rho, temp, dynamicThrustNet, dynamicThrustStats); %3 signifies mission 3 configuration
                        if (planes(index).performance.velocity3 == -1)
                            continue;
                        end

                        %Calculate scores
                        planes(index) = Mission2Score(planes(index));
                        planes(index) = Mission3Score(planes(index));
                        
                        %Calculate ground mission score
                        totalAssemblyTime = 200; %(s)
                        timePerPassenger = 5; %(s)
                        planes(index).performance.scoreGM = totalAssemblyTime + timePerPassenger*planes(index).performance.numPassengers;
                        
                        %Make sure all values are calculated
                        %if (SanityCheck(planes(index)))
                            index = index + 1;
                        %end

                    end %End m3Passengers Loop
                end %End m2PackageWeight Loop
            end %End powerSystem Loop
        end %End airfoil Loop
    end %End wingSpan Loop
end %End aspectRatio Loop

toc;

%friendlyFailCount = [failCountHeader; num2cell(spanFailCount)]; %Add headers to fail count table

scoresM2 = zeros(1, length(planes)); %Initilize arrays with 0s
scoresM3 = scoresM2; %Initilize arrays with 0s
scoresGM = scoresM2; %Initilize arrays with 0s

for (i = 1:index-1) %Load data into arrays that are easier to work with
    scoresM2(i) = planes(i).performance.score2;
    scoresM3(i) = planes(i).performance.score3;
    scoresGM(i) = planes(i).performance.scoreGM;
end

scoresM2 = scoresM2/max(scoresM2); %Normalize scores against best performers
scoresM3 = scoresM3/max(scoresM3);
scoresGM = scoresGM/max(scoresGM);

score = scoresM2 + scoresM3 + scoresGM;

[winners, indices] = maxk(score, 100); %Find the top planes
scores = score(indices);

%save("winners.mat", "winners"); %Save top results
%save("winners_scores.mat", "scores");
%save("planes.mat", "planes");

clear; %Clear variables to free RAM. RAM usage is the limiting factor in enabling the analysis to run and avoid crashing.