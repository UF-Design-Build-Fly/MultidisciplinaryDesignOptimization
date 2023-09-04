clc; %Clears the Command Window
clear; %Clears all variables
close all; %Closes are figure windows

tic %Starts stopwatch timer
%Using structs the way we do here generates warnings that slow matlab down. Comment this out (and restart matlab) when debugging.
warning('off','all')

%Define constants for aero equations
rho = 0.00235308; %Density air at Tuscon, Az (slug/ft^3)
%rho = 0.002391; %Density air at Whichita, Ks with average climate data from April 2021
temp = 293; %Temperature in kelvin at competition site

%Define plane properties to search
aspectRatios = [6 7 8 9 10]; %Wing aspect ratios
m2PackageWeight = 4:1:7; %Mission 2 Package Weight (pounds)
m3NumPassengers = 20:5:30; %Mission 3 Number of Passengers
wingSpans = 3:1:5; %Wing Spans (feet)
load("MotorSpreadsheet.mat");
numPowerSystems = height(MotorSpreadsheet); %Get number of prop/motor/battery configs
numPowerSystems = 50; %DEBUGGING: Only search first 50 to decrease runtime while redesigning

vertStabAspectRatio = 2; %From aero calculations done beforehand
horizStabAspectRatio = 4.5; %^^^

wheelWidth = 0.5; %(in)
wheelRadius = 1.5; %^^^
wheelSurfaceArea = (2*pi*(wheelRadius)^2 + 2*pi*wheelRadius*wheelWidth)/144;

[wings] = GenWingData(aspectRatios, wingSpans); %Call GenWingsData function to make airfoil data lookup table

maxSavedPlanes = 10000; %About 98% of aircraft will fail and be overwritten so maxSavedPlanes does not have to equal max iterations
planes(1:maxSavedPlanes) = struct(airplaneClass); %Create a matrix to hold all the computed planes
index = 1;
iteration = 1;
spanFailCount = zeros([length(wingSpans) 6]); %This creates a matrix to check which failure conditions are most prevailent at each span value
failCountHeader = {'span', 'space', 'ep weight', 'takeoff', 'moment', 'converge'};

for aspectRatioIndex = 1:size(aspectRatios)
    
    disp("Aspect Ratio: " + aspectRatioIndex + "/" + size(aspectRatios, 2));
    toc;

    for spanIndex = 1:length(wingSpans)

        spanFailCount(spanIndex, 1) = wingSpans(spanIndex);
        disp("Span: " + spanIndex + "/" + size(wingSpans, 2));
        toc;

        for airfoilIndex = 1:1:8
            for powerSystemIndex = 1:numPowerSystems
                for m2PackageWeightIndex = 1:length(m2PackageWeight)
                    for m3PassengersIndex = 1:length(m3NumPassengers)
                
                        %Start with a clean slate(overwrite failure flags) in case the previous plane failed
                        planes(index) = struct(AirplaneClass); 

                        planes(index).fuselage.wheelSA = wheelSurfaceArea;

                        %Set starting values for each plane
                        planes(index).wing.span = wingSpans(spanIndex);
                        planes(index).wing.aspectRatio = aspectRatios(aspectRatioIndex);
                        planes(index).performance.epWeight = m2PackageWeight(m2PackageWeightIndex);
                        planes(index).performance.numPassengers = m3NumPassengers(m3PassengersIndex);

                        %Set values from wings matrix into plane
                        planes(index).wing = SetWingData(planes(index).wing, wings, airfoilIndex, aspectRatioIndex, spanIndex);

                        %Set values from power system table into plane
                        planes(index).powerSystem = SetPowerSystemData(planes(index).powerSystem, MotorSpreadsheet, powerSystemIndex);

                        %Set payload and fuselage configuration
                        planes(index).fuselage = CalcFuselageData(planes(index).fuslage);
                        planes(index) = landingGear(planes(index), rho);
                        planes(index) = empennage(planes(index), horizStabAspectRatio, vertStabAspectRatio);

                        planes(index) = findTotalWeight(planes(index), m2PackageWeight(m2PackageWeightIndex), m3NumPassengers(m3PassengersIndex));


                        %Simulate mission 2
                        planes(index) = GenVelocityTest(planes(index), 2, rho, temp); %2 signifies mission 2 configuration
                        planes(index) = TakeoffChecker(planes(index), 2, rho);
                        planes(index) = mission2score(planes(index), m2PackageWeight(m2PackageWeightIndex));

                        %Simulate mission 3
                        planes(index) = GenVelocityTest(planes(index), 3, rho, temp); %3 signifies mission 3 configuration
                        planes(index) = TakeoffChecker(planes(index), 3, rho);
                        planes(index) = mission3score(planes(index), m3NumPassengers(m3PassengersIndex));
                        
                        %make sure all the calculated values make sense and meet
                        planes(index) = sanityCheck(planes(index), wingSpans(spanIndex), m3NumPassengers(m3PassengersIndex)); 
                        if planes(index).sanityFlag
                            index = index + 1;
                        else
                            if planes(index).takeoffFail
                                spanFailCount(1 + spanIndex, 4) = spanFailCount(1 + spanIndex, 4) + 1;
                            elseif planes(index).momentFail
                                spanFailCount(1 + spanIndex, 5) = spanFailCount(1 + spanIndex, 5) + 1;
                            elseif planes(index).convergeFail
                                spanFailCount(1 + spanIndex, 6) = spanFailCount(1 + spanIndex, 6) + 1;
                            end
                        end

                        iteration = iteration+1;
                    end
                end
            end
        end
    end
end

toc;

friendlyFailCount = [failCountHeader; num2cell(spanFailCount)]; %Add headers to fail count table

scoresM2 = zeros(1, length(planes)); %Initilize arrays with 0s
scoresM3 = scoresM2; %Initilize arrays with 0s
scoresGM = scoresM2; %Initilize arrays with 0s

for i = 1:index %Load data into arrays that are easier to work with
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

save("winners.mat", "winners"); %Save top results
save("winners_scores.mat", "scores");

clear; %Clear variables to free RAM. RAM usage is the limiting factor in enabling the analysis to run and avoid crashing.