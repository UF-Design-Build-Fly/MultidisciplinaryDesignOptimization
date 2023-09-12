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
dThrustNeuralNet = importKerasNetwork("Fitting_3_model.h5");
vmean = 49.708;
vstd = 45.637;
dmean = 10.523;
dstd = 4.487;
pmean = 7.308001;
pstd = 3.348207;
rpmmean = 12508.699;
rpmstd = 9050.014;
dThrustStats = [vmean, vstd, dmean, dstd, pmean, pstd, rpmmean, rpmstd];

%Define plane properties to search
aspectRatios = [6 7 8 9 10]; %Wing aspect ratios
m2PackageWeight = 4:1:5; %Mission 2 Package Weight (pounds) 4:1:8 for first run 2.5x
%m2PackageWeight*0.03108096 %Convert to slugs
m3NumPassengers = 20:5:25; %Mission 3 Number of Passengers  15:3:30 for first run 2.5x
wingSpans = 4:1:5; %Wing Spans (feet)                       2.5:1.25:5 for first run 1.33x
load("MotorSpreadsheet.mat");
numPowerSystems = height(MotorSpreadsheet); %Get number of prop/motor/battery configs 33x
numPowerSystems = 20; %DEBUGGING: Only search first 50 to decrease runtime while redesigning

vertStabAspectRatio = 2; %From aero calculations done beforehand
horizStabAspectRatio = 4.5; %^^^

[wings] = GenWingData(aspectRatios, wingSpans); %Call GenWingsData function to make airfoil data lookup table

maxSavedPlanes = 100; %About 98% of aircraft will fail and be overwritten so maxSavedPlanes does not have to equal max iterations
planes(1:maxSavedPlanes) = struct(AirplaneClass); %Create a matrix to hold all the computed planes
index = 1;
iteration = 1;
spanFailCount = zeros([length(wingSpans) 6]); %This creates a matrix to check which failure conditions are most prevailent at each span value
failCountHeader = {'span', 'space', 'ep weight', 'takeoff', 'moment', 'converge'};

for aspectRatioIndex = 1:size(aspectRatios, 2)
    
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

                        %Set starting values for each plane
                        planes(index).wing.span = wingSpans(spanIndex);
                        planes(index).wing.aspectRatio = aspectRatios(aspectRatioIndex);
                        planes(index).performance.m2Weight = m2PackageWeight(m2PackageWeightIndex);
                        planes(index).performance.numPassengers = m3NumPassengers(m3PassengersIndex);

                        %Set values from wings matrix into plane
                        planes(index).wing = WingClass.SetWingData(planes(index).wing, wings, airfoilIndex, aspectRatioIndex, spanIndex);

                        %Set values from power system table into plane
                        planes(index).powerSystem = PowerClass.SetPowerSystemData(planes(index).powerSystem, MotorSpreadsheet, powerSystemIndex);

                        %Set payload and fuselage configuration
                        planes(index).fuselage = FuselageClass.CalcFuselageData(planes(index));
                        planes(index).fuselage = FuselageClass.GenLandingGear(planes(index));
                        planes(index).empennage = EmpennageClass.GenEmpennage(planes(index), horizStabAspectRatio, vertStabAspectRatio);

                        planes(index) = FindTotalWeight(planes(index));


                        %Simulate mission 2
                        planes(index) = GenVelocityTest(planes(index), 2, rho, temp, dThrustNeuralNet, dThrustStats); %2 signifies mission 2 configuration
                        planes(index) = TakeoffChecker(planes(index), 2, rho);
                        planes(index) = Mission2Score(planes(index));

                        %Simulate mission 3
                        planes(index) = GenVelocityTest(planes(index), 3, rho, temp, dThrustNeuralNet, dThrustStats); %3 signifies mission 3 configuration
                        planes(index) = TakeoffChecker(planes(index), 3, rho);
                        planes(index) = Mission3Score(planes(index));
                        
                        %Simulate ground mission
                        totalAssemblyTime = 200; %(s)
                        timePerPassenger = 5; %(s)
                        planes(index).performance.scoreGM = totalAssemblyTime + timePerPassenger*planes(index).performance.numPassengers;
                        
                        %make sure all the calculated values make sense and meet
                        planes(index) = SanityCheck(planes(index)); 
                        if planes(index).sanityFlag
                            index = index + 1;
                        else
                            if planes(index).takeoffFail
                                spanFailCount(spanIndex, 4) = spanFailCount(spanIndex, 4) + 1;
                            %elseif planes(index).momentFail
                            %    spanFailCount(spanIndex, 5) = spanFailCount(spanIndex, 5) + 1;
                            elseif planes(index).convergeFail
                                spanFailCount(spanIndex, 6) = spanFailCount(spanIndex, 6) + 1;
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

%save("winners.mat", "winners"); %Save top results
%save("winners_scores.mat", "scores");

clear; %Clear variables to free RAM. RAM usage is the limiting factor in enabling the analysis to run and avoid crashing.