clc; %Clears the Command Window
clear; %Clears all variables
close all; %Closes are figure windows

tic %Starts stopwatch timer
%Using structs the way we do here generates warnings that slow matlab down. Comment this out (and restart matlab) when debugging.
%warning('off','all')

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
numPowerSystems = 50 %DEBUGGING: Only search first 50 to decrease runtime while redesigning

vertStabAspectRatio = 2; %From aero calculations done beforehand
horizStabAspectRatio = 4.5; %^^^

wheelWidth = 0.5; %(in)
wheelRadius = 1.5; %^^^
wheelSurfaceArea = (2*pi*(wheelRadius)^2 + 2*pi*wheelRadius*wheelWidth)/144;

[wings] = wingData(aspectRatios, wingSpans); %Call GenWingsData function to make airfoil data lookup table

maxSavedPlanes = 10000; %About 98% of aircraft will fail and be overwritten so maxSavedPlanes does not have to equal max iterations
planes(1:maxSavedPlanes) = struct(airplaneClass); %Create a matrix to hold all the computed planes
index = 1;
iteration = 1;
spanFailCount = zeros([length(wingSpans) 6]); %This creates a matrix to check which failure conditions are most prevailent at each span value
failCountHeader = {'span', 'space', 'ep weight', 'takeoff', 'moment', 'converge'};
friendlyFailCount = [failCountHeader; num2cell(spanFailCount)]; %May have to move this to after loop

for AR = 1:size(aspectRatios)
    disp("At Aspect ratio " + AR);
    toc;
    for spanIndex = 1:length(wingSpans) %NOTE: the wings function cant handle this iteration yet, all calls of plane(index) also need updated, waiting on wings update
        spanFailCount(1 + spanIndex, 1) = wingSpans(spanIndex);
        disp("At Span " + spanIndex);
        toc;
        for airfoil = 1:1:8
            for powerIndex = 1:numPowerSystems
                for electronicPackageIndex = 1:length(m2PackageWeight)
                    %keyboard
                    for antennaIndex = 1:length(m3NumPassengers) %Changed this for the rerun - use smarter logic otherwise %for every amount of syringes try up to the maximum number of vials

                        planes(index) = struct(airplaneClass); %make sure to start with a clean slate at this index as this code writes over the index of failed airplanes. Without cleanup things like failure flags stay set even when they shouldn't

                        planes(index).fuselage.wheelSA = wheelSurfaceArea;

                        %load starting values for each plane
                        planes(index).wing.span = wingSpans(spanIndex);
                        planes(index).performance.epWeight = m2PackageWeight(electronicPackageIndex);
                        planes(index).wing.aspectRatio = aspectRatios(AR); %ratio between length and width of wing.
                        planes(index).performance.antennaLength = m3NumPassengers(antennaIndex);

                        %read values from wings matrix into aircraft
                        %properties. See wingClass.m for property descriptions
                        planes(index).wing.clw = wings(airfoil, 1, AR, spanIndex);
                        planes(index).wing.clm = wings(airfoil, 2, AR, spanIndex);     %cl max
                        planes(index).wing.cd = wings(airfoil, 3, AR, spanIndex);     %cd i zero velocity coefficient of drag
                        planes(index).wing.clFlap = wings(airfoil, 4, AR, spanIndex);  %weight
                        planes(index).wing.weight = wings(airfoil, 5, AR, spanIndex);    %airfoil name
                        planes(index).wing.chord = wings(airfoil, 6, AR, spanIndex);    %airfoil name
                        planes(index).wing.planformArea = wings(airfoil, 7, AR, spanIndex);    %airfoil name
                        planes(index).wing.surfaceArea = wings(airfoil, 8, AR, spanIndex);    %airfoil name
                        planes(index).wing.name = wings(airfoil, 9, AR, spanIndex);    %airfoil name
                        planes(index).wing.thickness=wings(airfoil, 10, AR, spanIndex); %thickness of the wing (ft)


                        %load values from power system table into aircraft
                        planes(index) = powerSelections(planes(index), MotorSpreadsheet, powerIndex);

                        %set payload and fuselage configuration
                        %plane(index).fuselage.numSyringes = syringes(syringe_index);
                        %plane(index).fuselage.numVials = num_vials;
                        planes(index) = fuselage(planes(index)); %calculate fuselage size based on number of vials and syringes
                        planes(index) = landingGear(planes(index), rho);
                        planes(index) = empennage(planes(index), horizStabAspectRatio, vertStabAspectRatio);

                        planes(index) = findTotalWeight(planes(index), m2PackageWeight(electronicPackageIndex), m3NumPassengers(antennaIndex));

                        %Want to include sanity check here that throws
                        %planes out which don't meet space/weight requirements
                        planes(index) = volSanityCheck(planes(index), m2PackageWeight(electronicPackageIndex));
                        if planes(index).volSanityFlag == 0
                            %disp("Too Big!");
                            if planes(index).epFail
                                spanFailCount(1 + spanIndex, 3) = spanFailCount(1 + spanIndex, 3) + 1;
                            elseif planes(index).spaceFail
                                spanFailCount(1 + spanIndex, 2) = spanFailCount(1 + spanIndex, 2) + 1;
                            end
                            break;
                        end

                        %simulate mission 2
                        planes(index) = GenVelocityTest(planes(index), 2, rho, temp); %2 signifies mission 2 configuration
                        planes(index) = TakeoffChecker(planes(index), 2, rho);
                        planes(index) = mission2score(planes(index), m2PackageWeight(electronicPackageIndex));

                        %simulate mission 3
                        planes(index) = GenVelocityTest(planes(index), 3, rho, temp); %2 signifies mission 2 configuration
                        planes(index) = TakeoffChecker(planes(index), 3, rho);
                        planes(index) = mission3score(planes(index), m3NumPassengers(antennaIndex));

                        planes(index) = sanityCheck(planes(index), wingSpans(spanIndex), m3NumPassengers(antennaIndex)); %make sure all the calculated values make sense and meet
                        %competition requirements. needs to be updated for
                        %this year's competition

                        if planes(index).sanityFlag
                            index = index + 1;
                            %disp("Sane!");
                        else
                            if planes(index).takeoffFail
                                spanFailCount(1 + spanIndex, 4) = spanFailCount(1 + spanIndex, 4) + 1;
                            elseif planes(index).momentFail
                                spanFailCount(1 + spanIndex, 5) = spanFailCount(1 + spanIndex, 5) + 1;
                            elseif planes(index).convergeFail
                                spanFailCount(1 + spanIndex, 6) = spanFailCount(1 + spanIndex, 6) + 1;
                            end
                            %disp("insane!");
                        end
                        iteration = iteration+1;
                        
                    end
                    
                end
            end
        end
    end
end
toc;
%Code below is from post.m file to run post-processing on the data. We should probably move this into a function call later.
%THIS HAS NOT BEEN UPDATED

score2 = zeros(1, length(planes)); %mission 2 scores of all aircraft
score3 = score2; %mission 3 scores
%scoreg = score2; %ground scores
%vials = score2; %number of vials in given airplane
%syringes = score2; %number of syringes in given airplane
for i = 1:length(planes) %load data from structure into arrays that are easier to work with
    if (planes(i).sanityFlag == 1) || (planes(i).volSanityFlag == 1)  %all values remain 1 unless the airplane fails sanity check. If the last plane analyzed is not sane the original sanity check can't throw it out.
        score2(i) = (planes(i).performance.score2);
        score3(i) = (planes(i).performance.score3);
        %vials(i) = plane(i).fuselage.numVials;
        %syringes(i) = plane(i).fuselage.numSyringes;
        %scoreg(i) = 10 + 2*(3*syringes(i)/5) + 5*vials(i);%Baded on time to run, load and unload syringes, load vial
    end
end
%clear plane
[M2, I2] = max(score2); %find the airplanes with the best individual mission scores
[M3, I3] = max(score3);
%[Mg, Ig] = max(scoreg);
score2 = 1 + score2/M2; %normalize scores against best performers
score3 = 2 + score3/M3;
%scoreg = scoreg/Mg;
score = score2 + score3 + 1; %scoreg + 1;
[M, I] = maxk(score, 200); %find the top 200 airplanes
winners = planes(I);
scores = score(I);
%save("winnersAR" + Aspect_Ratios(1) + ".mat", "winners"); %save top 100 airplanes. It is impractical to save all airplanes checked.
%save("winnersAR" + Aspect_Ratios(1) + "_scores.mat", "scores");

%save("successes.mat", "plane");
%clear; %once finished don't keep hogging ram. Previous experience with running multiple instances on analysis on one computer shows that ram usage is the first limiting factor in enabling
%the analysis to run, so clearing it as often as possible is important to avoid hogging ram from other programs. This ram hogging is the leading cause of crashing for this code.
%scatter(vials,score);