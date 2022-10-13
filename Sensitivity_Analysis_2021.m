clear; clc;
%Ian-10/8/2022-Mostly done with debugging, still wamt to add more data processing, but can do that post analysis

warning('off','all') %using structs the way we do here generates a flood of warnings that slows matlab down. Comment this out (and restard matlab) when debugging.
tic

%Here we introduce constants for use in aero equations
rho = 0.00235308; %Density air at Tuscon, Az (slug/ft^3)
%rho = 0.002391; %Density air at Whichita, Ks with average climate data from April 2021
Temp = 293; %temperature in kelvin at competition site

%This section introduces values which are later iterated on in the for loop
Aspect_Ratios = 7;%[6, 7, 8, 9]; %, 10, 11, 12, 13]; %wing aspect ratios to consider
Electronic_Package_Weight = 4;%3:.5:5; %Electronic Package Weight for mission 2 in pounds
Antenna_Length = 30;%25:5:40; %Antenna Length in inches
span = 5;%3:1:6; %Span is a range of values, in feet
load("MotorSpreadsheet.mat");
Num_Power_Systems = height(MotorSpreadsheet);

VSAR = 2; %Vertical stabalizer aspect ratio - derived from aero calculations done beforehand
HSAR = 4.5; %Horizontal stab aspect ratio

width_wheel = 0.5; %width of wheels (in)
radius_wheel = 1.5; %wheel radius
sa_wheel = (2*pi*(radius_wheel)^2+ pi*2*radius_wheel*width_wheel)/144;

[wings] = wingData(Aspect_Ratios, span); %call wing function to make airfoil data lookup table
[wingrow, wingcol, wingpg, wingspan] = size(wings); %get indices to iterate over. Must also include size of wingspans this year

max_index = 100000; %roughly 15% of airplanes checked are succesful so preallocate enough memory for them. Dramatically speeds up computation.
plane(1:max_index) = struct(airplaneClass);%create a matrix to hold all the computed aircraft. Most aircraft will fail and be overwritten so max_index does not have to equal max iterations.
index = 1;
iterNum = 1;
spanfailcount = zeros([length(span) 6]); %This creates a matrix to check which failure conditions are most prevailent at each span value
% spanfailcount(1, 1) = "span";
% spanfailcount(1, 2) = 'space';
% spanfailcount(1, 3) = 'ep weight';
% spanfailcount(1, 4) = 'takeoff';
% spanfailcount(1, 5) = 'moment';
% spanfailcount(1, 6) = 'converge';
for AR = 1:wingpg
    disp("At Aspect ratio " + AR);
    toc;
    for spanIndex = 1:length(span) %NOTE: the wings function cant handle this iteration yet, all calls of plane(index) also need updated, waiting on wings update
        spanfailcount(1 + spanIndex, 1) = span(spanIndex);
        for airfoil = 1:1:8
            for powerIndex = 1:Num_Power_Systems
                for electronicPackageIndex = 1:length(Electronic_Package_Weight)
                    for antennaIndex = 1:length(Antenna_Length) %Changed this for the rerun - use smarter logic otherwise %for every amount of syringes try up to the maximum number of vials

                        plane(index) = struct(airplaneClass); %make sure to start with a clean slate at this index as this code writes over the index of failed airplanes. Without cleanup things like failure flags stay set even when they shouldn't

                        plane(index).fuselage.wheelSA = sa_wheel;

                        %load starting values for each plane
                        plane(index).wing.span = span(spanIndex);
                        plane(index).performance.epWeight = Electronic_Package_Weight(electronicPackageIndex);
                        plane(index).wing.aspectRatio = Aspect_Ratios(AR); %ratio between length and width of wing.
                        plane(index).performance.antennaLength = Antenna_Length(antennaIndex);

                        %read values from wings matrix into aircraft
                        %properties. See wingClass.m for property descriptions
                        plane(index).wing.clw = wings(airfoil, 1, AR, spanIndex);
                        plane(index).wing.clm = wings(airfoil, 2, AR, spanIndex);     %cl max
                        plane(index).wing.cd = wings(airfoil, 3, AR, spanIndex);     %cd i zero velocity coefficient of drag
                        plane(index).wing.clFlap = wings(airfoil, 4, AR, spanIndex);  %weight
                        plane(index).wing.weight = wings(airfoil, 5, AR, spanIndex);    %airfoil name
                        plane(index).wing.chord = wings(airfoil, 6, AR, spanIndex);    %airfoil name
                        plane(index).wing.planformArea = wings(airfoil, 7, AR, spanIndex);    %airfoil name
                        plane(index).wing.surfaceArea = wings(airfoil, 8, AR, spanIndex);    %airfoil name
                        plane(index).wing.name = wings(airfoil, 9, AR, spanIndex);    %airfoil name
                        plane(index).wing.thickness=wings(airfoil, 10, AR, spanIndex); %thickness of the wing (ft)


                        %load values from power system table into aircraft
                        plane(index) = powerSelections(plane(index), MotorSpreadsheet, powerIndex);

                        %set payload and fuselage configuration
                        %plane(index).fuselage.numSyringes = syringes(syringe_index);
                        %plane(index).fuselage.numVials = num_vials;
                        plane(index) = fuselage(plane(index)); %calculate fuselage size based on number of vials and syringes
                        plane(index) = landingGear(plane(index), rho);
                        plane(index) = empennage(plane(index), HSAR, VSAR);

                        plane(index) = findTotalWeight(plane(index), Electronic_Package_Weight(electronicPackageIndex), Antenna_Length(antennaIndex));

                        %Want to include sanity check here that throws
                        %planes out which don't meet space/weight requirements
                        plane(index) = volSanityCheck(plane(index), Electronic_Package_Weight(electronicPackageIndex));
                        if plane(index).volSanityFlag == 0
                            %disp("Too Big!");
                            if plane(index).epFail
                                spanfailcount(1 + spanIndex, 3) = spanfailcount(1 + spanIndex, 3) + 1;
                            elseif plane(index).spaceFail
                                spanfailcount(1 + spanIndex, 2) = spanfailcount(1 + spanIndex, 2) + 1;
                            end
                            break;
                        end

                        %simulate mission 2
                        plane(index) = GenVelocityTest(plane(index), 2, rho, Temp); %2 signifies mission 2 configuration
                        plane(index) = TakeoffChecker(plane(index), 2, rho);
                        plane(index) = mission2score(plane(index), Electronic_Package_Weight(electronicPackageIndex));

                        %simulate mission 3
                        plane(index) = GenVelocityTest(plane(index), 3, rho, Temp); %2 signifies mission 2 configuration
                        plane(index) = TakeoffChecker(plane(index), 3, rho);
                        plane(index) = mission3score(plane(index), Antenna_Length(antennaIndex));

                        plane(index) = sanityCheck(plane(index), span(spanIndex), Antenna_Length(antennaIndex)); %make sure all the calculated values make sense and meet
                        %competition requirements. needs to be updated for
                        %this year's competition

                        if plane(index).sanityFlag
                            index = index + 1;
                            %disp("Sane!");
                        else
                            if plane(index).takeoffFail
                                spanfailcount(1 + spanIndex, 4) = spanfailcount(1 + spanIndex, 4) + 1;
                            elseif plane(index).momentFail
                                spanfailcount(1 + spanIndex, 5) = spanfailcount(1 + spanIndex, 5) + 1;
                            elseif plane(index).convergeFail
                                spanfailcount(1 + spanIndex, 6) = spanfailcount(1 + spanIndex, 6) + 1;
                            end
                            %disp("insane!");
                            spanfailcount(spanIndex,3) = spanfailcount(spanIndex,3) + 1;
                        end
                        iterNum = iterNum+1;
                    end
                end
            end
        end
    end
end
toc;
%Code below is from post.m file to run post-processing on the data. We should probably move this into a function call later.
%THIS HAS NOT BEEN UPDATED

score2 = zeros(1, length(plane)); %mission 2 scores of all aircraft
score3 = score2; %mission 3 scores
%scoreg = score2; %ground scores
%vials = score2; %number of vials in given airplane
%syringes = score2; %number of syringes in given airplane
for i = 1:length(plane) %load data from structure into arrays that are easier to work with
    if (plane(i).sanityFlag == 1) || (plane(i).volSanityFlag == 1)  %all values remain 1 unless the airplane fails sanity check. If the last plane analyzed is not sane the original sanity check can't throw it out.
        score2(i) = (plane(i).performance.score2);
        score3(i) = (plane(i).performance.score3);
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
winners = plane(I);
scores = score(I);
%save("winnersAR" + Aspect_Ratios(1) + ".mat", "winners"); %save top 100 airplanes. It is impractical to save all airplanes checked.
%save("winnersAR" + Aspect_Ratios(1) + "_scores.mat", "scores");

%save("successes.mat", "plane");
%clear; %once finished don't keep hogging ram. Previous experience with running multiple instances on analysis on one computer shows that ram usage is the first limiting factor in enabling
%the analysis to run, so clearing it as often as possible is important to avoid hogging ram from other programs. This ram hogging is the leading cause of crashing for this code.
%scatter(vials,score);

