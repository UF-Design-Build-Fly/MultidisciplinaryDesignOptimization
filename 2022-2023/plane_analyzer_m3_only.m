clear; clc; clf; close all;
warning('off','all')
net = importKerasNetwork("Fitting_3_model.h5");
vmean = 49.708;
vstd = 45.637;
dmean = 10.523;
dstd = 4.487;
pmean = 7.308001;
pstd = 3.348207;
rpmmean = 12508.699;
rpmstd = 9050.014;

stats = [vmean, vstd, dmean, dstd, pmean, pstd, rpmmean, rpmstd];



%code to predict performance for one single aircraft
%Used to validate model against real-world performance

%Motor prop we used is index 2044

%selig7062 is first airfoil in wingData output
%span is 5 feet
%chord is 8.93 inches


rho = 0.00235308; %Density air at Tuscon, Az (slug/ft^3)
Temp = 293; %temperature in kelvin at competition site

Aspect_Ratios = 6;
AR = 1;

Electronic_Package_Weight = 0;%[0,2:.5:8]; %Electronic Package Weight for mission 2 in pounds

Antenna_Length = 0:5:40; %Antenna length, inches

span = 4.5;
spanIndex = 1;

load("MotorSpreadsheet8sComplete.mat");
Num_Power_Systems = height(MotorSpreadsheet);

VSAR = 3; %Vertical stabalizer aspect ratio - derived from aero calculations done beforehand
HSAR = 4.5; %Horizontal stab aspect ratio

width_wheel = 0.5; %width of wheels (in)
radius_wheel = 1.5; %wheel radius
sa_wheel = (2*pi*(radius_wheel)^2+ pi*2*radius_wheel*width_wheel)/144;

[wings] = wingData(Aspect_Ratios, span); %call wing function to make airfoil data lookup table
[wingrow, wingcol, wingpg, wingspan] = size(wings); %get indices to iterate over. Must also include size of wingspans this year

max_index = 10000;
plane(1:max_index) = struct(airplaneClass);


airfoil = 1;
index = 1;
iterNum = 0;
for powerIndex = 1:Num_Power_Systems
    for electronicPackageIndex = 1:length(Electronic_Package_Weight)
        for antennaIndex = 1:length(Antenna_Length)
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

            plane(index) = fuselage(plane(index)); %calculate fuselage size based on number of vials and syringes
            plane(index) = landingGear(plane(index), rho);
            plane(index) = empennage(plane(index), HSAR, VSAR);

            plane(index) = findTotalWeight(plane(index), Electronic_Package_Weight(electronicPackageIndex), Antenna_Length(antennaIndex));

            %Want to include sanity check here that throws
            %planes out which don't meet space/weight requirements
            plane(index) = volSanityCheck(plane(index), Electronic_Package_Weight(electronicPackageIndex));


            %simulate mission 2
            plane(index).performance.velocity2 = 0;
            %plane(index) = GenVelocityTest(plane(index), 2, rho, Temp, net, stats); %2 signifies mission 2 configuration
            %plane(index) = TakeoffChecker(plane(index), 2, rho);
            %plane(index) = mission2score(plane(index), Electronic_Package_Weight(electronicPackageIndex));

            %plane(index).performance.velocity3 = 0; %satisfy the sanity checker
            %simulate mission 3
            plane(index) = GenVelocityTest(plane(index), 3, rho, Temp, net, stats); %2 signifies mission 2 configuration
            plane(index) = TakeoffChecker(plane(index), 3, rho);
            plane(index) = mission3score(plane(index), Antenna_Length(antennaIndex));

            plane(index) = sanityCheck(plane(index), span(spanIndex), Antenna_Length(antennaIndex)); %make sure all the calculated values make sense and meet
            %competition requirements. needs to be updated for
            %this year's competition

            if plane(index).sanityFlag
                index = index + 1;
                %disp("Sane!");
            end
            iterNum = iterNum+1;
        end
    end
end

score2 = zeros(1, length(plane)); %mission 2 scores of all aircraft
score3 = score2; %mission 3 scores

for i = 1:length(plane) %load data from structure into arrays that are easier to work with
    if (plane(i).sanityFlag == 1) || (plane(i).volSanityFlag == 1)  %all values remain 1 unless the airplane fails sanity check. If the last plane analyzed is not sane the original sanity check can't throw it out.
        %score2(i) = (plane(i).performance.score2);
        score3(i) = (plane(i).performance.score3);
    end
end

%[M2, I2] = max(score2); %find the airplanes with the best individual mission scores
[M3, I3] = max(score3);

%score2 = 1 + score2/M2; %normalize scores against best performers
score3 = 2 + score3/M3;

score = score3;% + score3 + 1; %scoreg + 1;
[M, I] = maxk(score, 200); %find the top 200 airplanes
winners = plane(I);
save("winners_m3.mat", "winners"); %save top 100 airplanes. It is impractical to save all airplanes checked.
