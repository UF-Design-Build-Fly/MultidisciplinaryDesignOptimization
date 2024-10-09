%% =================================== Cleanup Environment ==================================
clc;					% Clears the Command Window
clear;					% Clears all variables
close all;				% Closes are figure windows
warning('off','all');	% Using structs the way we do here generates warnings that slow matlab down. Comment this out (and restart matlab) when debugging.

%% ==================================== Graph Parameters ====================================

wingClMax = 1.2;
wingSpan = 70/39.37; %m
wingAR = 6;

rho = 1.2127; %kg/m^3
windSpeed = 4.02336; %m/s
turnSpeedMultiplier = 0.8;

maxLift = @(airspeed) 1/2 * rho * (airspeed*turnSpeedMultiplier).^2 * (wingSpan^2/wingAR) * wingClMax; %N
turnAcc = @(airspeed, weight) maxLift(airspeed)./weight; %m/s^2
turnRad = @(airspeed, weight) airspeed.^2./turnAcc(airspeed, weight); %m

turnTime = @(airspeed, weight) 2*2*pi*turnRad(airspeed, weight)./(airspeed*turnSpeedMultiplier);
downWindTime = @(airspeed, windSpeed)(304.8./(airspeed+windSpeed));
upWindTime = @(airspeed, windSpeed)(304.8./(airspeed-windSpeed));
lapTime = @(airspeed, windSpeed, weight) (turnTime(airspeed, weight) + downWindTime(airspeed, windSpeed) + upWindTime(airspeed, windSpeed));

m2ScoreEqu = @(airspeed, emptyWeight, payloadWeight, windSpeed) (payloadWeight./(3*lapTime(airspeed, windSpeed, emptyWeight+payloadWeight)));

m3laps = @(airspeed, windSpeed, weight, glideTime) floor((300-glideTime)./lapTime(airspeed, windSpeed, weight));
m3bonus = @(bonusBox, gliderWeight) bonusBox./gliderWeight;
m3ScoreEqu = @(airspeed, windSpeed, emptyWeight, glideTime, bonusBox, gliderWeight) m3laps(airspeed, windSpeed, emptyWeight, glideTime) + m3bonus(bonusBox, gliderWeight);


%% ====================================== Gen M2 Graphs =====================================
% Parameters: Plane Speed, Empty Weight, Maximum Weight
basePlaneSpeed = 32; %m/s
planeSpeeds = linspace(0.7*basePlaneSpeed, 1.3*basePlaneSpeed, 10000);
psPercents = (planeSpeeds-basePlaneSpeed)/basePlaneSpeed;

baseEmptyWeight = 3.63; %kg
emptyWeights = linspace(0.7*baseEmptyWeight, 1.3*baseEmptyWeight, 10000);
ewPercents = (emptyWeights-baseEmptyWeight)/baseEmptyWeight;

basePayloadWeight = 3.63; %kg
payloadWeights = linspace(0.7*basePayloadWeight, 1.3*basePayloadWeight, 10000);
pwPercents = (payloadWeights-basePayloadWeight)/basePayloadWeight;

baseM2Score = m2ScoreEqu(basePlaneSpeed, baseEmptyWeight, basePayloadWeight, windSpeed);

%Vary Plane Speeds
psM2Scores = m2ScoreEqu(planeSpeeds, baseEmptyWeight, basePayloadWeight, windSpeed);
psM2ScoresPercents = (psM2Scores-baseM2Score)/baseM2Score;

%Vary Empty Weight
ewM2Scores = m2ScoreEqu(basePlaneSpeed, emptyWeights, basePayloadWeight, windSpeed);
ewM2ScoresPercents = (ewM2Scores-baseM2Score)/baseM2Score;

%Vary Payload Weight
pwM2Scores = m2ScoreEqu(basePlaneSpeed, baseEmptyWeight, payloadWeights, windSpeed);
pwM2ScoresPercents = (pwM2Scores-baseM2Score)/baseM2Score;

figure;
hold on;
grid on;
plot(psPercents*100, psM2ScoresPercents*100);
plot(ewPercents*100, ewM2ScoresPercents*100);
plot(pwPercents*100, pwM2ScoresPercents*100);
xlim([-30, 30]);
legend("Plane Speed", "Empty Weight", "Payload Weight", "Location", "best");
xlabel("% Change in Parameter");
ylabel("% Change in M2 Score");
title("M2 Score Sensitivity Analysis");


%% ====================================== Gen M3 Graphs =====================================
% Parameters: Plane Speed, Empty Weight, Glide Time, Glider Weight
baseGlideTime = 60; %s
glideTimes = linspace(0.7*baseGlideTime, 1.3*baseGlideTime, 10000);
gtPercents = (glideTimes-baseGlideTime)/baseGlideTime;

baseGliderWeight = 0.423; %lbs (chosen so that upper bound is 0.55 lbs)
gliderWeights = linspace(0.7*baseGliderWeight, 1.3*baseGliderWeight, 10000);
gwPercents = (gliderWeights-baseGliderWeight)/baseGliderWeight;

bonusBox = 2.5;

baseM3Score = m3ScoreEqu(basePlaneSpeed, windSpeed, baseEmptyWeight, baseGlideTime, bonusBox, baseGliderWeight);

%Vary Plane Speeds
psM3Scores = m3ScoreEqu(planeSpeeds, windSpeed, baseEmptyWeight, baseGlideTime, bonusBox, baseGliderWeight);
psM3ScoresPercents = (psM3Scores-baseM3Score)/baseM3Score;

%Vary Empty Weight
ewM3Scores = m3ScoreEqu(basePlaneSpeed, windSpeed, emptyWeights, baseGlideTime, bonusBox, baseGliderWeight);
ewM3ScoresPercents = (ewM3Scores-baseM3Score)/baseM3Score;

%Vary Glide Time
gtM3Scores = m3ScoreEqu(basePlaneSpeed, windSpeed, baseEmptyWeight, glideTimes, bonusBox, baseGliderWeight);
gtM3ScoresPercents = (gtM3Scores-baseM3Score)/baseM3Score;

%Vary Glider Weight
gwM3Scores = m3ScoreEqu(basePlaneSpeed, windSpeed, baseEmptyWeight, baseGlideTime, bonusBox, gliderWeights);
gwM3ScoresPercents = (gwM3Scores-baseM3Score)/baseM3Score;

figure;
hold on;
grid on;
plot(psPercents*100, psM3ScoresPercents*100);
plot(ewPercents*100, ewM3ScoresPercents*100);
plot(gtPercents*100, gtM3ScoresPercents*100);
plot(gwPercents*100, gwM3ScoresPercents*100);
xlim([-30, 30]);
legend("Plane Speed", "Empty Weight", "Glide Time", "Glider Weight", "Location", "best");
xlabel("% Change in Parameter");
ylabel("% Change in M3 Score");
tit = sprintf("M3 Score Sensitivity Analysis - Bonus Box: %g", bonusBox);
title(tit);
