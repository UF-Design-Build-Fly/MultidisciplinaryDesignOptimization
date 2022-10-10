function [plane] = mission3score(plane, Antenna_Length)
%M3 score = 2 + [N_(antenna length / mission time) / Max_(antenna length / mission time)]
%If mission (3 laps) takes longer than 5 minutes, sanity check will throw
%out plane
lapdist = (500*4)+(500*pi)+(250*pi); %estimates the overall lap distance in ft
time3 = (3*(lapdist/plane.performance.velocity3))/60; %Calculates the time taken to fly 3 laps 
plane.performance.score3 = Antenna_Length/time3; %Calculates mission 3 score

end