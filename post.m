warning('off', 'all');

plane1 = load("plane1.mat");
plane1 = plane1.plane;

plane2 = load("plane2.mat");
plane2 = plane2.plane;

plane3 = load("plane3.mat");
plane3 = plane3.plane;

plane4 = load("plane4.mat");
plane4 = plane4.plane;

allPlanes = [plane1, plane2, plane3, plane4];
score2 = zeros(1, length(allPlanes)); %mission 2 score
score3 = score2; %mission 3 score
scoreg = score2; %ground score
vials = score2; %number of vials in given airplane
syringes = score2; %number of syringes in given airplane
for i = 1:length(allPlanes) %load data from structure into arrays that are easier to work with
    if(allPlanes(i).sanityFlag == 1)
        score2(i) = (allPlanes(i).performance.score2);
        score3(i) = floor((allPlanes(i).performance.score3));
        vials(i) = allPlanes(i).fuselage.numVials;
        syringes(i) = allPlanes(i).fuselage.numSyringes;
        scoreg(i) = 10 + 2*(3*syringes(i)/5) + vials(i)*5;%run, load and unload syringes, load vial       
    end
end
[M2, I2] = max(score2);
[M3, I3] = max(score3);
[Mg, Ig] = max(scoreg);
score2 = 1 + score2/M2;
score3 = 2 + score3/M3;
scoreg = scoreg/Mg;
score = score2 + score3 + scoreg + 1;
[M, I] = max(score);
scatter(vials,score);