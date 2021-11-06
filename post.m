score2 = zeros(1, length(plane)); %mission 2 score
score3 = score2; %mission 3 score
scoreg = score2; %ground score
vials = score2; %number of vials in given airplane
syringes = score2; %number of syringes in given airplane
for i = 1:length(plane) %load data from structure into arrays that are easier to work with
    if(plane(i).sanityFlag == 1) %all values remain zero unless the airplane passes sanity check. If the last plane analyzed is not sane the original sanity check can't throw it out.
        score2(i) = (plane(i).performance.score2);
        score3(i) = floor((plane(i).performance.score3));
        vials(i) = plane(i).fuselage.numVials;
        syringes(i) = plane(i).fuselage.numSyringes;
        scoreg(i) = 10 + 2*(3*syringes(i)/5) + 5*vials(i);%run, load and unload syringes, load vial
    end
end
clear plane; %reduce ram usage asap
[M2, I2] = max(score2); %find the airplanes with the best individual mission scores
[M3, I3] = max(score3);
[Mg, Ig] = max(scoreg);
score2 = 1 + score2/M2; %normalize scores against best performers
score3 = 2 + score3/M3;
scoreg = scoreg/Mg;
score = score2 + score3 + scoreg + 1;
[M, I] = maxk(score, 100); %find the top 100 airplanes
winners = plane(I);
save("winnersAR" + Aspect_Ratios(1) + ".mat", "winners"); %save top 100 airplanes. It is impractical to save all airplanes checked.
clear; %once finished don't keep hogging ram. Previous experience with running multiple instances on analysis on one computer shows that ram usage is the first limiting factor in enabling
       %the analysis to run, so clearing it as often as possible is important to avoid hogging ram from other programs. This ram hogging is the leading cause of crashing for this code.
%scatter(vials,score);