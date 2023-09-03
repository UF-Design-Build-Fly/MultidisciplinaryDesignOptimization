%Score visualizer
% speed = zeros(1,200);
% for i = 1:length(winners)
%     speed(i) = winners(i).performance.score3;
% end
% 
% plot(1:length(winners), speed)
% 
% %spreadsheet data cleaner
% load('MotorSpreadsheet8sComplete.mat')
% indices = [0];
% j = 1;
% for i = 1:height(MotorSpreadsheet)
% %     if(table2array(MotorSpreadsheet(i,'cells')) > max)
% %         max = table2array(MotorSpreadsheet(i,'cells'));
% %     end
%     if(table2array(MotorSpreadsheet(i,'cells')) >= 8)
%           indices(j) = i;
%           j = j + 1;
%     end
% 
% end

% 
score2 = zeros(length(winners),1);
score3 = zeros(length(winners),1);
for i = 1:length(winners)
    score2(i) = winners(i).performance.score2;
    score3(i) = winners(i).performance.score3;
end
[M2, I2] = max(score2);
[M3, I3] = max(score3);
score2 = 1 + score2/M2;
score3 = 2 + score3/M3;
scores = score2 + score3 + 1;

% for i = 1:length(winners)
%     winners(i).performance.score1 = score(i);
% end
% 
% [M, I] = maxk(score,200);
% newScore = score(I);
% 
% close all;
% plot(score)
% figure
% plot(newScore)
close all;
% scores = zeros(1,length(winners));
% for i = 1:length(winners)
%     scores(i) = winners(i).performance.score3;
%     %disp(winners(i).power.motorName)
% end
plot(maxk(scores,200),LineWidth=1)
title("Concave Property of Aircraft Configuration vs. Estimated Total Mission Performance", FontSize=14)
xlabel("Aircraft Rank", FontSize=16)
ylabel("Scaled Score Estimate", FontSize=16)