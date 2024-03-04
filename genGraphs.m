for (i = 1:length(planes)) %Load data into arrays that are easier to work with
    x(i) = i;
    y(i) = planes(i).performance.score2 + planes(i).performance.score3 + planes(i).performance.scoreGM;
end

plot(x, y, '.')
xlabel('Number of Passengers')
ylabel('Score')
%zlabel('Total Score')
grid on