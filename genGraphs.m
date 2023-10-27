x = zeros(1, index-1); %Initilize arrays with 0s
y = scoresM2; %Initilize arrays with 0s

for (i = 1:index-1) %Load data into arrays that are easier to work with
    x(i) = planes(i).performance.numPassengers;
    y(i) = planes(i).wing.span;
end

plot(x, score, '.')
xlabel('Number of Passengers')
ylabel('Score')
%zlabel('Total Score')
grid on