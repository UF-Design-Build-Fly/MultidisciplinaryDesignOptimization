x = zeros(1, length(planes));
y = zeros(1, length(planes));

for (i = 1:length(planes)) %Load data into arrays that are easier to work with
	x(i) = i;
	y(i) = planes(i).performance.scoreTotal;
end

plot(x, y, '.')
xlabel('Number of Passengers')
ylabel('Score')
%zlabel('Total Score')
grid on