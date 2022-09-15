[wings] = wing;
[wingrow, wingcol, wingpg] = size(wings);

for AR = 1:wingpg
    for af = 1:wingrow
        plane(index).wing.cl = wings(af, 1, AR);
        plane(index).wing.clm = wings(af, 2, AR);     %cl max
        plane(index).wing.cd0 = wings(af, 3, AR);     %cd i zero velocity coefficient of drag
        plane(index).wing.weight = wings(af, 4, AR);       %weight
        plane(index).wing.name = wings(af, 5, AR);    %airfoil name
		plane(index).wings.aspectRatio = Aspect_Ratios(AR); %ratio between length and width of wing.
    end
end

       
   