function [thrust] = dynamicThrust(diameter, pitch, rpm, velocity, net, stats)
%dynamicThrust: Use neural network trained on apc dynamic thrust data to
%effectively predict thrust output given a propeller and forward velocity
%Expect to take v in feet per second. Model is trained in miles per hour
velocity = velocity/1.467; %convert to mph

velocity = (velocity - stats(1))/stats(2); %normalize by mean and std. Neural networks require normalized data
diameter = (diameter - stats(3))/stats(4);
pitch = (pitch - stats(5))/stats(6);
rpm = (rpm - stats(7))/stats(8);

thrust = predict(net, [velocity, diameter, pitch, rpm]);

end