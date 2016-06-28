clear,clc, close all

load('angle_results/angle_data_7')
angles(1,1,1) = angles(1,1,2);

for i = 1:length(angles(1,3,:))
    a(i) = angles(1,1,i);
    b(i) = angles(1,2,i);
    c(i) = angles(1,3,i);
    if i > 50
        stepp(i) = 10.2;
    else
        stepp(i) = 6.4;
    end
end

figure
plot(a)
hold on
% plot(b)
% plot(c)
plot(stepp)
grid on
title('Step input, no buffer')
xlabel('Frames')
ylabel('Rotation (degrees)')
legend('Output','Step input')






