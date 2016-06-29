clear,clc, close all

load('angle_results/angle_data_8_3')

angles(1,2,1) = angles(1,2,2);

for i = 1:length(angles(1,3,:))
    a(i) = angles(1,1,i);
    b(i) = angles(1,2,i);
    c(i) = angles(1,3,i);
    if i > 61
        stepp(i) = 10.7;
    else
        stepp(i) = 7.75;
    end
end

figure
% plot(a)
hold on
plot(b)
% plot(c)
plot(stepp)
grid on
title('Step input, buffer size 3')
xlabel('Frames')
ylabel('Rotation (degrees)')
legend('Output','Step input')


% Y = fft(a,251);
% 
% Pyy = Y.*conj(Y)/251;
% f = 1000/251*(0:127);
% hold off
% figure
% plot(f,Pyy(1:128))
% title('Power spectral density')
% xlabel('Frequency (Hz)')

load('save_times_1.mat')

time(1) = sum(pic_time(6:40))/(40-6);
time(2) = sum(suf_time(6:40))/(40-6);
time(3) = sum(ore_time(6:40))/(40-6);
time(4) = sum(matches_time(6:40))/(40-6);
time(5) = sum(dd_time(6:40))/(40-6);
time(6) = sum(match_time(6:40))/(40-6);
time(7) = sum(orginal_time(6:40))/(40-6);
time(8) = sum(angle_cal_time(6:40))/(40-6);
time(9) = sum(filter_time(6:40))/(40-6);
time(10) = sum(filter_angle_time(6:40))/(40-6);

time

overall = sum(time)
freq = 1/overall




