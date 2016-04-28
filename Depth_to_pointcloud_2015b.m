clear, clc, close all

%colorDevice = imaq.VideoDevice('kinect',1)
depthDevice = imaq.VideoDevice('kinect',2);


%step(colorDevice);
step(depthDevice);

%colorImage = step(colorDevice);
depthImage = step(depthDevice);

ptCloud = pcfromkinect(depthDevice,depthImage);

player = pcplayer([-1.5    1.5],[-1.5    1.5],[0.5 2.5],...
    'VerticalAxis','y','VerticalAxisDir','down');

xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');


for i = 1:500
    %colorImage = step(colorDevice);
    depthImage = step(depthDevice);
%     depthImage(depthImage<800) = 0; %trimming minium and maxium length
    depthImage(depthImage>800) = 0;
    ptCloud = pcfromkinect(depthDevice,depthImage);

    view(player,ptCloud);
end


%release(colorDevice);
release(depthDevice);

















