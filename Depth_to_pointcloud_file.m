clear, clc, close all

depthDevice = imaq.VideoDevice('kinect',2);

step(depthDevice);

depthImage = step(depthDevice);

ptCloud = pcfromkinect(depthDevice,depthImage);

depthImage = step(depthDevice);

for i = 1:424;
    for j = 1:512;
        if depthImage(i,j) > 1300
            depthImage(i,j) = 0;
        end
        if depthImage(i,j) < 800
            depthImage(i,j) = 0;
        end
    end
end
ptCloud = pcfromkinect(depthDevice,depthImage);

release(depthDevice);

pcwrite(ptCloud,'box_moved','PLYFormat','binary');
