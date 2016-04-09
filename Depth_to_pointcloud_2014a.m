clear, clc, close all

depthDevice = imaq.VideoDevice('kinect',2)

step(depthDevice);


depthImage = step(depthDevice);


xyzPoints = depthToPointCloud(depthImage,depthDevice);


pcshow(xyzPoints,'VerticalAxis','y','VerticalAxisDir','down');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

release(depthDevice);