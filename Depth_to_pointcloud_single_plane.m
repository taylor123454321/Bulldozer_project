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

%pcshow(ptCloud)

maxDistance = 0.01;
referenceVector = [0,1,0.3];
maxAngularDistance = 0.8;

[model1,inlierIndices,outlierIndices] = pcfitplane(ptCloud,maxDistance,referenceVector,maxAngularDistance);
plane1 = select(ptCloud,inlierIndices);
remainPtCloud = select(ptCloud,outlierIndices);

% roi = [-inf,inf;0.4,inf;-inf,inf];
% sampleIndices = findPointsInROI(ptCloud,roi);

% [model2,inlierIndices,outlierIndices] = pcfitplane(remainPtCloud,maxDistance,'SampleIndices',sampleIndices);
% plane2 = select(remainPtCloud,inlierIndices);
% remainPtCloud = select(remainPtCloud,outlierIndices);

figure
pcshow(plane1)
title('First Plane')

% figure
% pcshow(plane2)
% title('Second Plane')

figure
pcshow(remainPtCloud)
title('Remaining Point Cloud')





