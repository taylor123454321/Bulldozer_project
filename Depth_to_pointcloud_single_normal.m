clear, clc, close all

depthDevice = imaq.VideoDevice('kinect',2);

step(depthDevice);

for i = 1:20
    depthImage = step(depthDevice);
end
tic
taken = 31;
filter = 10;

clear('depthImage')
for i = 1:taken
    image = step(depthDevice);
    %     image(image<600) = 0; %trimming minium and maxium length
    image(image>800) = 0;
    image(:,1:20) = 0;
    image(:,424-20:424) = 0;
    image(1:20,:) = 0;
    image(:,512-20:512) = 0;
    cat = image;
end


% sum_ = uint16(zeros(424,512));
% elements = uint16(zeros(424,512));
% dImage_filtered = uint16(zeros(424,512,taken-filter));
% 
% for b = 1:filter
%     sum_(:,:) = sum_(:,:) + depthImage(:,:,b);
%     for i = 1:424
%         for j = 1:512
%             if depthImage(i,j,b) ~= 0
%                 elements(i,j) = elements(i,j) + 1;
%             end
%         end
%     end
% end
% 
% for b = filter+1:taken
%     sum_(:,:) = sum_(:,:) - depthImage(:,:,b-filter) + depthImage(:,:,b);
%     for i = 1:424
%         for j = 1:512
%             if depthImage(i,j,b-filter) ~= 0
%                 elements(i,j) = elements(i,j) - 1;
%             end
%             if depthImage(i,j,b) ~= 0
%                 elements(i,j) = elements(i,j) + 1;
%             end
%         end
%     end
%     dImage_filtered(:,:,b-filter) = sum_(:,:)./elements(:,:);
% end
ptCloud = pcfromkinect(depthDevice,cat);%dImage_filtered(:,:,1));
pcshow(ptCloud)

release(depthDevice);
tic
normals = pcnormals(ptCloud,30);

toc

x = ptCloud.Location(1:10:end,1:10:end,1);
y = ptCloud.Location(1:10:end,1:10:end,2);
z = ptCloud.Location(1:10:end,1:10:end,3);
u = normals(1:10:end,1:10:end,1);
v = normals(1:10:end,1:10:end,2);
w = normals(1:10:end,1:10:end,3);

sensorCenter = [0,0,0.3];
for k = 1 : numel(x)
    p1 = sensorCenter - [x(k),y(k),z(k)];
    p2 = [u(k),v(k),w(k)];
    % Flip the normal vector if it is not pointing towards the sensor.
    angle = atan2(norm(cross(p1,p2)),p1*p2');
    if angle > pi/2 || angle < -pi/2
        u(k) = -u(k);
        v(k) = -v(k);
        w(k) = -w(k);
    end
end


xx = x;yy = y;zz = z;
for i = 1:43
    for j = 1:52
        xx(i,j) = 0;
        yy(i,j) = 0;
        zz(i,j) = 0;
    end
end

pcshow(ptCloud)
title('Adjusted Normals of Point Cloud')
hold on
quiver3(xx, yy, zz, u, v, w);

A = zeros(132,3);
F = zeros(132,3);

k = 1;
for i = 1:43
    for j = 1:52
        if ~isnan(u(i,j))
            if ~isnan(v(i,j))
                if ~isnan(w(i,j))
                    A(k,1) = u(i,j);
                    A(k,2) = v(i,j);
                    A(k,3) = w(i,j);
                    F(k,1) = x(i,j);
                    F(k,2) = y(i,j);
                    F(k,3) = z(i,j);
                    k = k + 1;
                end
            end
        end
    end
end

leng = length(A);

uu = A(:,1);
vv = A(:,2);
ww = A(:,3);
% xx = zeros(leng,1);
% yy = zeros(leng,1);
% zz = zeros(leng,1);
xx = F(:,1);
yy = F(:,2);
zz = F(:,3);

B = eye(leng + 2);
k = 2;

for i = 1:leng
    B(k,k) = A(i,2);
    B(k,k-1) = A(i,1);
    B(k,k+1) = A(i,3);
    k = k + 1;
end

%quiver3(xx, yy, zz, uu, vv, ww);
%quiver(xx, yy, uu, vv);

grid on












