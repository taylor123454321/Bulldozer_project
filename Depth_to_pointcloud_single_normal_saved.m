clear, clc, close all

depthDevice = imaq.VideoDevice('kinect',2);

step(depthDevice);

for i = 1:20
    depthImage = step(depthDevice);
end

ptCloud = pcfromkinect(depthDevice,depthImage);

player = pcplayer([-1.5    1.5],[-1.5    1.5],[0.5 2.5],...
    'VerticalAxis','y','VerticalAxisDir','down');

xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');


%for i = 1:500
depthImage = step(depthDevice);

for i = 1:424;
    for j = 1:512;
        if depthImage(i,j) > 1700
            depthImage(i,j) = 0;
        end
        if depthImage(i,j) < 800
            depthImage(i,j) = 0;
        end
    end
end
ptCloud = pcfromkinect(depthDevice,depthImage);


release(depthDevice);

normals = pcnormals(ptCloud);
%view(player,ptCloud);
%end


x = ptCloud.Location(1:10:end,1:10:end,1);
y = ptCloud.Location(1:10:end,1:10:end,2);
z = ptCloud.Location(1:10:end,1:10:end,3);
u = normals(1:10:end,1:10:end,1);
v = normals(1:10:end,1:10:end,2);
w = normals(1:10:end,1:10:end,3);

sensorCenter = [0,-0.3,0.3];
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

% for i = 1:43
%     for j = 1:52
%         if x(i,j) < -0.4 || x(i,j) > 0.4
%             x(i,j) = NaN;
%         end
%         if y(i,j) > 0.33
%             y(i,j) = NaN;
%         end
% %         if z(i,j) > 1.5
% %             z(i,j) = NaN;
% %         end
%     end
% end

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
%surf(xx,yy,zz)
%view(0,-45)


save('normal_vectors','xx','yy','zz','u','v','w');




