clear, clc, close all


ptCloud = pcread('box_moved.ply');

normals = pcnormals(ptCloud);
pcshow(ptCloud);
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
%         if x(i,j) > 1.5
%             x(i,j) = NaN;
%         end
%         if y(i,j) > 1.5
%             y(i,j) = NaN;
%         end
%         if z(i,j) > 1.5
%             z(i,j) = NaN;
%         end
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





