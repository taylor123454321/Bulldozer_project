clear, clc, close all

depthDevice = imaq.VideoDevice('kinect',2);

step(depthDevice);

for i = 1:20
    depthImage = step(depthDevice);
end

%ptCloud = pcfromkinect(depthDevice,depthImage);

player = pcplayer([-1.5    1.5],[-1.5    1.5],[0.5 2.5],...
    'VerticalAxis','y','VerticalAxisDir','down');

xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');



for i = 1:5
    depthImage(:,:,i) = step(depthDevice);
end


for i = 1:424
    for j = 1:512
        for b = 1:5
            if depthImage(i,j,b) == 0;
                depthImage(i,j,:) = 0;
            end
        end
    end
end

for b = 1:5
    for i = 1:424
        for j = 1:512
            if depthImage(i,j,b) > 1300
                depthImage(i,j,b) = 0;
            end
            if depthImage(i,j,b) < 800
                depthImage(i,j,b) = 0;
            end
        end
    end
end

dImage = uint16(zeros(424,512));

for b = 1:5
    for i = 1:424
        for j = 1:512
            dImage(i,j) = dImage(i,j) + depthImage(i,j,b);
        end
    end
end

ptCloud = pcfromkinect(depthDevice,dImage);

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

A = zeros(132,3);

k = 1;
for i = 1:43
    for j = 1:52
        if ~isnan(u(i,j))
            if ~isnan(v(i,j))
                if ~isnan(w(i,j))
                    A(k,1) = u(i,j);
                    A(k,2) = v(i,j);
                    A(k,3) = w(i,j);
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
xx = zeros(leng,1);
yy = zeros(leng,1);
zz = zeros(leng,1);

B = eye(leng + 2);
k = 2;

for i = 1:leng
    B(k,k) = A(i,2);
    B(k,k-1) = A(i,1);
    B(k,k+1) = A(i,3);
    k = k + 1;
end

B_ = eig(B);

%quiver3(xx, yy, zz, uu, vv, ww);
%quiver(xx, yy, uu, vv);

theta_C = zeros(1,leng);
theta_D = zeros(1,leng);
theta_E = zeros(1,leng);

for i = 1:leng
    theta_C(i) = atan(A(i,1)/A(i,2));
    theta_D(i) = atan(A(i,2)/A(i,3));
    theta_E(i) = atan(A(i,1)/A(i,3));
end

vector_count = zeros(leng,3);
trim_c = 0.15;
trim_d = 0.15;
trim_e = 0.15;


for i = 1:leng
    for j = 1:leng
        
        if theta_C(j)-trim_c < theta_C(i) && theta_C(i) < theta_C(j)+trim_c
            vector_count(i,1) = vector_count(i,1) + 1;
        end
        if theta_D(j)-trim_d < theta_D(i) && theta_D(i) < theta_D(j)+trim_d
            vector_count(i,2) = vector_count(i,2) + 1;
        end
        if theta_E(j)-trim_e < theta_E(i) && theta_E(i) < theta_E(j)+trim_e
            vector_count(i,3) = vector_count(i,3) + 1;
        end
    end
    %     theta_C(i) = NaN;
    %     theta_D(i) = NaN;
    %     theta_E(i) = NaN;
end

% for i = 1:length
%     for j = 1:3
%         if vector_count(i,j) < 18
%             vector_count(i,j) = NaN;
%         end
%
%     end
% end

% scatter(1:length,vector_count(:,1))
% hold on
% scatter(1:length,vector_count(:,2))
% hold on
% scatter(1:length,vector_count(:,3))


% close all
%
% plot(vector_count(:,1))
% hold on
% plot(vector_count(:,2))
% hold on
% plot(vector_count(:,3))

grid on












