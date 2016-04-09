clear, clc, close all

depthDevice = imaq.VideoDevice('kinect',2);

step(depthDevice);

for i = 1:20
    depthImage = step(depthDevice);
end

player = pcplayer([-1.5    1.5],[-1.5    1.5],[0.5 2.5],...
    'VerticalAxis','y','VerticalAxisDir','down');

xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');



frames = 5;
max_matches = 500;
vector_time = zeros(max_matches,frames);

for p = 1:frames
    x = 0;
    y = 0;
    z = 0;
    xx = 0;
    yy = 0;
    zz = 0;
    u = 0;
    v = 0;
    w = 0;
    uu = 0;
    vv = 0;
    ww = 0;
    normals = 0;
    
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
    
    %ptCloud = pcdownsample(ptCloud,'random',0.80);
    
    normals = pcnormals(ptCloud);
    
    
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
    %     for i = 1:43
    %         for j = 1:52
    %             xx(i,j) = 0;
    %             yy(i,j) = 0;
    %             zz(i,j) = 0;
    %         end
    %     end
    
    %pcshow(ptCloud)
    %title('Adjusted Normals of Point Cloud')
    %hold on
    quiver3(xx, yy, zz, u, v, w);
    %surf(xx,yy,zz)
    %view(0,-45)
    
    A = 0;
    A = zeros(132,3);
    F = 0;
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
    xx = zeros(leng,1);
    yy = zeros(leng,1);
    zz = zeros(leng,1);
    
    B = 0;
    B = eye(leng + 2);
    k = 2;
    
    for i = 1:leng
        B(k,k) = A(i,2);
        B(k,k-1) = A(i,1);
        B(k,k+1) = A(i,3);
        k = k + 1;
    end
    
    pcshow(ptCloud)
    hold on
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
    trim_c = 0.07;
    trim_d = trim_c;
    trim_e = trim_c;
    
    
    for i = 1:leng
        for j = 1:leng
            
            if (theta_C(j)-trim_c > theta_C(i)) || (theta_C(i) < theta_C(j)+trim_c)
                vector_count(i,1) = vector_count(i,1) + 1;
            end
            if theta_D(j)-trim_d > theta_D(i) || theta_D(i) < theta_D(j)+trim_d
                vector_count(i,2) = vector_count(i,2) + 1;
            end
            if theta_E(j)-trim_e > theta_E(i) || theta_E(i) < theta_E(j)+trim_e
                vector_count(i,3) = vector_count(i,3) + 1;
            end
        end
        %     theta_C(i) = NaN;
        %     theta_D(i) = NaN;
        %     theta_E(i) = NaN;
    end
    
    for i = 1:leng
        for j = 1:3
            if vector_count(i,j) < 20
                vector_count(i,j) = NaN;
            end
        end
    end
    
    
    k = 1;
    for i = 1:leng
        if ~isnan(vector_count(i,2))
            vector_time(k,p) = vector_count(i,1);
            k = k + 1;
        end
    end
end
release(depthDevice);

trim = 0.001;
n = 1;
for i = 1:length(F)
    for j = 1:length(F)
        if (F(i,1)-trim) < F(j,1) && F(i,1)+trim > F(j,1)
            if F(i,2)-trim < F(j,2) && F(i,2)+trim > F(j,2)
                if F(i,3)-trim < F(j,3) && F(i,3)+trim > F(j,3)
                    n = n + 1;
                end
            end
        end
    end
end

n
% total = zeros(1,frames);
%
% for i = 1:frames
%     for j = 1:max_matches
%         total(i) = total(i) + vector_time(j,i);
%         if vector_time(j,i)
%             break
%         end
%     end
%     vector_avg(i) = total(i)/j;
% end
%
% vector_avg
% mean(vector_avg)

k = 1;
for i = 1:leng
    if isnan(vector_count(i,1)) || isnan(vector_count(i,2)) || isnan(vector_count(i,3))
        A_(k,:) = A(i,:);
        k = k + 1;
    end
end
xxx = zeros(k-1,1);
yyy = zeros(k-1,1);
zzz = zeros(k-1,1);



% close all


% scatter(1:leng,vector_count(:,1))
% hold on
% scatter(1:leng,vector_count(:,2))
% hold on
% scatter(1:leng,vector_count(:,3))

% plot(vector_count(:,1))
% hold on
% plot(vector_count(:,2))
% hold on
% plot(vector_count(:,3))

% for i = 1:frames
%     hold on
%     scatter(1:max_matches,vector_time(:,i))
% end

% for i = 1:frames
%     hold on
%     scatter3(ones(1,max_matches)*i,1:max_matches,vector_time(:,i))
% end


% for i = 1:frames
%     hold on
%     plot(vector_time(:,i))
% end

% quiver3(xxx, yyy, zzz, A_(:,1), A_(:,3), A_(:,2));

grid on








