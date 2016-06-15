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
normals = pcnormals(ptCloud,80);

toc

x = ptCloud.Location(1:10:end,1:10:end,1);
y = ptCloud.Location(1:10:end,1:10:end,2);
z = ptCloud.Location(1:10:end,1:10:end,3);
u = normals(1:10:end,1:10:end,1);
v = normals(1:10:end,1:10:end,2);
w = normals(1:10:end,1:10:end,3);

sensorCenter = [0,0,0];
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

A = zeros(20,3);
F = zeros(20,3);

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
trim_1 = 0.07;
trim_2 = trim_1;
trim_3 = trim_2;
directions = 3;
p = 2;

if leng > 20
    vector_count = zeros(leng,1);
    m = 1;
    
    for i = 1:leng % searching for matching vectors based on angle
        for j = 1:leng
            if A(j,1) > (A(i,1)-trim_1) && A(j,1) < (A(i,1)+trim_1)
                if A(j,2) > (A(i,2)-trim_2) && A(j,2) < (A(i,2)+trim_2)
                    if A(j,3) > (A(i,3)-trim_3) && A(j,3) < (A(i,3)+trim_3)
                        if i ~= j
                            vector_count(i,m) = j;
                            m = m + 1;
                        end
                    end
                end
            end
        end
        m = 1;
    end
    
    [i,width] = size(vector_count);
    
    m = 1;
    n = 0;
    matches = zeros(leng,1);
    for i = 1:leng % clearing zero lines
        if vector_count(i,1) ~= 0
            for j = 1:width
                vector_count_2(m,j) = vector_count(i,j);
                if vector_count(i,j) ~= 0
                    n = n + 1;
                end
            end
            vector_match(m,1) = i;
            matches(i,1) = n;
            n = 0;
            m = m + 1;
        end
    end
    
    [leng_2,width] = size(vector_count_2);
    matches_cut = matches; % trimming the matches
    biggest_match = max(matches);
    matching = 4; % ratio
    corr_value = 0.5;
    
    cut_off = 0.15*biggest_match;
    % 0.75 depends on how many matches/vectors there are
    % cut of the amount of surface area covered with vectors wanting to read to
    % noise. A calculation for this is needed.
    for i = 1:length(matches);
        if matches(i,1) < cut_off
            matches_cut(i,1) = 0;
        end
    end
    
    clear('main_vector')
    m = 1;
    for i = 1:length(matches_cut)
        if matches_cut(i,1) ~= 0
            main_vector(m,:) = A(i,:);
            m = m + 1;
        end
    end
    if main_vector ~= 0
        count = zeros(1,length(main_vector(:,1)));
        if p == 2 || vector_overall(1,1) == 0  || vector_overall(2,2) == 0  || vector_overall(3,3) == 0
            [vectors_from_matches,check] = find_n(main_vector, directions, corr_value); % finds independent vectors from the matches
        else
            vectors_from_matches = vector_overall;
        end
        v = vectors_from_matches;
        leng_mid = length(main_vector);
        vectors_to_sum = zeros(1,3,directions);
        correrlation = zeros(1,3);
        
        for i = 1:leng_mid % searching for matching vectors based on correclation
            for j = 1:directions
                correrlation(j) = corr(vectors_from_matches(j,:)',main_vector(i,:)');
            end
            [k,corr_max] = max(correrlation);
            vectors_to_sum(i,:,corr_max) = main_vector(i,:);
        end
    end
end

% 
% for i = 1:length(main_vector(:,1))
%     for j = 1:3
%         u(i,j) = main_vector(i,j)*1.25;
%     end
% end


mm = zeros(length(u(:,1)),1);
% quiver3(mm, mm, mm, u(:,1), u(:,2), u(:,3),'black');





