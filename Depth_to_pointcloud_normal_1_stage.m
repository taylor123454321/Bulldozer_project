clear, clc, close all

depthDevice = imaq.VideoDevice('kinect',2);
load('base_vector')
vector_base = vector_total;
clear('vector_total')

step(depthDevice);

for i = 1:30
    depthImage = step(depthDevice);
end

tic
frames = 30;
filter = 4;

clear('depthImage')
for i = 1:filter
    image = step(depthDevice);
    %     image(image<600) = 0; %trimming minium and maxium length
    image(image>800) = 0;
    image(:,1:20) = 0;
    image(:,424-20:424) = 0;
    image(1:20,:) = 0;
    image(:,512-20:512) = 0;
    depthImage(:,:,i) = image;
end

sum_ = uint16(zeros(424,512));
elements = uint16(zeros(424,512));
dImage_filtered = uint16(zeros(424,512,frames-filter));

for b = 1:filter
    sum_(:,:) = sum_(:,:) + depthImage(:,:,b);
    for i = 1:424
        for j = 1:512
            if depthImage(i,j,b) ~= 0
                elements(i,j) = elements(i,j) + 1;
            end
        end
    end
end


% ptCloud = pcfromkinect(depthDevice,dImage_filtered(:,:,1));
% pcshow(ptCloud)
% release(depthDevice);

directions = 3;
filter_vector = 8;
filter_angle = 2;
match_angle = 20;
vectors_to_sum_overall = zeros(directions,3,filter_vector);
angle_to_sum_overall = zeros(directions,directions,filter_angle);
vector_total = zeros(directions,3,frames);
angle_total = zeros(directions,directions,frames);
vector_overall = zeros(directions,3);
max_matches = 500;
vector_time = zeros(max_matches,frames);
angle_ = zeros(directions,directions);
angles = zeros(1,directions,frames);
angle_overall = zeros(directions,directions);
trim_1 = 0.125;
trim_2 = trim_1;
trim_3 = trim_2;


for p = filter+1:frames
    p
    clear('depthImage')
    for i = 1:3
        depthImage = step(depthDevice);
    end
    image = step(depthDevice);
    %     image(image<600) = 0; %trimming minium and maxium length
    image(image>800) = 0;
    image(:,1:20) = 0;
    image(:,424-20:424) = 0;
    image(1:20,:) = 0;
    image(:,512-20:512) = 0;
    depthImage(:,:,p) = image;
    
    sum_(:,:) = sum_(:,:) - depthImage(:,:,p-filter) + depthImage(:,:,p);
    for i = 1:424
        for j = 1:512
            if depthImage(i,j,p-filter) ~= 0
                elements(i,j) = elements(i,j) - 1;
            end
            if depthImage(i,j,p) ~= 0
                elements(i,j) = elements(i,j) + 1;
            end
        end
    end
    dImage_filtered(:,:,p) = sum_(:,:)./elements(:,:);
    
    ptCloud = pcfromkinect(depthDevice,dImage_filtered(:,:,p));
    
    normals = pcnormals(ptCloud);
    
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
    
    clear('A')
    %     xx = x;yy = y;zz = z;
    A = zeros(20,3);
    F = zeros(20,3);
    
    k = 1;
    for i = 1:43
        for j = 1:52
            if ~isnan(u(i,j))
                if ~isnan(v(i,j))
                    if ~isnan(w(i,j))
                        A(k,1) = u(i,j); % vectors that values
                        A(k,2) = v(i,j); % A point
                        A(k,3) = w(i,j); % F location
%                         F(k,1) = x(i,j);
%                         F(k,2) = y(i,j);
%                         F(k,3) = z(i,j);
                        k = k + 1;
                    end
                end
            end
        end
    end

    %     AA(:,:,p) = A;
    uu = A(:,1);
    vv = A(:,2);
    ww = A(:,3);
    leng = length(A);
    xx = zeros(leng,1);
    yy = zeros(leng,1);
    zz = zeros(leng,1);
    
    if p == filter+1
        pcshow(ptCloud)
        hold on
    end
%     quiver3(xx, yy, zz, uu, vv, ww);
    
    main_vector = 0;
    vectors = 0;
    v = 0;
    vectors_from_matches = 0;

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
            if p == filter+1 || vector_overall(1,1) == 0  || vector_overall(2,2) == 0  || vector_overall(3,3) == 0
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
            
            % averaging the matches around the independent vectors
            vectors_from_matches = zeros(directions,3);
            for i = 1:directions
                sums = remove_zeros(vectors_to_sum(:,:,i));
                vectors_from_matches(i,:) = sum(sums)/length(sums(:,1));
            end
            
            m = 1;
            for i = 1:directions %takes the best fitting vectors and finds vectors around it from orginal data then averages them
                for j = 1:leng
                    if angle_betweend(A(j,:),vectors_from_matches(i,:)) < match_angle
                        vectors_to_average(m,:,i) = A(j,:);
                        m = m + 1;
                    end
                end
            end
            
            clear('vectors')
            vectors = [0 0 0];
            for i = 1:length(vectors_to_average(1,1,:))
                sums = remove_zeros(vectors_to_average(:,:,i));
                vectors(i,:) = sum(sums)/length(sums(:,1));
            end
            vectors;
            if  length(vectors(:,1)) < directions
                vectors(directions,:) = NaN;
            end
            
            vector_total(:,:,p-1) = vectors;
            filter_index = p-filter_vector:p;
            
            for i = 1:filter_vector
                if filter_index(i) < 1
                    filter_index(i) = 1;
                end
            end
            
            for i = 1:filter_vector
                vectors_to_sum_overall(:,:,i) = vector_total(:,:,filter_index(i));
            end
            
            for i = 1:directions
                for j = 1:3
                    vector_overall(i,j) = sum(vectors_to_sum_overall(i,j,:),'omitnan')/filter_vector;
                end
            end
            
            for i = 1:length(vector_overall(:,1))
                for j = 1:length(vector_base(:,1))
                    angle_(i,j) = angle_betweend(vector_overall(i,:),vector_base(j,:));
                end
            end
            
            vector_total(:,:,p-1) = vector_overall;
            
            angle_total(:,:,p-1) = angle_;
            angle_filter_index = p-filter_angle:p;
            
            for i = 1:filter_angle
                if angle_filter_index(i) < 1
                    angle_filter_index(i) = 1;
                end
            end
            
            for i = 1:filter_angle
                angle_to_sum_overall(:,:,i) = angle_total(:,:,angle_filter_index(i));
            end
            
            for i = 1:directions
                for j = 1:directions
                    angle_overall(i,j) = sum(angle_to_sum_overall(i,j,:),'omitnan')/filter_angle;
                end
            end
            angle_total(:,:,p-1) = angle_overall;
            angles(:,:,p) = min(angle_overall);
            angles(:,:,p)
            
            for i = 1:directions
                for j = 1:3
                    w(i,j) = vectors(i,j)*2;
                end
            end
            for i = 1:length(main_vector(:,1))
                for j = 1:3
                    u(i,j) = main_vector(i,j)*1.25;
                end
            end
            for i = 1:directions
                for j = 1:3
                    z(i,j) = v(i,j)*1.5;
                end
            end
            for i = 1:directions
                for j = 1:3
                    y(i,j) = vectors_from_matches(i,j)*1.75;
                end
            end
            
            mm = zeros(length(u(:,1)),1);
            ll = zeros(length(z(:,1)),1);
            nn = zeros(length(y(:,1)),1);
            oo = zeros(length(w(:,1)),1);
            
            quiver3(mm, mm, mm, u(:,1), u(:,2), u(:,3),'black');
            quiver3(ll, ll, ll, z(:,1), z(:,2), z(:,3),'green');
            quiver3(nn, nn, nn, y(:,1), y(:,2), y(:,3),'cyan');
            quiver3(oo, oo, oo, w(:,1), w(:,2), w(:,3),'red');
            
        end
    end
end
pcshow(ptCloud)
release(depthDevice);


% quiver3(xxx, yyy, zzz, A_(:,1), A_(:,3), A_(:,2));

% grid on
% % u = main_vector_p
% % for i = 1:length(main_vector_p(:,1))
% %     for j = 1:3
% %         main_vector_p(i,j) = main_vector_p(i,j)*1.5;
% %     end
% % end
% %
% %
% % mm = zeros(length(main_vector_p(:,1)),1);
% %
% % quiver3(mm, mm, mm, main_vector_p(:,1), main_vector_p(:,2), main_vector_p(:,3),'black');

% hold off
% plot(matches)
% matches
%
% vector_count_2
%
% for i = 1:length(main_vector(:,1))
%     for j = 1:3
%         u(i,j) = main_vector(i,j)*1.25;
%     end
% end
% for i = 1:directions
%     for j = 1:3
%         z(i,j) = v(i,j)*1.5;
%     end
% end
% for i = 1:directions
%     for j = 1:3
%         y(i,j) = vectors_from_matches(i,j)*1.75;
%     end
% end
for i = 1:directions
    for j = 1:3
        w(i,j) = vector_overall(i,j)*3;
    end
end
%
% hold on
% quiver3(xx, yy, zz, uu, vv, ww);
%
% mm = zeros(length(u(:,1)),1);
% ll = zeros(length(z(:,1)),1);
% nn = zeros(length(y(:,1)),1);
clear('oo')
oo = zeros(length(w(:,1)),1);
%
% quiver3(mm, mm, mm, u(:,1), u(:,2), u(:,3),'black');
% quiver3(ll, ll, ll, z(:,1), z(:,2), z(:,3),'green');
% quiver3(nn, nn, nn, y(:,1), y(:,2), y(:,3),'red');
quiver3(oo, oo, oo, w(:,1), w(:,2), w(:,3),'green');

for i = 1:length(angles(1,3,:))
    a(i) = angles(1,1,i);
    b(i) = angles(1,2,i);
    c(i) = angles(1,3,i);
end

figure
plot(a)
hold on
plot(b)
plot(c)
grid on









