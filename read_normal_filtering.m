clear, clc, close all

load('A_for_trip_2.mat');
frames = 5;
max_matches = 500;
vector_time = zeros(max_matches,frames);

leng = length(A);

uu = A(:,1);
vv = A(:,2);
ww = A(:,3);
%         xx = F(:,1); yy = F(:,2); zz = F(:,3);
xx = zeros(leng,1);
yy = zeros(leng,1);
zz = zeros(leng,1);

hold on
quiver3(xx, yy, zz, uu, vv, ww);

%     theta_C = zeros(1,leng); theta_D = zeros(1,leng); theta_E =
%     zeros(1,leng);

%     for i = 1:leng
%         theta_C(i) = atan(A(i,1)/A(i,2)); theta_D(i) =
%         atan(A(i,2)/A(i,3)); theta_E(i) = atan(A(i,1)/A(i,3));
%     end

vector_count = zeros(leng,1);
trim_1 = 0.06;
trim_2 = trim_1;
trim_3 = trim_2;
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
cut_off = 0.4*biggest_match; 
% 0.75 depends on how many matches/vectors there are
% cut of the amount of surface area covered with vectors wanting to read to
% noise. A calculation for this is needed.
for i = 1:length(matches);
    if matches(i,1) < cut_off
        matches_cut(i,1) = 0;
    end
end

m = 1;
for i = 1:length(matches_cut)
    if matches_cut(i,1) ~= 0
        main_vector(m,:) = A(i,:);
        m = m + 1;
    end
end

leng_mid = length(main_vector); 
directions = 3;
vectors = find_n(main_vector, directions, 0.5) % finds independent vectors from the matches
v = vectors;
vectors_to_sum = zeros(1,3,directions);

for i = 1:leng_mid % searching for matching vectors based on correclation
    for j = 1:directions
        correrlation(j) = corr(vectors(j,:)',main_vector(i,:)');
    end
    correrlation;
    [k,corr_max] = max(correrlation);
    vectors_to_sum(i,:,corr_max) = main_vector(i,:);
end
vectors_to_sum;

% NON ZERO NON ZERO NON ZERO NON ZERO NON ZERO 
m = 1;
for i = 1:length(vectors_to_sum(:,1,1))
    if vectors_to_sum(i,1,1) ~= 0
        vectors_to_sum_1(m,:) = vectors_to_sum(i,:,1);
        m = m + 1;
    end
end
m = 1;
for i = 1:length(vectors_to_sum(:,1,2))
    if vectors_to_sum(i,1,2) ~= 0
        vectors_to_sum_2(m,:) = vectors_to_sum(i,:,2);
        m = m + 1;
    end
end
m = 1;
for i = 1:length(vectors_to_sum(:,1,3))
    if vectors_to_sum(i,1,3) ~= 0
        vectors_to_sum_3(m,:) = vectors_to_sum(i,:,3);
        m = m + 1;
    end
end
% NON ZERO NON ZERO NON ZERO NON ZERO NON ZERO 

vectors = zeros(directions,3);
vectors(1,:) = sum(vectors_to_sum_1)/length(vectors_to_sum_1(:,1));
vectors(2,:) = sum(vectors_to_sum_2)/length(vectors_to_sum_2(:,1));
vectors(3,:) = sum(vectors_to_sum_3)/length(vectors_to_sum_3(:,1));

vectors

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
        y(i,j) = vectors(i,j)*1.75;
    end
end

mm = zeros(length(u(:,1)),1);
ll = zeros(length(z(:,1)),1);
nn = zeros(length(y(:,1)),1);

quiver3(mm, mm, mm, u(:,1), u(:,2), u(:,3),'black');
quiver3(ll, ll, ll, z(:,1), z(:,2), z(:,3),'green');
quiver3(nn, nn, nn, y(:,1), y(:,2), y(:,3),'red');

grid on

%
% vector_count_2
% close all
%
% matches_sorted = sort(matches);
% deltax = 1;
% for i = 2:deltax:length(matches_sorted)-1
%     f(i) = [matches_sorted(i+deltax) - matches_sorted(i-deltax)] / (2*deltax);
% end
%
%
% plot(f)




