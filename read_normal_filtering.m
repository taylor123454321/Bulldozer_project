clear, clc, close all

load('A_for_trip_2.mat');

frames = 5
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
for i = 1:leng
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

m = 1;
vector = zeros(1,5);

matches_cut = matches;
biggest_match = max(matches);
cut_off = 0.6*biggest_match;
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
vector_count_3 = zeros(leng_mid,1);
trim_1 = 0.1;
trim_2 = trim_1;
trim_3 = trim_2;
m = 1;
vectors = zeros(3,3);


for i = 1:leng_mid % searching for matching vectors based on correclation
    for j = 1:3
        if i ~= j
            if corr(vectors(j,:)',main_vector(i,:)') > 0.8
                vectors_to_sum(m,:) = main_vector(j,:);
                m = m + 1;
            end
        end
    end
    vectors_to_sum
    m = 1;
end
vector_count_3

% [i,width_2] = size(vector_count_2);
%
% m = 1;
% n = 0;
% for i = 1:leng_mid
%     if vector_count_2(i,1) ~= 0
%         for j = 1:width_2
%             vector_count_2(m,j) = vector_count_2(i,j);
%             if vector_count_2(i,j) ~= 0
%                 n = n + 1;
%             end
%         end
%         vector_match_2(m,1) = i;
%         matches_2(i,1) = n;
%         n = 0;
%         m = m + 1;
%     end
% end
%
% matches_cut_2 = matches_2;
% biggest_match_2 = max(matches_2);
% cut_off_2 = 0.8*biggest_match;
% for i = 1:length(matches_2);
%     if matches_2(i,1) < cut_off_2
%         matches_cut_2(i,1) = 0;
%     end
% end

% m = 1;
% for i = 1:length(matches_cut_2)
%     if matches_cut_2(i,1) ~= 0
%         main_vector_mid(m,:) = main_vector(i,:);
%         m = m + 1;
%     end
% end




for i = 1:length(main_vector(:,1))
    for j = 1:3
        u(i,j) = main_vector(i,j)*1.2;
    end
end
% for i = 1:length(main_vector_mid(:,1))
%     for j = 1:3
%         y(i,j) = main_vector_mid(i,j)*1.4;
%     end
% end

mm = zeros(length(u(:,1)),1);
% nn = zeros(length(y(:,1)),1);

quiver3(mm, mm, mm, u(:,1), u(:,2), u(:,3),'black');
% quiver3(nn, nn, nn, y(:,1), y(:,2), y(:,3),'red');
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




