clear, clc, close all

load('A_for_trip_4.mat');

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
trim_1 = 0.015;
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

%     for i = 1:leng_2
%         if vector_count_2(i,width) ~= 0
%             main_vector_p(m,:) = A(vector_match(i,1),:); vector(1) = i; m
%             = m + 1;
%         end if vector_count_2(i,width-1) ~= 0 && vector(1) ~= i
%             main_vector_p(m,:) = A(vector_match(i,1),:); vector(2) = i; m
%             = m + 1;
%         end if vector_count_2(i,width-2) ~= 0 && vector(1) ~= i  &&
%         vector(2) ~= i
%             main_vector_p(m,:) = A(vector_match(i,1),:); m = m + 1;
%         end
%     end

biggest_match = max(matches);
cut_off = round(0.7*biggest_match);
for i = 1:length(matches);
    if matches(i,1) < cut_off
        matches(i,1) = 0;
    end
end

m = 1;
for i = 1:length(matches)
    if matches(i,1) ~= 0
        main_vector_p(m,:) = A(i,:);
        m = m + 1;
    end
end


% total = zeros(1,frames);
%
% for i = 1:frames
%     for j = 1:max_matches
%         total(i) = total(i) + vector_time(j,i); if vector_time(j,i)
%             break
%         end
%     end vector_avg(i) = total(i)/j;
% end
%
% vector_avg mean(vector_avg)

% k = 1; for i = 1:leng
%     if isnan(vector_count(i,1)) || isnan(vector_count(i,2)) ||
%     isnan(vector_count(i,3))
%         A_(k,:) = A(i,:); k = k + 1;
%     end
% end xxx = zeros(k-1,1); yyy = zeros(k-1,1); zzz = zeros(k-1,1);

% close all



% scatter(1:leng,vector_count(:,1)) hold on
% scatter(1:leng,vector_count(:,2)) hold on
% scatter(1:leng,vector_count(:,3))

% plot(vector_count(:,1)) hold on plot(vector_count(:,2)) hold on
% plot(vector_count(:,3))

% for i = 1:frames
%     hold on scatter(1:max_matches,vector_time(:,i))
% end

% for i = 1:frames
%     hold on
%     scatter3(ones(1,max_matches)*i,1:max_matches,vector_time(:,i))
% end


% for i = 1:frames
%     hold on plot(vector_time(:,i))
% end

% quiver3(xxx, yyy, zzz, A_(:,1), A_(:,3), A_(:,2));

% grid on
u = main_vector_p;
for i = 1:length(main_vector_p(:,1))
    for j = 1:3
        main_vector_p(i,j) = main_vector_p(i,j)*1.5;
    end
end


mm = zeros(length(main_vector_p(:,1)),1);

quiver3(mm, mm, mm, main_vector_p(:,1), main_vector_p(:,2), main_vector_p(:,3),'black');

matches
%
% vector_count_2


