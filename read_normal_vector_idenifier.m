clear, clc, close all

load('FF.mat')

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

length = length(A);

uu = A(:,1);
vv = A(:,2);
ww = A(:,3);
xx = zeros(length,1);
yy = zeros(length,1);
zz = zeros(length,1);

B = eye(length + 2);
k = 2;

for i = 1:length
    B(k,k) = A(i,2);
    B(k,k-1) = A(i,1);
    B(k,k+1) = A(i,3);
    k = k + 1;
end

B_ = eig(B);

quiver3(xx, yy, zz, uu, vv, ww);
%quiver(xx, yy, uu, vv);

theta_C = zeros(1,length);
theta_D = zeros(1,length);
theta_E = zeros(1,length);

for i = 1:length
    theta_C(i) = atan(A(i,1)/A(i,2));
    theta_D(i) = atan(A(i,2)/A(i,3));
    theta_E(i) = atan(A(i,1)/A(i,3));
end

vector_count = zeros(length,3);
trim_c = 0.15;
trim_d = 0.15;
trim_e = 0.15;


for i = 1:length
    for j = 1:length
        
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


plot(vector_count(:,1))
hold on
plot(vector_count(:,2))
hold on
plot(vector_count(:,3))

grid on







