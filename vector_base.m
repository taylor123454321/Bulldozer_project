function [vector_location, vector_number] = vector_base(F)

trim = 0.05;
n = 1;
m = 1;

[row_len, col_len, frames] = size(F);
match = zeros(row_len,10);

left_x = -trim;
right_x = trim;
left_y = -trim;
right_y = trim;
left_z = -trim;
right_z = trim;

for i = 1:row_len
    for o = 2:frames
        for j = 1:row_len
            if (F(j,1,o) > (F(i,1,1)+left_x)) && (F(j,1,o) < (F(i,1,1)+right_x))          % X
                if (F(j,2,o) > (F(i,2,1)+left_y)) && (F(j,2,o) < (F(i,2,1)+right_y))      % Y
                    if (F(j,3,o) > (F(i,3,1)+left_z)) && (F(j,3,o) < (F(i,3,1)+right_z))  % Z
                        left_x = F(j,1,o) - trim;
                        right_x = F(j,1,o) + trim;
                        left_y = F(j,2,o) - trim;
                        right_y = F(j,2,o) + trim;
                        left_z = F(j,3,o) - trim;
                        right_z = F(j,3,o) + trim;
                        
                        match(i,m,o) = j;
                        n = n + 1;
                        m = m + 1;
                    end
                end
            end
        end
    end
    left_x = -trim;
    right_x = trim;
    left_y = -trim;
    right_y = trim;
    left_z = -trim;
    right_z = trim;
    m = 1;
end

[l, width] = size(match);

% for i = 1:row_len
%     for j = 1:width
%         if i == match(i,j);
%             match(i,j) = 0;
%         end
%     end
% end

vector_spot= zeros(row_len,3);

sum1 = 0;
sum2 = 0;
sum3 = 0;
l = 0;

for i = 1:row_len
    for j = 1:col_len
        for k = 1:width
            if match(i,k) ~= 0
                sum1 = F(match(i,k),1) + sum1;
                sum2 = F(match(i,k),2) + sum2;
                sum3 = F(match(i,k),3) + sum3;
                l = l + 1;
            end
        end
        vector_spot(i,1) = sum(sum1)/l;
        vector_spot(i,2) = sum(sum2)/l;
        vector_spot(i,3) = sum(sum3)/l;
        sum1 = 0;
        sum2 = 0;
        sum3 = 0;
        l = 0;
    end
end

% for i = 1:row_len
%     fprintf('%d, %d %d %d %d %d %d %d\n',i,match(i,1),match(i,2),match(i,3),match(i,4),match(i,5),match(i,6),match(i,6));
% end


scatter3(F(:,1),F(:,2),F(:,3),'b')
hold on
scatter3(vector_spot(:,1),vector_spot(:,2),vector_spot(:,3))


end













