function [vectors, rot, axis, uncert] = match_vectors(base_vector, searching_vector, directions, round)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(searching_vector(:,1))
    for j = 1:length(base_vector(:,1))
        angle(i,j) = angle_betweend(searching_vector(i,:),base_vector(j,:));
    end
end

history_1 = mean(min(angle));

matches = zeros(length(angle));

for i = 1:length(base_vector(:,1))
    [value,col] = min(angle(:,i));
    matches(i,col) = 1;
end

[closest_vector,index] = min(min(angle));
if round == 2
    for i = 1:3
        if angle(index,i) == closest_vector
            angle(index,i) = angle(3,3);
        end
    end
    [closest_vector,index] = min(min(angle));
end

axis = searching_vector(index,:);

rot = max(min(angle));

vectors = zeros(directions,3);
for i = 1: directions
    vectors(i,:) = rot_3space(searching_vector(i,:),axis,rot);
end


for i = 1:length(vectors(:,1))
    for j = 1:length(base_vector(:,1))
        angle_2(i,j) = angle_betweend(vectors(i,:),base_vector(j,:));
    end
end
 
history_2 = mean(min(angle_2));
if history_2 > history_1
    rot = -rot;
    vectors_2 = zeros(directions,3);
    for i = 1: directions
        vectors_2(i,:) = rot_3space(searching_vector(i,:),axis,rot);
    end
    for i = 1:length(vectors(:,1))
        for j = 1:length(base_vector(:,1))
            angle_3(i,j) = angle_betweend(vectors_2(i,:),base_vector(j,:));
        end
    end
    history_3 = mean(min(angle_3));
    if history_3 < history_2 && history_3 < history_1
        vectors = vectors_2;
    end
end

for i = 1:length(vectors(:,1))
    for j = 1:length(base_vector(:,1))
        angle_4(i,j) = angle_betweend(vectors(i,:),base_vector(j,:));
    end
end

mean(min(angle_4))
uncert = min(angle_4);


end

