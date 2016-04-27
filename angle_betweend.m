function [ angle ] = angle_betweend( vector_1, vector_2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i = 1:length(vector_1(:,1))
    angle(i,:) = acosd(dot(vector_1(i,:),vector_2(i,:))/(norm(vector_1(i,:))*norm(vector_2(i,:))));
end
end

