function [ output ] = better_match(vector_base, vector_test )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
output = zeros(length(vector_base(:,1)),1);

for i = 1:length(vector_base(:,1))
    output(i) = angle_betweend(vector_base(i,:),vector_test(i,:);
end

end

