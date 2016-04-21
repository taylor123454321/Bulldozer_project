function [out] = check_vec_corr(vectors,test)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 0 to for match
% 1 for no match

out = 0;

for i = 1:length(vectors(:,1))
    if corr(vectors(i,:)',test') < 0.6
        out = 1;
        break;
    end
end


end

