function [out] = check_vec_corr(vectors,test, corr_value)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 0 for match
% 1 for no match

out = 0;

for i = 1:length(vectors(:,1))
    co = corr(vectors(i,:)',test');
    if co > corr_value %&& isnan(co) ~= 1;
        out = 1;
        break;
    end
end


end

