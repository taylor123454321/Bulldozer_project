function [out] = check_vec_corr(vectors,test, corr_value)
% This function checks for correletion between one incoming vector and the
% set of vectors already picked
% 0 for match
% 1 for no match

out = 0;

for i = 1:length(vectors(:,1))
    co = corr(vectors(i,:)',test');
    if co > corr_value
        out = 1;
        break;
    end
end


end

