function [vectors] = find_n(main_vector, directions, corr_value)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

vectors = zeros(directions,3);
vectors(1,:) = main_vector(1,:);
m = 2;
for i = 1:length(main_vector(:,1))
    if check_vec_corr(vectors, main_vector(i,:), corr_value) == 0
        % 0 for match
        % 1 for no match
        vectors(m,:) = main_vector(i,:);
        m = m + 1;
    end
    if m > directions
        break;
    end
end

end

