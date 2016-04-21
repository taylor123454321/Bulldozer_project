function [vectors] = find_3(main_vector, directions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

vectors = zeros(directions,3);
vectors(1,:) = main_vector(1,:);
m = 2;
for i = 1:length(main_vector(:,1))
    if check_vec_corr(vectors, main_vector(i,:)) == 1
        vectors(m,:) = main_vector(i,:);
        m = m + 1;
    end
    if m > directions
        break;
    end
end

end

