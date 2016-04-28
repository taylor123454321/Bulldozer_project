function [ output ] = remove_zeros( input )
% this functions takes all the rows of zeros and removes them

m = 1;
output = [0 0 0];
for i = 1:length(input(:,1))
    if input(i,1) ~= 0 && input(i,2) ~= 0
        output(m,:) = input(i,:);
        m = m + 1;
    end
end

end

