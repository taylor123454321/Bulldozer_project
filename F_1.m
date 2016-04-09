clear, clc, close all

load('FF.mat');

trim = 0.0296;
n = 1;
m = 1;

[row_len, col_len] = size(F);
for i = 1:row_len
    for j = 1:row_len
        if (F(j,1) > (F(i,1)-trim)) && (F(j,1) < (F(i,1)+trim))
            if (F(j,2) > (F(i,2)-trim)) && (F(j,2) < (F(i,2)+trim))
                if (F(j,3) > (F(i,3)-trim)) && (F(j,3) < (F(i,3)+trim))
                    match(i,m) = j;
                    n = n + 1;
                    m = m + 1;
                end
            end
        end
    end
    m = 1;
end


% for i = 1:length(match)
%     for j = 1:col_len
%         for k = 1:row_len
%             spots_vectors = 
%             
%         end
%     end
% end


n

n/length(F)

