clear,clc,close all

load('base_vector_to_avg')
frames = 5;
directions = 3;

vector_total = zeros(directions,3);
for i = 1:directions
    for j = 1:3
        vector_total(i,j) = sum(vector_overall(i,j,:))/frames;
    end
end
for i = 1:5
    vector_total - vector_overall(:,:,i);
end

save('base_vector', 'vector_total')