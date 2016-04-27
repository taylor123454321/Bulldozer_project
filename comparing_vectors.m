clear, clc, close all

load('vectors_to_compare_1')
vectors_1 = vectors;
load('vectors_to_compare_2')
vectors_2 = vectors;
load('vectors_to_compare_3')
vectors_3 = vectors;
load('vectors_to_compare_4')
vectors_4 = vectors;
load('vectors_to_compare_5')
vectors_5 = vectors;

clear('vectors')

quiver3([0;0;0], [0;0;0], [0;0;0], vectors_1(:,1), vectors_1(:,2), vectors_1(:,3),'blue');
hold on
% quiver3([0;0;0], [0;0;0], [0;0;0], vectors_2(:,1), vectors_2(:,2), vectors_2(:,3),'green');
% quiver3([0;0], [0;0], [0;0], vectors_3(:,1), vectors_3(:,2), vectors_3(:,3),'green');
% quiver3([0;0;0], [0;0;0], [0;0;0], vectors_4(:,1), vectors_4(:,2), vectors_4(:,3),'cyan');
% quiver3([0;0], [0;0], [0;0], vectors_5(:,1), vectors_5(:,2), vectors_5(:,3),'red');
rotate3d on
directions = 3;
history = 0;
vector_matched = vectors_2;

for i = 1:20
    [vector_matched, angle, axis, uncert(i,:)] = match_vectors(vectors_1, vector_matched, directions, i);
    uncert
end


rotation = angle_betweend(vector_matched,vectors_2)


oo = zeros(length(vector_matched(:,1)),1);
quiver3(oo, oo, oo, vector_matched(:,1), vector_matched(:,2), vector_matched(:,3),'black');
% quiver3(oo, oo, oo, vector_matched_2(:,1), vector_matched_2(:,2), vector_matched_2(:,3),'cyan');
% quiver3(oo, oo, oo, vector_matched_3(:,1), vector_matched_3(:,2), vector_matched_3(:,3),'red');

% oo = [0;0;0];
% i = [1 0 0; 0 0 0; 0 0 0];
% quiver3(oo, oo, oo, i(:,1), i(:,2), i(:,3),'cyan');


rotate3d on


