function [ vector_out ] = rot_3space(vector, axis, angle)
% This rotates a vector about an axis


vector_out = vector*cosd(angle) + cross(axis,vector)*sind(angle) + axis*dot(axis,vector)*(1-cosd(angle));


end

