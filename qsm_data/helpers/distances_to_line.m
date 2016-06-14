function [d,V,h,B] = distances_to_line(Q,LDirec,LPoint)
%
% USAGE: [d,V,h,B] = distances_to_line(Q,LDirec,LPoint)
%
% Calculates the distances of the points, given in the rows of the
% matrix Q, to the line defined by one of its point and its direction.

if size(LDirec,1) == 1
    LDirec = LDirec';
end
LDirec = LDirec/norm(LDirec);% normalize the vector

A = mat_vec_subtraction(Q,LPoint);
h = A*LDirec;

B = repmat(LDirec',length(Q(:,1)),1);
B = [h.*B(:,1) h.*B(:,2) h.*B(:,3)];
V = A-B;

d = sqrt(sum(V.*V,2));