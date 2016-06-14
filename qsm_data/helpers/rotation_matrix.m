function R = rotation_matrix(u,k)
%Rotation matrix.
%When multiplied to a vector performs the rotation of it around the axis
%defined by <u> and by the angle specified in <k>.
R = zeros(3,3);
c = cos(k); 
s = sin(k);
R(1,:) = [u(1)^2+(1-u(1)^2)*c  u(1)*u(2)*(1-c)-u(3)*s  u(1)*u(3)*(1-c)+u(2)*s];
R(2,:) = [u(1)*u(2)*(1-c)+u(3)*s  u(2)^2+(1-u(2)^2)*c  u(2)*u(3)*(1-c)-u(1)*s];
R(3,:) = [u(1)*u(3)*(1-c)-u(2)*s  u(2)*u(3)*(1-c)+u(1)*s  u(3)^2+(1-u(3)^2)*c];
