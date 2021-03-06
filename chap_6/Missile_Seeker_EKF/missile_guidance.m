function a_cmd = missile_guidance(r, v, vM)
%
% the missile guidance law
%

N = 3;
R2 = r' * r;
OmL = cross(r,v)/R2;
a_cmd = N * cross(OmL, vM);
end