function z = reentry_meas(x, Rd, params, status)
%
%   reentry measurement model
%

% radar cooridates
xR = params.xR;
yR = params.yR;

% x(k)
x1 = x(1);
x2 = x(2);

% measurement noise
if status == 'sy'
    v = sqrt(Rd) * randn(2,1);
else
    v = zeros(2,1);
end

% the radar measurements : a distance and an elevation angle
r = sqrt((x1-xR)^2 + (x2-yR)^2) + v(1);
the = atan2((x2-yR), (x1-xR)) + v(2);

z = [r; the];

end