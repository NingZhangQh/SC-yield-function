% 20200702 Smoothed Classic yield function (2D)
% Author: Ning Zhang
%   This function provides a sample programming for the new yield function
%   In this code, the yield value, gradient, and Hessian for a given stress
% are sovled
% 
% [Output]
%   f1: the yield value
%   a1: the gradient
%   H1: the Hessian
%   fres: some status variables
% [Input]
%   sig  = [sigx, sigy, txy, sigz] is the stress vector 
%   pars = [M, K, t, m, alpha, beta, gama] is the parameter vector
%   where 
%       M K alpha beta gama are model parameters for GC yied function 
%       m <<0 is the comrpessive cap strength (if m==-inf, then m is not valid)
%       0 <= t <= K/M is the tensile strength
%   needaH is an option
%       true: calculate f1, a1, and H1
%       false: only calculate f1
% [Note]
%   the code is designed for 2d plane stress/strain, where the stresses is ordered by [sigx, sigy txy, sigz] 
%   
% [Example]
%   sig =   [5,2,0.1,4]; % f12
%   c = 14e6; phi = 0.4; beta  = 0.995; t = 0; m = -inf;
% 
%   sinf = sin(phi); cosf = cos(phi); sqrt3 = sqrt(3);
%   gbar  = 6/pi*atan(sinf/sqrt3);
%   alpha = 1/cos(pi/6*(gbar+1));
%   gama  = 1-gbar;
% 
%   tt= 6/(sqrt3*(3-sinf));
%   M = sinf*tt; 
%   K = c*cosf*tt;
%  
%   pars  = [M, K, t, m, alpha, beta gama];
%   needaH = 1; sig2 = sig;
% 
%   [f1, a1y, H1y, fres] = sub_SC(sig, pars, needaH);

function [f1, a1, H1, fres] = sub_SC(sig, pars, needaH)
% prepare parameters
sqrt3 = sqrt(3); pi_6  = pi/6;
M = pars(1); K = pars(2);
t = pars(3); m = pars(4);
alpha = pars(5); beta  = pars(6); gama  = pars(7);

%% get basic variables
sigz = sig(4);
p  = (sig(1) + sig(2) + sigz)/3;
sx = sig(1) - p; sy = sig(2) - p; sz = sigz - p;
txy= sig(3);

Jsqr = 0.5*(sx*sx + sy*sy + sz*sz) + txy*txy;
if Jsqr > 0
    J    = sqrt(Jsqr);
    sbar3= J * Jsqr;
    J3 = sx*sy*sz - sz*txy*txy;
    sin3t = -1.5*sqrt3*J3/sbar3;
    if sin3t < -1
        sin3t = -1;
    elseif sin3t > 1
        sin3t = 1;
    end
else % J2 =0
    sin3t = 0; J3 = 0;
    J = 0; sbar3 = 0;
end

%% cal f
Gama = alpha*cos(1/3*acos(-beta*sin3t) - gama*pi/6);
JG   = J*Gama;

a1 = [0;0;0;0];
H1 = [0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0];
fres= sub_calF(p, JG, K, M, t, m); % f, num, a12,a13, p12,p13
f1  = fres(1);
if ~needaH
    return
end

%% pD, J3D, JD, thetaDcos3,  JDD, J3DD
% prepare
pD  = [1/3;1/3; 1/3; 0];
JD  = [0;0;0;0];
J3D = [0;0;0;0];
JDD = [0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0];
J3DD= [0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0];

% sbarD
tt = 0.5/J;
JD(1) = sx*tt;
JD(2) = sy*tt;
JD(3) = sz*tt;
JD(4) = 2*txy*tt; % ok

% J3D
tt = Jsqr/3;
J3D(1) = sy*sz + tt;
J3D(2) = sz*sx + tt;
J3D(3) = sx*sy-txy*txy + tt;
J3D(4) = -2*sz*txy;   % ok

% lodaD
thetaDcos3t = -sqrt3/2/sbar3 * (J3D - 3*J3/J*JD);

% sbarDD
JDD(1,1) = (1/3 - sx*sx/4/Jsqr)/J;
JDD(1,2) = (-1/6 - sx*sy/4/Jsqr)/J;   JDD(2,1) = JDD(1,2);
JDD(1,3) = (-1/6 - sx*sz/4/Jsqr)/J;   JDD(3,1) = JDD(1,3);
JDD(1,4) = (-txy*sx/2/Jsqr)/J;        JDD(4,1) = JDD(1,4);

JDD(2,2) = (1/3 - sy*sy/4/Jsqr)/J;
JDD(2,3) = (-1/6 - sy*sz/4/Jsqr)/J;   JDD(3,2) = JDD(2,3);
JDD(2,4) = (-txy*sy/2/Jsqr)/J;        JDD(4,2) = JDD(2,4);

JDD(3,3) = (1/3 - sz*sz/4/Jsqr)/J;
JDD(3,4) = (-txy*sz/2/Jsqr)/J;        JDD(4,3) = JDD(3,4);

JDD(4,4) = (1 - txy*txy/Jsqr)/J;

% J3DD
J3DD(1,1) = (sx-sy-sz)/3;
J3DD(1,2) = 2/3*sz;            J3DD(2,1) = J3DD(1,2);
J3DD(1,3) = 2/3*sy;            J3DD(3,1) = J3DD(1,3);
J3DD(1,4) = 2/3*txy;           J3DD(4,1) = J3DD(1,4);

J3DD(2,2) = (sy-sx-sz)/3;
J3DD(2,3) = 2/3*sx;            J3DD(3,2) = J3DD(2,3);
J3DD(2,4) = 2/3*txy;           J3DD(4,2) = J3DD(2,4);

J3DD(3,3) = (sz-sx-sy)/3;
J3DD(3,4) = -4/3*txy;          J3DD(4,3) = J3DD(3,4);

J3DD(4,4) = -2*sz;

%% get c1, c2, c3 for f1
d  = -alpha*beta/( sqrt(1 - beta^2*sin3t^2) )...
    *sin(acos(-beta*sin3t)/3 - gama*pi_6);
c1  = M;
c2  = Gama - d*sin3t;
c3  = -sqrt3/2/Jsqr *d;

% the H1
betasin3t = beta*sin3t;
dDtheta_cos3t = beta*beta / (1-betasin3t*betasin3t)...
    *(3*d*sin3t-Gama);

c2D  = (-2*d-dDtheta_cos3t*sin3t) * thetaDcos3t;
c3D  = -sqrt3/Jsqr/2 * (dDtheta_cos3t*thetaDcos3t - 2*d/J * JD);

order = [1,2,4,3];
for ii = 1: 4
    mm = order(ii);
    a1(mm) = c1 * pD(ii) + c2 * JD(ii) + c3 * J3D(ii);
    a1(mm) = c1 * pD(ii) + c2 * JD(ii) + c3 * J3D(ii);
    H1(mm,1)= c2D(ii) * JD(1) + c2 * JDD(ii,1) + c3D(ii)*J3D(1) + c3 * J3DD(ii,1);
    H1(mm,2)= c2D(ii) * JD(2) + c2 * JDD(ii,2) + c3D(ii)*J3D(2) + c3 * J3DD(ii,2);
    H1(mm,4)= c2D(ii) * JD(3) + c2 * JDD(ii,3) + c3D(ii)*J3D(3) + c3 * J3DD(ii,3);
    H1(mm,3)= c2D(ii) * JD(4) + c2 * JDD(ii,4) + c3D(ii)*J3D(4) + c3 * J3DD(ii,4);
end

%% get gradient a and Hessian H
a2 = [1/3;  1/3; 0;  1/3];
a3 = [-1/3;-1/3; 0; -1/3];

% prepare data
num = fres(2); a12 = fres(3); a13 = fres(4); d12 = fres(5); d13 = fres(6);

% combined with f2
if num == 1
elseif num == 2
    a1 = a2;
elseif num == 3
    a1 = a3;
elseif num == 12 || num == 123
    da = a1 - a2; pi_a = pi/a12;
    gx = 0.5 *sin(0.5*pi_a*d12);
    gxx= 0.25*cos(0.5*pi_a*d12)*pi_a;
    % modify a1 and H1
    a1 = 0.5*(a1 + a2) + gx*da;
    H1 = (0.5 + gx) * H1 + gxx*(da*da');
end

% combined with f3 
if num == 13 || num == 123
    da = a1 - a3; pi_a = pi/a13;
    gx = 0.5 *sin(0.5*pi_a*d13);
    gxx= 0.25*cos(0.5*pi_a*d13)*pi_a;
    % modify a1 and H1
    a1 = 0.5*(a1 + a3) + gx*da;
    H1 = (0.5 + gx) * H1 + gxx*(da*da');
end

end

% get data
function res = sub_calF(p, JG, K, M, t, m)
res = [0,0,0,0,0,0];
a12 = K - M*t;
a13 = K - M*m;

f1 = M*p-K+JG;
f2 = p - t;
f3 = m - p;

fnum= 1;
d12 = f1-f2;
if t~=inf &&-a12 <= d12 && d12 <= a12
    f12 = 0.5*(f1+f2+a12) - a12/pi*cos(0.5*pi*d12/a12);
    f1  = f12;  fnum = fnum*10+2;
elseif t~=inf && f2 > f1
    f1  = f2;   fnum = 2;
end

d13 = f1-f3;
if m~=-inf && -a13 <= d13 && d13 <= a13
    f13 = 0.5*(f1+f3+a13) - a13/pi*cos(0.5*pi*d13/a13);
    f1 = f13;    fnum = fnum*10+3;
elseif m~=-inf && f3>f1
    f1 = f3;     fnum = 3;
end
res(1)=f1; res(2)=fnum; res(3)=a12; res(4)=a13; res(5) = d12; res(6) = d13;
end
