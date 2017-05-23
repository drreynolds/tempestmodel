function [c] = max_stable_step(A,b,dir)
% Usage: [c] = max_stable_step(A,b,dir)
% 
% Determines the magnitude of the maximum linearly stable step size
% in the direction indicated by the complex number dir.  This is
% computed by first creating the linear stability function R(z)
% for the Dahlquist test problem
%      y' = lambda*y,  y(0) = 1.
% We then ensure that
%      |R(0)| \le 1
% and we determine a value 
%      zmax = u*cmax, cmax \in \R, cmax > 0,  u = dir/|dir|
% such that 
%      |R(zmax)| > 1.
% If we cannot find a value zmax, then we return with z = inf.
% Otherwise we then test a large set of values in the interval
% (0,cmax) to determine a reasonable c_ub > 0 such that it is the
% minimum of all tested values still satisfying
%      |R(u*c_ub)| > 1.
% With that in place, we then perform a bisection algorithm within
% the interval [0,c_ub] to determine the maximum value of c \in \R,
% c>0, such that 
%      |R(u*d)| \le 1  for all d \le c,
% to within 4 significant digits.
%
% We then return the value  c.
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2017, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% verify that |dir| ~= 0, and set unit direction vector
if (abs(dir) == 0)
   error('Unsuable direction (cannot have zero magnitude)')
end
u = dir/abs(dir);

% set tolerance for bisection algorithm
tol = 1e-6;

% generate stability function handle
s = length(b);
b = reshape(b,1,s);
e = ones(s,1);
R = @(z) abs(1 + z*b*((eye(s)-z*A)\e));

% verify that R(0) \le 1
if (R(0) > 1)
   error('Method is unstable at z=0')
end

% scan to find cmax
cmax = eps;
maxtries=512;
for i=1:maxtries
   if (R(cmax*u) > 1)
      break;
   end
   cmax = 2*cmax;
end
if (i == maxtries)  % no max found, return with c=inf
   c = inf;
   return
end
   

% perform bisection algorithm
c_lb = 0;
c_ub = cmax;
c = 0.5*(c_lb + c_ub);
Ra = R(c_lb*u)-1;
Rb = R(c_ub*u)-1;
Rc = R(c*u)-1;
if (Ra*Rb > 0)
   error('Unusable bisection interval');
end
for i=1:100000
   if (Ra*Rc < 0)
      c_ub = c;
      Rb = Rc;
   else
      c_lb = c;
      Ra = Rc;
   end
   c = 0.5*(c_lb + c_ub);
   Rc = R(c*u)-1;
   if ( 0.5*(c_ub-c_lb) < tol )
      break;
   end
end

% end of function
