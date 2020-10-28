function out = TreesNoSavanna
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,phi0,phi1,s,theta,alpha)
dydt=[(1-kmrgd(1))*(phi0+(phi1-phi0)/(1+exp(-(kmrgd(1)-theta)/s))-alpha*kmrgd(1));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(TreesNoSavanna);
y0=[0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,phi0,phi1,s,theta,alpha)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,phi0,phi1,s,theta,alpha)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,phi0,phi1,s,theta,alpha)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,phi0,phi1,s,theta,alpha)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,phi0,phi1,s,theta,alpha)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,phi0,phi1,s,theta,alpha)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,phi0,phi1,s,theta,alpha)
