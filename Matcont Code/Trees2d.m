function out = Trees2d
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
function dydt = fun_eval(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
dydt=[mu*kmrgd(2) + nu*(1-kmrgd(1)-kmrgd(2)) - beta*kmrgd(1)*(1-kmrgd(1)-kmrgd(2));
beta*kmrgd(1)*(1-kmrgd(1)-kmrgd(2))-(w0 + (w1-w0)/(1+exp(-(kmrgd(1)-t1)/s1)))*kmrgd(2)-mu*kmrgd(2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Trees2d);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,mu,nu,beta,s1,t1,w0,w1)
