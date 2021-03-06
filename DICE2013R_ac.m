% -------------------------------------------------------------------------
%  This file is part of DICE2013R-mc.
% 
% DICE2013R-mc -- Dynamic Integrated model of Climate and Economy 
%                   2013R - Matlab and Casadi.
% Copyright (C) 2016 Timm Faulwasser, Christopher M. Kellett, and 
%       Steven R. Weller. 
% 
% DICE2013R-mc can be downloaded from https://github.com/cmkellett/DICE2013R-mc.
% 
% DICE2013R-mc is free software; you can redistribute it and/or modify it 
% under the terms of the GNU Lesser General Public License (LGPL) as 
% published by the Free Software Foundation; either version 3 of the 
% License, or (at your option) any later version.
% 
% DICE2013R-mc is distributed in the hope that it will be useful, but 
% without any warranty; without even the implied warranty of merchantability 
% or fitness for a particular purpose.  See the GNU Lesser General Public 
% License for more details.                        
%
% -------------------------------------------------------------------------
%
%   This code provides a bare bones replication of the GAMS DICE2013R
%   (http://aida.wss.yale.edu/~nordhaus/homepage/Web-DICE-2013-April.htm)
%   for use in MATLAB with CasADi.  See the accompanying documentation
%   for more details.
%
%
% -------------------------------------------------------------------------
% Release notes
%
%  - 1.0 initial release requires casadi-matlabR2014b-v3.1.1 or v.3.0.0 
%  - 1.1 fixes compatability with casadi-matlabR2014b-v3.2.1
%  - 2.0 adds parameters for DICE2016R-091916ap, with switch to choose
%       between 2013R (110513_vanilla) and 2016R parameter values
%  - 2.1 fixes minor parameter errors in 2016R and time indexing errors
%       in sigma and the radiative forcing
% -------------------------------------------------------------------------

clc
clear all
close all

% import casadi name space
import casadi.*
test      = SX.sym('test', 1, 2)


addpath('/Users/alexandercampbell/casadi-osx-matlabR2015a-v3.5.1')
import casadi.*
x = MX.sym('x')
disp(jacobian(sin(x),x))

%% ========================================================================
% define data of the DICE Optimal Control Problem (OCP)
% =========================================================================
Param_set = 2013; % Choose DICE 2013, 2016, or 1 for custom values 
Params1 = set_DICE_parameters(Param_set); 
Params2 = I_DICE_params(Param_set); % get DICE parameters
%%
Params = Params1;
N = Params.N;
p = Params2.p;
q = Params2.q;
s=4;

switch Param_set
    case 1
        x0 = [0.8 0.0068 830.4 1527 10010 135]'; % initial condition 
    % IT IS RECOMMENDED THAT CHANGES BE LIMITED TO THE ABOVE!!!
    % Initial conditions below correspond to DICE2013R and DICE2016R.
    case 2013
        x0 = [0.8 0.0068 830.4 1527 10010 135 0]'; % initial condition 
    case 2016
        x0 = [0.85 0.0068 851 460 1740 223 0]'; % initial condition
end
x0 = [x0;0]; % extra state for objective computation

nx = length(x0); % # of states
nu = 2; % # of inputs
nv = 2; % # of variables of interest (consumption & emission) needed for efficient computation of SCC 

% ========================================================================
% construct guess for decision variables of Nonlinear Program (NLP)
% =========================================================================
u_guess(1,:) = [0 0 0 0.5*ones(1,N-5) 0 0 0];
u_guess(2,:) = [0.5*ones(1,N) 0];

x_guess = zeros(nx, N+1);

for i = 1:N
    if i == 1
        x_guess(:,1) = x0(1:nx);
    end
  [fi, hi] = dice_dynamics(x_guess(:,i),u_guess(:,i),i,Params);
  x_guess(:, i+1) = fi;
  v_guess(:, i)   = hi;
end
v_guess(:, N+1)   = hi;

xuv_guess = [x_guess(:);u_guess(:);v_guess(:)];

% ========================================================================
% construct bounds on decision variables of NLP
% =========================================================================
u_LB = [zeros(1,N+1); [zeros(1,N-10) Params.optlsrv*ones(1,10) 0]];
%u_LB = [zeros(1,N+1); [zeros(1,N) 0]];
%u_UB = [[Params.miu0 1.0*ones(1,28) 1*ones(1,N-28)]; ones(1,N+1)];
u_UB = [zeros(1,3) ones(1,N-2); ones(1,N+1)];
u_LB = u_LB(:,1:N+1);
u_UB = u_UB(:,1:N+1);

x_LB = zeros(nx, N+1);
x_LB(end,:) = -inf; % no lower bound on objective 
x_UB = [inf*ones(nx, N+1)]; % no upper bound on states 
if Params2.p == 2
    %x_UB(1,:) = 2*ones(1,N+1); % 2 degree upper bound
    x_UB(3,:) = 550*2.13*ones(1,N+1); % 2 degree upper bound
end 

v_LB = -inf*ones(nv, N+1);
v_UB =  inf*ones(nv, N+1);

xuv_LB = [x_LB(:); u_LB(:); v_LB(:)];
xuv_UB = [x_UB(:); u_UB(:); v_UB(:)];

%% ========================================================================
% define NLP
% =========================================================================
% allocate CASADI variables
x      = SX.sym('x', nx, N+1); % states + value 
u      = SX.sym('u', nu, N+1); % inputs
v      = SX.sym('v', nv, N+1); % variables of interest, i.e. emissions & consumption
eq_con = SX.sym('eq_con', nx+nv, N+1); % (nx+nv) * N+1 constraints for the dynamics and variables of interest

% loop over dynamics 
for  i = 1:N
    if  i == 1
       eq_con(1:nx,1) = x(1:nx,1) - x0(1:nx); % x(:,1) = x0;    
    end
    % equality constraints for dynamics
    [fi, hi] = dice_dynamics(x(1:nx,i),u(:,i),i,Params,v(:,i));
    eq_con(1:nx,i+1) = x(1:nx,i+1) - fi; 
    
    % equality constraints assigning emmissions and consumption 
    eq_con(nx+1:end, i) = (v(:,i) - hi);
    if i == N
     eq_con(nx+1:end, i+1) = v(:,i+1);    % dummy constraint at i = N+1
    end
    
end

% define the objective (Mayer term)
obj = ((5 * Params.scale1 * x(end, N+1)) - Params.scale2);

% define nlp variables
nlp = struct('x', [x(:);u(:);v(:)], 'f', obj, 'g', eq_con(:));

% options for IPOPT
opts = struct;
opts.ipopt.max_iter    = 3000;
opts.ipopt.print_level = 3;%0,3
opts.print_time        = 0;
opts.ipopt.acceptable_tol =1e-12;
opts.ipopt.acceptable_obj_change_tol = 1e-12;
opts.ipopt.mu_strategy                 = 'adaptive';

% create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp);

% solve the NLP

lbg = zeros(nx+nv, N+1);
lbg(end-1:end, end) =  0;% % bound for dummy constraint
ubg = zeros(nx+nv, N+1); 
ubg(end-1:end, end) =  0; % bound for dummy constraint

res = solver( 'x0' , xuv_guess,...      % solution guess
              'lbx', xuv_LB,...         % lower bound on x
              'ubx', xuv_UB,...         % upper bound on x
              'lbg', 0,...              % lower bound on g
              'ubg', 0);                % upper bound on g


%% ========================================================================
% assign and plot results
% =========================================================================
s = 3

x_res = full(res.x);
% reshape solution to matrix variables
v_opt4 = x_res(end-length(v(:))+1:end);
v_opt4 = reshape(v_opt4, nu, N+1);

u_opt4 = x_res(end-length(v(:))-length(u(:))+1:end-length(v(:)));
u_opt4 = reshape(u_opt4, nu, N+1);
u_opt4 = u_opt4(:,1:end-1); % last value is not meaningful

x_opt4 = x_res(1:end-length(v(:))-length(u(:)));
x_opt4 = reshape(x_opt4, nx, N+1);

% Marginals and Social Cost of Carbon
lam  = full(res.lam_g);
lam  = reshape(lam, nx+nv, N+1);
lamE = lam(end-1,1:end);
lamC = lam(end,1:1:end);
if Param_set == 2013
    SCC(s,:) = -1000*[lamE(1:end-1)]./[.00001+lamC(1:end-1)];
else
    SCC(s,:) = -1000*[lamE(1:end-1)]./[lamC(1:end-1)];
end

%% Default plot results
s = 3

x_res = full(res.x);
% reshape solution to matrix variables
v_opt4 = x_res(end-length(v(:))+1:end);
v_opt4 = reshape(v_opt4, nu, N+1);

u_opt4 = x_res(end-length(v(:))-length(u(:))+1:end-length(v(:)));
u_opt4 = reshape(u_opt4, nu, N+1);
u_opt4 = u_opt4(:,1:end-1); % last value is not meaningful

x_opt4 = x_res(1:end-length(v(:))-length(u(:)));
x_opt4 = reshape(x_opt4, nx, N+1);

% Marginals and Social Cost of Carbon
lam  = full(res.lam_g);
lam  = reshape(lam, nx+nv, N+1);
lamE = lam(end-1,1:end);
lamC = lam(end,1:1:end);
if Param_set == 2013
    SCC4(1,:) = -1000*[lamE(1:end-1)]./[.00001+lamC(1:end-1)];
else
    SCC4(1,:) = -1000*[lamE(1:end-1)]./[lamC(1:end-1)];
end