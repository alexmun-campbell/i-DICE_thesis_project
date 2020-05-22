% -------------------------------------------------------------------------
%  This file is based on the DICE. I've started to adapt it to use an
%  investment-based abatement decision. 

clc
clear all
close all

% import casadi name space
import casadi.*

%% ========================================================================
% define data of the DICE Optimal Control Problem (OCP)
% =========================================================================
Param_set = 2013; % Choose DICE 2013, 2016, or 1 for custom values 
Params2 = I_DICE_params(Param_set); % get DICE parameters
Params3 = I_DICE_params(Param_set);

%%

Params = Params3;
N = Params.N;

switch Param_set
    case 1
        x0 = [0.8 0.0068 830.4 1527 10010 135 0]'; % initial condition 
    % IT IS RECOMMENDED THAT CHANGES BE LIMITED TO THE ABOVE!!!
    % Initial conditions below correspond to DICE2013R and DICE2016R.
    case 2013
        x0 = [0.8 0.0068 830.4 1527 10010 135 0.039 0]'; % initial condition 
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
u_guess(1,:) = [0 zeros(1,N-3) 0 0 0];
u_guess(2,:) = [0.2*ones(1,N) 0];
%u_guess(2,:) = [0*ones(1,N) 0];

x_guess = zeros(nx, N+1);

for i = 1:N
    if i == 1
        x_guess(:,1) = x0(1:nx);
    end
  [fi, hi] = I_DICE_dynamics(x_guess(:,i),u_guess(:,i),i,Params);
  x_guess(:, i+1) = fi;
  v_guess(:, i)   = hi;
end
v_guess(:, N+1)   = hi;

xuv_guess = [x_guess(:);u_guess(:);v_guess(:)];

% ========================================================================
% construct bounds on decision variables of NLP
% =========================================================================
%u_LB = [zeros(1,N+1); [zeros(1,N-10) Params.optlsrv*ones(1,10) 0]];
u_LB = [zeros(1,N+1); [zeros(1,N) 0]];
%u_LB = [zeros(1,N+1); [zeros(1,N) 0]];
u_UB = [0.5*ones(1,N+1); ones(1,N+1)];
%u_UB = [ones(1,N+1); ones(1,N+1)];
u_LB = u_LB(:,1:N+1);
u_UB = u_UB(:,1:N+1);

x_LB = zeros(nx, N+1);
x_LB(end,:) = -inf; % no lower bound on objective 
x_UB = [inf*ones(nx, N+1)]; % no upper bound on states 
if Params2.p == 2
    %x_UB(1,:) = 2*ones(1,N+1); % 2 degree upper bound
    x_UB(3,:) = 550*2.13*ones(1,N+1); % 450 degree upper bound
end 
x_UB(7,:) = 1; % mu state is limited 

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

tol = 0.0001;
theta = Params.theta1(1);
cost_aim = C(7,15);
cost = cost_aim+10
z = 1;
while abs(cost - cost_aim) > tol
% loop over dynamics
Params.theta1 = ones(1,N+1).*theta;
for  i = 1:N
    if  i == 1
       eq_con(1:nx,1) = x(1:nx,1) - x0(1:nx); % x(:,1) = x0;    
    end
    % equality constraints for dynamics
    [fi, hi] = I_DICE_dynamics(x(1:nx,i),u(:,i),i,Params,v(:,i));
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
          
x_res = full(res.x);
xtemp = x_res(1:end-length(v(:))-length(u(:)));
xtemp = reshape(xtemp, nx, N+1);
cost = xtemp(8,15);
theta = theta*(cost_aim/cost);
cost_aim/cost
end 