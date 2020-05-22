function [Gross_Economic_Output, E_Industrial, Net_Output, Per_Cap_Consumption, Damages_Fraction, ...
    Atm_Carbon_ppm, Marg_Cost_Abatement, Consumption] = ...
    i_compute_auxiliary_quantities(x_opt,u_opt, Params,nx);
% Computes additional outputs based on the solution of the DICE
%   optimal control problem.
% Model type is indicated by p. 1 = endogenous policy. 2 = 2 degrees C
% specification is indicated by q. 1 = only transition costs. 2 = linear
% level costs
% N2 : period over which costs are equalised 


% DICE inputs: mitigation rate and savings rate
psi = u_opt(1,:);
s  = u_opt(2,:);

% DICE endogenous states
TATM = x_opt(1,:);
TLO = x_opt(2,:);
MATM = x_opt(3,:);
MUP = x_opt(4,:);
MLO = x_opt(5,:);
K = x_opt(6,:);
mu = x_opt(7,:);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                           Unpack Parameters
% -------------------------------------------------------------------------

N           = Params.N;
sigma       = Params.sigma;
A_TFP       = Params.A_TFP;
L           = Params.L;
theta1      = Params.theta1;
F_EX        = Params.F_EX;
E_Land      = Params.E_Land;
eta         = Params.eta;
M_AT_Base   = Params.M_AT_Base;
delta       = Params.delta;
gamma       = Params.gamma;
theta2      = Params.theta2;
a2          = Params.a2;
a3          = Params.a3;
pb          = Params.pb;
deltaPB     = Params.deltaPB;
alpha       = Params.alpha;
rho         = Params.rho;
xi1         = Params.xi1;
xi2         = Params.xi2;
Phi_T       = Params.Phi_T;
Phi_M       = Params.Phi_M;
zeta11      = Phi_M(1,1);
zeta12      = Phi_M(1,2);
deltaEA     = Params.deltaEA;

p           = Params.p;
q           = Params.q;

if p == 1
    T = 23; % this is for endogenous policy 
end 

if p == 2
    T = 11; % this is for 2 degrees C policy 
end 

%cut = N; % this is when LCCC is completely foreward looking
cut = T; % this is when LCCC is only worked out across the cost period




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                           Compute quantities of interest
% -------------------------------------------------------------------------
for i=1:N
    Gross_Economic_Output(i) = A_TFP(i)*(K(i)^gamma)*((L(i)/1000)^(1-gamma));
end 

% GDP growth
Yg = ones(1,1);
for i = 1:(N-1)
    Yg(:,i) = (Gross_Economic_Output(:,i+1)-Gross_Economic_Output(:,i))./Gross_Economic_Output(:,i);
end 
Yg = [0, Yg, ones(1,N+1)*Yg(end)];
Yg1 = 1 + Yg;
Y0 = Gross_Economic_Output(1);
for i = 1:(2*N)
    Y(i)=Y0*prod(Yg1(1:i));
end 
for i = 1:N
    deflate(i) = ((1-rho)*(1-deltaEA))^(i-1);
end 
for i = 1:N
    deflate2(i) = ((1-rho)*(1-deltaEA)*(1-Yg(i)))^(i-1);
end 
for i = 1:N
    denom(i) = sum(deflate.*sigma(i).*Y(i:(i+N-1)));
end 
for i = 1:N
    denom2(i) = sum(deflate.*sigma(i));
end 
for i = 1:N
    denom3(i) = sum(deflate2.*sigma(i));
end 


for i=1:N
    Damages(i) = 1 - a2*(TATM(i)^a3); 
    Damages_Fraction(i) = 1 - Damages(i);
    if q == 1
        Costs(i) = theta1(i)*(psi(i)^theta2);
    end 
    if q == 2
        Costs(i) = mu(i)*theta1(i)*(psi(i)^theta2);
    end 
    if q == 3
        Costs(i) = theta1(i)*(psi^theta2)+theta3*(mu^theta4);
    end 
    Net_Output(i) = (Damages(i) - Costs(i))*Gross_Economic_Output(i);
    E_Industrial(i) = sigma(i)*(1-mu(i))*Gross_Economic_Output(i); 
    Per_Cap_Consumption(i) = 1000*(1-s(i)) * Net_Output(i) / L(i);
    Atm_Carbon_ppm(i) = MATM(i)/2.13;
    Consumption(i) = Net_Output(i) * (1-s(i));
    
    %Marg_Cost_Abatement(i) = 1000 * (theta2*Costs(i)*Y(i))/(psi(i)*denom(i));
    Marg_Cost_Abatement(i) = 1000 * (theta2*Costs(i))/(psi(i)*denom2(i));
    %Marg_Cost_Abatement(i) = 1000 * (theta2*Costs(i)/(psi(i)*denom3(i)));
end