function plot_results2(u, x, u2, x2, SCC, SCC2, Params, Params2, Param_set,type)
% Plot the solution of DICE2013R_mc

N           = Params.N;
BaseYear    = Params.BaseYear;
sigma       = Params.sigma;
A_TFP       = Params.A_TFP;
L           = Params.L;
theta1      = Params.theta1;
F_EX        = Params.F_EX;
E_Land      = Params.E_Land;
p           = Params2.p;
q           = Params2.q;


if p == 1
    T = 23; % this is for endogenous policy 
end 

if p == 2
    T = 11; % this is for 2 degrees C policy 
end 

T = 19

nx = length(x(:,1));
nx2 = length(x2(:,1));

% Translate time indices into years
end_year = BaseYear + (N+1)*5;
years = BaseYear:5:end_year;
EndYear = BaseYear + N*5;   % Used for plotting axes
EndYear2 = BaseYear + T*5;

% Create cost flows variable
C(1,:) = x(7,:);
C(2,:) = x2(8,:);
c_flow(:,1) = C(:,1);
for i = 2:length(C)
    c_flow(:,i)=C(:,i)-C(:,i-1);
end 

% Compute auxiliary outputs
[Gross_Economic_Output(1,:), E_Industrial(1,:), Net_Output(1,:), Per_Cap_Consumption(1,:), ...
Damages_Fraction(1,:), Atm_Carbon_ppm(1,:), ... 
Marg_Cost_Abatement(1,:), Consumption(1,:)] = ...
compute_auxiliary_quantities(x,u,Params,nx);
[Gross_Economic_Output(2,:), E_Industrial(2,:), Net_Output(2,:), Per_Cap_Consumption(2,:), ...
Damages_Fraction(2,:), Atm_Carbon_ppm(2,:), ... 
Marg_Cost_Abatement(2,:), Consumption(2,:)] = ...
i_compute_auxiliary_quantities(x2,u2,Params2,nx2);

%load('/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/Cons_baseline.mat')
%Consumption_loss = (Consumption(1:2,:)./Cons_baseline(1:N))-1

figure(1), title('Environmental variables')

subplot(2,3,1) 
plot(years(1:end-1),x(1,:))
hold on
plot(years(1:end-1),x2(1,:))
ylabel('TAT')
xlabel('years') 
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,3,2) 
plot(years(1:end-1),x(2,:))
hold on
plot(years(1:end-1),x2(2,:))
ylabel('TLO')
xlabel('years'),
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,3,3) 
plot(years(1:end-1),x(3,:))
hold on
plot(years(1:end-1),x2(3,:))
ylabel('MAT')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,3,4) 
plot(years(1:end-1),x(4,:))
hold on
plot(years(1:end-1),x2(4,:))
ylabel('MUP')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,3,5) 
plot(years(1:end-1),x(5,:))
hold on
plot(years(1:end-1),x2(5,:))
ylabel('MLO')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,3,6) 
plot(years(1:end-2),Atm_Carbon_ppm(1,:))
hold on
plot(years(1:end-2),Atm_Carbon_ppm(2,:))
ylabel('Atmospheric Carbon ppm')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

if type == 1
figure(2), title('Economic variables') 

subplot(2,5,1) 
plot(years(1:end-1),x(6,:))
hold on
plot(years(1:end-1),x2(6,:))
title('Capital')
ylabel('USD')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,2) 
plot(years(1:end-2),u(1,:))
hold on
plot(years(1:end-2),u2(1,:))
title('Emissions Control Rate')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,3) 
plot(years(1:end-2),u(2,:))
hold on
plot(years(1:end-2),u2(2,:))
title('Savings Rate')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,4) 
plot(years(1:end-2),u(1,:))
hold on
plot(years(1:end-2),x2(7,(1:(end-1))))
title('Mu')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,5) 
plot(years(1:end-2),SCC)
hold on 
plot(years(1:end-2),SCC2)
xlabel('years')
if Param_set == 2013
    ylabel('2005 USD')
else
    ylabel('2010 USD')
end
title('Social Cost of Carbon')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,6)
plot(years(1:11),c_flow(1,1:11))
hold on 
plot(years(1:11),c_flow(2,1:11))
title('Financial Effort')
ylabel('USD')
xlabel('years')

subplot(2,5,7) 
plot(years(1:end-2),Net_Output(1,:))
hold on
plot(years(1:end-2),Net_Output(2,:))
title('Net Output')
ylabel('USD Trillions?')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,8) 
plot(years(1:end-2),Per_Cap_Consumption(1,:))
hold on
plot(years(1:end-2),Per_Cap_Consumption(2,:))
title('Per Capita Consumption')
ylabel('USD')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,9) 
plot(years(1:end-2),Damages_Fraction(1,:))
hold on
plot(years(1:end-2),Damages_Fraction(2,:))
title('Damages Fraction')
ylabel('% Output')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

subplot(2,5,10) 
plot(years(1:end-2),Marg_Cost_Abatement(1,:))
hold on
plot(years(1:end-2),Marg_Cost_Abatement(2,:))
title('Marginal Cost of Abatement')
ylabel('$/tCO2')
xlabel('years')
temp = axis; axis([BaseYear EndYear temp(3) temp(4)]);

else 
figure(3), title('Economic variables') 

subplot(2,5,1) 
plot(years(1:T),x(6,1:T))
hold on
plot(years(1:T),x2(6,1:T))
title('Capital')
ylabel('USD')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,2) 
plot(years(1:T),u(1,1:T))
hold on
plot(years(1:T),u2(1,1:T))
title('Emissions Control Rate')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,3) 
plot(years(1:T),u(2,1:T))
hold on
plot(years(1:T),u2(2,1:T))
title('Savings Rate')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,4) 
plot(years(1:T),u(1,1:T))
hold on
plot(years(1:T),x2(7,1:T))
title('Mu')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,5) 
plot(years(1:T),SCC(1:T))
hold on 
plot(years(1:T),SCC2(1:T))
xlabel('years')
if Param_set == 2013
    ylabel('2005 USD')
else
    ylabel('2010 USD')
end
title('Social Cost of Carbon')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,6)
plot(years(1:T),c_flow(1,1:T))
hold on 
plot(years(1:T),c_flow(2,1:T))
title('Financial Effort')
ylabel('USD')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,7) 
plot(years(1:T),Net_Output(1,1:T))
hold on
plot(years(1:T),Net_Output(2,1:T))
title('Net Output')
ylabel('USD Trillions?')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

%subplot(2,5,8) 
%plot(years(1:T),Consumption_loss(2,1:T))
%title('Per Capita Consumption')
%ylabel('USD')
%xlabel('years')
%temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,9) 
plot(years(1:T),Damages_Fraction(1,1:T))
hold on
plot(years(1:T),Damages_Fraction(2,1:T))
title('Damages Fraction')
ylabel('% Output')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(2,5,10) 
plot(years(1:T),Marg_Cost_Abatement(1,1:T))
hold on
plot(years(1:T),Marg_Cost_Abatement(2,1:T))
title('Marginal Cost of Abatement')
ylabel('$/tCO2')
xlabel('years')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

end 