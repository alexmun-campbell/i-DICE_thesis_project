function plot_results3(u, x, u2, x2, u3, x3, SCC, SCC2, SCC3, Params, Params2, Params3, Param_set,type)

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

T = 16;

nx = length(x(:,1));
nx2 = length(x2(:,1));
nx3 = length(x3(:,1));

% Translate time indices into years
end_year = BaseYear + (N+1)*5;
years = BaseYear:5:end_year;
EndYear = BaseYear + N*5;   % Used for plotting axes
EndYear2 = BaseYear + T*5;

% Create cost flows variable
C(1,:) = x(7,:);
C(2,:) = x2(8,:);
C(3,:) = x3(8,:);
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
[Gross_Economic_Output(3,:), E_Industrial(3,:), Net_Output(3,:), Per_Cap_Consumption(3,:), ...
Damages_Fraction(3,:), Atm_Carbon_ppm(3,:), ... 
Marg_Cost_Abatement(3,:), Consumption(3,:)] = ...
i_compute_auxiliary_quantities(x3,u3,Params3,nx3);


load('/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/Cons_baseline.mat')
Consumption_loss = 100*((Consumption(1:3,:)./Cons_baseline(1:N))-1)

figure(1)

subplot(4,2,1) 
plot(years(1:T),x(1,1:T))
hold on
plot(years(1:T),x2(1,1:T))
plot(years(1:T),x3(1,1:T))
title('Average Air Temperature ')
ylabel('Degrees Centigrade')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);
leg = legend('DICE','I-DICE central','I-DICE transition','Location','northoutside')
newPosition = [0.35 0.03 0.3 0.025];
newUnits = 'normalized';
set(leg,'Position', newPosition,'Units', newUnits);



subplot(4,2,2) 
plot(years(1:T),Atm_Carbon_ppm(1,1:T))
hold on
plot(years(1:T),Atm_Carbon_ppm(2,1:T))
plot(years(1:T),Atm_Carbon_ppm(3,1:T))
title('Atmospheric Emissions')
ylabel('ppm')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(4,2,3) 
plot(years(1:T),SCC(1:T))
hold on 
plot(years(1:T),SCC2(1:T))
plot(years(1:T),SCC3(1:T))
if Param_set == 2013
    ylabel('2005 USD')
else
    ylabel('2010 USD')
end
title('Social Cost of Emissions')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(4,2,4) 
plot(years(1:T),Marg_Cost_Abatement(1,1:T))
hold on
plot(years(1:T),Marg_Cost_Abatement(2,1:T))
plot(years(1:T),Marg_Cost_Abatement(3,1:T))
title('Marginal Cost of Abatement / LCCC')
ylabel('$/tCO2')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(4,2,5)
plot(years(1:T),c_flow(1,1:T))
hold on 
plot(years(1:T),c_flow(2,1:T))
plot(years(1:T),c_flow(3,1:T))
title('Financial Effort')
ylabel('$ trillions')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(4,2,6) 
plot(years(1:T),u(1,1:T))
hold on
plot(years(1:T),x2(7,1:T))
plot(years(1:T),x3(7,1:T))
title('Mu')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(4,2,7) 
plot(years(1:T),u(2,1:T))
hold on
plot(years(1:T),u2(2,1:T))
plot(years(1:T),u3(2,1:T))
title('Savings Rate')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);

subplot(4,2,8) 
plot(years(1:T),Consumption_loss(1,1:T))
hold on
plot(years(1:T),Consumption_loss(2,1:T))
plot(years(1:T),Consumption_loss(3,1:T))
title('% Consumption loss')
ylabel('%')
temp = axis; axis([BaseYear EndYear2 temp(3) temp(4)]);


set(gcf, 'Position',  [0,0, 600, 900])
%annotation('textbox',[0, 0, 1, 1],'String','epsilon = 0.01','EdgeColor','None')


%saveas(gcf,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Document/figures/final_p2_450.png')

end 