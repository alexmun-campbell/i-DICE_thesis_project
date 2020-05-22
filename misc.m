%% MISC 
% This script contains odds and ends, such as calibration workings and
% random diagnostics

% Calculate output 
Params = set_DICE_parameters(Param_set); % get DICE parameters
Params2 = I_DICE_params(Param_set); % get i-DI

%%
J = zeros(1,1);
mu = zeros(1,N);
S = zeros(1,N);
K = zeros(1,N+1);
C = zeros(1,N+1);
psi = zeros(1,N);
%%
z=3
p = Params2.p;
q = Params2.q;
s=3
a = (s-1)*3+z; % this is an indexing device 

% DICE inputs: mitigation rate and savings rate
if z==1 
    J(a,1) = -x_opt(8,end);
    mu(a,:) = u_opt(1,:);
    S(a,:)  = u_opt(2,:);
    K(a,:) = x_opt(6,:);
    C(a,:) = x_opt(7,:);
    psi(a,:) = zeros(1,length(u_opt));
    
    % Compute auxiliary outputs
    [Gross_Output(a,:), E_Industrial(a,:), Net_Output(a,:),...
        Per_Cap_Consumption(a,:), Damages_Fraction(a,:), ...
    Atm_Carbon_ppm(a,:), Marg_Cost_Abatement(a,:), Consumption(a,:)] = ...
    compute_auxiliary_quantities(x_opt,u_opt,Params1,nx);
end 

if z==2
    J(a,:) = -x_opt2(nx,end);
    psi(a,:) = u_opt2(1,:);
    S(a,:)  = u_opt2(2,:);
    K(a,:) = x_opt2(6,:);
    mu(a,:) = x_opt2(7,1:(end-1));
    C(a,:) = x_opt2(8,:);
    
    % Compute auxiliary outputs
    [Gross_Output(a,:), E_Industrial(a,:), Net_Output(a,:), Per_Cap_Consumption(a,:), ...
    Damages_Fraction(a,:), Atm_Carbon_ppm(a,:), ... 
    Marg_Cost_Abatement(a,:), Consumption(a,:)] = ...
    i_compute_auxiliary_quantities(x_opt2,u_opt2,Params2,nx);
end 

if z==3
    J(a,:) = -x_opt3(nx,end);
    psi(a,:) = u_opt3(1,:);
    S(a,:)  = u_opt3(2,:);
    K(a,:) = x_opt3(6,:);
    mu(a,:) = x_opt3(7,1:(end-1));
    C(a,:) = x_opt3(8,:);
    
        % Compute auxiliary outputs
    [Gross_Output(a,:), E_Industrial(a,:), Net_Output(a,:), Per_Cap_Consumption(a,:), ...
    Damages_Fraction(a,:), Atm_Carbon_ppm(a,:), ... 
    Marg_Cost_Abatement(a,:), Consumption(a,:)] = ...
    i_compute_auxiliary_quantities(x_opt3,u_opt3,Params3,nx);
end 
%% Calculate effort 

theta1 = Params.theta1
theta2 = Params.theta2
sigma = Params.sigma

plot(C(1,:))
hold on 
plot(C(2,:))
plot(C(3,:))


%% 
mu_temp = 0.01:0.01:1;

% Total abatement costs chart
TAC = theta1(1:end-1)'*(mu_temp.^theta2);
plot(mu_temp,TAC(1,:))
hold on 
plot(mu_temp,TAC(20,:))
plot(mu_temp,TAC(40,:))
plot(mu_temp,TAC(60,:))
%scatter(test(2,:),test(1,:))
legend({'T=1','T=20','T=40','T=60'},'location','northwest')
ylabel('Abatement costs as % output')
xlabel('Level of abatement')
saveas(gcf,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Document/figures/DICE_TAC.png')

% Generalised MAC for DICE 
for i = 1:60
    MAC_DICE(i,:) = Params.pb*(1-Params.deltaPB)^(i-1).*(mu_temp.^(theta2-1));
end 
plot(mu_temp,MAC_DICE(1,:))
hold on 
plot(mu_temp,MAC_DICE(20,:))
plot(mu_temp,MAC_DICE(40,:))
plot(mu_temp,MAC_DICE(60,:))
legend({'T=1','T=20','T=40','T=60'},'location','northwest')
title('DICE: Marginal Abatement Costs')
ylabel('$/tCO2')
xlabel('Level of abatement')
saveas(gcf,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Document/figures/DICE_MAC.png')

% Specific MAC for DICE
MAC1 = 1000*(theta1(1:(end-1)).*theta2./sigma(1:(end-2))).*(mu(1,:).^(theta2-1));
plot(MAC1)

% Manual MAC for DICE

Gross_Emissions = sigma(1:60).*Gross_Output(1:60)./1000

z=1
for i = 1:20
    MAC11(i)=((C(z,i+1)-C(z,i))/(mu(z,i+1)-mu(z,i)))./sigma(i);
end 

plot(MAC1(1:11))
hold on
plot(MAC11)


plot(MAC11)
hold on 
plot(C(1,1:11))
plot(mu(1,1:11))

% Marginal abatement costs chart for i-DICE 
Params = I_DICE_params(Param_set);
for i = 1:60
    MAC2(i,:) = Params.pb*(1-Params.deltaPB)^(i-1).*(mu_temp.^(theta2-1));
end


%% Calculate discounted costs 

beta = Params.rho;
s=2

for i = 1:60
    C_flow(:,i)=C(:,i+1)-C(:,i);
end 

cut = 11
plot(C_flow(1,1:cut))
hold on 
plot(C_flow(2,1:cut))
plot(C_flow(3,1:cut))


% Sum of the costs of the first 11 periods, which is the legnth of time
% that the economy is transitioning to ful decarbonisation


plot(c_flow(1,1:11))
hold on 
plot(c_flow(2,1:11))
legend('DICE','i-DICE')
title('Endogenous abatement comparison')
saveas(gcf,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Document/figures/effort_compare4.png')

% And now we calculate the flat cost required in i-DICE in order to equal
% the discounted sum of DICE costs 

%% THIS IS TO EXPORT TO R

writematrix(mu,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/mu.csv')

writematrix(c_flow,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/cflow.csv')

Cons_baseline = Consumption
save('/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/Cons_baseline.mat','Cons_baseline')

importdata('/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/Cons_baseline.mat','Cons_baseline')

writematrix(Consumption,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/Consumption2.csv')

writematrix(ans,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Code/Data_analysis/Cons_baseline.csv')
Cons_baseline = ans

asdf = 1

%%
rho = Params.rho
delta = 0.005;
sig = Params.sigma(1)

Yg = ones(1,1);
for i = 1:(N-1)
    Yg(:,i) = (Gross_Output(:,i+1)-Gross_Output(:,i))./Gross_Output(:,i);
end 

Yg = [0, Yg, ones(1,60)*Yg(end)]; 
Yg1 = 1 + Yg

Y0 = Gross_Output(1);
for i = 1:(2*N)
    Y(i)=Y0*prod(Yg1(1:i));
end 
for i = 1:N
    deflate(i) = ((1-rho)*(1-delta))^(i-1);
end 

for j = 1:N
    denom(j) = sum(deflate.*sigma(1:N).*Y(j:(j+N-1)));
end 

Sig = Params1.sigma;
Sig_g = ones(1,1);
for i = 1:(N)
    Sig_g(:,i) = (Sig(:,i+1)-Sig(:,i))./Sig(:,i);
end 



plot(Yg)
plot(Y)
hold on 
plot(Gross_Output)


(1-rho)*(1-delta)
(1-rho-delta)



for i = 1:60
    test1(i) = (1-rho-delta)^(i-1)*sig*Gross_Output(i)
end 
plot(test1)

%% 
 
for i = 1:N
    SCCdef(i) = (1+rho)^(-(i-1)*5)*SCC(i)
end
for i = 1:N
    SCCdef(i) = (1+Yg(i))^(-(i-1))*SCCdef(i)
end
for i = 1:N
    SCCdef(i) = (1+Sig_g(i))^(-(i-1))*SCCdef(i)
end
BaseYear = Params.BaseYear
end_year = BaseYear+38*5
years = BaseYear:5:end_year;
plot(years,SCC(1:39))
hold on 
plot(years,SCCdef(1:39))
legend({'Carbon price','Deflated carbon price'},...
    'location','northwest')
title('DICE carbon price, with endogenous policy')
ylabel('$/tCO2')
temp = axis; axis([BaseYear  end_year temp(3) temp(4)]);
%saveas(gcf,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Document/figures/Carbon_price_deflated.png')

%%

T = 16;
BaseYear = Params.BaseYear;
end_year = BaseYear+(T-1)*5;
years = BaseYear:5:end_year;

T2 = 8;
BaseYear2 = Params.BaseYear;
end_year2 = BaseYear+(T2-1)*5;
years2 = BaseYear2:5:end_year2;

figure
subplot(2,2,1)
plot(years,SCC2(1,1:T))
hold on 
plot(years,Marg_Cost_Abatement(2,1:T))
title('Endogenous policy, central costs')
ylim([0 150])
ylabel('$')
temp = axis; axis([BaseYear  end_year temp(3) temp(4)]);
leg = legend('SCC','LCCC','location','northwest')
newPosition = [0.35 0.03 0.3 0.025];
newUnits = 'normalized';
set(leg,'Position', newPosition,'Units', newUnits);



subplot(2,2,2)
plot(years,SCC3(1,1:T))
hold on 
plot(years,Marg_Cost_Abatement(3,1:T))
title('Endogenous policy, transitional costs')
ylim([0 150])
ylabel('$')
temp = axis; axis([BaseYear  end_year temp(3) temp(4)]);

subplot(2,2,3)
plot(years,SCC2(2,1:T))
hold on 
plot(years,Marg_Cost_Abatement(5,1:T))
title('450ppm policy, central costs')
ylim([0 300])
ylabel('$')
temp = axis; axis([BaseYear2  end_year2 temp(3) temp(4)]);

subplot(2,2,4)
plot(years,SCC3(2,1:T))
hold on 
plot(years,Marg_Cost_Abatement(6,1:T))
title('450ppm policy, transitional costs')
ylim([0 300])
ylabel('$')
temp = axis; axis([BaseYear2  end_year2 temp(3) temp(4)]);

saveas(gcf,'/Users/alexandercampbell/Documents/SSE_MECON/Thesis/Document/figures/SCC_LCCC.png')


%%

beta = 0.015;
ep = 0.02;
mu0 = 1;
for i = 1:200
    tester(i)=((1-beta)*(1-ep))^i;
end 
for i = 1:200
    tester2(i) = sum(tester(1:i));
end 

tester1 = ones(1,200)*(beta+ep)
plot(tester)
hold on 
plot(tester1)

