% define initial conditions
total_start_conf = 4.82;
G10 = 0.36 * total_start_conf;
S0 = 0.14 * total_start_conf;
SD0 = 0.04 * total_start_conf;
G20 = 0.46 * total_start_conf;
death0 = 0;

y0 = [G10, S0, SD0, G20, death0]';
time = [0:1:310]';

%define the drug doses for the ATRi and the PARPi
ATRidose =  [0,0,0,0.1,0.3,1,0.1,0.3,0.3,1,1];
PARPidose = [0,0.2,1,0,0,0,1,0.2,1,1,0.2 ];

model_output=zeros(length(time),length(ATRidose));

for i = 1:numel(ATRidose)
    [tout, yout] = ode45(@(t, y) models_ddr(t, y, ATRidose(i), PARPidose(i)), time, y0);
    model_output(:,i) = sum(yout,2);
end


figure;

subplot(2, 3, 1);
plot(time,model_output(:,1),'k','LineWidth',2)

xlabel('Time (hours)')
ylabel('Cell Confluency')
legend('DMSO Sim','Location','northwest');
set(gca,'FontSize',15)
title('(I) Cell Confluency - DMSO')
xlim([0, 310])
ylim([0,100])


subplot(2, 3, 2);

plot(time,model_output(:,2),'r','LineWidth',2)
hold on
plot(time,model_output(:,3),'color', [1, 0.5, 0],'LineWidth',2)

xlabel('Time (hours)')
ylabel('Cell Confluency')
legend('0.2\muM PARPi Sim','1\muM PARPi Sim','Location','northwest');
set(gca,'FontSize',15)
title('(II) Cell Confluency - PARPi')
xlim([0, 310])


subplot(2, 3, 3);

plot(time,model_output(:,4),'color',[0.49, 0.18, 0.56],'LineWidth',2)
hold on
plot(time,model_output(:,5),'g','LineWidth',2)
hold on
plot(time,model_output(:,6),'c','LineWidth',2)

xlabel('Time (hours)')
ylabel('Cell Confluency')
legend('0.1\muM ATRi Sim','0.3\muM ATRi Sim','1\muM ATRi Sim','Location','northwest');
set(gca,'FontSize',15)
title('(III) Cell Confluency - ATRi')
xlim([0, 310])
ylim([0,100])


subplot(2, 3, 4);

plot(time,model_output(:,7),'b','LineWidth',2)

xlabel('Time (hours)')
ylabel('Cell Confluency')
legend('0.1\muM ATRi + 1\muM PARPi Sim','Location','northwest');
set(gca,'FontSize',15)
title('(IV) Cell Confluency - Combinations')
xlim([0, 310])


subplot(2, 3, 5);

plot(time,model_output(:,8),'color',[0.93, 0.69, 0.13],'LineWidth',2)
hold on
plot(time,model_output(:,9),'color',[0.72, 0.27, 1.00],'LineWidth',2)

xlabel('Time (hours)')
ylabel('Cell Confluency')
legend('0.3\muM ATRi + 0.2\muM PARPi Sim','0.3\muM ATRi + 1\muM PARPi Sim','Location','northwest');
set(gca,'FontSize',15)
title('(V) Cell Confluency - Combinations')
xlim([0, 310])


subplot(2, 3, 6);

plot(time,model_output(:,10),'color',[0.64, 0.08, 0.18],'LineWidth',2)
hold on
plot(time,model_output(:,11),'b','LineWidth',2)

xlabel('Time (hours)')
ylabel('Cell Confluency')
legend('1\muM ATRi + 0.2\muM PARPi Sim','1\muM ATRi + 1\muM PARPi Sim','Location','northwest');
set(gca,'FontSize',15)
title('(VI) Cell Confluency - Combinations')
xlim([0, 310])



function dydt = models_ddr(t, y, ATRidose, PARPidose)
G1 = y(1);
S = y(2);
SD= y(3);
G2 = y(4);
death = y(5);

p0=0.3558;
k3=0.1420;
scale_par=1.824;
h_O=1.0962;
EC50_O=0.1275;
Emax_O=0.5609;
h_C=1.5187;
EC50_C=0.2579;
Emax_C=1;
cp=100;
q0=1;

k1_scaled=scale_par*0.0678;  %rate at which cells leave y(1)
k2_scaled=scale_par*0.1742;  %rate at which cells leave y(2)
k4_scaled=scale_par*0.0530;  %rate at which cells leave y(4)

%Scale all rate parameters to logistic growth
pop_size = G1+S+SD+G2+death;
k1_scaled=k1_scaled*(1-pop_size/cp);
k2_scaled=k2_scaled*(1-pop_size/cp);
k3_scaled=k3*(1-pop_size/cp);
k4_scaled=k4_scaled*(1-pop_size/cp);

% Calculate the drug effects
eff_PAPRi = Emax_O * PARPidose ^ h_O / (PARPidose ^ h_O + EC50_O ^ h_O );
eff_ATRi =  Emax_C * ATRidose  ^ h_C / (ATRidose  ^ h_C + EC50_C ^ h_C );
eff_comb = eff_PAPRi + eff_ATRi - eff_PAPRi * eff_ATRi;

% Define the equations
dydt = zeros(5, 1);
dydt(1) = 2*k4_scaled*y(4) - k1_scaled*y(1);
dydt(2) = k1_scaled*p0*(1-eff_PAPRi)*y(1) - k2_scaled*y(2) + k3*(q0*(1-eff_comb))*y(3);
dydt(3) = k1_scaled*(1-p0*((1-eff_PAPRi)))*y(1) - k3*y(3);
dydt(4) = k2_scaled*y(2) - k4_scaled*y(4);
dydt(5) = k3_scaled*(1-q0*(1-eff_comb))*y(3);

end

