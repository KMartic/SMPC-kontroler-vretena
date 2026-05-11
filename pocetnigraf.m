clear all; 
close all; 
clc;

Tspan = 20;
T = 0.01;
t_s  = 0:0.01:Tspan;
Td = 0.5;
t_d = 0:Td:Tspan;


J1 = 10;
J2 = 640; 
c = 400; 
d = 20; 
ip = 12.5; 
FS = 18;


A = [-d/(J1*ip^2), -c/(J1*ip), d/(J1*ip);
     1/ip, 0, -1;
     d/(J2*ip), c/J2, -d/J2 ];
B = [1/J1; 0; 0];
C = [0 0 1];
D = [0];


%modeli
sys_c = ss(A,B,C,D);  
sys_d  = c2d(sys_c,Td,'zoh');  

num = [ip*d ip*c];
den = [(J1*J2*ip^2) (J1*ip^2*d + J2*d) (J1*ip^2*c + J2*c) 0];
G = tf(num,den);


% ulaz koji kasni
u_cont = zeros(length(t_s),1);
u_cont(t_s>=1) = 1;

u_disc = zeros(length(t_d),1);
u_disc(t_d>=1) = 1;

% odzivi
y_ss_cont = lsim(sys_c, u_cont, t_s);
y_ss_disc = lsim(sys_d,  u_disc, t_d);
y_tf = lsim(G, u_cont(:,1), t_s);

%% simulink
sim("busilica.slx")

%% GRAFOVI
figure (1)

subplot(3,1,1)
plot(t_s, y_ss_cont, 'LineWidth',2); 
hold on
stairs(t_d, y_ss_disc, 'LineWidth',2)
grid on
title('Odziv na jednadžbe prostora stanja','FontSize', FS)
xlabel('t [s]','FontSize', FS)
ylabel('ω_2','FontSize', FS)
legend('ω_2 - kontinuirani odziv p.s.','ω_2 - diskretni odziv p.s.','Location','best')


subplot(3,1,2)
plot(t_s, y_tf, 'LineWidth',2)
grid on
title('Odziv na prijenosnu funkciju','FontSize', FS)
xlabel('t [s]','FontSize', FS)
ylabel('ω_2','FontSize', FS)
legend('ω_2 - prijenosna funkcija','Location','best')


subplot(3,1,3)
plot(t, w2, 'LineWidth',2)
grid on
title('Odziv Simulink modela','FontSize', FS)
xlabel('t [s]','FontSize', FS)
ylabel('ω_2','FontSize', FS)
legend('ω_2 - Simulink','Location','best')

p = pole(G)

figure(2)
pzmap(G); 
grid on
ax = gca;
h = findobj(ax,'Type','line');          
set(h,'MarkerSize',10,'LineWidth',1.5);