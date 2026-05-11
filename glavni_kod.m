clc
clear all

global N w2_ref YR Ycmax Ycmin Umax Umin q qN r n m K Acl Sox Sou AR BRa BRb

N = 15;  % pr. horizont
Td = 0.5;   % vrijeme diskretizacije
Tspan = 120;
Nmc = 1;  % broj mc simulacija
mu_m2 = 0;    % ocekivanje suma
sigma_m2 = 10;      % devijacija suma
w2_ref = pi/30*100; % zeljena brzina vrtnje
epsilon = 0.05;   % dizajnerski parametar
brojac = 1;
broj_prekrsenih = 0;
FS = 18;

Sx = [];
YR = [];
Ycmax = [];
Ycmin = [];
Umax = [];
Umin = [];

% ogranicenja

ymax = 0.2; 
ymin = -0.2;
umax = 200; 
umin = -200; 

% tezine
q = 1e3;
qN = 1e3;
r = 1;

J1 = 10;
J2 = 640; 
c = 400; 
d = 20; 
ip = 12.5; 


A = [-d/(J1*ip^2), -c/(J1*ip), d/(J1*ip);
     1/ip, 0, -1;
     d/(J2*ip), c/J2, -d/J2 ];
B = [1/J1; 0; 0];
C = [0 0 1];
D = [0 0];
G = [0; 0; -1/J2];


%diskretizacjia
sys_c = ss(A, [B G], C, D);
sys_d = c2d(sys_c, Td, 'zoh');
Ad = sys_d.A;
Bd = sys_d.B(:,1);
Gd = sys_d.B(:,2);
Cd = sys_d.C;
Dd = sys_d.D;
[n, ~] = size(Ad);
[~, m] = size(Bd);


%% lqr

Qx = blkdiag(0, 0, q);
Ru = r*eye(m);

% K i PN (Riccati)
try
    [Kdlqr, PN_lqr] = dlqr(Ad, Bd, Qx, Ru);
    K = -Kdlqr;                                     
catch
    [K, PN_lqr] = LQR_regulator(Ad, Bd, Qx, Ru); 
end

Acl = Ad + Bd*K;



% dopr pocetnog stanja Sx i bu ulaza Su
Su = zeros(n*N, N*m);
for i = 1:N
    Sx = [Sx; Ad^i];
    for j = 1:i
        Su((i-1)*n+1:i*n, (j-1)*m+1:j*m) = Ad^(i-j) * Bd;
    end
end

% ciljane matrice
for i = 1:N
    YR = [YR; w2_ref];
    Umax = [Umax; umax];
    Umin = [Umin; umin];
end

% matrice izl varijable
C_izl = repmat({Cd}, 1, N);
C_izl_dia = blkdiag(C_izl{:});
Sox = C_izl_dia * Sx;
Sou = C_izl_dia * Su;

% matrica ogranicenja
C_ogr = repmat({[0 1 0]}, 1, N);
C_ogr_dia = blkdiag(C_ogr{:});
Scx = C_ogr_dia * Sx;
Scu = C_ogr_dia * Su;

%% vjerojatnosno ogranicenje

var0 = zeros(n);  % pocetna kovarijanca pogreske
Qw  = sigma_m2^2; % varijanca poremecaja

% vektori ogranicenja stanja
xmax    = [NaN; ymax(1); NaN];
xmin    = [NaN; ymin(1); NaN];
epsilon_m = [NaN; epsilon; NaN];

[nmax, nmin, numcc] = vjerojatnosna_ogranicenja(Acl, Gd, Qw, N, xmax, xmin, var0, epsilon_m);

alfa_redovi = (2:n:n*N).';  

% suzavanje ogranicenja alfe
Ycmax = nmax(alfa_redovi); 
Ycmin = nmin(alfa_redovi);  

% Ycmax = ymax * ones(N,1);
% Ycmin = ymin * ones(N,1); 

%% graf suzavanja vjerojatnosnih ogranicenja za delta alfa kroz horizont
% (prikazuje originalne granice i suzene (tightened) granice koje MPC koristi)
t_pred = (1:N)' * Td; % vrijeme predikcije: 1..N koraka u buducnost

figure(3); clf; hold on; grid on;
plot(t_pred, Ycmax, 'LineWidth', 2);
plot(t_pred, Ycmin, 'LineWidth', 2);
yline(ymax,'k--','LineWidth',1);
yline(ymin,'k--','LineWidth',1);
xlabel('N[-]','FontSize',FS);
ylabel('Δα [rad]','FontSize',FS);
title('Suzavanje vjerojatnosnih ogranicenja za Δα','FontSize',FS);
legend('gornja suzena granica', ...
       'donja suzena granica', ...
       'Location','best');
hold off;

% formiranje lin ogranicenja
Iu = eye(m*N);

AR_u  = [Iu; -Iu];
BRa_u = zeros(2*m*N, n);
BRb_u = [Umax; -Umin];

AR_x  = [Scu; -Scu];
BRa_x = [-Scx; Scx];
BRb_x = [Ycmax; -Ycmin];

AR  = [AR_u; AR_x];
BRa = [BRa_u; BRa_x];
BRb = [BRb_u; BRb_x];

%ukljucenost ogranicejna

active = ~isnan(BRb);
AR  = AR(active, :);
BRa = BRa(active, :);
BRb = BRb(active, :);

% AR  = zeros(0, m*N);
% BRa = zeros(0, n);
% BRb = zeros(0, 1);


%% monte carlo

for mc = 1:Nmc


    %% sum
    t_m2 = (0:Td:Tspan)'; 
    m2_values = mu_m2 + sigma_m2 * randn(size(t_m2));
    m2 = timeseries(m2_values, t_m2);

    res = sim('busilica.slx');

    % krsenje ogranicenja
    krsenje = d_a > ymax | d_a < ymin;

    brojac_krsenja = nnz(krsenje);

    if brojac_krsenja > 0
        broj_prekrsenih = broj_prekrsenih + 1;
    end

    posto = (broj_prekrsenih/Nmc)*100;
   


%% grafovi

legenda = sprintf('ime');   

fig = figure(1);

subplot(3,1,1);
hold on; grid on;
plot(t, w1, 'LineWidth', 1.5, 'DisplayName', legenda);
title('Brzina na strani motora','FontSize', FS)
xlabel('t [s]','FontSize', FS);
ylabel('ω_1 [rad/s]','FontSize', FS);
legend('show','Location','best');
lgd.FontSize = 9;

subplot(3,1,2);
hold on; grid on;
    if isempty(findobj(gca,'Tag','wgranice'))
        yline(ymax,'k--','LineWidth',1,'HandleVisibility','off');
        yline(ymin,'k--','LineWidth',1,'HandleVisibility','off');
end
plot(t, d_a, 'LineWidth', 1.5, 'DisplayName', legenda);
title('Kut torzije vratila','FontSize', FS)
xlabel('t [s]','FontSize', FS);
ylabel('Δα [rad]','FontSize', FS);
legend('show','Location','best');


subplot(3,1,3);
hold on; grid on;
if isempty(findobj(gca,'Tag','w2_ref'))
         yline(w2_ref,'k--','LineWidth',1,'HandleVisibility','off');
end
plot(t, w2, 'LineWidth', 1.5, 'DisplayName', legenda);
title('Brzina na strani alata - izlaz sustava','FontSize', FS)
xlabel('t [s]','FontSize', FS);
ylabel('ω_2 [rad/s]','FontSize', FS);
legend('show','Location','best');


figure(2);

subplot(2,1,1);
hold on; grid on;
    if isempty(findobj(gca,'Tag','ugranice'))
        yline(umax,'k--','LineWidth',1,'HandleVisibility','off');
        yline(umin,'k--','LineWidth',1,'HandleVisibility','off');
end
stairs(t, u, 'LineWidth', 1, 'DisplayName', legenda);
title('Ulazni moment - upravljačka varijabla','FontSize', FS)
xlabel('t [s]','FontSize', FS);
ylabel('m1 [Nm] - u','FontSize', FS);
legend('show','Location','best');

subplot(2,1,2);
hold on; grid on;
    if isempty(findobj(gca,'Tag','w2_ref'))
        yline(w2_ref,'k--','LineWidth',1,'HandleVisibility','off');
end
plot(t, w2, 'LineWidth', 1.5, 'DisplayName', legenda);
title('Brzina na strani alata - izlaz sustava','FontSize', FS)
xlabel('t [s]','FontSize', FS);
ylabel('ω_2 [rad/s]','FontSize', FS);
legend('show','Location','best');

    fprintf('broj simulacije: %d\n ', brojac);
    brojac = brojac + 1;

end

fprintf('broj izvrsenih simulacija: %d ', Nmc);
fprintf('broj simulacija koje su prekršile granicu za deltaalfa: %d / %d\n', broj_prekrsenih, Nmc);
fprintf('postotak simulacija koje su prekršile granicu za deltaalfa: %.2f%%\n', posto);