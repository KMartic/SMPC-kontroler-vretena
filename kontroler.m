function u = kontroler(w1, d_a, w2)

import casadi.*

global N w2_ref YR Ycmax Ycmin Umax Umin q qN r n m K Acl Sox Sou AR BRa BRb

% pocetno stanje
x0 = [w1; d_a; w2];

% Odlucivacke varijable i parametar (stanje)
U_sym = SX.sym('U', m*N, 1);
p_sym = SX.sym('p', n, 1);

% predikcija omega2
Yf = Sox*p_sym + Sou*U_sym;
y_tren = [p_sym(3); Yf];  % omega 2

% ciljna funkcija

z_kum1 = y_tren(N+1) - w2_ref;
J_kum1 = qN * (z_kum1^2);

z_kum2 = y_tren(1:N) - w2_ref; 
J_kum2 = q * sum1(z_kum2.^2);

J_u = r * sum1(U_sym.^2);

J = J_kum2 + J_kum1 + J_u;

% lin ogranicenja
BR = BRa*p_sym + BRb;   
g_sym = AR*U_sym - BR;     

% NLP 
nlp = struct('x', U_sym, 'p', p_sym, 'f', J, 'g', g_sym);

opts = struct;
opts.ipopt.print_level = 0;
opts.ipopt.max_iter = 50;
opts.ipopt.tol = 1e-6;
opts.print_time = 0;

solver = nlpsol('solver', 'ipopt', nlp, opts);

% ogranicenja
lbg = -inf(size(g_sym));
ubg = zeros(size(g_sym));

lbx = -inf(m*N,1);
ubx = inf(m*N,1);

% pocetna pretpostavnka i izracun u
U0 = zeros(m*N, 1);

sol = solver('x0', U0, 'p', x0, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
Uopt = full(sol.x);

% samo prvi korak uzmi
u = Uopt(1:m);

end
