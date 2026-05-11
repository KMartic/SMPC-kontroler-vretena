function [K, PN] = LQR_regulator(A, B, Qx, Ru)


PNprev = ones(size(Qx));
PN = zeros(size(Qx));

while ~isequal(PN, PNprev)
    PNprev = PN;
    for i = 1:100
        PN = Qx + A'*PN*A - (A'*PN*B) * ((Ru + B'*PN*B) \ (B'*PN*A));
    end
end

K = -(Ru + B'*PN*B) \ (B'*PN*A);
end
