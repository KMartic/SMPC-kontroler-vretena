function [nmax, nmin, numcc] = vjerojatnosna_ogranicenja(A, G, Qw, N, xmax, xmin, var0, epsilon)


nx = size(A,1);

SIGMA = zeros(N, nx);
var = var0;
for i = 1:N
    var = A*var*A' + G*Qw*G';

    SIGMA(i,:) = sqrt(max(diag(var), 0)).';
end

nmax_mat = zeros(nx, N);
nmin_mat = zeros(nx, N);
numcc = 0;

for i = 1:N
    for j = 1:nx

        if (isnan(xmax(j))) == 0
            nmax_mat(j,i) = -norminv(1 - epsilon(j), -xmax(j), SIGMA(i,j));
            numcc = numcc + 1;
        else
            nmax_mat(j,i) = NaN;
        end


        if (isnan(xmin(j))) == 0
            nmin_mat(j,i) =  norminv(1 - epsilon(j),  xmin(j), SIGMA(i,j));
            numcc = numcc + 1;
        else
            nmin_mat(j,i) = NaN;
        end
    end
end


auxnmax = zeros(N*nx, 1);
auxnmin = zeros(N*nx, 1);
for i = 1:N
    auxnmax((i-1)*nx+1:i*nx, 1) = nmax_mat(:,i);
    auxnmin((i-1)*nx+1:i*nx, 1) = nmin_mat(:,i);
end

nmax = auxnmax;
nmin = auxnmin;

end
