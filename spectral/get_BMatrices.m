function [Bsc,Bcs] = get_BMatrices(n)

Bsc = zeros(n);
Bcs = zeros(n);

for i = 1:n
    for j = 1:n
        if mod(i-j,2) == 1
            Bsc(i,j) = 2/pi * (1/(i+j) + 1/(i-j));
            Bcs(i,j) = 2/pi * (1/(i+j) - 1/(i-j));
        end
    end
end