function [c1] = get_cVectors(n)

c1 = zeros(n,1);

% c1
for i = 1:n
    if mod(i,2) == 1
        c1(i) = -4/(i^2*pi^2);
    end
end
