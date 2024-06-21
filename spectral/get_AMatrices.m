function [A0,A1s,A2s,A3s,A1c] = get_AMatrices(n)

A0 = eye(n);
A1s = zeros(n);
A2s = A1s;
A3s = A1s;
A1c = A1s;

% A1s
for i = 1:n
    for j = 1:n
        if mod(i-j,2) == 1
            A1s(i,j) = 2/pi^2*(1/(i+j)^2 - 1/(i-j)^2);
        end
    end
end
A1s = A1s + 1/2*eye(n);

% A2s
for i = 1:n
    for j = 1:n
        if i == j
            A2s(i,j) = 1/3 - 2/(pi^2*(i+j)^2);
        else
            A2s(i,j) = (-1)^(i-j) * 2/pi^2*(1/(i-j)^2 - 1/(i+j)^2);
        end
    end
end

% A3s
for i = 1:n
    for j = 1:n
        if i == j
            A3s(i,j) = 1/4 - 3/(pi^2*(i+j)^2);
        elseif mod(i-j,2) == 1
            A3s(i,j) = 3/pi^2*(1/(i+j)^2 - 1/(i-j)^2) ...
                    + 12/pi^4*(1/(i-j)^4 - 1/(i+j)^4);
        else
            A3s(i,j) = 3/pi^2*(1/(i-j)^2 - 1/(i+j)^2);
        end
    end
end

% -------------------------------------------------------------------------
% A1c
for i = 1:n
    for j = 1:n
        if mod(i-j,2) == 1
            A1c(i,j) = 2/pi^2*(-1/(i+j)^2 - 1/(i-j)^2);
        end
    end
end
A1c = A1c + 1/2*eye(n);
