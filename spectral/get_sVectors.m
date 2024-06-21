function [s0,s1,s2,s3,s4,s5,s6] = get_sVectors(n)

s0 = zeros(n,1);
s1 = s0;
s2 = s0;
s3 = s0;
s4 = s0;
s5 = s0;
s6 = s0;

% s0
for i = 1:n
    if mod(i,2) == 1
        s0(i) = 4/(i*pi);
    end   
end

% s1
for i = 1:n
    s1(i) = (-1)^(i+1)*2/(i*pi);
end

% s2
for i = 1:n
    if mod(i,2) == 1
        s2(i) = 2/(i*pi) - 8/(i^3*pi^3);
    else
        s2(i) = -2/(i*pi);
    end   
end

% s3
for i = 1:n
    s3(i) = (-1)^(i+1)*(2/(i*pi) - 12/(i^3*pi^3));
end

% s4
for i = 1:n
    if mod(i,2) == 1
        s4(i) = 2/(i*pi) - 24/(i^3*pi^3) + 96/(i^5*pi^5);
    else
        s4(i) = -2/(i*pi) + 24/(i^3*pi^3);
    end
end

% s5
for i = 1:n
    s5(i) = (-1)^(i+1)*(2/(i*pi) - 40/(i^3*pi^3) + 240/(i^5*pi^5));
end

% s6
for i = 1:n
    if mod(i,2) == 1
        s6(i) = 2/(i*pi) - 60/(i^3*pi^3) + 720/(i^5*pi^5) - 2880/(i^7*pi^7);
    else
        s6(i) = -2/(i*pi) + 60/(i^3*pi^3) - 720/(i^5*pi^5);
    end
end
