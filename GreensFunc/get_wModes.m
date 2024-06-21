function [W_,W_Int3_] = get_wModes(dn_,dp_,x0,J,Jinv,M4invRw,sizeM,nMax,x,num_x,num_x0)

% stacked vector of wn, wn_x,..., wn_x3
W_ = zeros(sizeM,num_x);
W_Int3_ = zeros(sizeM,num_x);

M4invRw = reshape(M4invRw, [nMax,1,num_x0]);

for j = 1:num_x % because of memory limits, do this outer loop over x and integrate in x0 immediately
    X = x(j);

    W_(:,j) = get_W_j(dn_,dp_,X,x0,J,Jinv,M4invRw,sizeM,nMax);
    W_Int3_(:,j) = get_W_Int3_j(dn_,dp_,X,x0,J,Jinv,M4invRw,sizeM,nMax);

    if mod(j,20) == 0
        disp("Done up to x = " + X); 
    end
end