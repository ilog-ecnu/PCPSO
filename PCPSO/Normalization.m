function [PopObj,z,znad] = Normalization(PopObj,z,znad)


    [N,M] = size(PopObj);
    %disp(size(PopObj,1));
    %% Update the ideal point
    z = min(z,min(PopObj,[],1));
    
    %% Update the nadir point
    % Identify the extreme points
    W = zeros(M) + 1e-6;
    W(logical(eye(M))) = 1;
    ASF = zeros(N,M);   
    for i = 1 : M
        ASF(:,i) = max(abs((PopObj-repmat(z,N,1))./(repmat(znad-z,N,1)))./repmat(W(i,:),N,1),[],2);
    end
    [~,extreme] = min(ASF,[],1);
    % Calculate the intercepts
    Hyperplane = (PopObj(extreme,:)-repmat(z,M,1))\ones(M,1);
    a = (1./Hyperplane)' + z;
    if any(isnan(a)) || any(a<=z)
        a = max(PopObj,[],1);
    end
    znad = a;
    
    %% Normalize the population
    PopObj = (PopObj-repmat(z,N,1))./(repmat(znad-z,N,1));
end