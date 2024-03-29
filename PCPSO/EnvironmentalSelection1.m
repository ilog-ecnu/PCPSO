function [Population,z,znad] = EnvironmentalSelection1(Population,W,N,z,znad)


    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    St = find(FrontNo<=MaxFNo);

    %% Normalization
    [PopObj,z,znad] = Normalization(Population(St).objs,z,znad);
    
    %% theta-non-dominated sorting
    tFrontNo = tNDSort(PopObj,W);
    
    %% Selection
    if size(PopObj,1) < N
        Next      = St(tFrontNo<=MaxFNo);
        Population = Population(Next);
    else
        Next      = St(tFrontNo<=2);
        %disp(size(Next,2));
        if size(Next,2) < N/3
        % Population for next generation
            Next      = St(tFrontNo<=2);
            Population = Population(Next);
        else
            Population = Population(Next);
        end
    end
end

function tFrontNo = tNDSort(PopObj,W)
% Do theta-non-dominated sorting

    N  = size(PopObj,1);
    NW = size(W,1);
    
    %% Calculate the d1 and d2 values for each solution to each weight
    normP  = sqrt(sum(PopObj.^2,2));
    Cosine = 1 - pdist2(PopObj,W,'cosine');
    d1     = repmat(normP,1,size(W,1)).*Cosine;
    d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);
    %% Clustering
    [~,class] = min(d2,[],2);  
    %% Sort
    theta = zeros(1,NW) + 5;
    theta(sum(W>1e-4,2)==1) = 1e6;
    tFrontNo = zeros(1,N);
    for i = 1 : NW
        C = find(class==i);
        [~,rank] = sort(d1(C,i)+theta(i)*d2(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
        
    end
end