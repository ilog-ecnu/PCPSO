function [FNmean,FNbest] = SubPopRank(Populations)


    Population = [Populations{:}];
    K          = length(Populations);
    Flag       = [];
    j          = 1;
    for i = 1 : K
        Flag(j:j+length(Populations{i})-1) = i;
        j = j + length(Populations{i});
    end
    FrontNoAll = NDSort(Population.objs,inf);
    FNmean     = zeros(1,K);
    FNbest     = zeros(1,K);
    for i = 1 : K
        FNmean(i) = mean(FrontNoAll(Flag==i));
        FNbest(i) = min(FrontNoAll(Flag==i));
    end
end