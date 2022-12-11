function SubPopulation = Divide(Population,SubCount)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    skills = [];
    for i = 1:size(Population,2)
        skills = [skills,Population(i).add(end)];
    end

    for i = 1 : SubCount
        SubPopulation{i} = Population(skills==i);
    end
end