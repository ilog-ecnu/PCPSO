classdef MOEMT < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <large/none>
% Multi-objective multifactorial evolutionary algorithm
% rmp --- 1 --- Random mating probability

%------------------------------- Reference --------------------------------
% A. Gupta, Y. Ong, L. Feng, and K. C. Tan, Multiobjective multifactorial
% optimization in evolutionary multitasking, IEEE Transactions on
% Cybernetics, 2017, 47(7): 1652-1665.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            rmp = Algorithm.ParameterSet(0.8);
            optimiser = 3;
            psi = 2;
            
            %% Initialize population
            Population    = Problem.Initialization();
            groups = 1;
            methodToSelectxPrimeSolutions = 2;
            xPrimeList = WOF_selectxPrimes(Population, 4, methodToSelectxPrimeSolutions); 
            transformedProblemPopulationSize = 10;
            w ={};
            fitness = [];
            Dummy = [];
            for i = 1:4
                gamma = i*5;
                xPrime = xPrimeList(i);
                res = [];
                [G,gamma]                   = WOF_createGroups(Problem,gamma,xPrime,groups);
                GlobalDummy         = WOFcreateGlobalDummy(gamma, xPrime, G, Problem, transformedProblemPopulationSize, psi,optimiser);
                Dummy = [Dummy,GlobalDummy];
                for j = 1:size(Population,2)
                    ans = Weight(Population(j).dec,gamma,GlobalDummy);
                    solution = WeightIndividual(ans,GlobalDummy,Problem);
                    res = [res,solution];
                end
                fitness = [fitness;calrank(res)];
                w{i} = res;
            end
            
            [~,I] = max(fitness);
            for i = 1:size(Population,2)
                Population(i).adds([w{I(i)}(i).dec,I(i)]);
            end

            %% Optimization
            count = 0;
            while Algorithm.NotTerminated(Population)
                Offspring = CreateOff(Problem,Population,rmp,Dummy,xPrimeList);
                Population = EnviSelect([Population,Offspring],Problem.N);
                SubPopulation = Divide(Population,4);
                
                count = count + 1;
                if(mod(count,5)==0)
                    tpop = [];
                    sum = 1;
                    for i = 1:4
                        if (size(SubPopulation{i},2))>Problem.N/4
                            tpop = [tpop,SubPopulation{i}(Problem.N/4+1:end)];
                            SubPopulation{i} = SubPopulation{i}(1:Problem.N/4);
                        end
                    end
                    for i = 1:4
                        gamma = i*5;
                        if (size(SubPopulation{i},2))<Problem.N/4
                            tmp = Problem.N/4-size(SubPopulation{i},2);
                            for j = sum:sum+tmp-1
                                ans = Weight(tpop(j).dec,gamma,Dummy(i));
                                tpop(j).add = [ans,i];
                                SubPopulation{i} = [SubPopulation{i},tpop(j)];
                            end
                        end
                        sum = sum+Problem.N/4-size(SubPopulation{i},2);
                    end
                    xPrimeList = [];
                for i = 1:size(SubPopulation,2)
                    xPrime = WOF_selectxPrimes(SubPopulation{i}, 1, methodToSelectxPrimeSolutions);
                    xPrimeList = [xPrimeList, xPrime];
                end
                Dummy = [];
                for i = 1:4
                    gamma = i*5;
                    xPrime = xPrimeList(i);
                    [G,gamma]                   = WOF_createGroups(Problem,gamma,xPrime,groups);
                    GlobalDummy         = WOFcreateGlobalDummy(gamma, xPrime, G, Problem, transformedProblemPopulationSize, psi,optimiser);
                    GlobalDummy.xPrime = xPrime;
                    Dummy = [Dummy,GlobalDummy];
                end
                end
                
                Population = [SubPopulation{:}];
            end
        end
    end
end

function GlobalDummy = WOFcreateGlobalDummy(gamma, xPrime, G, Global, populationSize, psi,optimiser)
    % Creates a dummy object. Needed to simulate the global class. Its
    % necessary to include this method into the Platemo
    % framework. 
    GlobalDummy = {};
    GlobalDummy.lower       = zeros(1,gamma);
    GlobalDummy.upper       = ones(1,gamma).*2.0;
    if or(optimiser == 2,optimiser == 4)
        [uniW,GlobalDummy.N]    = UniformPoint(populationSize,Global.M);
        GlobalDummy.uniW        = uniW;
    else
        GlobalDummy.N           = populationSize;
    end
    GlobalDummy.xPrime      = xPrime;
    GlobalDummy.G           = G;
    GlobalDummy.psi         = psi;
    GlobalDummy.xPrimelower = Global.lower;
    GlobalDummy.xPrimeupper = Global.upper;
    GlobalDummy.isDummy     = true;
    GlobalDummy.Global      = Global;
end

function WeightPopulation = WOFcreateInitialWeightPopulation(N, gamma, GlobalDummy,Problem,Population)
    %creates an initial population for the transformed problem
    decs = rand(N,gamma).*2.0;
    WeightPopulation = [];
    for i = 1:N
        solution = decs(i,:);
     %   solution = Weight(decs(i,:),GlobalDummy,Problem,Population(i));
        WeightPopulation = [WeightPopulation; solution];
    end
end