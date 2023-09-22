classdef PCPSO < ALGORITHM
% <multi/many> <real> <large/none> <constrained/none>

%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            [V,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Population2    = Problem.Initialization();
            K = 2;
            M = cell(1,K);
            index = 0;
            Populations = cell(1,K);
            
            [z,znad]      = deal(min(Population.objs),max(Population.objs));
            [z2,znad2]      = deal(min([Population.objs;Population2.objs]),max([Population.objs;Population2.objs]));
            tag = 0;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                
                if Problem.FE<Problem.maxFE/3
                    Fitness1 = calFitness(Population.objs);
                    
                    if length(Population) >= 2
                        Rank = randperm(length(Population),floor(length(Population)/2)*2);
                    else
                        Rank = [1,1];
                    end 
                    Loser  = Rank(1:end/2);
                    Winner = Rank(end/2+1:end);
                    Change = Fitness1(Loser) >= Fitness1(Winner);
                    Temp   = Winner(Change);
                    Winner(Change) = Loser(Change);
                    Loser(Change)  = Temp;
                    Offspring      = Operator(Population(Loser),Population(Winner));
                    [Population2,z2,znad2]  = EnvironmentalSelection1([Population2,Offspring],V,Problem.N,z2,znad2);
                    
                    temp = Offspring;
                    
                    Pop = repmat(Population2,1,ceil(size(temp.objs,1)/size(Population2.objs,1)));
                    Pop = Pop(:,1:size(temp.objs,1));
                    
                    Offspring1  = OperatorGA([Pop,temp]);
                    [Population2,z2,znad2]  = EnvironmentalSelection2([Population2,Offspring,Offspring1],V,Problem.N,z2,znad2);
                    
                    [Population,z2,znad2]  = EnvironmentalSelection2([Population,Offspring,Offspring1],V,Problem.N,z2,znad2);
                    
                else
                    tag = tag + 1;
                    if index == 0
                        index       = randperm(floor(size(Population.objs,1)/K)*K); 
                        temp        = reshape(index,K,floor(size(Population.objs,1)/K));
                        for i = 1:K
                            Populations{i} = Population(temp(i,:));
                        end
                        index = 1;
                    end
                
                    [~,rank] = sort(SubPopRank(Populations));
                    for i = 1:K
                        Population = Populations{i};
                        Fitness = calFitness(Population.objs);
                        
                        if length(Populations{i}) >= 2
                            Rank = randperm(length(Populations{i}),floor(length(Populations{i})/2)*2);
                        else
                            Rank = [1,1];
                        end
                        
                        Loser  = Rank(1:end/2);
                        Winner = Rank(end/2+1:end);
                        Change = Fitness(Loser) >= Fitness(Winner);
                        Temp   = Winner(Change);
                        Winner(Change) = Loser(Change);
                        Loser(Change)  = Temp;
                        Offspring      = Operator(Populations{i}(Loser),Populations{i}(Winner));
                        [Pop,FrontNo,~] = ES(Populations{i},floor(Problem.N/K));
                        Pop = Pop(FrontNo==1);
                        if size(Pop.objs,1) > floor(Problem.N/K)/2
                            Pop = Populations{i}(Winner);
                        end
                        M{rank(i)} = mean(Pop.objs,1);

                        if i>1
                            Mean = M{rank(i)};
                            if i>2
                                for j = 2:i-1
                                    Mean = Mean+M{rank(j)};
                                end
                            end
                            Mean = Mean/(i-1);
                            gv = Mean;
                            Mean = repmat(gv,size(Offspring.objs,1)+size(Populations{i}.objs,1),1);
                            tmp = [Populations{i},Offspring];
                            dis = sqrt(sum((tmp.objs-Mean).^2,2));
                            mmax = max(dis);
                            dis = dis/mmax;  
                            mmax = max(dis);
                            dis = dis*(-1)+mmax;
                            [Populations{i},z2,znad2]  = EnvironmentalSelection3([Populations{i},Offspring],V,floor(Problem.N/K),z2,znad2,dis);


                        else
                            [Populations{i},z2,znad2]  = EnvironmentalSelection2([Populations{i},Offspring],V,floor(Problem.N/K),z2,znad2);

                        end
                    end
                    %merge and divide
                    if tag==50 &&  K > 1
                        tag = tag-50;
                        if K > 1
                            [ss,index] = SubPopSimility(Populations);
                        end
                        if ss < 0.1
                            K = K-1;
                            i = index(1);
                            j = index(2);
                            [Populations{i},z2,znad2]  = EnvironmentalSelection4([Populations{i},Populations{j}],V,floor(Problem.N/K),z2,znad2);
                            Populations(j) = [];
                            M(j) = [];
                        elseif ss > 0.3
                            for i =1:K
                                [Populations{i},z2,znad2]  = EnvironmentalSelection4(Populations{i},V,floor(Problem.N/(K+1)),z2,znad2);
                            end
                            K = K + 1;
                            Populations{K}   = Problem.Initialization(floor(Problem.N/(K+1)));
                        end
                    end
                    Population = [Populations{:}];
                end
            end

        end
    end
end

function Fitness = calFitness(PopObj)
% Calculate the fitness by shift-based density

    N      = size(PopObj,1);
    fmax   = max(PopObj,[],1);
    fmin   = min(PopObj,[],1);
    PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    Dis    = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    Fitness = min(Dis,[],2);
end

function [Population,FrontNo,CrowdDis] = ES(Population,N)
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    if size(Population.objs,1)>=N
        Next(Last(Rank(1:N-sum(Next)))) = true;
    else
        Next(Last(Rank(1:size(Population.objs,1)-sum(Next)))) = true;
    end
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function [ss,index] = SubPopSimility(Pop)
    K  = length(Pop);
    ss = 1000000;
    for i = 1 : K-1
        for j = i+1 : K
          %  score_KL = sum(sum(Pop{i}.objs.*log(eps+Pop{i}.objs./(Pop{j}.objs+eps))));
            vec1 = Pop{i}.objs;
            vec2 = Pop{j}.objs;
            disp([size(vec1,1),size(vec2,1)]);
            if any(vec1(:))
                vec1 = vec1/sum(vec1(:));
            end
 
            if any(vec2)
                vec2 = vec2/sum(vec2(:));
            end
            score_JS = (sum(sum(vec1.* log(eps + vec1./((vec1+vec2)./2+eps))))+sum(sum(vec2.* log(eps + vec2./((vec1+vec2)./2+eps)))))./2;
            disp(score_JS);
            if score_JS < ss
                ss = score_JS;
                index = [i,j];
            end
        end
    end
end
            