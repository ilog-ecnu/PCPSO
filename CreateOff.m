function Offspring = CreateOff(Problem,Population,rmp,Dummy,xPrimeList)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parent selection
    Offspring = [];
    for i = 1 : floor(length(Population)/2)
        P1 = Population(i);
        P2 = Population(i+floor(length(Population)/2));
        I1 = P1.add(end);
        I2 = P2.add(end);
        if (I1 == I2)
            OffDec1   = cross_mu([P1.add(1:end-1);P2.add(1:end-1)],size(P1.obj,2));
            GlobalDummy = Dummy(I1);
            xPrimeVars = xPrimeList(I1);
            for i = 1:size(OffDec1,1)
                W = OffDec1(1,:);
                x1 = WOF_transformationFunctionMatrixForm(xPrimeVars.dec,W(GlobalDummy.G),GlobalDummy.xPrimeupper,GlobalDummy.xPrimelower, GlobalDummy.psi);
                x1 = Problem.Evaluation(x1);
                x1.adds([W,I1]);
                Offspring = [Offspring,x1];
            end
            
            
        elseif (rand<rmp)
            GlobalDummy1 = Dummy(I1);
            xPrimeVars = xPrimeList(I1);
            P1 = P1.adds(1:end-1);
            x1 = WOF_transformationFunctionMatrixForm(xPrimeVars.dec,P1(GlobalDummy1.G),GlobalDummy1.xPrimeupper,GlobalDummy1.xPrimelower, GlobalDummy1.psi);
            GlobalDummy2 = Dummy(I2);
            xPrimeVars = xPrimeList(I2);
            P2 = P2.adds(1:end-1);
            x2 = WOF_transformationFunctionMatrixForm(xPrimeVars.dec,P2(GlobalDummy2.G),GlobalDummy2.xPrimeupper,GlobalDummy2.xPrimelower, GlobalDummy2.psi);
            OffDec2   = OperatorGA(Problem,[x1;x2]);
            Off = Problem.Evaluation(OffDec2);
            ans1 = Weight(x1,size(GlobalDummy1.lower,2),GlobalDummy1);
            ans1 = [ans1,I1];
            ans2 = Weight(x1,size(GlobalDummy2.lower,2),GlobalDummy2);
            ans2 = [ans2,I2];
            Off(1).adds(ans1);
            Off(2).adds(ans2);
            Offspring = [Offspring,Off];
        else
            OffDec1        = mutation(Problem,P1.dec);
            OffDec2        = mutation(Problem,P2.dec);
            Off = Problem.Evaluation([OffDec1;OffDec2]);
            Off(1).adds(P1.adds);
            Off(2).adds(P2.adds);
            Offspring = [Offspring,Off];
        end
    end


end
function Offspring = mutation(Problem,Offspring)
            N=1;
            [proC,disC,proM,disM] = deal(1,20,1,20);
            D=100;
            Lower = repmat(Problem.lower,2*N,1);
            Upper = repmat(Problem.upper,2*N,1);
            Site  = rand(2*N,D) < proM/D;
            mu    = rand(2*N,D);
            temp  = Site & mu<=0.5;
            Offspring       = min(max(Offspring,Lower),Upper);
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
            temp = Site & mu>0.5; 
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
            Offspring = Offspring(1,:);
end



function Offspring = cross_mu(Parent,M,Parameter)

    if nargin > 2
        [proC,disC,proM,disM] = deal(Parameter{:});
    else
        [proC,disC,proM,disM] = deal(1,20,1,20);
    end
    if isa(Parent(1),'SOLUTION')
        calObj = true;
        Parent = Parent.decs;
    else
        calObj = false;
    end
    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
    lower = zeros(1,D);
    upper = zeros(1,D)+2;
            %% Genetic operators for real encoding
            % Simulated binary crossover
            beta = zeros(N,D);
            mu   = rand(N,D);
            beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
            beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
            beta = beta.*(-1).^randi([0,1],N,D);
            beta(rand(N,D)<0.5) = 1;
            beta(repmat(rand(N,1)>proC,1,D)) = 1;
            Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                         (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
            % Polynomial mutation
            Lower = repmat(lower,2*N,1);
            Upper = repmat(upper,2*N,1);
            Site  = rand(2*N,D) < proM/D;
            mu    = rand(2*N,D);
            temp  = Site & mu<=0.5;
            Offspring       = min(max(Offspring,Lower),Upper);
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
            temp = Site & mu>0.5; 
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    if calObj
        Offspring = SOLUTION(Offspring);
    end
end