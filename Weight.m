function ans = Weight(Dec,gamma,GlobalDummy)
            
            if nargin > 0
                xPrimeVars = GlobalDummy.xPrime.dec;
                xPrimeSize = size(xPrimeVars,2);
                pWert = 0.2;
  %              Dec  = max(min(Dec,GlobalDummy.upper),GlobalDummy.lower);
                u = (Dec-xPrimeVars)./(pWert.*(GlobalDummy.xPrimeupper-GlobalDummy.xPrimelower))+1;
                ans = [];
                for i =1:gamma
                    t = median(u(GlobalDummy.G==i));
                    ans = [ans,t];
                end
                ans(ans<0) = 0;
                ans(ans>2) = 2;
%                 if gamma == 15
%                     disp(ans);
%                     
%                 end
          %      obj = WOF_WeightIndividual;
                
                % Set the infeasible decision variables to boundary values


         %       Popdec = Pop.dec;
%                 x = WOF_transformationFunctionMatrixForm(xPrimeVars,Popdec(GlobalDummy.G),GlobalDummy.xPrimeupper,GlobalDummy.xPrimelower, GlobalDummy.psi);
%                 x = x(1:end-1);
%                 ans = [];
% 
%                 ans = [ans,Popdec(end)];

            end
            
 end
