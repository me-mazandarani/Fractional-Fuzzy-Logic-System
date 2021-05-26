%% FFIS201210
% FFIS.m is a function for the implementation of fractional fuzzy inference
% System
%% The reference
% Mehran Mazandarani, Li Xiu, Fractional Fuzzy Inference System:
% The New Generation of Fuzzy Inference Systems,
% IEEE Access, Vol. 8, pp. 126066-126082, 2020. 
% DOI: 10.1109/ACCESS.2020.3008064 
%  Link: https://ieeexplore.ieee.org/document/9136681

% If there is any question contact the authors of the above paper.  
% Date of issue 2020-12-10
% Before using FFIS.m read carefully the help file, help_FFIS.pdf.
% This program is the first version of the function FFIS.m (FFIS_201210)
% The algorith of this function has not been written in an optimal manner.
% Thus, naturally, it may take more time than that one may expect for 
% getting the output.




function output=FFIS(fis,Inputs,Fids)


Nrules=max(size(fis.Rules));
No_Inputs=size(fis.Inputs,2);

% Initialization for the aggregation process
h_output=0.01;
output_range=fis.Output.Range;
% Discretizing the Range of output
y=linspace(output_range(1),output_range(2),diff(output_range)/h_output);

out=zeros(1,numel(y));

% % % % % % % % % %


for i=1:Nrules
      
    TD=zeros(1,No_Inputs);
    Input_MF_No=fis.Rules(i).Antecedent;
    
    for j=1:No_Inputs
        Input_MF_Param=fis.Inputs(1,j).MembershipFunctions(1,Input_MF_No(j)).Parameters;
        Input_MF_Type=fis.Inputs(1,j).MembershipFunctions(1,Input_MF_No(j)).Type;
        Input_MF= fismf(Input_MF_Type,Input_MF_Param);
        TD(j) = evalmf(Input_MF,Inputs(j));
    end
    switch fis.AndMethod
        case 'min'
            Mu_star=min(TD);
        otherwise
            Mu_star=prod(TD);
    end
    
    Output_MF_No=fis.Rules(i).Consequent;
    
    % % %     The max operator for the Aggregation Method
  
    out=max(out,FTR(fis,Output_MF_No,Mu_star,y,Fids));
    
end
defuz=fis.DefuzzificationMethod;
% % % Defuzzification of the output
output = defuzz(y,out,defuz);

%% Fractional Translation Rule

function out=FTR(fis,Output_MF_No,Mu_star,y,Fids)


Output_MF_Param=fis.Outputs(1,1).MembershipFunctions(1,Output_MF_No).Parameters;
out=FMF(y,Output_MF_Param,Mu_star,Fids(2*Output_MF_No-1:2*Output_MF_No));

%%

%% Fractional Membership Function
% This function works only for the triangular and trapezoidal membership 
% functions.

function out=FMF(y,Output_MF_Param,Mu_star,Fidx)

Fidx_value=str2func(['@(x)',num2str(Fidx{1})]);
switch Fidx{2}
    case {'b',1}
        betta=Fidx_value(Mu_star);
        cs=HMFunction_betta(Output_MF_Param,Mu_star,0);
        bs=HMFunction_betta(Output_MF_Param,Mu_star,betta);
        as=HMFunction_betta(Output_MF_Param,0,betta);
        as=min(bs,as);
        
        if   abs(cs-Output_MF_Param(end))<=10^(-8)
            cs=Output_MF_Param(end);
        end
        
        out=Mu_star*trapmf(y,[as,bs,cs,Output_MF_Param(end)]);
        
    case {'a',0}
        alpha=Fidx_value(Mu_star);
        bs=HMFunction_alpha(Output_MF_Param,Mu_star,0);
        cs=HMFunction_alpha(Output_MF_Param,Mu_star,alpha);
        ds=HMFunction_alpha(Output_MF_Param,0,alpha);
        ds=max(cs,ds);
        
        if   abs(bs-Output_MF_Param(1))<=10^(-8)
            bs=Output_MF_Param(1);
        end
        
        out=Mu_star*trapmf(y,[Output_MF_Param(1),bs,cs,ds]);
        
end

%% Fractional Horizontal MF (beta form)
function outgr=HMFunction_betta(Output_MF_Param,Mu_star,betta)

s=numel(Output_MF_Param);
switch s
    case 3
        a=Output_MF_Param(1);
        b=Output_MF_Param(2);
        c=Output_MF_Param(3);
        
        outgr=(a+(b-a)*Mu_star)+(1-Mu_star)*(c-a)*(1-betta);
        
    case 4
        a=Output_MF_Param(1);
        b=Output_MF_Param(2);
        c=Output_MF_Param(3);
        d=Output_MF_Param(4);
        
        outgr=(a+(b-a)*Mu_star)+((d-a)-Mu_star*(d-a+b-c))*(1-betta);
        
end



end

%% Fractional Horizontal MF (alpha form)
function outgr=HMFunction_alpha(Output_MF_Param,Mu_star,alpha)

s=numel(Output_MF_Param);
switch s
    case 3
        a=Output_MF_Param(1);
        b=Output_MF_Param(2);
        c=Output_MF_Param(3);
        
        outgr=(a+(b-a)*Mu_star)+(1-Mu_star)*(c-a)*alpha;
        
    case 4
        a=Output_MF_Param(1);
        b=Output_MF_Param(2);
        c=Output_MF_Param(3);
        d=Output_MF_Param(4);
        
        outgr=(a+(b-a)*Mu_star)+((d-a)-Mu_star*(d-a+b-c))*alpha;
        
end



end


end






end







end
