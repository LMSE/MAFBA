function solution = get_solution_ec_model(model)
    
    solution = optimizeCbModel(model);
    solution.phiPmax = model.ub(end-1);
    solution.phiMmax = model.ub(end);
    
    if solution.stat == 1
        biomassRxn = find(model.c);
        solution.f = solution.x(biomassRxn);
        fprintf('The maximum growth rate is %1.4f \n', solution.f);
        solution.x = solution.x(1:length(model.rxns));
        Mrxns = [model.protGroup(1).rxns; model.protGroup(4).rxns; model.protGroup(5).rxns];
       
        %Adding the sum of each sector in the overall proteome
        costP_m = model.w_c(Mrxns).* abs(solution.x(Mrxns));
        solution.phiP_m = sum(costP_m); % phiP_m
        costP_c = model.w_c(model.protGroup(1).rxns).* abs(solution.x(model.protGroup(1).rxns));
        solution.phiCm = sum(costP_c); %phiCm
        costP_r = model.w_c(model.protGroup(3).rxns).* abs(solution.x(model.protGroup(3).rxns));
        solution.phiR = sum(costP_r); %phiR
        costP_ec = model.w_c(model.protGroup(2).rxns).* abs(solution.x(model.protGroup(2).rxns));
        solution.phiEc = sum(costP_ec); %phiEc
        costP_er = model.w_c(model.protGroup(4).rxns).* abs(solution.x(model.protGroup(4).rxns));
        solution.phiEr = sum(costP_er); %phiEr
        costP_em = model.w_c(model.protGroup(5).rxns).* abs(solution.x(model.protGroup(5).rxns));
        solution.phiEm = sum(costP_em); %phiEm
        subset = solution.x(1:length(solution.x)-2);     % Access elements from startIndex to endIndex
        costP_p = model.w_c(1:length(solution.x)-2).* abs(subset);
        solution.phiP = sum(costP_p); %phiP used
        solution.phiCc = solution.phiPmax - solution.phiP ; %phiCc

        %Adding the sum of each sector in the membrane proteome
        costM_m = model.w_m(Mrxns).* abs(solution.x(Mrxns));
        solution.phiM = sum(costM_m); % phiM 
        costM_c = model.w_m(model.protGroup(1).rxns).* abs(solution.x(model.protGroup(1).rxns));
        solution.phiCm_m = sum(costM_c); %phiCm
        costM_em = model.w_m(model.protGroup(5).rxns).* abs(solution.x(model.protGroup(5).rxns));
        solution.phiEm_m = sum(costM_em); %phiEm
        costM_er = model.w_m(model.protGroup(4).rxns).* abs(solution.x(model.protGroup(4).rxns));
        solution.phiEr_m = sum(costM_er); %phiEr

    else
        solution.x = zeros( length(model.rxns), 1);
        solution.f = 0;
        
        solution.phiP_m = 0; solution.phiCm = 0; solution.phiCc = 0;solution.phiR = 0;...
        solution.phiEc = 0; solution.phiEr = 0; solution.phiEm = 0; solution.phiP = 0;
        solution.phiM =0; solution.phiCm_m =0; solution.phiEm_m =0; solution.phiEr_m =0;
        
    end
    
end