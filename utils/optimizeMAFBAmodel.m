function [solution, model] = optimizeMAFBAmodel(model, varargin)
    % Set default values
    options = varargin;
    [phiM_max, gamma, Sol] = parseOptionalArguments(options);

    % Technical setup and checks
    [model, nRxnsExt, nMix, nRxns, Mrxns, phi_max] = setupModelAndChecks(model, phiM_max, gamma);

    % Solve the optimization problem
    solution = getSolution(model, Sol, nMix, nRxns, phi_max, Mrxns);
end

function [phiM_max, gamma, Sol] = parseOptionalArguments(options)
    phiM_max = 0.20;
    gamma = 1;
    Sol = 1;

    optargin = size(options, 2);
    if mod(optargin, 2) ~= 0
        error('Error, optargin is not an even number.');
    end
    if optargin > 0
        for i = 1:2:optargin  % Corrected loop definition
            switch lower(options{i})  % Convert to lowercase for case-insensitivity
                case 'phimmax'
                    phiM_max = options{i + 1};
                case 'gamma'
                    gamma = options{i + 1};
                case 'incl_sol'
                    Sol = options{i + 1};
                otherwise
                    error('Optional argument not recognized.');
            end
        end
    end
end


function [model, nRxnsExt, nMix, nRxns, Mrxns, phi_max,phiM_max] = setupModelAndChecks(model, phiM_max, gamma)
    checkModelFields(model);
    [nMets,nRxns,prevS]= getModelSmatrix(model);
    
    [posVec, negVec, mixVec, posVec_m, negVec_m, mixVec_m, Smixed, ...
      nRxnsExt, nMix,Mrxns] = categorizeReactions(model, nRxns, prevS);
  
    [model,phi_max,phiM_max] = splitAndBuildMatrix(model,  nMix, mixVec, nMets, nRxns, ...
      nRxnsExt, posVec, negVec, posVec_m, negVec_m, mixVec_m, gamma, Smixed,prevS);

end

function checkModelFields(model)
    if ~isfield(model, 'w_c') || ~isfield(model, 'w_m') || ~isfield(model, 'protGroup')
        error('Invalid MAFBA model: Required fields missing.');
    end
end

function [nMets,nRxns,prevS]= getModelSmatrix(model)
    [nMets,nRxns]=size(model.S); 
    if strcmp(model.rxns{nRxns},'Memb_Const')
%         fprintf('model already with both protein constraints\n')
        nMets = nMets-2;
        nRxns = nRxns-2;
        prevS = model.S(1:nMets,1:nRxns); 
    elseif strcmp(model.rxns{nRxns},'Prot_Const')
%         fprintf('model already with one protein constraint\n')
        nMets = nMets-1;
        nRxns = nRxns-1;
        prevS = model.S(1:nMets,1:nRxns);
    
    else           
%         fprintf('Integrating both contraints to the model\n') 
        prevS = model.S;
    end 
end

function [posVec, negVec, mixVec, posVec_m, negVec_m, mixVec_m, Smixed, ...
    nRxnsExt, nMix,Mrxns] = categorizeReactions(model, nRxns, prevS)
    posVec = zeros(nRxns, 1); nPos = 0;
    negVec = zeros(nRxns, 1); nNeg = 0;
    mixVec = zeros(nRxns, 1); nMix = 0;
    Mrxns = [model.protGroup(4).rxns; model.protGroup(1).rxns; model.protGroup(5).rxns]; % Membrane Reactions
    posVec_m = zeros(nRxns, 1); nPos_m = 0;
    negVec_m = zeros(nRxns, 1); nNeg_m = 0;
    mixVec_m = zeros(nRxns, 1); nMix_m = 0;

    for r = 1:nRxns
        if model.w_c(r) ~= 0
            if model.ub(r) >= 0 && model.lb(r) >= 0 % irreversible positive reactions
                nPos = nPos + 1;
                posVec(nPos) = r;
                if ismember(r, Mrxns, 'rows')
                    nPos_m = nPos_m + 1;
                    posVec_m(nPos_m) = r;
                end
            elseif model.ub(r) <= 0 && model.lb(r) < 0 % irreversible negative reactions
                nNeg = nNeg + 1;
                negVec(nNeg) = r;
                if ismember(r, Mrxns, 'rows')
                    nNeg_m = nNeg_m + 1;
                    negVec_m(nNeg_m) = r;
                end
            else % reversible reactions
                nMix = nMix + 1;
                mixVec(nMix) = r;
                if ismember(r, Mrxns, 'rows')
                    nMix_m = nMix_m + 1;
                    mixVec_m(nMix_m) = r;
                end
            end
        end
    end

    % Position vectors of each classified reaction
    posVec = posVec(1:nPos);
    negVec = negVec(1:nNeg);
    mixVec = mixVec(1:nMix);
    posVec_m = posVec_m(1:nPos_m);
    negVec_m = negVec_m(1:nNeg_m);
    mixVec_m = mixVec_m(1:nMix_m);

    Smixed = sparse(-prevS(:, mixVec)); % stoich. mat. for extra Reactions

    nRxnsExt = nRxns + nMix; % extend the number of reactions
end

function [model,phi_max,phiM_max] = splitAndBuildMatrix(model,  nMix, mixVec, nMets, nRxns, ...
    nRxnsExt, posVec, negVec, posVec_m, negVec_m, mixVec_m, gamma, Smixed,prevS)
    lb = model.lb;
    ub = model.ub;
    c = model.c;
    names = model.rxnNames;
    ids = model.rxns;

    for n = 1:nMix
        lb(mixVec(n)) = 0; % the original one is made irreversible, set value of v(+)
        lb(nRxns + n) = 0; % a new, positive, variable is introduced
        ub(nRxns + n) = -model.lb(mixVec(n)); % set value of v(-)
        c(nRxns + n) = -model.c(mixVec(n)); % obj. vec. minimize or maximize the opposite of the actual obj func
        names(nRxns + n) = model.rxnNames(mixVec(n));
        ids(nRxns + n,:) = model.rxns(mixVec(n),:);
    end

    rrxn = (nRxns + 1):nRxnsExt;
    prxn_v = zeros(nMets, 1); % Protein constraint added as a reaction
    mrxn_v = zeros(nMets, 1); % Membrane constraint added as a reaction
    pos_prxn = nRxnsExt + 1;
    pos_mrxn = nRxnsExt + 2;

    wv_c = model.w_c(posVec);
    wv_c = [wv_c; (-model.w_c(negVec))]; % minus sign!, makes protein cost +ve 
    wv_c = [wv_c; (model.w_c(mixVec))];
    wv_c = [wv_c; model.w_c(mixVec)]; % here a + is needed, since both x+ and x- are >0
    wv_c = [wv_c; -1]; % We will use these coefficients to activate the proteome constraints  

    % Proteome constraints
    phi_Q = 0;
    for n = 1:length(model.protGroup)
        phi_Q = phi_Q + model.protGroup(n).phi0;
    end

    % adjusting phiM_max according to the gamma ratio: phiM_max / phi_max
    if gamma > 0
        phiM_max = gamma * (1 - phi_Q);
    elseif gamma < 0
        error('gamma value cannot be negative.');
    end

    phi_max = 1 - phi_Q;

    % Updating boundaries to limit the flux of the proteome reactions constraints
    lb(pos_prxn) = 0; % phi_max
    ub(pos_prxn) = phi_max; % phi_max
    ub(pos_mrxn) = phiM_max; % phiM_max
    lb(pos_mrxn) = 0; % we don't want negative values, for that we set the w coefficient to -1 / less than phiM_max
    c(pos_prxn) = 0; c(pos_mrxn) = 0; 
    names(pos_prxn) = cellstr('Proteome Constraint');
    names(pos_mrxn) = cellstr('Membrane Constraint');
    ids(pos_prxn,:) = cellstr('Prot_Const');
    ids(pos_mrxn,:) = cellstr('Memb_Const');


%   Common error is to use single comma for an horizontal concatenation
%   instead of using semicolon for a vertical concatenation
    Sac = sparse(1, [posVec; negVec; mixVec; rrxn'; pos_prxn], wv_c, 1, nRxnsExt + 2);

    wv_m = [model.w_m(posVec_m); (-model.w_m(negVec_m)); (model.w_m(mixVec_m)); (model.w_m(mixVec_m)); -1];
    mxr = find(ismember(mixVec, mixVec_m));
    mxr = nRxns + mxr;

    Sac_m = sparse(1, [posVec_m; negVec_m; mixVec_m; mxr; pos_mrxn], wv_m, 1, nRxnsExt + 2);        

    % ------ TEST -------------------------
%     checkConcatenationDimensions(prevS , Smixed, prxn_v, mrxn_v , Sac,Sac_m);
    % ----------------------------------
    
    % building the extended S matrix, as well as the constraints parameters
    A=[[prevS , Smixed, prxn_v, mrxn_v] ; Sac;Sac_m];
    % Update the model fields with the new constraints
    model = updateModelFields(model, nMets, nRxnsExt, A, lb, ub, c, names, ids);
    
end

function model = updateModelFields(model, nMets, nRxnsExt, A, lb, ub, c, names, ids)
    model.S = A;
    model.lb = lb;
    model.ub = ub;
    model.c = c;
    model.rxnNames = names;
    model.rxns = ids;
    [model.nmets, model.nrxns] = size(model.S);
    model.mets(nMets + 1) = cellstr('Protein_Allocation');
    model.mets(nMets + 2) = cellstr('Membrane_Allocation');
    model.b(nMets + 1) = 0; % rhs for the equalities
    model.b(nMets + 2) = 0; % rhs for the equalities
    model.csense(nMets + 1) = 'E'; % equality
    model.csense(nMets + 2) = 'E'; % equality
    model.subSystems{nRxnsExt + 1} = 'Unassigned'; model.subSystems{nRxnsExt + 2} = 'Unassigned';   

    model.grRules{nRxnsExt + 1, 1} = ''; model.grRules{nRxnsExt + 2, 1} = ''; 
    model.rules{nRxnsExt + 1, 1} = ''; model.rules{nRxnsExt + 2, 1} = ''; 
    model.metFormulas{nMets + 1, 1} = 'E'; % equality
    model.metFormulas{nMets + 2, 1} = 'E'; % equality
    model.metNames{nMets + 1} = 'P'; % equality
    model.metNames{nMets + 2} = 'M'; % equality
    model.metCharges(nMets + 1) = 0; model.metCharges(nMets + 2) = 0; % fix metCharges
    model.rxnGeneMat(nRxnsExt + 1, 1:length(model.genes)) = zeros(1, length(model.genes)); % fix rxnGeneMat
    model.rxnGeneMat(nRxnsExt + 2, 1:length(model.genes)) = zeros(1, length(model.genes));

    % Extract field names as a list
    fieldList = fieldnames(model);
    for f = 1:length(fieldList)
        field = fieldList{f};
        if startsWith(field, 'met')
            if length(model.(field)) < length(model.mets) 
                model.(field){nMets + 1, 1} = 'nan' ; model.(field){nMets + 2, 1} = 'nan' ;
            end
        elseif startsWith(field, 'rxns')
            if length(model.(field)) < length(model.rxns)
                model.(field){nRxnsExt + 1, 1} = 'nan' ; model.(field){nRxnsExt + 2, 1} = 'nan' ;
            end
        elseif length(model.(field)) < length(model.rxns) && length(model.(field)) > length(model.mets)
                if isnumeric(model.(field))
                    model.(field)(nRxnsExt + 1, 1) = nan ; model.(field)(nRxnsExt + 2, 1) = nan ;
                else
                    model.(field){nRxnsExt + 1, 1} = 'nan' ; model.(field){nRxnsExt + 2, 1} = 'nan' ;
                end
        end
    end
    model.rxnisprot_wCoef_cID = cellfun(@num2str, num2cell(model.w_c), 'UniformOutput', false); % annotate proteome allocation coefficients (W_c)
    model.rxnisprot_wCoef_mID = cellfun(@num2str, num2cell(model.w_m), 'UniformOutput', false); % annotate membrane allocation coefficients (W_m)
    
end

function solution = getSolution(model, Sol, nMix, nRxns, phi_max, Mrxns)
    if Sol == 1
        solution = optimizeCbModel(model);
    else
        solution = 'not included';
    end

    if solution.stat == 1
        % Put together fluxes after splitting
        for i = 1:nMix
            solution.x(mixVec(i)) = solution.x(mixVec(i)) - solution.x(nRxns + i);
        end

        % Adding the sum of each sector in the solution
        wx_m = model.w_c(Mrxns) .* abs(solution.x(Mrxns));
        solution.phiM = sum(wx_m); % phiM
        wx_c = model.w_c(model.protGroup(1).rxns) .* abs(solution.x(model.protGroup(1).rxns));
        solution.phiCm = sum(wx_c); % phiCm
        wx_cc = model.w_c(model.protGroup(7).rxns) .* abs(solution.x(model.protGroup(7).rxns));
        solution.phiCc = (sum(wx_cc) * -1) + (phi_max - solution.x(end - 1)); % phiCc
        wx_r = model.w_c(model.protGroup(3).rxns) .* abs(solution.x(model.protGroup(3).rxns));
        solution.phiR = sum(wx_r); % phiR
        wx_ec = model.w_c(model.protGroup(2).rxns) .* abs(solution.x(model.protGroup(2).rxns));
        solution.phiEc = sum(wx_ec); % phiEc
        wx_er = model.w_c(model.protGroup(4).rxns) .* abs(solution.x(model.protGroup(4).rxns));
        solution.phiEr = sum(wx_er); % phiEr
        wx_em = model.w_c(model.protGroup(5).rxns) .* abs(solution.x(model.protGroup(5).rxns));
        solution.phiEm = sum(wx_em); % phiEm
        wx_p = model.w_c(1:nRxns) .* abs(solution.x(1:nRxns));
        solution.phi_max = sum(wx_p); % phi max
    else
        solution.x = zeros(nRxns, 1);
        solution.f = 0;
        solution.phiM = 0; solution.phiCm = 0; solution.phiCc = 0;
        solution.phiR = 0; solution.phiEc = 0; solution.phiEr = 0;
        solution.phiEm = 0; solution.phi_max = 0;
    end
end

function checkConcatenationDimensions(prevS , Smixed, prxn_v, mrxn_v , Sac,Sac_m)
    disp('Dimensions of vectors:');
    disp(['prevS: ', num2str(size(prevS))]);
    disp(['Smixed: ', num2str(size(Smixed))]);
    disp(['prxn_v: ', num2str(size(prxn_v))]);
    disp(['mrxn_v: ', num2str(size(mrxn_v))]);
    disp(['Sac: ', num2str(size(Sac))]);
    disp(['Sac_m: ', num2str(size(Sac_m))]);
    V1 = [prevS , Smixed, prxn_v, mrxn_v ];
    disp(['V1: ', num2str(size(V1))]);
%     A=[[prevS , Smixed, prxn_v, mrxn_v] ; Sac;Sac_m];
%     disp(['A: ', num2str(size(A))]);

   
end
