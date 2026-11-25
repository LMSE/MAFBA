function model = addProteinGroups_v1(model, optC, phiEopt)
% addProteinGroups_v1 - Adds the 'protGroup' field to a COBRA model with membrane allocation groups.
%
% Inputs:
%    model : COBRA metabolic model
%    optC  : Chooses the content of the phiC group (see setPhiC subfunction).
%
% Output:
%    model : Modified model with the additional 'protGroup' field.

% Default value for optC
if ~exist('optC', 'var')
    optC = 'glc';
end
if ~exist('phiEopt', 'var')
    phiEopt = 'subSys';
end

% Reactions
Biomass = find(model.c);

protGroup = initializeDataStructure(7,{'name','phi0','rxns'});

% - Allocating reactions to each protein group according to the subsystem in the model
if ~isfield(model, 'protGroup')
    % Define phiR group
    protGroup(3) = createPhiRGroup(Biomass);

    % Define phiCm and phiCc groups
    [model, protGroup(1), protGroup(7)] = setPhiC(model, optC);

    % Remove overlapping reactions
    protGroup(7).rxns = removeRxnsOverlap(protGroup(7).rxns, protGroup(1).rxns);

    % Define phiEc, phiEr, and phiEm groups
    [protGroup(2), protGroup(4), protGroup(5)] = setPhiEList(model, phiEopt);

    % Remove overlapping reactions
    protGroup(2).rxns = removeRxnsOverlap(protGroup(2).rxns, [protGroup(1).rxns; protGroup(7).rxns; protGroup(3).rxns]);
    protGroup(4).rxns = removeRxnsOverlap(protGroup(4).rxns, [protGroup(1).rxns; protGroup(2).rxns; protGroup(3).rxns; protGroup(7).rxns]);
    protGroup(5).rxns = removeRxnsOverlap(protGroup(5).rxns, [protGroup(1).rxns; protGroup(2).rxns; protGroup(3).rxns; protGroup(7).rxns; protGroup(4).rxns]);

    % Define phiQ and others group
    protGroup(6) = createPhiQGroup(model, protGroup);

    % Assign protGroup to the model
    model.protGroup = protGroup;
else
    % Update phiCm and phiCc groups
    [model, protGroup(1), protGroup(7)] = setPhiC(model, optC);
    model.protGroup(1) = protGroup(1);
    model.protGroup(7) = protGroup(7);
end

end

function group = createPhiRGroup(Biomass)
    % Create phiR group
    group.name = 'phiR';
    group.phi0 = 0.066;
    group.rxns = Biomass;
end

function [model, phiCm, phiCc] = setPhiC(model, optC)
    % Set phiCm and phiCc groups
    phiCm.name = 'phiCm';
    phiCm.phi0 = 0;

    phiCc.name = 'phiCc';
    phiCc.phi0 = 0;

    % Implement the logic to set rxns_Cm and rxns_Cc
    [model, rxns_Cm, rxns_Cc] = setPhiCLogic(model, optC);

    phiCm.rxns = rxns_Cm;
    phiCc.rxns = rxns_Cc;
end

function group = createPhiQGroup(model, protGroup)
    % Create phiQ and others group
    group.name = 'phiQ and others';
    group.phi0 = 0;
    group.rxns = (1:length(model.rxns))';

    % Remove overlapping reactions
    for i = 1:numel(protGroup)
        group.rxns = removeRxnsOverlap(group.rxns, protGroup(i).rxns);
    end
end

function rxns = removeRxnsOverlap(rxns1, rxns2)
    % Remove overlapping reactions between two sets of reactions
    rxns = setdiff(rxns1, rxns2);
end

function [phiEc, phiEr, phiEm] = setPhiEList(model, phiEopt)
    % Set phiEc, phiEr, and phiEm groups
    phiEc.name = 'phiEc';
    phiEc.phi0 = 0;
    phiEr.name = 'phiEr';
    phiEr.phi0 = 0;
    phiEm.name = 'phiEm';
    phiEm.phi0 = 0;

    if isfield(model, 'rxnislocationID')
        phiEopt = 'location';
        [phiEc.rxns, phiEr.rxns, phiEm.rxns] = set_PhiE_list_wLocation(model, phiEopt);
    else
        [phiEc.rxns, phiEr.rxns, phiEm.rxns] = set_PhiE_list(model, phiEopt);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiC rxns vector and bounds %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,rxns_Cm, rxns_Cc] = setPhiCLogic(model,optC)
% This function also updates the bound for the exchange reactions

EX_glc=find(strcmp(model.rxns,'EX_glc__D_e_b')); % glc transport

% Upper glycolysis
if strcmp(optC,'default') || strcmp(optC,'glc')
    rxns_Cm=find(ismember(model.rxns,{'GLCptspp','GLCt2pp'} ));
    rxns_Cc=EX_glc;model.ub(EX_glc)= 1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'g6p')
    EX_g6p=find(strcmp(model.rxns,'EX_g6p_e')); % g6p transport
    t_g6p=find(strcmp(model.rxns,'G6Pt6_2pp')); % g6p transport
    rxns_Cm=t_g6p;
    rxns_Cc=EX_g6p;
    model.lb(EX_g6p)=-1000;
    model.ub(EX_glc)=0; 
elseif strcmp(optC,'man6p')
    EX_man6p=find(strcmp(model.rxns,'EX_man6p_e')); % man6p transport
    t_man6p=find(strcmp(model.rxns,'MAN6Pt6_2pp')); % man6p transport
    rxns_Cc=EX_man6p;
    rxns_Cm=t_man6p;
    model.lb(EX_man6p)=-1000;
    model.ub(EX_glc)=0; 
elseif strcmp(optC,'mnl')
    EX_mnl=find(strcmp(model.rxns,'EX_mnl_e')); % mannitol transport
    t_mnl=find(strcmp(model.rxns,'MNLptspp')); % mannitol transport
    rxns_Cm=t_mnl;
    rxns_Cc=EX_mnl;
    model.lb(EX_mnl)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'lcts')
    EX_lcts=find(strcmp(model.rxns,'EX_lcts_e_') ); % lactose transport
    t_lcts=find(strcmp(model.rxns,'LCTSt') ); % lactose transport
    rxns_Cc=EX_lcts;
    rxns_Cm=t_lcts;
    model.lb(EX_lcts)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'sucr')
    EX_sucr=find(strcmp(model.rxns,'EX_sucr_e_') ); % sucrose transport
    t_sucr=find(strcmp(model.rxns,'SUCpts') ); % lactose transport
    rxns_Cc=EX_sucr;
    rxns_Cm=t_sucr;
    model.lb(EX_sucr)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'malt')
    EX_malt=find(strcmp(model.rxns,'EX_malt_e_') ); % maltose transport
    t_malt=find(ismember(model.rxns,{'MALTabc','MALTpts'} )); % maltose transport
    rxns_Cc=EX_malt;
    rxns_Cm=t_malt;
    model.lb(EX_malt)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'fru')
    EX_fru=find(strcmp(model.rxns,'EX_fru_e_')); % fructose transport
    t_fru=find(ismember(model.rxns,{'FRUpts','FRUpts2'} )); % fructose transport
    rxns_Cc=EX_fru;     rxns_Cm=t_fru;
    model.lb(EX_fru)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'sbt')
    EX_sbt=find(strcmp(model.rxns,'EX_sbt_D_e_')); % sorbitol
    t_sbt=find(strcmp(model.rxns,'SBTpts')); % sorbitol
    rxns_Cc=EX_sbt;
    rxns_Cm=t_sbt;
    model.lb(EX_sbt)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'man')
    EX_man=find(strcmp(model.rxns,'EX_man_e_')); % mannose
    t_man=find(strcmp(model.rxns,'MANpts')); % mannose
    rxns_Cc=EX_man;
    rxns_Cm=t_man;
    model.lb(EX_man)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'glcn')
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    t_glcn=find(strcmp(model.rxns,'GLCNt2r')); % gluconate
    rxns_Cc=EX_glcn;
    rxns_Cm=t_glcn;
    model.lb(EX_glcn)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'arab')
    EX_arab=find(strcmp(model.rxns,'EX_arab_L_e_')); % arabinose
    t_arab=find(ismember(model.rxns,{'ARBabc','ARBt2r'} )); % arabinose
    rxns_Cc=EX_arab;
    rxns_Cm=t_arab;
    model.lb(EX_arab)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'rib')
    EX_rib=find(strcmp(model.rxns,'EX_rib_D_e_')); % ribose
    t_rib=find(strcmp(model.rxns,'RIBabc')); % ribose
    rxns_Cc=EX_rib;
    rxns_Cm=t_rib;
    model.lb(EX_rib)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'xyl')
    EX_xyl=find(strcmp(model.rxns,'EX_xyl_D_e_')); % xylose
    t_xyl=find(ismember(model.rxns,{'XYLabc','XYLt2'} )); % xylose
    rxns_Cc=EX_xyl;
    rxns_Cm=t_xyl;
    model.lb(EX_xyl)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'gal')
    EX_gal=find(strcmp(model.rxns,'EX_gal_e_')); % galactose
    t_gal=find(ismember(model.rxns,{'GALabc','GALt2'} )); % galactose
    rxns_Cc=EX_gal;
    rxns_Cm=t_gal;
    model.lb(EX_gal)=-1000;
    model.ub(EX_glc)=0;


elseif strcmp(optC,'g6p+glcn')
    EX_g6p=find(strcmp(model.rxns,'EX_g6p_e_')); % g6p transport
    t_g6p=find(strcmp(model.rxns,'G6Pt6_2')); % g6p transport
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    t_glcn=find(strcmp(model.rxns,'GLCNt2r')); % gluconate
    rxns_Cc=[EX_glcn,EX_g6p];
    rxns_Cm=[t_g6p,t_glcn];
    model.lb(EX_glcn)=-1000;
    model.lb(EX_g6p)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'man6p+glcn')
    EX_man6p=find(strcmp(model.rxns,'EX_man6p_e_')); % man6p transport
    t_man6p=find(strcmp(model.rxns,'MAN6Pt6_2')); % man6p transport
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    t_glcn=find(strcmp(model.rxns,'GLCNt2r')); % gluconate
    rxns_Cc=[EX_glcn,EX_man6p];
    rxns_Cm=[t_man6p,t_glcn];
    model.lb(EX_glcn)=-1000;
    model.lb(EX_man6p)=-1000;
    model.ub(EX_glc)=0;
    
% Mid glycolysis
elseif strcmp(optC,'glcr')
    EX_glcr=find(strcmp(model.rxns,'EX_glcr_e_')); % glucorate
    t_glcr=find(strcmp(model.rxns,'GLCRt2r')); % glucorate
    rxns_Cc=EX_glcr;
    rxns_Cm=t_glcr;
    model.lb(EX_glcr)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'dha')
    EX_dha=find(strcmp(model.rxns,'EX_dha_e_')); % dihydroxyacetone
    t_dha=find(strcmp(model.rxns,'DHAt')); % dihydroxyacetone
    rxns_Cc=EX_dha;
    rxns_Cm=t_dha;
    model.lb(EX_dha)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'glyc')
    EX_glyc=find(strcmp(model.rxns,'EX_glyc_e_')); % glycerol
    t_glyc=find(strcmp(model.rxns,'GLYCt')); % glycerol
    rxns_Cc=EX_glyc;
    rxns_Cm=t_glyc;
    model.lb(EX_glyc)=-1000;
    model.ub(EX_glc)=0;

%TODO: Continue adding more reactions
    
    % TCA C-lim
elseif strcmp(optC,'pyr')
    EX_pyr=find(strcmp(model.rxns,'EX_pyr_e_')); % pyruvate
    t_pyr=find(strcmp(model.rxns,'PYRt2r')); % pyruvate
    rxns_Cc=EX_pyr;
    rxns_Cm=t_pyr;
    model.lb(EX_pyr)=-1000;
    model.ub(EX_glc)=0; 
elseif strcmp(optC,'akg')
    EX_akg=find(strcmp(model.rxns,'EX_akg_e_')); % akg
    t_akg=find(strcmp(model.rxns,'AKGt2r')); % akg
    rxns_Cc=EX_akg;
    rxns_Cm=t_akg;
    model.lb(EX_akg)=-1000;
    model.ub(EX_glc)=0; 
elseif strcmp(optC,'cit')
    EX_cit=find(strcmp(model.rxns,'EX_cit_e_')); % citrate
    t_cit=find(strcmp(model.rxns,'CITt7')); % citrate
    rxns_Cc=EX_cit;
    rxns_Cm=t_cit;
    model.lb(EX_cit)=-1000;
    model.ub(EX_glc)=0; 
elseif strcmp(optC,'succ')
    EX_succ=find(strcmp(model.rxns,'EX_succ_e_')); % succinate
    t_succ=find(ismember(model.rxns,{'CITt7','SUCCabc','SUCCt2_2','SUCCt2_3','SUCCt2b','SUCFUMt','TARTRt7','SUCCabc'} )); % succinate
    rxns_Cc=EX_succ;
    rxns_Cm=t_succ;
    model.lb(EX_succ)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'fum')
    EX_fum=find(strcmp(model.rxns,'EX_fum_e_')); % fumarate
    t_fum=find(ismember(model.rxns,{'FUMt2_2','FUMt2_3', 'SUCFUMt' } )); % fumarate
    rxns_Cc=EX_fum;
    rxns_Cm=t_fum;
    model.lb(EX_fum)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'mal')
    EX_mal=find(strcmp(model.rxns,'EX_mal_L_e_') | strcmp(model.rxns,'EX_mal-L_e_')); % malate
    t_mal=find(ismember(model.rxns,{'MALt2_2','MALt2_3'} ));  % malate
    rxns_Cc=EX_mal;
    rxns_Cm=t_mal;
    model.lb(EX_mal)=-1000;
    model.ub(EX_glc)=0;
    
% other C-lim
elseif strcmp(optC,'glcn')
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % lactate transport
    t_glcn=find(strcmp(model.rxns,'GLCNt2r'));
    rxns_Cc=EX_glcn;
    rxns_Cm=t_glcn;
    model.lb(EX_glcn)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'lac')
    EX_lac=find(strcmp(model.rxns,'EX_lac_D_e_') | strcmp(model.rxns,'EX_lac_D_e_')); % lactate transport
    t_lac=find(strcmp(model.rxns,'D_LACt2'));
    rxns_Cc=EX_lac;
    rxns_Cm=t_lac;
    model.lb(EX_lac)=-1000;
    model.ub(EX_glc)=0; 
elseif strcmp(optC,'ac')
    EX_ac=find(strcmp(model.rxns,'EX_ac_e_')); % acetate
    t_ac=find(strcmp(model.rxns,'ACt2r')); % acetate
    rxns_Cc=EX_ac;
    rxns_Cm=t_ac;
    model.lb(EX_ac)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'for')
    EX_for=find(strcmp(model.rxns,'EX_for_e_')); % formate
    t_for=find(strcmp(model.rxns,'FORt')); % formate
    rxns_Cc=EX_for;
    rxns_Cm=t_for;
    model.lb(EX_for)=-1000;
    model.ub(EX_glc)=0;
elseif strcmp(optC,'etoh')
    EX_etoh=find(strcmp(model.rxns,'EX_etoh_e_')); % ethanol
    t_etoh=find(strcmp(model.rxns,'ETOHt2r')); % ethanol
    rxns_Cc=EX_etoh;    
    rxns_Cm=t_etoh;
    model.lb(EX_etoh)=-1000;
    model.ub(EX_glc)=0;
else
    disp('Error in ''addProteinGroupsToModel'': ''optC'' is not valid. Exiting.');
    return;
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiE rxns vector  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phiEcVec, phiErVec, phiEmVec]=set_PhiE_list(model,phiEopt)

nRxns=length(model.rxns);
phiEcVec=[];
[subSysL_q, subSysL_m, tprt_diff, Resp_rxnlist] = getSubsysList(model.modelID);

% find the constraint vector
if strcmpi(phiEopt, 'All')    
    for n=1:nRxns
        if n~=nGrowth
            phiEcVec=cat(1,phiEcVec,n);
        end
    end

elseif strcmpi(phiEopt, 'subSys')
    M_subSysList = cell([length(model.subSystems) 1]);
    for k1 = 1:length(model.subSystems)
        M_subSysList{k1} = cellstr(model.subSystems{k1,1});
    end
    M_subSysList = M_subSysList{1:k1};
    % Convert M_subSysList to a cell array of character vectors
    M_subSysList = cellstr(M_subSysList);
    % Convert subSysL_q to a cell array of character vectors
    subSysL_q = cellstr(subSysL_q);
    
    phiErVec = find(ismember(model.rxns, Resp_rxnlist));
    phiEcVec=find(~ismember(M_subSysList,subSysL_q));
    phiEcVec= setdiff(phiEcVec,phiErVec);
    phiEmVec= find(ismember(M_subSysList,subSysL_m));
    difVec = find(ismember(model.rxns, tprt_diff)); 
    phiEmVec= setdiff(phiEmVec,difVec);   % Remove difussion transport
%     fprintf('%i reactions inserted in phiE.\n',length(phiEcVec));
else
    disp('No valid option in set_PhiE_list.');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiE rxns vector  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phiEcVec, phiErVec, phiEmVec]=set_PhiE_list_wLocation(model,phiEopt)

    nRxns=length(model.rxns);
    phiEcVec=[];

    % find the constraint vector
    if strcmpi(phiEopt, 'All')    
        for n=1:nRxns
            if n~=nGrowth
                phiEcVec=cat(1,phiEcVec,n);
            end
        end
    elseif strcmpi(phiEopt, 'location')
        % This option works with iJR904, iAF1260 and iJO1366

        M_loclist = model.location;
        M_subSysList = cell([length(model.subSystems) 1]);
        for k1 = 1:length(model.subSystems)
            M_subSysList(k1) =  cellstr(model.subSystems{k1});
        end
        M_subSysList = M_subSysList(1:k1);

        phiErVec= find(contains(M_loclist , 'Inner Membrane') & ...
                      contains(M_subSysList , 'Oxidative Phosphorylation'));
        phiEcVec= find(~contains(M_loclist , 'Inner Membrane') & ...
                       ~contains(M_loclist , 'Protein_constraint'));
        phiEmVec= find(contains(M_loclist , 'Inner Membrane'));
        phiEmVec= setdiff(phiEmVec,phiErVec);

        % Diffusion Transport will be removed since the cost would be neglected
    %     since kcat very high and mw 1 = 1/10000000;
    else
        disp('No valid option in set_PhiE_list.');
    end
end

% function vec3=removeRxnsOverlap(vec1,vec2)
% % Subtracts from vec1 all entries  which also appear in vec2
% vec3=vec1(~ismember(vec1,vec2));
% end

function [subSysL_q, subSysL_m, tprt_diff, Resp_rxnlist]=getSubsysList(model_id)

% Function to choice the list to use to classified the reactions using the
% subsystem attribute of the reaction
% 

% TODO> CHANGE THIS FUNCTION TO IF FORM TO SELECT LIST FOR EACH MODEL

% Create a nested structure
Smapobj = struct();
Smapobj.iJR904.subSysL_q ={'Transport Extracellular'; ...
        'Transport Inner Membrane'; 'Transport Outer Membrane'; ...
        'Transport Outer Membrane Porin'; 'Putative Transporters'; ...
        'Exchange'; 'Others'; 'Alternate Carbon Metabolism'; ''};
Smapobj.iJR904.subSysL_m = {'Transport Extracellular'; ...
        'Transport Inner Membrane'; 'Putative Transporters'};
Smapobj.iJR904.tprt_diff= {'CO2t', 'DHAt', 'FORt', 'GLYALDt', ...
        'H2Ot', 'O2t', 'UREAt'};
Smapobj.iJR904.Resp_rxnlist = { 'ATPS4r'; 'CYTBD'; 'CYTBO3'; 'FDH2';'FDH3'; 
            'NADH5'; 'NADH6'; 'NADH7'; 'NADH8'; 'NADH9'};
    
    
Smapobj.iJO1366.subSysL_q={'Transport Extracellular'; ...
        'Transport Inner Membrane'; 'Transport Outer Membrane'; ...
        'Transport Outer Membrane Porin'; 'Putative Transporters'; ...
        'Intracellular demand'; 'Extracellular exchange'; ...
        'Unassigned'; 'Alternate Carbon Metabolism'; ''};
Smapobj.iJO1366.subSysL_m={'Transport Extracellular'; ...
        'Transport Inner Membrane'; 'Putative Transporters'};
Smapobj.iJO1366.tprt_diff={'CO2tpp'; 'DHAtpp'; 'DMSOtpp'; ...
        'ETOHtrpp'; 'FORtppi'; 'FALDtpp'; 'GLYALDtpp'; 'GUAtpp'; ...
        'H2Otpp'; 'H2tpp'; 'MEOHtrpp'; 'UREAtpp'; 'O2t'};
Smapobj.iJO1366.Resp_rxnlist={ 'ATPS4rpp'; 'CYTBD2pp'; 'CYTBDpp';...
            'CYTBO3_4pp'; 'FDH4pp';'FDH5pp'; 'NADH16pp';'NADH5';...
            'NADH17pp'; 'NADH18pp';'NADH9'; 'GLCDpp'; 'HYD1pp'; ...
            'HYD2pp'; 'HYD3pp'; 'NADH16pp'; 'NADH17pp'; 'NADH18pp'; ...
            'NO3R1pp'; 'NO3R2pp'; 'SUCDi'; 'THD2pp'; 'TMAOR1'; ...
            'TMAOR1pp'; 'TMAOR2'; 'TMAOR2pp'}; 

Smapobj.iML1515.subSysL_q ={'Transport Extracellular'; ...
        'Transport Inner Membrane'; 'Transport Outer Membrane'; ...
        'Transport Outer Membrane Porin'; 'Putative Transporters'; ...
        'Intracellular demand'; 'Extracellular exchange'; ...
        'Unassigned'; 'Alternate Carbon Metabolism'; ''};
Smapobj.iML1515.subSysL_m ={'Transport Extracellular'; ...
        'Transport Inner Membrane'; 'Putative Transporters'};
Smapobj.iML1515.tprt_diff ={'CO2tpp'; 'DHAtpp'; 'DMSOtpp'; ...
        'ETOHtrpp'; 'FORtppi'; 'FALDtpp'; 'GLYALDtpp'; 'GUAtpp'; ...
        'H2Otpp'; 'H2tpp'; 'MEOHtrpp'; 'UREAtpp'; 'O2tpp';'NOtpp';...
        'N2Otpp';'SO2tpp'; 'MG2tpp';'CLtipp'; 'ALAtpp'; 'DMSOtpp';...
        'URAtpp'; 'ASPtpp'; 'FUMt1pp'; 'SUCCt1pp' };
    
Smapobj.iML1515.Resp_rxnlist={ 'ATPS4rpp'; 'CYTBD2pp'; 'CYTBDpp';...
            'CYTBO3_4pp'; 'FDH4pp';'FDH5pp'; 'NADH16pp';'NADH5';...
            'NADH17pp'; 'NADH18pp';'NADH9'; 'GLCDpp'; 'HYD1pp'; ...
            'HYD2pp'; 'HYD3pp'; 'NADH16pp'; 'NADH17pp'; 'NADH18pp'; ...
            'NO3R1pp'; 'NO3R2pp'; 'SUCDi'; 'THD2pp'; 'TMAOR1'; ...
            'TMAOR1pp'; 'TMAOR2'; 'TMAOR2pp';'DMSOR1'; 'DMSOR1pp';...
            'DMSOR2'; 'DMSOR2pp'}; 
        
structureModel = getfield(Smapobj, model_id);        
% Access values using a string as the key
subSysL_q = cellstr(structureModel.subSysL_q); 
subSysL_m = cellstr(structureModel.subSysL_m);  
tprt_diff = cellstr(structureModel.tprt_diff);  
Resp_rxnlist = cellstr(structureModel.Resp_rxnlist);  

end

function dataStructure = initializeDataStructure(numStructures, fields)
    % INITIALIZEDATASTRUCTURE initializes a cell array to store structures


    % Create and assign structures with the specified fields
    for i = 1:numStructures
        dataStructure(i) = createStructure(fields);
    end
end

function structure = createStructure(fields)
    % Create a structure with specified fields
    
    % Initialize the structure
    structure = struct();

    % Initialize fields with default values or leave them empty
    for j = 1:length(fields)
        structure.(fields{j}) = [];
    end
end

