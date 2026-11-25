function [solMAFBA,ec_model] = setUpModel_v1(pre_model,sigma_mode,addZYprot,sigma_param,gamma, testPrint)
%Function to set up a enzyme constrained model from a SBML model using both
%constraints the overall proteome allocation constraint and the membrane 
%allocation constrain. 
%   Detailed explanation goes here

if ~exist('sigma_param','var')
    sigma_param = 6.977; %optimized value for parametrized w_mode  
end
if ~exist('sigma_mode','var')
    sigma_mode = 'parametrized_totalweight';
end
if ~exist('addZYprot','var')
    addZYprot = 0;
end
if ~exist('gamma','var')
    gamma = 1;
end
if ~exist('testPrint','var')
    testPrint = 1;
end

% List of metabolites (names or IDs)
metaboliteList = {'D-Glucose'}; % Metabolite names or specific IDs

% Call the function
pre_model = set_CarbonUptakeReactions(pre_model, metaboliteList);

% Remove default bounds for glucose, oxygen and CO2
% Get the indexes for each reaction forward or reverse
ex_glc_r_f=find(strcmp(pre_model.rxns,'EX_glc__D_e_f'));
ex_glc_r_b=find(strcmp(pre_model.rxns,'EX_glc__D_e_b'));
ex_o2_r_f=find(strcmp(pre_model.rxns,'EX_o2_e_f'));
ex_o2_r_b=find(strcmp(pre_model.rxns,'EX_o2_e_b'));
ex_co2_r_f=find(strcmp(pre_model.rxns,'EX_co2_e_f'));
ex_co2_r_b=find(strcmp(pre_model.rxns,'EX_co2_e_b'));
% oxygen
pre_model.lb(ex_o2_r_b)=0; pre_model.ub(ex_o2_r_b)=20; %reverse
pre_model.ub(ex_o2_r_f)=0; pre_model.ub(ex_o2_r_f)=1000; %forward
% CO2
pre_model.lb(ex_co2_r_b)=0; pre_model.ub(ex_co2_r_b)=0; %reverse
pre_model.ub(ex_co2_r_f)=0; pre_model.ub(ex_co2_r_f)=1000; %forward

% Knocking out the GLCDe flux, as a correction to the model.
glcd_r=find(strcmp(pre_model.rxns,'GLCDpp'));
pre_model.lb(glcd_r)=0;
pre_model.ub(glcd_r)=0;

% % Printing the upper boundary of the glucose exchange
% ex_glc_r_b = find(strcmp(pre_model.rxns, 'EX_glc__D_e_b'));
% fprintf('1 The upper boundary of the glucose exchange is: %f\n', pre_model.ub(ex_glc_r_b));

% Knocking out loops
% knock_loops = {'ACALDtex_b', 'SUCCtex_b', 'CO2tex_b', 'FE2tex_b', 'PItex_b', 'NH4tex_b'};
% pre_model.ub(ismember(pre_model.rxns, knock_loops))=1000;

% Add protein expression simulation reactions
if addZYprot == 1 
    [pre_model, reactionIndices] = addZYProteins(pre_model);
end

%% 

% Add protein groups to the model.
% Categorize the reactions in a protein group requires specific list for 
% diferent models
pre_model=addProteinGroups_v1(pre_model);

% Glucose
pre_model.lb(ex_glc_r_b)=0; pre_model.ub(ex_glc_r_b)=30; %reverse
pre_model.lb(ex_glc_r_f)=0; pre_model.ub(ex_glc_r_f)=0; %forward


% % Printing the upper boundary of the glucose exchange
% ex_glc_r_b = find(strcmp(pre_model.rxns, 'EX_glc__D_e_b'));
% fprintf('2.1 The upper boundary of the glucose exchange is: %f\n', pre_model.ub(ex_glc_r_b));
%% 

% Set the offsets for the different groups.
pre_model.protGroup(1).phi0 = 0; % phiCm_0
pre_model.protGroup(2).phi0 = 0; % phiEc_0
pre_model.protGroup(3).phi0 = 0.066; % phiR_0
pre_model.protGroup(4).phi0 = 0; % phiEr_0
pre_model.protGroup(5).phi0 = 0; % phiEm_0
pre_model.protGroup(6).phi0 = 0.45; % phiQ_0
pre_model.protGroup(7).phi0 = 0; % phiCc_0


% Select the weight according with the w_mode given
if strcmp(sigma_mode,'homogeneous')
%     weight =  0.000354; % model_iML1515 homogeneous
    pre_model=setConstraintWeights(pre_model,0,sigma_param, 'mode', 'homogenous', 'rescaleFlag', 0);
%     disp(pre_model.w_c)
elseif strcmp(sigma_mode,'heterogeneous')
%     weight =  0.000354; % model_iML1515 homogeneous
    pre_model=setConstraintWeights(pre_model,0,sigma_param, 'mode', 'heterogeneous', 'rescaleFlag', 0);
%     disp(pre_model.w_c)
elseif strcmp(sigma_mode,'parametrized_totalweight')
%     weight =  7.12; % model_iML1515 parametrized
    pre_model=setConstraintWeights(pre_model,0,sigma_param, 'mode', 'parametrized', ...
            'mbw_mode', 'total', 'pMTP', 0.32, 'testPrint', testPrint ,'rescaleFlag', 2);
elseif strcmp(sigma_mode,'parametrized_TMbw')
%     weight =  7.12; % model_iML1515 parametrized
    pre_model=setConstraintWeights(pre_model,0,sigma_param, 'mode', 'parametrized', ...
    'mbw_mode', 'transmembrane', 'pMTP', 0.32, 'testPrint', testPrint ,'rescaleFlag', 2);
elseif strcmp(sigma_mode,'parametrized_MSA') || strcmp(sigma_mode,'parametrized_MSA_only')
    pre_model=setConstraintWeights(pre_model,0,sigma_param, 'mode', 'parametrized', ...
    'mbw_mode', 'MSA', 'pMTP', 1, 'testPrint', testPrint ,'rescaleFlag', 2);
    
else
    error('Sigma mode (sigma_mode) not recognized, use either homogenous, parametrized_totalweight, parametrized_TMbw or parametrized_MSA')
end

% 0.0007596 , 0.2683 ,  0.000354 , 0.4152 , 0.7150, 6.9
% Now we have to set the weights for the different groups.
% (note that the function creates the field model.weights if needed)
% % All the reactions with the fine grained approach
% pre_model=setWeights_iML(pre_model,0,weight, 'mode', 'fineg', ...
%             'mbw_mode', 'total', 'pMTP', 0.32, 'rescaleFlag', 1);
% R_rxn = 73;
% disp(pre_model.w_c(R_rxn));
% disp(pre_model.w_m(R_rxn));

% pre_model=setWeights_iML(pre_model,2,weight ); % 0.00083
% pre_model=setWeights_iML(pre_model,4,weight );
% pre_model=setWeights_iML(pre_model,5,weight ); % 0.00062

pre_model=setConstraintWeights(pre_model,3,0.169, 'mode', 'homogenous', 'rescaleFlag', 0);
pre_model=setConstraintWeights(pre_model,6,0,  'mode', 'homogenous', 'rescaleFlag', 0);
pre_model=setConstraintWeights(pre_model,7,0,  'mode', 'homogenous', 'rescaleFlag', 0); % wCc

if addZYprot == 1 
    %Correct RecP reactions
    pre_model.w_c(reactionIndices(1)) = 1; pre_model.w_c(reactionIndices(2)) = 1; 
    pre_model.w_m(reactionIndices(1)) = 0; pre_model.w_m(reactionIndices(2)) = 1; 
end

% % Printing the upper boundary of the glucose exchange
% ex_glc_r_b = find(strcmp(pre_model.rxns, 'EX_glc__D_e_b'));
% fprintf('2.2 The upper boundary of the glucose exchange is: %f\n', pre_model.ub(ex_glc_r_b));

% disp(pre_model.w_c)

%Run the FBA
[solMAFBA,ec_model] = optimizeMAFBAmodel(pre_model, 'gamma', gamma);
% % Printing the upper boundary of the glucose exchange
% ex_glc_r_b = find(strcmp(ec_model.rxns, 'EX_glc__D_e_b'));
% fprintf('3 The upper boundary of the glucose exchange is: %f\n', ec_model.ub(ex_glc_r_b));

if strcmp(sigma_mode,'parametrized_MSA')
    fMSA = gamma;
    ec_model = implementMSAModel(ec_model, fMSA);
    solMAFBA = get_solution_ec_model(ec_model);
    fprintf('The maximum growth rate is %1.4f (mode: %s)\n', solMAFBA.f, sigma_mode);
elseif strcmp(sigma_mode,'parametrized_MSA_only')
    fMSA = gamma;
    ec_model = implementMSAModel(ec_model, fMSA);
    ec_model.S(end-1,1:end-1) = 0; % Set the last row (except last column) of S matrix to zeros    
    solMAFBA = get_solution_ec_model(ec_model);
    fprintf('The maximum growth rate is %1.4f (mode: %s)\n', solMAFBA.f, sigma_mode);
else
    fprintf('The maximum growth rate is %1.4f \n', solMAFBA.f);
end

% Clear useless variables
% clear('ex_glc_r','ex_o2_r', 'ex_co2_r', 'glcd_r', 'weight');
end

function [pre_model, reactionIndices] = addZYProteins(pre_model)
    % Add proteins Z and Y to the model
    
    if ~isfield(pre_model, 'rxns') 
        error('Invalid model format. Please provide a valid COBRA model structure.');
    end
    
    if nargin < 1
        error('Missing input argument: pre_model');
    end
    
    % prot_Z is a cytosolic protein
    % prot Y is a membrane protein
    pre_model = addReaction(pre_model, 'EX_pre_z','reactionFormula', 'pre_z <=> ');
    pre_model = addReaction(pre_model, 'recP_Z','reactionFormula', 'pre_z -> prot_Z');
    pre_model = addReaction(pre_model, 'EX_Z','reactionFormula', 'prot_Z <=> ');

    pre_model = addReaction(pre_model, 'EX_pre_y','reactionFormula', 'pre_y <=> ');
    pre_model = addReaction(pre_model, 'recP_Y','reactionFormula', 'pre_y -> prot_Y');
    pre_model = addReaction(pre_model, 'EX_Y','reactionFormula', 'prot_Y <=> ');

    reactionNamesToSearch = {'recP_Z', 'recP_Y'};
    reactionIndices = zeros(size(reactionNamesToSearch));

    for i = 1:length(reactionNamesToSearch)
        reactionIndex = find(strcmp(pre_model.rxns, reactionNamesToSearch{i}));
        if ~isempty(reactionIndex)
            reactionIndices(i) = reactionIndex;
        else
            warning(['Reaction ' reactionNamesToSearch{i} ' not found in the model.']);
        end
%         Location specific annotations
        if strcmp(reactionNamesToSearch{i}, 'recP_Y')
            pre_model.location{reactionIndex} = 'Inner Membrane';
            pre_model.subSystems{reactionIndex} = 'Putative Transporter';
            pre_model.mw_t(reactionIndex) = 1000;
        else
            pre_model.location{reactionIndex} = 'Cytoplasmic';
            pre_model.subSystems{reactionIndex} = 'Unassigned';
            pre_model.mw_t(reactionIndex) = 0;
        end
%         General annotations main reaction
        pre_model.kcat(reactionIndex) = 1/3600; 
        pre_model.mw_c(reactionIndex) = 1000;
        pre_model.lb(reactionIndex) = 0;
        pre_model.ub(reactionIndex) = 0;
%         General annotations exchange side reactions
        pre_model.kcat(reactionIndex-1) = 1e12; 
        pre_model.mw_c(reactionIndex-1) = 0;
        pre_model.mw_t(reactionIndex-1) = 0;
        pre_model.lb(reactionIndex-1) = -1000;
        pre_model.ub(reactionIndex-1) = 0; 
        pre_model.location{reactionIndex-1} = 'Exchange';
        pre_model.subSystems{reactionIndex-1} = 'Unassigned';
        
        pre_model.kcat(reactionIndex+1) = 1e12; 
        pre_model.mw_c(reactionIndex+1) = 0;
        pre_model.mw_t(reactionIndex+1) = 0;
        pre_model.lb(reactionIndex+1) = 0;
        pre_model.ub(reactionIndex+1) = 1000; 
        pre_model.location{reactionIndex+1} = 'Exchange';
        pre_model.subSystems{reactionIndex+1} = 'Unassigned';
    end
    
end

