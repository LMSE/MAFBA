function ec_model = setUp_GlucoseGrowth(ec_model,dateVersion, modelV)
    % Check if attributes are provided, otherwise set default
    if nargin < 1 || isempty(ec_model)
        ec_model = 'EMPTY'; % Default 
    end
    if nargin < 2 || isempty(dateVersion)
        dateVersion = 'AlterDfv2'; % Default 
    end
    if nargin < 3 || isempty(modelV)
        modelV = 0; % Default 
    end

    % Check if ec_model is provided, otherwise load the default model
    if isstruct(ec_model)
        disp('ec_model exists.');
        % Your code here for when ec_model exists

    else
        disp('ec_model does not exist.');
        % Your code here for when ec_model does not exist
        gammas = [1];
        saveMATLAB = 0;
        saveSBML = 0;
        opt = 1;
        alearn = 1;
        [~, ec_models] = loadEcModel(gammas, 'modelV',modelV, 'dateVersion', dateVersion, ...
            'saveMATLAB', saveMATLAB , 'saveSBML', saveSBML, 'opt' , opt, 'alearn', alearn);
        ec_model = ec_models{1};

    end

% Remove default bounds for glucose, oxygen, and CO2
    % Get the indexes for each reaction forward or reverse
    ex_glc_r_f = find(strcmp(ec_model.rxns, 'EX_glc__D_e_f'));
    ex_glc_r_b = find(strcmp(ec_model.rxns, 'EX_glc__D_e_b'));
    ex_o2_r_f = find(strcmp(ec_model.rxns, 'EX_o2_e_f'));
    ex_o2_r_b = find(strcmp(ec_model.rxns, 'EX_o2_e_b'));
    ex_co2_r_f = find(strcmp(ec_model.rxns, 'EX_co2_e_f'));
    ex_co2_r_b = find(strcmp(ec_model.rxns, 'EX_co2_e_b'));
     
    % Set bounds for oxygen
    ec_model.lb(ex_o2_r_b) = 0; 
    ec_model.ub(ex_o2_r_b) = 20; % Reverse
    ec_model.lb(ex_o2_r_f) = 0; 
    ec_model.ub(ex_o2_r_f) = 1000; % Forward
    
    % Set bounds for CO2
    ec_model.lb(ex_co2_r_b) = 0; 
    ec_model.ub(ex_co2_r_b) = 0; % Reverse
    ec_model.lb(ex_co2_r_f) = 0; 
    ec_model.ub(ex_co2_r_f) = 1000; % Forward
    
    % Set bounds for Glucose
    ec_model.lb(ex_glc_r_b) = 0; 
    ec_model.ub(ex_glc_r_b) = 30; % Reverse
    ec_model.lb(ex_glc_r_f) = 0; 
    ec_model.ub(ex_glc_r_f) = 1000; % Forward

    % Knock out the GLCDe flux as a correction to the model
    glcd_r = find(strcmp(ec_model.rxns, 'GLCDpp'));
    ec_model.lb(glcd_r) = 0;
    ec_model.ub(glcd_r) = 0;

    % Knock out loops (if needed)
    % knock_loops = {'ACALDtex_b', 'SUCCtex_b', 'CO2tex_b', 'FE2tex_b', 'PItex_b', 'NH4tex_b'};
    knock_loops = {'LDH_D2'};
    ec_model.lb(ismember(ec_model.rxns, knock_loops)) = 0;
    ec_model.ub(ismember(ec_model.rxns, knock_loops)) = 0;
     

end
