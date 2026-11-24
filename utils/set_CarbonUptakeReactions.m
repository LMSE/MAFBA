function modifiedModel = set_CarbonUptakeReactions(inputModel, metaboliteList)
    % MODIFYTRANSPORTREACTIONS Modify COBRA model to restrict metabolite transport
    % INPUTS:
    %   inputModel - COBRA model structure
    %   metaboliteList - Cell array of metabolite IDs or names to modify
    % OUTPUT:
    %   modifiedModel - Modified COBRA model with updated transport reactions

    % Initialize the modified model as a copy of the input model
    modifiedModel = inputModel;

    % Convert metabolite names to metabolite IDs (if provided as names)
    metaboliteIDList = {};
    for i = 1:length(metaboliteList)
        metabolite = metaboliteList{i};

        % First check if the input looks like a metabolite ID (contains '_e' suffix)
        if ismember(metabolite, inputModel.mets) || contains(metabolite, '_e')
            fprintf('Using metabolite ID: %s\n', metabolite);
            if ~ismember(metabolite, inputModel.mets)
                warning('Metabolite ID "%s" not found in model mets. Will attempt to find it with variations.', metabolite);
                % Try common variations (with brackets, etc.)
                candidates = {
                    metabolite,
                    [metabolite '[e]'],
                    strrep(metabolite, '_e', '[e]'),
                    strrep(metabolite, '_e', '(e)')
                };
                
                for j = 1:length(candidates)
                    if ismember(candidates{j}, inputModel.mets)
                        metabolite = candidates{j};
                        fprintf('Found variation: %s\n', metabolite);
                        break;
                    end
                end
            end
            metaboliteIDList{end+1} = metabolite;
        elseif ismember(metabolite, inputModel.metNames) % Metabolite name
            % Find metabolite IDs for `[e]` and `[c]` compartments only
            fprintf('Using metabolite name: %s\n', metabolite);
            metIndexes = find(strcmp(inputModel.metNames, metabolite));
            if isempty(metIndexes)
                warning('Metabolite name "%s" not found. Checking for partial matches.', metabolite);
                % Try partial string matching
                for k = 1:length(inputModel.metNames)
                    if contains(inputModel.metNames{k}, metabolite, 'IgnoreCase', true)
                        metIndexes(end+1) = k;
                    end
                end
            end
            
            for idx = 1:length(metIndexes)
                metID = inputModel.mets{metIndexes(idx)};
                if endsWith(metID, '[e]') || endsWith(metID, '[c]') || ...
                   endsWith(metID, '(e)') || endsWith(metID, '(c)') || ...
                   contains(metID, '_e') || contains(metID, '_c')
                    metaboliteIDList{end+1} = metID;
                    fprintf('Added metabolite ID: %s\n', metID);
                end
            end
        else
            % If not found by exact match, try to find metabolites that contain this string
            potentialMatches = {};
            for j = 1:length(inputModel.metNames)
                if contains(inputModel.metNames{j}, metabolite, 'IgnoreCase', true)
                    potentialMatches{end+1} = inputModel.metNames{j};
                end
            end
            
            if ~isempty(potentialMatches)
                warning('Metabolite "%s" not found exactly. Found similar metabolites: %s', ...
                    metabolite, strjoin(potentialMatches(1:min(3,length(potentialMatches))), ', '));
                if length(potentialMatches) > 3
                    fprintf('... and %d more potential matches\n', length(potentialMatches)-3);
                end
            end
            
            error('Metabolite "%s" not found in the model.', metabolite);
        end
    end

    % Validate that all metabolites in the final list exist in the model
    metaboliteIDList = unique(metaboliteIDList); % Remove duplicates
    missingMets = setdiff(metaboliteIDList, modifiedModel.mets);
    if ~isempty(missingMets)
        error('The following metabolites were not found in the model: %s', strjoin(missingMets, ', '));
    end

    % Process each metabolite ID
    for j = 1:length(metaboliteIDList)
        metaboliteID = metaboliteIDList{j};
        metaboliteIndex = find(strcmp(modifiedModel.mets, metaboliteID));
        if isempty(metaboliteIndex)
            warning('Metabolite ID %s not found in model. Skipping.', metaboliteID);
            continue;
        end

        % Find reactions involving the specified metabolite
        relatedReactions = findRxnsFromMets(modifiedModel, metaboliteID);
        fprintf('Found %d reactions for metabolite %s\n', length(relatedReactions), metaboliteID);

        for i = 1:length(relatedReactions)
            rxnID = relatedReactions{i};
            rxnIndex = find(strcmp(modifiedModel.rxns, rxnID));

            % Check if the reaction subsystem includes Transport or Exchange
            if isfield(modifiedModel, 'subSystems') && ~isempty(modifiedModel.subSystems{rxnIndex})
                if contains(modifiedModel.subSystems{rxnIndex}, 'Transport', 'IgnoreCase', true)
                    % Logic based on metabolite compartment
                    if (endsWith(metaboliteID, '[e]') || contains(metaboliteID, '_e')) && modifiedModel.S(metaboliteIndex, rxnIndex) < 0
                        % Metabolite in [e] compartment and negative coefficient
                        modifiedModel.ub(rxnIndex) = 1000; % Allow forward direction
                        modifiedModel.lb(rxnIndex) = 0;   % Block reverse reaction
                        fprintf('  Set transport reaction %s (forward, external uptake)\n', rxnID);
                    elseif (endsWith(metaboliteID, '[c]') || contains(metaboliteID, '_c')) && modifiedModel.S(metaboliteIndex, rxnIndex) > 0
                        % Metabolite in [c] compartment and positive coefficient
                        modifiedModel.ub(rxnIndex) = 1000; % Allow forward direction
                        modifiedModel.lb(rxnIndex) = 0;   % Block reverse reaction
                        fprintf('  Set transport reaction %s (forward, internal production)\n', rxnID);
                    else
                        % Block reactions that don't match the criteria
                        modifiedModel.ub(rxnIndex) = 0;
                        modifiedModel.lb(rxnIndex) = 0;
                        fprintf('  Blocked transport reaction %s (not matching criteria)\n', rxnID);
                    end
                elseif contains(modifiedModel.subSystems{rxnIndex}, 'Exchange', 'IgnoreCase', true)
                    % Logic based on metabolite compartment
                    if (endsWith(metaboliteID, '[e]') || contains(metaboliteID, '_e')) && modifiedModel.S(metaboliteIndex, rxnIndex) > 0
                        % Metabolite in [e] compartment and positive coefficient
                        modifiedModel.ub(rxnIndex) = 1000; % Allow forward direction
                        modifiedModel.lb(rxnIndex) = 0;   % Block reverse reaction
                        fprintf('  Set exchange reaction %s (forward)\n', rxnID);
                    else
                        % Block reactions that don't match the criteria
                        modifiedModel.ub(rxnIndex) = 0;
                        modifiedModel.lb(rxnIndex) = 0;
                        fprintf('  Blocked exchange reaction %s (not matching criteria)\n', rxnID);
                    end
                end
            end
        end
    end

    fprintf('Modified transport reactions for the specified metabolites based on their coefficients and compartments.\n');
end
