function [costsP, costsM, solution] = getReactionCosts(model)
    solution = get_solution_ec_model(model);

    % Check if the last values of model.rxns is 'Memb_const'
    if strcmp(model.rxns{end}, 'Memb_Const')
        coeffsP = model.S(end-1, :);
        coeffsM = model.S(end, :);
    % Check if the last values of model.rxns is 'Prot_const'
    elseif strcmp(model.rxns{end}, 'Prot_Const')
        coeffsP = model.S(end, :);    
        coeffsM = zeros(1, length(model.rxns));
    else
        error('Unexpected format for model.rxns. Unable to determine coefficients.');
    end

    % Calculate the cost of each reaction
    costsP = abs(solution.x) .* coeffsP';
    costsM = abs(solution.x) .* coeffsM';
    
end