function model = setConstraintWeights(model,group,sigma_parameter,varargin)
% function model = setConstraintWeights(model,group,myW,rescaleFlag)
% This function can be used to set the values of the weights of the
% different proteome sectors, instead of manually editing the 'w' field in
% the model.
%
% Inputs:
% 
%   model : a valid CAFBA model (including the protGroup field)
%   myW : value for the weights in the selected group (default: 0)
%   group: proteome group to be varied (integer number from 1 to 4)
%   rescaleFlag: this is useful in the case of randomly generated weights.
%          If 0 (default), the function sets all the weights in the group
%          to the value specified by myW. If 1, instead, the function
%          rescales all weights in the proteome sector proportionally,
%          so as to fix their average value to myW without affecting the
%          ratios of different weights, w_i/w_j for any i,j in the group.
%
% Output:
%
%   model: a valid CAFBA model (including the protGroup and w fields).
%
% See also the setWeights_rand() function.

% Check the inputs
if ~exist('model','var')
    error('No model was provided. Exiting from setWeights().');
end
if ~isfield(model,'protGroup') % extra check on the model
    error('protGroup field missing. Exiting from setWeights().');
end
if ~isfield(model,'w_c') % extra check on the model
    model.w_c =zeros(length(model.rxns),1);
end
if ~isfield(model,'w_m') % extra check on the model
    model.w_m=zeros(length(model.rxns),1);
end
if ~isfield(model,'sa_m') % extra check on the model
    model.sa_m=zeros(length(model.rxns),1);
end

optargin=size(varargin,2);
if floor(optargin/2)~=optargin/2
    error('Error, optargin is not an even number.');
elseif optargin>0
    for i=1:2:(optargin-1)
        switch varargin{i}
            % mode,pMTP,rescaleFlag, mbw_mode
            case 'mode'
                mode=varargin{i+1};
            case 'pMTP'
                pMTP=varargin{i+1};
            case 'rescaleFlag'
                rescaleFlag=varargin{i+1};
            case 'mbw_mode'
                mbw_mode=varargin{i+1};
            case 'testPrint'
                testPrint=varargin{i+1};
            case 'numIssue_t'
                numIssue_t=varargin{i+1};
            case 'shape'
                shape = varargin{i+1};
            otherwise
                error('Optional argument not recognized.');
        end
    end
end

% Defaults
if ~exist('sigma_parameter','var') || isempty(sigma_parameter)
    sigma_parameter=0;
end
if ~exist('group','var') || isempty(group)
    group=0; % all reactions
end
if ~exist('rescaleFlag','var')
    rescaleFlag=0;
end
if ~exist('mode','var')
    mode='homogenous';
end
if ~exist('mbw_mode','var')
    mbw_mode='total';
end
if ~exist('pMTP','var')
    pMTP=0.32; % g_prot/g_dcw
end
if ~exist('testPrint','var')
    testPrint=0; 
end
if ~exist('numIssue_t','var')
    numIssue_t = 1e-6; 
end
if ~exist('shape','var')
    shape = 1 ;
end
if strcmp(mode,'homogenous')
    model=mode_homogenous(model,group,sigma_parameter,rescaleFlag,numIssue_t);
    return
end
if strcmp(mode,'heterogeneous')
    model=mode_heterogeneous(model,group,sigma_parameter,rescaleFlag,numIssue_t,shape);
    return
end

if strcmp(mode,'parametrized')
    model=mode_parametrized(model,group,sigma_parameter,pMTP,rescaleFlag,mbw_mode, testPrint, numIssue_t);
%     model=mode_parametrized(model,group,myW,TPC,rescaleFlag,mbw_mode,testPrint,numIssue_t)
    return
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION- Set all the weights according to the %
% homogeneous approach.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model=mode_homogenous(model,group,sigma_parameter,rescaleFlag,numIssue_t)

% All rxns in group
if group>0
    reactions=model.protGroup(group).rxns; 
elseif group==0
    reactions = 1:length(model.rxns); % All the reactions
else
    error('group>length(model.protGroup) or group<0 in setWeights()');
end

% Just the group
if rescaleFlag==0
    % All weights are uniformly set to 'myW'
    for ri=1:length(reactions)
        r = reactions(ri);
%         disp(r);
        model.w_c(r)=sigma_parameter;
        if strcmp(model.location{r}, 'Inner Membrane')
            model.w_m(r)=sigma_parameter;
        end
    end
elseif rescaleFlag==1 
    % All weights in the group are rescaled such that
    % the new average is 'myW'
    average=0;
    for ri=1:length(reactions)
        r=reactions(ri);
        average=average+model.w_c (r);
    end
    average=average/length(reactions);
%     disp(average);
    for ri=1:length(reactions)
        r=reactions(ri);
        model.w_c(r)=(sigma_parameter/average)*model.w_c(r);
        model.w_m(r)=(sigma_parameter/average)*model.w_m(r);
    end       
else
    % All weights in the group are rescaled by a factor 'myW'
    for ri=1:length(reactions)            
        model.w_c(ri)=sigma_parameter*model.w_c(ri);
        if strcmp(model.location{ri}, 'Inner Membrane')
            model.w_m(ri)=sigma_parameter*model.w_m(ri);
        end
    end
end

% Check for potential numeric issues
% Find the indices of values that exceed the threshold
indices = abs(model.w_c) < numIssue_t;
% Update the values that exceed the threshold
model.w_c (indices) = 0; % For example, setting them to zero
% Find the indices of values that exceed the threshold
indices = abs(model.w_m) < numIssue_t;
% Update the values that exceed the threshold
model.w_m(indices) = 0; % For example, setting them to zero

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION- Set all the weights according to the %
% heterogeneous approach.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model=mode_heterogeneous(model,group,sigma_parameter,rescaleFlag,numIssue_t, shape)

% All rxns in group
if group>0
    reactions=model.protGroup(group).rxns; 
elseif group==0
    reactions = 1:length(model.rxns); % All the reactions
else
    error('group>length(model.protGroup) or group<0 in setWeights()');
end

% Just the group
if rescaleFlag==0
    % Option 0: extract the weights from a x=10^unif. pdf, that is, p(x) ~ 1/x
    Nrr = length(reactions);
    valuesRand = sigma_parameter*generateExpBoxRand(Nrr,shape);
    
    for ri=1:length(reactions)
        r = reactions(ri);
        val = valuesRand(ri);
        model.w_c(r)= val;
        if strcmp(model.location{r}, 'Inner Membrane')
            model.w_m(r)=val;
        end
    end
elseif rescaleFlag==1 
    % Option 1: extract the weights from a gamma distribution
    Nrr = length(reactions);
    valuesRand = sigma_parameter*generateGammaRnd(Nrr,shape);
    
    for ri=1:length(reactions)
        r = reactions(ri);
        val = valuesRand(ri);
        model.w_c(r)= val;
        if strcmp(model.location{r}, 'Inner Membrane')
            model.w_m(r)=val;
        end
    end
    
else
    % Option 3: extract the weights from a lognormal distribution
    Nrr = length(reactions);
    valuesRand = sigma_parameter*generateLogNormalRand(Nrr,shape);
    for ri=1:length(reactions)
        r = reactions(ri);
        val = valuesRand(ri);
        model.w_c(r)= val;
        if strcmp(model.location{r}, 'Inner Membrane')
            model.w_m(r)=val;
        end
    end
end

% Check for potential numeric issues
% Find the indices of values that exceed the threshold
indices = abs(model.w_c) < numIssue_t;
% Update the values that exceed the threshold
model.w_c (indices) = 0; % For example, setting them to zero
% Find the indices of values that exceed the threshold
indices = abs(model.w_m) < numIssue_t;
% Update the values that exceed the threshold
model.w_m(indices) = 0; % For example, setting them to zero

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION- Set all the weights according to the %
% fine grained approach.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model=mode_parametrized(model,group,sigma_parameter,TPC,rescaleFlag,mbw_mode,testPrint,numIssue_t)
    if testPrint
        fprintf('Rescaling the w cost by a factor of %1.4f \n', sigma_parameter);
    end
    
    % All rxns in group
    if group > 0
        reactions = model.protGroup(group).rxns; 
    elseif group == 0
        reactions = 1:length(model.rxns); % All the reactions
    else
        error('group>length(model.protGroup) or group<0 in setConstraintWeights()');
    end

    % Calculate weights based on different scenarios
    for ri = 1:length(reactions)
        r = reactions(ri);
        mw_c = model.mw_c(r) / 1000;
        mw_m = model.mw_m(r) / 1000;
        sa_m = model.sa_m(r);
        kcat = model.kcat(r) * 3600;
        model.w_c(r) = (mw_c / TPC) / kcat;
        if strcmp(mbw_mode, 'transmembrane')
            model.w_m(r) = (mw_m / TPC) / kcat;
        elseif strcmp(mbw_mode, 'total')
            indexes_PMb = mw_m > 0;
            mw_c_filtered = mw_c * indexes_PMb;
            model.w_m(r) = (mw_c_filtered / TPC) / kcat;
        elseif strcmp(mbw_mode, 'MSA')
            kcat = model.kcat(r);
            model.w_m(r) = sa_m / kcat;
        end
    end

    % Rescale weights based on the rescaleFlag
    if rescaleFlag == 0
        % No rescaling needed, weights are calculated using only 'TPC'
    elseif rescaleFlag == 1 
        % Rescale weights to achieve the new average 'myW'
        average_c = mean(model.w_c(reactions));
        average_m = mean(model.w_m(reactions));
        for ri = 1:length(reactions)
            r = reactions(ri);
            model.w_c(r) = (sigma_parameter / average_c) * model.w_c(r);
            model.w_m(r) = (sigma_parameter / average_m) * model.w_m(r);
        end
    elseif rescaleFlag == 2 
        % Rescale weights by a factor 'sigma'
        model.w_c(reactions) = sigma_parameter * model.w_c(reactions);
        model.w_m(reactions) = sigma_parameter * model.w_m(reactions);
    else
        % Default: Rescale weights by a factor 'sigma'
        model.w_c(reactions) = sigma_parameter * model.w_c(reactions);
        model.w_m(reactions) = sigma_parameter * model.w_m(reactions);
    end
  
    % Check for potential numeric issues
    % Find the indices of values that exceed the threshold
    indices = abs(model.w_c) < numIssue_t;
    % Update the values that exceed the threshold
    model.w_c(indices) = 0; % For example, setting them to zero
    % Find the indices of values that exceed the threshold
    indices = abs(model.w_m) < numIssue_t;
    % Update the values that exceed the threshold
    model.w_m(indices) = 0; % For example, setting them to zero

end 

   function x=generateGammaRnd(N,k)
    %generateGammaRnd(N,k)
    % Generates N random numbers according to a Gamma distribution
    %  f(x|k,theta) ~ x^(k-1) e^(-x/theta) ~ gampdf(x,k,theta)
    %
    % The shape parameter k>0 is chosen by the user, while
    % the scale parameter b is set such that the expected value
    % of the random variabile is 1 for any a:
    %
    % E(x)=1 
    % V(x)= 1/k
    %
    % (tends to a gaussian for k->\infty)
    %
    if ~exist('k','var') || k <= 0
        k=1; % default value
    end
    theta=gamma(k)/gamma(k+1);

    x= gamrnd(k,theta,[N,1]);
end

function x=generateExpBoxRand(N,w)
    % Generates a random variable distributed as the exponential of an
    % uniform distribuited number. If p(x)=1/w for x in [a,a+w], and calling
    % y=c^x, we get the distribution:
    % p(y)=1/(ln(c)*w*y) in the interval [c^a,c^(a+w)].
    % The expected value is:
    %
    % <y>=c^a/(ln(c)*w)*(c^w-1)
    %
    % Choosing <y>=1, we find a in terms of c and w:
    %
    % a=ln(ln(c)*w/(c^w-1))/ln(c)
    %
    if w==0
        x=ones(N,1);
    else
        c=10; % If c=10, w is the number of decades spanned by y (delta)
        a=log(log(c)*w/(c^w-1))/log(c);
        x= c.^(a + w*rand(N,1));
    end
end

function x=generateLogNormalRand(N,sigma)
    % Extracts N numbers from a lognormal distribution

    if sigma==0
        x=ones(N,1);
    else
        x=zeros(N,1);
        for i=1:N
            x(i)=10^( normrnd(0,sigma));
        end
    end
end