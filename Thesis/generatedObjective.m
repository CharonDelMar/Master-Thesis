function [obj, grad] = generatedObjective(inputVariables, extraParams)
%generatedObjective Compute objective function value and gradient
%
%   OBJ = generatedObjective(INPUTVARIABLES, EXTRAPARAMS) computes the
%   objective value OBJ at the point INPUTVARIABLES, using the extra
%   parameters in EXTRAPARAMS.
%
%   [OBJ, GRAD] = generatedObjective(INPUTVARIABLES, EXTRAPARAMS)
%   additionally computes the objective gradient value GRAD at the current
%   point.
%
%   Auto-generated by prob2struct on 28-Sep-2023 14:27:04

%% Map solver-based variables to problem-based.
a = inputVariables(:);

%% Compute objective function.
eqnsQuad = zeros(80, 1);
HeqmvecCell = cell(80, 1);

for i = 0:79
    Heq = extraParams{3*i+1};
    feq = extraParams{3*i+2}(:);
    rhseq = extraParams{3*i+3}(:);
    Heqmvec = Heq * inputVariables(:);
    HeqmvecCell{i+1} = Heqmvec;
    eqnsQuad(i+1) = 0.5*dot(inputVariables(:), Heqmvec) + dot(feq, inputVariables(:)) + rhseq;
end

obj = eqnsQuad;

if nargout > 1
    %% Compute objective gradient.
    % To call the gradient code, notify the solver by setting the
    % SpecifyObjectiveGradient option to true.
    jacQuad = zeros(numel(inputVariables(:)), 80);
    for i = 0:79
        feq = extraParams{3*i+2}(:);
        Heqmvec = HeqmvecCell{i+1};
        jacQuad(:,i+1) = Heqmvec + feq;
    end
    
    grad = jacQuad;
end

end