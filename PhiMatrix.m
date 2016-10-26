
function Phi = PhiMatrix(basisFunctions, X)
    Phi = ones(length(X),length(basisFunctions)+1);
    
    for i=1:numel(basisFunctions)
        baseFunc = basisFunctions{i};
        Phi(:,i+1) = baseFunc(X);
    end
end