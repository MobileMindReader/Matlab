function y = phi(basisFunctions, weightParameters, X)
    if length(basisFunctions) ~= length(weightParameters)-1
        disp('Error');
        return
    end
    y = weightParameters(1)*ones(1,length(X));
    for i=1:length(basisFunctions)
        baseFunc = basisFunctions{i};
        y = y + weightParameters(i+1)*baseFunc(X);
    end
end