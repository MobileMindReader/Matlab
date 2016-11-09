function y = phi(functions, weightParameters, X)
    if length(functions) ~= length(weightParameters)
        disp('Error');
        return
    end
    
    y = zeros(1,length(X));
    
    for i=1:length(functions)
        func = functions{i};
        y = y + weightParameters(i)*func(X);
    end
end