
function Phi = PhiMatrix(functions, X)
    Phi = ones(length(X),length(functions));
    
    for i=1:numel(functions)
        func = functions{i};
        Phi(:,i) = func(X);
    end
        
%     Phi2 = ones(length(X),length(functions));
%     for i=2:numel(functions)
%         func = functions{i};
%         Phi2(:,i) = func(X);
%     end
end