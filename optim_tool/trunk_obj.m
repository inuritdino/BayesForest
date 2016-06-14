function z = trunk_obj(X,data,model)
% Objective function for the trunk-data. Most commonly the trunk-data
% contains diameter at base, breast height diameter, and height.

[~,A] = model(X);
if(isempty(A))
    disp('Simulation had errors. Set max distance.');
    z = 10.0;
else
    % DiamBase
    tf = strcmp('DiamBase',A.colheaders);
    DiamBase = A.data(tf);
    % Breast height diameter
    tf = strcmp('BreastDiam',A.colheaders);
    BreastDiam = A.data(tf);
    % Height
    tf = strcmp('Height',A.colheaders);
    Height = A.data(tf);
    
    % Calculate score
    z1 = abs( (DiamBase - data.DiamBase)/data.DiamBase );
    z2 = abs( (BreastDiam - data.BreastDiam)/data.BreastDiam );
    z3 = abs( (Height - data.Height)/data.Height );
    
    %z = (z1 + z2 + z3) / 3;
    % Weighted sum
    z = 0.1*z1 + 0.1*z2 + 0.8*z3;
    
    % Plot
    figure(10);
    stem([z1 z2 z3]);
    title(['Dist = ' num2str(z)]);
end
end