function result = KWWFit(TAU, G2)

% Set up fittype and options.
ft = fittype( 'b*exp(-(2.*a1*x)^a2) + 1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 1];
opts.Upper = [Inf Inf Inf];
opts.StartPoint = [0.5 1.5 0.9];

% StartPointMat = ...
%     [   0.2238    0.7513    1.0000
%     0.7431    0.3922    1.0000
%     0.2769    0.0462    1.0000
%     0.9502    0.0344    1.0000
%     0.1869    0.4898    1.0000
%     0.9593    0.5472    1.0000
%     0.2543    0.8143    1.0000
%     0.2511    0.6160    1.0000];
% for Idx = 1:length(StartPointMat)
%     opts.StartPoint = StartPointMat(Idx,:);
%     
%     
    % Fit model to data.
    Idx = 1;
    [G2Fitresult, gof,output] = fit( TAU, G2, ft, opts );
    G2Fit{Idx} = feval(G2Fitresult,TAU);
    CoefArray{Idx} = coeffvalues(G2Fitresult);
    Rsd{Idx} = output.residuals;
    Rsquare(Idx) = gof.rsquare;
% end
[~,minIdx] = min(Rsquare);

result = {G2Fit{minIdx},CoefArray{minIdx},Rsd(minIdx)};
    
end


