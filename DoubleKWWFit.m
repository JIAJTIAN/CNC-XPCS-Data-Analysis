function result = DoubleKWWFit(TAU, G2)

% Set up fittype and options.
ft = fittype( 'b0*(b1*exp(-(2.*a1*x)^a2) + (1-b1)*exp(-(2.*a3*x)^a4)) + 1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0];
% opts.StartPoint = [0.5 0.5 0.5 0.5 0.5 0.5];

opts.StartPoint = [0.913375856139019 0.63235924622541 0.0975404049994095 0.278498218867048 0.546881519204984 0.957506835434298];
opts.Upper = [Inf Inf Inf Inf Inf 1];
% Fit model to data.

[G2Fitresult, gof, output] = fit( TAU, G2, ft, opts );
G2Fit = feval(G2Fitresult,TAU);
CoefArray = coeffvalues(G2Fitresult);
Rsd = output.residuals;

% end

result = {G2Fit,CoefArray,Rsd};

end


