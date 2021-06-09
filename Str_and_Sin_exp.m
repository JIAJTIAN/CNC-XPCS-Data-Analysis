function result = StrExpFit(TAU, G2) % Kohlrausch model

% Set up fittype and options.
ft = fittype( 'exp(-2 * (a1*x)^a2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [1 0.6868];
opts.Upper = [Inf Inf];

% Fit model to data.
[G2Fitresult, gof,output] = fit( TAU, G2, ft, opts );
G2FitKWW = feval(G2Fitresult,TAU);
CoefArrayKWW = coeffvalues(G2Fitresult);
RsdKWW = output.residuals;


ftsingle = fittype('exp(-(2*a1*x))', 'independent', 'x', 'dependent', 'y');
 opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
 opts.StartPoint = [1];
opts.Upper = [Inf];

[G2Fitresult, gof,output] = fit( TAU, G2, ftsingle, opts);
G2FitSingle = feval(G2Fitresult, TAU);
CoefArraySingle = coeffvalues(G2Fitresult);
RsdSingle = output.residuals;

 result = {G2FitKWW, CoefArrayKWW, RsdKWW, G2FitSingle, CoefArraySingle, RsdSingle};
end


