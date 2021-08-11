% KWW fit Function

function kww= KWWFunc(xdata,params)
bg = params(1);

gamma = params(2);

x = xdata;

kww = exp(-2 .* gamma .* xdata) +bg;

end