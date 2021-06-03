% Equilibrum Check for g2 and g2b
function result = EquilChk(TAU,G2,G2b,Diff,PassRate)
DiffIndex = abs(G2 - G2b) > G2.*Diff; % Equilibruim check
TotalNum = length(TAU);
PassNum = length(find(DiffIndex==1));

if PassNum/TotalNum<PassRate
    result = 1;
else
    result = 0;
end

end