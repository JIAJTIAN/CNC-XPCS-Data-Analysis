%% Segment Data Packaging
function result = SgPackaging(TauCellRp,G2CellRp)
TauCat = [];
G2Cat = [];

for Index = 1:length(TauCellRp) %adds all of the TauCellRp data cells together
    TauCat = cat(1,TauCellRp{Index},TauCat);
    G2Cat = cat(1,G2CellRp{Index},G2Cat);
end
[TAU,SortIndex] = sort(TauCat); 
G2 = G2Cat(SortIndex);

result = [TAU,G2];

end