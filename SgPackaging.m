%% Segment Data Packaging
function result = SgPackaging(TauCellRp,G2CellRp)
% TauPack = TauCellRp{1};
% G2Pack = G2CellRp{1};
% 
% for Index = 2:length(TauCellRp)
%     TAU2 = TauCellRp{Index};
%     Data2 = G2CellRp{Index};
%     
%     Result = Stitching(TauPack,TAU2,G2Pack,Data2);
%     TauPack = Result{1};
%     G2Pack = Result{2};
% end

TauPack = [];
G2Pack = [];
for Index = 1:length(TauCellRp)
TauPack = cat(1,TauCellRp{Index},TauPack);
G2Pack = cat(1,G2CellRp{Index},G2Pack);

end
[TAU,SortIndex] = sort(TauPack);
G2 = G2Pack(SortIndex);

TAUUnique = unique(TAU);

for UniqIdx = 1:length(TAUUnique)
    G2Idx = find(TAU == TAUUnique(UniqIdx));
    G2Unique(UniqIdx) = mean(G2(G2Idx));
end


result = [TAUUnique,G2Unique'];

end