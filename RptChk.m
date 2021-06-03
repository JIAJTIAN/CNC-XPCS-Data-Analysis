% Repetition Check 

function result = RptChk(TAU2,R1,R2,Diff)

for Index = 1:min(length(R1),length(R2))
    
    if abs(R1(Index) - R2(Index)) > R1(Index)*Diff % Equilibruim check
        R(Index) = R1(Index);
        TAU(Index) = TAU2(Index);
    else
        R(Index) = (R1(Index)+R2(Index))/2;
        TAU(Index) = TAU2(Index);
    end
    
end
result = [transpose(TAU), transpose(R)];

end