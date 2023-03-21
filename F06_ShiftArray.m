function [Ashift] = F06_ShiftArray(A,n)
%F06_SHIFTARRAY Shift an array and append zeros

n=round(n);
for ctr=1:length(n)
    nT = n(ctr);
if nT>=0 %
    nT = min(nT,length(A));
    Ashift(ctr,:) = [A(nT+1:end) zeros(1,nT)];    
else
    nT=round(nT);
    nT = max(nT,-length(A));
    Ashift(ctr,:) = [zeros(1,-nT) A(1:end+nT)];
end

end

