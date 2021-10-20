function mXout = subtractSNR(mX,snrmin, snrmax)

if snrmin==0 | snrmax==0
    warning('no surrounding DFT bin subtraction performed'); 
    mXout = mX; 
    return
end

baseline = zeros(size(mX));
for xi=1:length(mX)
    ind1 = max(1,xi-snrmax); 
    ind2 = max(1,xi-snrmin); 
    ind3 = min(length(mX),xi+snrmin); 
    ind4 = min(length(mX),xi+snrmax); 

    baseline(:,xi) = mean(mX(:,[ind1:ind2,ind3:ind4]), 2); 
end
mXout = mX-baseline; 

