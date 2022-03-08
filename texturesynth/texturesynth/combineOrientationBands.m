function ch = combineOrientationBands(scale, Nor, rpyr0, pind0)
  bandNums = [1:Nor] + (scale-1)*Nor+1;  %ori bands only
  ind1 = pyrBandIndices(pind0, bandNums(1));
  indN = pyrBandIndices(pind0, bandNums(Nor));
  bandInds = [ind1(1):indN(length(indN))];
  
  % Make fake pyramid, containing dummy hi, ori, lo
  fakePind = [pind0(bandNums(1),:);pind0(bandNums(1):bandNums(Nor)+1,:)];
  fakePyr = [zeros(prod(fakePind(1,:)),1);...
	 rpyr0(bandInds); zeros(prod(fakePind(size(fakePind,1),:)),1);];
  ch = reconSFpyr(fakePyr, fakePind, [1]);     % recon ori bands only
