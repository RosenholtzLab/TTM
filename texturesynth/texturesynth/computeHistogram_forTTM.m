function ht = computeHistogram_forTTM(X, bins)
[hn, hx] = histo(X, bins);
ht.n = hn;
ht.x = hx;
[mu sigmasq skw krt] = moments_forTTM(X);
ht.moments = [mu sigmasq skw krt];