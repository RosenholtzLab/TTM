function [mu sigmasq skw krt] = marginalStats_forTTM(x)
    mu = mean2(x);
    sigmasq = var_forTTM(x);
    skw = skew2_forTTM(x, mu, sigmasq);
    krt = kurt2_forTTM(x, mu, sigmasq);
return