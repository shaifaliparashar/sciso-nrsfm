function ex = coeff_order(c,fc,thc,th)
or = ceil(log10(abs(c.*thc')));
or_max = max(or);
ex = fc(or>=(or_max-th))*c(or>=(or_max-th));