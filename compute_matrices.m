function  [res, de, sum_de, ff, f_opt] = compute_matrices(c9,c8,c7,c6,c5,c4,c3,c2,c1,c0,th,th2)

syms f real

fc = [f^12,f^11,f^10,f^9,f^8,f^7,f^6,f^5,f^4,f^3,f^2,f,1];

thc = subs(fc,th^2);

sum_de = 0;
p1=size(c9{1},1);
p2 =size(c9{1},2);
% create equations

parfor t = 1 : p2*(p1-1) % locally at each point
    t
    [i, j]   = find_indices(t,p2,p1);
    c = [c9{1}(1,i);c9{2}(1,i);c9{3}(1,i);c9{4}(1,i);c9{5}(1,i);c9{6}(1,i);c9{7}(1,i);c9{8}(1,i);c9{9}(1,i);c9{10}(1,i);c9{11}(1,i);c9{12}(1,i);c9{13}(1,i)];
    e19 = coeff_order(c,fc,thc,th2);
    c = [c8{1}(1,i);c8{2}(1,i);c8{3}(1,i);c8{4}(1,i);c8{5}(1,i);c8{6}(1,i);c8{7}(1,i);c8{8}(1,i);c8{9}(1,i);c8{10}(1,i);c8{11}(1,i);c8{12}(1,i);c8{13}(1,i)];
    e18 = coeff_order(c,fc,thc,th2);
    c = [c7{1}(1,i);c7{2}(1,i);c7{3}(1,i);c7{4}(1,i);c7{5}(1,i);c7{6}(1,i);c7{7}(1,i);c7{8}(1,i);c7{9}(1,i);c7{10}(1,i);c7{11}(1,i);c7{12}(1,i);c7{13}(1,i)];
    e17 = coeff_order(c,fc,thc,th2);
    c = [c6{1}(1,i);c6{2}(1,i);c6{3}(1,i);c6{4}(1,i);c6{5}(1,i);c6{6}(1,i);c6{7}(1,i);c6{8}(1,i);c6{9}(1,i);c6{10}(1,i);c6{11}(1,i);c6{12}(1,i);c6{13}(1,i)];
    e16 = coeff_order(c,fc,thc,th2);
    c = [c5{1}(1,i);c5{2}(1,i);c5{3}(1,i);c5{4}(1,i);c5{5}(1,i);c5{6}(1,i);c5{7}(1,i);c5{8}(1,i);c5{9}(1,i);c5{10}(1,i);c5{11}(1,i);c5{12}(1,i);c5{13}(1,i)];
    e15 = coeff_order(c,fc,thc,th2);
    c = [c4{1}(1,i);c4{2}(1,i);c4{3}(1,i);c4{4}(1,i);c4{5}(1,i);c4{6}(1,i);c4{7}(1,i);c4{8}(1,i);c4{9}(1,i);c4{10}(1,i);c4{11}(1,i);c4{12}(1,i);c4{13}(1,i)];
    e14 = coeff_order(c,fc,thc,th2);
    c = [0;c3{1}(1,i);c3{2}(1,i);c3{3}(1,i);c3{4}(1,i);c3{5}(1,i);c3{6}(1,i);c3{7}(1,i);c3{8}(1,i);c3{9}(1,i);c3{10}(1,i);c3{11}(1,i);c3{12}(1,i)];
    e13 = coeff_order(c,fc,thc,th2);
    c = [0;c2{1}(1,i);c2{2}(1,i);c2{3}(1,i);c2{4}(1,i);c2{5}(1,i);c2{6}(1,i);c2{7}(1,i);c2{8}(1,i);c2{9}(1,i);c2{10}(1,i);c2{11}(1,i);c2{12}(1,i)];
    e12 = coeff_order(c,fc,thc,th2);
    c = [0;0;c1{1}(1,i);c1{2}(1,i);c1{3}(1,i);c1{4}(1,i);c1{5}(1,i);c1{6}(1,i);c1{7}(1,i);c1{8}(1,i);c1{9}(1,i);c1{10}(1,i);c1{11}(1,i)];
    e11 = coeff_order(c,fc,thc,th2);
    c = [0;0;c0{1}(1,i);c0{2}(1,i);c0{3}(1,i);c0{4}(1,i);c0{5}(1,i);c0{6}(1,i);c0{7}(1,i);c0{8}(1,i);c0{9}(1,i);c0{10}(1,i);c0{11}(1,i)];
    e10 = coeff_order(c,fc,thc,th2);
    
    c = [c9{1}(j,i);c9{2}(j,i);c9{3}(j,i);c9{4}(j,i);c9{5}(j,i);c9{6}(j,i);c9{7}(j,i);c9{8}(j,i);c9{9}(j,i);c9{10}(j,i);c9{11}(j,i);c9{12}(j,i);c9{13}(j,i)];
    e29 = coeff_order(c,fc,thc,th2);
    c = [c8{1}(j,i);c8{2}(j,i);c8{3}(j,i);c8{4}(j,i);c8{5}(j,i);c8{6}(j,i);c8{7}(j,i);c8{8}(j,i);c8{9}(j,i);c8{10}(j,i);c8{11}(j,i);c8{12}(j,i);c8{13}(j,i)];
    e28 = coeff_order(c,fc,thc,th2);
    c = [c7{1}(j,i);c7{2}(j,i);c7{3}(j,i);c7{4}(j,i);c7{5}(j,i);c7{6}(j,i);c7{7}(j,i);c7{8}(j,i);c7{9}(j,i);c7{10}(j,i);c7{11}(j,i);c7{12}(j,i);c7{13}(j,i)];
    e27 = coeff_order(c,fc,thc,th2);
    c = [c6{1}(j,i);c6{2}(j,i);c6{3}(j,i);c6{4}(j,i);c6{5}(j,i);c6{6}(j,i);c6{7}(j,i);c6{8}(j,i);c6{9}(j,i);c6{10}(j,i);c6{11}(j,i);c6{12}(j,i);c6{13}(j,i)];
    e26 = coeff_order(c,fc,thc,th2);
    c = [c5{1}(j,i);c5{2}(j,i);c5{3}(j,i);c5{4}(j,i);c5{5}(j,i);c5{6}(j,i);c5{7}(j,i);c5{8}(j,i);c5{9}(j,i);c5{10}(j,i);c5{11}(j,i);c5{12}(j,i);c5{13}(j,i)];
    e25 = coeff_order(c,fc,thc,th2);
    c = [c4{1}(j,i);c4{2}(j,i);c4{3}(j,i);c4{4}(j,i);c4{5}(j,i);c4{6}(j,i);c4{7}(j,i);c4{8}(j,i);c4{9}(j,i);c4{10}(j,i);c4{11}(j,i);c4{12}(j,i);c4{13}(j,i)];
    e24 = coeff_order(c,fc,thc,th2);
    c = [0;c3{1}(j,i);c3{2}(j,i);c3{3}(j,i);c3{4}(j,i);c3{5}(j,i);c3{6}(j,i);c3{7}(j,i);c3{8}(j,i);c3{9}(j,i);c3{10}(j,i);c3{11}(j,i);c3{12}(j,i)];
    e23 = coeff_order(c,fc,thc,th2);
    c = [0;c2{1}(j,i);c2{2}(j,i);c2{3}(j,i);c2{4}(j,i);c2{5}(j,i);c2{6}(j,i);c2{7}(j,i);c2{8}(j,i);c2{9}(j,i);c2{10}(j,i);c2{11}(j,i);c2{12}(j,i)];
    e22 = coeff_order(c,fc,thc,th2);
    c = [0;0;c1{1}(j,i);c1{2}(j,i);c1{3}(j,i);c1{4}(j,i);c1{5}(j,i);c1{6}(j,i);c1{7}(j,i);c1{8}(j,i);c1{9}(j,i);c1{10}(j,i);c1{11}(j,i)];
    e21 = coeff_order(c,fc,thc,th2);
    c = [0;0;c0{1}(j,i);c0{2}(j,i);c0{3}(j,i);c0{4}(j,i);c0{5}(j,i);c0{6}(j,i);c0{7}(j,i);c0{8}(j,i);c0{9}(j,i);c0{10}(j,i);c0{11}(j,i)];
    e20 = coeff_order(c,fc,thc,th2);
    res{t} =  [e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0,0,0,0,0,0,0,0;...
        0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0,0,0,0,0,0,0;...
        0,0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0,0,0,0,0,0;...
        0,0,0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0,0,0,0,0;...
        0,0,0,0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0,0,0,0;...
        0,0,0,0,0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0,0,0;...
        0,0,0,0,0,0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0,0;...
        0,0,0,0,0,0,0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10,0;...
        0,0,0,0,0,0,0,0,e19,e18,e17,e16,e15,e14,e13,e12,e11,e10;...
        e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0,0,0,0,0,0,0,0;...
        0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0,0,0,0,0,0,0;...
        0,0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0,0,0,0,0,0;...
        0,0,0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0,0,0,0,0;...
        0,0,0,0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0,0,0,0;...
        0,0,0,0,0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0,0,0;...
        0,0,0,0,0,0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0,0;...
        0,0,0,0,0,0,0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20,0;...
        0,0,0,0,0,0,0,0,e29,e28,e27,e26,e25,e24,e23,e22,e21,e20];
    
    de{t} = det(res{t});
    sum_de = sum_de + de{t}^2;
    
end

disp('calculate derivatives')
% first and second order derivatives
df = diff(sum_de,f); df2 = diff(df,f);

disp('calculate f')
% solve first order derivative
ff = vpasolve(df==0,f);
ff = real(ff);
ff = ff(ff>=0);
ff = sqrt(ff);
ff = ff(ff> 0.5*th);
ff = ff(ff< th);

disp('pick local minima')
% pick local minimas
if ~isempty(ff)
    for i = 1: length(ff)
        f_loc(i) = subs(df2,f,ff(i)^2);
    end
    ff = ff(f_loc>0);
end
f_opt = median(ff);
%save('f_estimate.mat','de','sum_de','ff','f','f_opt')


