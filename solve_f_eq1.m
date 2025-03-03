function [ff, res,de,sum_de] = solve_f_eq1(c9,c8,c7,c6,c5,c4,c3,c2,c1,c0,th,th2)

syms f real
sc = 1;
fc = [f^12,f^11,f^10,f^9,f^8,f^7,f^6,f^5,f^4,f^3,f^2,f,1];

thc = subs(fc,th^2);
tic
t = 0;
% create equations
for i = 1 : size(c9{1},2) % locally at each point
    i
    for k = 1: 1%size(c9{1},1)
        for j = 1: min(9,size(c9{1},1))
            if k ~= j
                c = [c9{1}(k,i);c9{2}(k,i);c9{3}(k,i);c9{4}(k,i);c9{5}(k,i);c9{6}(k,i);c9{7}(k,i);c9{8}(k,i);c9{9}(k,i);c9{10}(k,i);c9{11}(k,i);c9{12}(k,i);c9{13}(k,i)];
                e19 = coeff_order(c,fc,thc,th2);
                c = [c8{1}(k,i);c8{2}(k,i);c8{3}(k,i);c8{4}(k,i);c8{5}(k,i);c8{6}(k,i);c8{7}(k,i);c8{8}(k,i);c8{9}(k,i);c8{10}(k,i);c8{11}(k,i);c8{12}(k,i);c8{13}(k,i)];
                e18 = coeff_order(c,fc,thc,th2);
                c = [c7{1}(k,i);c7{2}(k,i);c7{3}(k,i);c7{4}(k,i);c7{5}(k,i);c7{6}(k,i);c7{7}(k,i);c7{8}(k,i);c7{9}(k,i);c7{10}(k,i);c7{11}(k,i);c7{12}(k,i);c7{13}(k,i)];
                e17 = coeff_order(c,fc,thc,th2);
                c = [c6{1}(k,i);c6{2}(k,i);c6{3}(k,i);c6{4}(k,i);c6{5}(k,i);c6{6}(k,i);c6{7}(k,i);c6{8}(k,i);c6{9}(k,i);c6{10}(k,i);c6{11}(k,i);c6{12}(k,i);c6{13}(k,i)];
                e16 = coeff_order(c,fc,thc,th2);
                c = [c5{1}(k,i);c5{2}(k,i);c5{3}(k,i);c5{4}(k,i);c5{5}(k,i);c5{6}(k,i);c5{7}(k,i);c5{8}(k,i);c5{9}(k,i);c5{10}(k,i);c5{11}(k,i);c5{12}(k,i);c5{13}(k,i)];
                e15 = coeff_order(c,fc,thc,th2);
                c = [c4{1}(k,i);c4{2}(k,i);c4{3}(k,i);c4{4}(k,i);c4{5}(k,i);c4{6}(k,i);c4{7}(k,i);c4{8}(k,i);c4{9}(k,i);c4{10}(k,i);c4{11}(k,i);c4{12}(k,i);c4{13}(k,i)];
                e14 = coeff_order(c,fc,thc,th2);
                c = [0;c3{1}(k,i);c3{2}(k,i);c3{3}(k,i);c3{4}(k,i);c3{5}(k,i);c3{6}(k,i);c3{7}(k,i);c3{8}(k,i);c3{9}(k,i);c3{10}(k,i);c3{11}(k,i);c3{12}(k,i)];
                e13 = coeff_order(c,fc,thc,th2);
                c = [0;c2{1}(k,i);c2{2}(k,i);c2{3}(k,i);c2{4}(k,i);c2{5}(k,i);c2{6}(k,i);c2{7}(k,i);c2{8}(k,i);c2{9}(k,i);c2{10}(k,i);c2{11}(k,i);c2{12}(k,i)];
                e12 = coeff_order(c,fc,thc,th2);
                c = [0;0;c1{1}(k,i);c1{2}(k,i);c1{3}(k,i);c1{4}(k,i);c1{5}(k,i);c1{6}(k,i);c1{7}(k,i);c1{8}(k,i);c1{9}(k,i);c1{10}(k,i);c1{11}(k,i)];
                e11 = coeff_order(c,fc,thc,th2);
                c = [0;0;c0{1}(k,i);c0{2}(k,i);c0{3}(k,i);c0{4}(k,i);c0{5}(k,i);c0{6}(k,i);c0{7}(k,i);c0{8}(k,i);c0{9}(k,i);c0{10}(k,i);c0{11}(k,i)];
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
                t= t+ 1;
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
                %         tic
                %         de{i} = det(res2);
                %         toc
            end
        end
    end
end
toc


% evaluate determinants
tic
len = length(res);
sum_de = 0;
parfor i = 1:len
   de{i} = det(res{i});
   sum_de = sum_de + de{i}^2;
   i
end
toc
delete(gcp('nocreate'))



% first and second order derivatives
df = diff(sum_de,f); df2 = diff(df,f);


% solve first order derivative
tic
ff = vpasolve(df==0,f);
ff = real(ff);
ff = ff(ff>=0);
ff = sqrt(ff);
ff = ff(ff> 0.6*th);
ff = ff(ff< th);
toc

% pick local minimas
tic
if ~isempty(ff)
for i = 1: length(ff)
    f_loc(i) = subs(df2,f,ff(i)^2);
end
ff = ff(f_loc>0);
end
% % compute residual
% for i = 1: length(ff)
%     %for j = 1:length(de)
%         f_res(i) = subs(sum_de,f,ff(i)^2);
%     %end
% end
% 
% % find solution that minimizes sum of deviation
% f_res2 = sum(f_res');
% [m idx]= sort(f_res2);
% 
% 
% 
% mx = max(max(f_res2));
% f_res2 = double(f_res2./mx);
% toc
% fres = median(f_res2');
% f_s = double(ff);
% f_th = median(f_res(:));
% for i= 30:len
%  ff = vpasolve(de{i}==0);
% ff = real(ff);
% ff = ff(ff>=0);
% f_s{i} = sqrt(ff);
% end