% example script
clear all;close all;

% add libraries
addpath(genpath('BBS'));
%addpath('SfTv0_3');
addpath(genpath('gloptipoly3'));
addpath(genpath('SeDuMi_1_3'));
addpath(genpath('schwarps'));
addpath(genpath('l1magic'));

% load dataset
 %load tshirt.mat
 load tshirt_data.mat
% K = scene.K;
% 
% idx = 1:10; % the images the algorithm needs to be evaluated on (in order)
% par = 2e-3; % schwarzian parameter.. needs to be tuned (usually its something close to 1e-3)
% [Pgth,Ngth,q_n,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_tshirt_dataset(idx,10, scene, par);
% 
% 
% % Christoffel Symbols (see equation 15 in the paper)
% %T1 = [-2*k1 -k2;-k2 0];
% %T2 = [0 -k1;-k1 -2*k2];
% % Christoffel Symbols change of variable  (see equation 10 in the paper) written in terms of k1b and k2b
% 
% % coeff of k1       % coeff of k2        % constant term
% T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
% T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;
% 
% % k1b = -T2_12 = a*k1 + b*k2 + t1;
% % k2b = -T1_12 = c*k1 + d*k2 + t2;
% a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;
% 
% % write metric tensor with k1b and k2b: g11, g12, g22
% u = I2u;
% v = I2v;
% 
% % Metric tensor (see equation 26 in the paper)
% % pullback metric tensor (see equation 8 in the paper): g11b,g12b,g22b
% u1 = repmat(I1u ,length(idx)-1,1);
% v1 = repmat(I1v ,length(idx)-1,1);


% % primary equations
% [a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,b214,b124,b034,b204,b114,...
%     b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000] = create_equations_view(a,b,c,d,t1,t2,u,u1,v,v1);
% 
% % resultant equations in terms of f and x1
% eq1_c9 = coeff_eq1_c9(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c8 = coeff_eq1_c8(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c7 = coeff_eq1_c7(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c6 = coeff_eq1_c6(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c5 = coeff_eq1_c5(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c4 = coeff_eq1_c4(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c3 = coeff_eq1_c3(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c2 = coeff_eq1_c2(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c1 = coeff_eq1_c1(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);
% 
% eq1_c0 = coeff_eq1_c0(a304,a214,a124,a034,a204,a114,a024,a302,a212,a122,a032,a202,a112,a022,a102,a012,a002,a300,a210,a120,a030,a200,a110,a020,a100,a010,a000,b304,...
%     b214,b124,b034,b204,b114,b024,b302,b212,b122,b032,b202,b112,b022,b102,b012,b002,b300,b210,b120,b030,b200,b110,b020,b100,b010,b000);

tic
res = solve_f_eq1(eq1_c9,eq1_c8,eq1_c7,eq1_c6,eq1_c5,eq1_c4,eq1_c3,eq1_c2,eq1_c1,eq1_c0,5e3,3);
toc

tic
len = length(res);
sum_de2 = 0;
parfor i = 1:len
   de{i} = det(res{i});
   de2{i}= de{i}^2;
   sum_de2 = sum_de2 + de2{i};
end
toc

tic
ff = vpasolve(sum_de2==0);
ff = real(ff);
ff = ff(ff>=0);
f_s = sqrt(ff);
toc

tic
[c9,c8,c7,c6,c5,c4,c3,c2,c1,c0,d9,d8,d7,d6,d5,d4,d3,d2,d1,d0] = create_polynomial_coefficients_uncalib_res(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
toc
tic
res = solve_polynomial(eq);
toc
% recover first order derivatives on rest of the surfaces
k1_all = [res(:,1)';a.*repmat(res(:,1)',length(idx)-1,1) + b.*repmat(res(:,2)',length(idx)-1,1) + t1];
k2_all = [res(:,2)';c.*repmat(res(:,1)',length(idx)-1,1) + d.*repmat(res(:,2)',length(idx)-1,1) + t2];

u_all = [I1u;I2u]; v_all = [I1v;I2v];

% find normals on all surfaces N= [N1;N2;N3]
N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
n = sqrt(N1.^2+N2.^2+N3.^2);
N1 = N1./n ; N2 = N2./n; N3 = N3./n;

N = [N1(:),N2(:),N3(:)]';
N_res = reshape(N(:),3*length(idx),length(u_all));

% find indices with no solution
idx = find(res(:,1)==0);
N_res(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = [];

% Integrate normals to find depth
P_grid=calculate_depth(N_res,u_all,v_all,1e0);


% compare with ground truth
[P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,qgth,Pgth);
[N,err_n] = compare_with_Ngth(P2,qgth,Ngth);

% plot results
for i=1:size(u_all,1)
     figure(i)
    plot3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),'go');
    hold on;
      plot3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),'ro');
      %quiver3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),N(3*(i-1)+1,:),N(3*(i-1)+2,:),N(3*(i-1)+3,:));
      %quiver3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),Ngth(3*(i-1)+1,:),Ngth(3*(i-1)+2,:),Ngth(3*(i-1)+3,:));
    hold off;
    axis equal;
end
mean(err_p')
mean(err_n')
