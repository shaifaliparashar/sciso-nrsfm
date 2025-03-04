function [eq130,eq121,eq112,eq103,eq120,eq111,eq102,eq110,eq101,eq100,eq230,eq221,eq212,eq203,eq220,eq211,eq202,eq210,eq201,eq200] = create_polynomial_coefficients(a,b,c,d,t1,t2,e,e1,u,u1,v,v1)

% metric tensor at image i written in terms of CS: g11, g12, g22 

g1120 = a.^2.*e;
g1111 = a.*b.*e.*2;
g1102 = b.^2.*e;
g1110 = (a.*e.*t1.*2) -(a.*u.*2);
g1101 = (b.*e.*t1.*2) -(b.*u.*2);
g1100 = (e.*t1.^2) -(t1.*u.*2) + 1;

g1220 = a.*c.*e;
g1211 = (a.*d.*e) + (b.*c.*e);
g1202 = b.*d.*e;
g1210 = -(a.*v)-(c.*u)+(a.*e.*t2)+(c.*e.*t1);
g1201 = -(b.*v)-(d.*u)+(b.*e.*t2)+(d.*e.*t1);
g1200 = -(t2.*u)-(t1.*v)+(e.*t1.*t2);

g2220 = c.^2.*e;
g2211 = c.*d.*e.*2;
g2202 = d.^2.*e;
g2210 = (c.*e.*t2.*2) -(c.*v.*2);
g2201 = (d.*e.*t2.*2) -(d.*v.*2);
g2200 = (e.*t2.^2) -(t2.*v.*2) + 1;

% pullback metric tensor at image i: g11b, g12b, g22b 

g11b20 = a.^2.*e1;
g11b11 = a.*b.*e1.*2;
g11b02 = b.^2.*e1;
g11b10 = -(a.^2.*u1.*2)-(a.*b.*v1.*2);
g11b01 = -(b.^2.*v1.*2)-(a.*b.*u1.*2);
g11b00 = a.^2+b.^2;

g12b20 = a.*c.*e1;
g12b11 = (a.*d.*e1) + (b.*c.*e1);
g12b02 = b.*d.*e1;
g12b10 = -(v1.*((a.*d) + (b.*c)))-(a.*c.*u1.*2);
g12b01 = -(u1.*((a.*d) + (b.*c)))-(b.*d.*v1.*2);
g12b00 = a.*c+b.*d;

g22b20 = c.^2.*e1;
g22b11 = c.*d.*e1.*2;
g22b02 = d.^2.*e1;
g22b10 = -(c.^2.*u1.*2)-(c.*d.*v1.*2);
g22b01 = -(d.^2.*v1.*2)-(c.*d.*u1.*2);
g22b00 = c.^2+d.^2;

%eq1 = g11*g22b - g11b*g22
eq130 = g1110.*g22b20 + g1120.*g22b10 - g11b10.*g2220 - g11b20.*g2210;
eq121 = g1101.*g22b20 + g1110.*g22b11 + g1111.*g22b10 + g1120.*g22b01 - g11b01.*g2220 - g11b10.*g2211 - g11b11.*g2210 - g11b20.*g2201;
eq112 = g1101.*g22b11 + g1102.*g22b10 + g1110.*g22b02 + g1111.*g22b01 - g11b01.*g2211 - g11b02.*g2210 - g11b10.*g2202 - g11b11.*g2201;
eq103 = g1101.*g22b02 + g1102.*g22b01 - g11b01.*g2202 - g11b02.*g2201;

eq120 = g1100.*g22b20 + g1110.*g22b10 + g1120.*g22b00 - g11b00.*g2220 - g11b10.*g2210 - g11b20.*g2200;
eq111 = g1100.*g22b11 + g1101.*g22b10 + g1110.*g22b01 + g1111.*g22b00 - g11b00.*g2211 - g11b01.*g2210 - g11b10.*g2201 - g11b11.*g2200;
eq102 = g1100.*g22b02 + g1101.*g22b01 + g1102.*g22b00 - g11b00.*g2202 - g11b01.*g2201 - g11b02.*g2200;

eq110 = g1100.*g22b10 + g1110.*g22b00 - g11b00.*g2210 - g11b10.*g2200;
eq101 = g1100.*g22b01 + g1101.*g22b00 - g11b00.*g2201 - g11b01.*g2200;

eq100 = g1100.*g22b00 - g11b00.*g2200;

%eq2 = g12*g22b - g12b*g22
eq230 = g1210.*g22b20 + g1220.*g22b10 - g12b10.*g2220 - g12b20.*g2210;
eq221 = g1201.*g22b20 + g1210.*g22b11 + g1211.*g22b10 + g1220.*g22b01 - g12b01.*g2220 - g12b10.*g2211 - g12b11.*g2210 - g12b20.*g2201;
eq212 = g1201.*g22b11 + g1202.*g22b10 + g1210.*g22b02 + g1211.*g22b01 - g12b01.*g2211 - g12b02.*g2210 - g12b10.*g2202 - g12b11.*g2201;
eq203 = g1201.*g22b02 + g1202.*g22b01 - g12b01.*g2202 - g12b02.*g2201;

eq220 = g1200.*g22b20 + g1210.*g22b10 + g1220.*g22b00 - g12b00.*g2220 - g12b10.*g2210 - g12b20.*g2200;
eq211 = g1200.*g22b11 + g1201.*g22b10 + g1210.*g22b01 + g1211.*g22b00 - g12b00.*g2211 - g12b01.*g2210 - g12b10.*g2201 - g12b11.*g2200;
eq202 = g1200.*g22b02 + g1201.*g22b01 + g1202.*g22b00 - g12b00.*g2202 - g12b01.*g2201 - g12b02.*g2200;

eq210 = g1200.*g22b10 + g1210.*g22b00 - g12b00.*g2210 - g12b10.*g2200;
eq201 = g1200.*g22b01 + g1201.*g22b00 - g12b00.*g2201 - g12b01.*g2200;

eq200 = g1200.*g22b00 - g12b00.*g2200;



