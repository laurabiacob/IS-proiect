%% Incarcare date primite
t = scope18(:,1);
u = scope18(:,2);
y1 = scope18(:,3);
y2 = scope18(:,4); %sistem cu zero
%subplot(211)
%plot(t,[u,y1-3,y2-6]), grid, legend('u','y1','y2') %, title('y_1'), legend('u','y1')
%subplot(212)
plot(t,u,t,y1), title('y_2'), legend('u','y2')
%% Determinarea amplitudinilor si a amplificarii maxime
% pentru y1
i1 = 377; i2 = 396; %u
i3 = 383; i4 = 403; %y

Ain_y1 = (u(i2)-u(i1))/2;
Aout_y1 = (y1(i4)-y1(i3))/2;
Mr1 = Aout_y1/Ain_y1;

roots([7.654 0 -7.654 0 1]); %4*Mr^2 0 -4*Mr^2 1
z1 = 0.393;
Tr = t(420)-t(383);
fr = 1/Tr;
wr1 = 2*pi*fr;
wn = wr1/sqrt(1-2*(z1^2));
k1 = mean(y1)/mean(u);
Ph_r1 = (t(377)-t(384))*wr1;

H_y1 = tf(k1*wn^2,[1 2*z1*wn wn^2]);
A = [0, 1; -wn^2, -2*z1*wn];
B = [0; k1*wn^2];
C = [1, 0];
D = 0;
yc = lsim(A,B,C,D,u,t,[y1(1),(y1(2)-y1(1))/(t(2)-t(1))]);
plot(t, [yc y1]), legend('yc','y1')
Empn1 = norm(y1-yc)/norm(y1-mean(y1));
%% pentru y2
i5 = 377; i6 = 396; %u
i7 = 383; i8 = 400; %y

Ain_y2 = (u(i6)-u(i5))/2;
Aout_y2 = (y2(i8)-y2(i7))/2;
Mr2 = Aout_y2/Ain_y2;
roots([8.8011 0 -8.8011 0 1]);
%z2 = 0.4995;
z2 = 0.3915;
Tr2 = t(418)-t(383);
fr2 = 1/Tr2;
wr2 = 2*pi*fr;
wn2 = wr2/sqrt(1-2*z2^2);

Tzw = sqrt(Mr2^2*4*z2^2*(1-z2^2)-1);
Tz = Tzw/wr2;
Ph_r2 = (t(377)-t(381))*wr2;
k2 = mean(y2)/mean(u);

H_y2 = tf(k2*wn2^2*[Tz 1],[1 2*z2*wn2 wn2^2]);
A = [0, 1; -wn2^2, -2*z2*wn2];
B = [k2*wn2^2*Tz; k2*wn2^2-2*z2*wn2^3*k2*Tz];
C = [1 0];
D = 0;
yc1 = lsim(A,B,C,D,u,t,[y2(1),(y2(2)-y2(1))/(t(2)-t(1))-k2*wn2^2*Tz*u(1)]);
plot(t, [yc1 y2]), legend('yc','y2')
Empn2 = norm(y2-yc1)/norm(y2-mean(y2))
%% Estimarea Locului de Transfer
plot(t,u,t,y2)

u1 = 206; yo1 = 208; 
u2 = 239; yo2 = 240; 
u3 = 267; yo3 = 269; % w1, w2 si M1, M2

u5 = 317; yo5 = 320; 
u6 = 338; yo6 = 342; %w4 M4
 
%u9 = 632; yo9 = 639;
%u10 = 644; yo10 = 650;
%u11 = 654; yo11 = 660;

% u7 = 940; yo7 = 946;
% u8 = 948; yo8 = 952;

u9 = 736; yo9 = 742;
u10 = 746; yo10 = 752;
u11 = 755; yo11 = 761;

u12 = 507; yo12 = 514;
u13 = 521; yo13 = 527;
u14 = 534; yo14 = 540;

u7 = 675; yo7 = 681;
%u7 = 861; yo7 = 867;

u8 = 686; yo8 = 692;
%u8 = 869; yo8 = 875;

w1 = pi/(t(u2)-t(u1));
w2 = pi/(t(u3)-t(u2));
w4 = pi/(t(u6)-t(u5));
w6 = pi/(t(u10)-t(u9));
w7 = pi/(t(u13)-t(u12));
w5 = pi/(t(u8)-t(u7));

M1 = (y2(yo1)-y2(yo2))/(u(u1)-u(u2));
M2 = (y2(yo2)-y2(yo3))/(u(u2)-u(u3));
M4 = (y2(yo5)-y2(yo6))/(u(u5)-u(u6));
M5 = (y2(yo7)-y2(yo8))/(u(u7)-u(u8));
M6 = (y2(yo9)-y2(yo10))/(u(u9)-u(u10));
M7 = (y2(yo13)-y2(yo14))/(u(u13)-u(u14));

w = [w1 w4 wr2 w7 w5 w6];
M = [M1 M4 Mr2 M7 M5 M6];

semilogx(w,20*log10(M),'*'), grid
hold on
bode(H_y2), grid 

ph1 = (t(u1)-t(yo1))*w1;
ph2 = (t(u2)-t(yo2))*w2;
ph4 = (t(u6)-t(yo6))*w4;
ph5 = (t(u7)-t(yo7))*w5;
ph6 = (t(u9)-t(yo9))*w6;
ph7 = (t(u12)-t(yo12))*w7;

ph1d = rad2deg(ph1);
ph2d = rad2deg(ph2);
ph4d = rad2deg(ph4);
phrd = rad2deg(Ph_r2);
ph5d = rad2deg(ph5);
ph6d = rad2deg(ph6);
ph7d = rad2deg(ph7);

phd = [ph1d ph4d phrd ph7d ph5d ph6d];

semilogx(w,phd,'*'), grid
%%
subplot(211), semilogx(w,20*log10(M),'*-')
subplot(212), semilogx(w,phd,'*-')
%% Modele parametrice sistem y1
dt = t(2)-t(1);
dy1 = iddata(y1,u,dt);
%% Validare prin autocorelatie cu ARMAX
Marmax_y1 = armax(dy1,[2 1 2 1]);
Hz_y1 = tf(Marmax_y1.B,Marmax_y1.A,dt,'variable','z^-1');
Hc = d2c(Hz_y1);
[Num,Den] = tfdata(Marmax_y1,"v");
[A,B,C,D] = tf2ss(Num,Den); 
ya = dlsim(A',C',B',D,u);
plot(t,[y1 ya]), legend('ya','y1')
empn9 = norm(y1-ya)/norm(y1-mean(y1));

% figure
% resid(Marmax_y1,dy1,'corr',5)
% figure
% compare(dy1,Marmax_y1), shg
%% Validare prin intercorelatie cu OE
Moe_y1 = oe(dy1,[1 2 1]); %2 1 1
Hz_oe1 = tf(Moe_y1.B,Moe_y1.F,dt,'variable','z^-1');
Hc_oe1 = d2c(Hz_oe1);
[num,den] =tfdata(Moe_y1,'v')
[A,B,C,D] = tf2ss(num,den);
yo = dlsim(num,den,u);
plot(t,[y1 yo]), legend('y1','yo')
empn = norm(y1-yo)/norm(y1-mean(y1));
% resid(Moe_y1,dy1,5)
% figure
% compare(dy1,Moe_y1), shg
%% Modele parametrice sistem y2
dt = t(2)-t(1);
dy2 = iddata(y2,u,dt);
%% Validare prin autocorelatie cu ARMAX
Marmax_y2 = armax(dy2,[2 2 3 0])
Hz_y2 = tf(Marmax_y2.B,Marmax_y2.A,dt,'variable','z^-1');
Hc2 = d2c(Hz_y2);
% figure
% resid(Marmax_y2,dy2,'corr',5)
% figure
% compare(dy2,Marmax_y2), shg
%% Validare prin intercorelatie cu OE
Moe_y2 = oe(dy2,[2 2 0]);
Hz_oe2 = tf(Moe_y2.B,Moe_y2.F,dt,'variable','z^-1');
Hc_oe2 = d2c(Hz_oe2)
resid(Moe_y2,dy2,5)
figure
compare(dy2,Moe_y2), shg
%% MEP
d_id = iddata(y1,u,dt);
Mpem = pem(d_id);
%resid(Mpem,d_id)
figure
subplot(211)
compare(d_id,Mpem), shg
subplot(212)
compare(d_id,Marmax_y1), shg
%%
d_id2 = iddata(y2,u,dt);
Mpem = pem(d_id2);
%resid(Mpem,d_id)
figure
subplot(211)
compare(d_id2,Mpem), shg
subplot(212)
compare(d_id2,Marmax_y2), shg
