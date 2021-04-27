%%
%{
W-M函数人工模拟粗糙表面，并进行接触面的载荷分析。参考文章如下：
[1]YAN W, KOMVOPOULOS K. Contact analysis of elastic-plastic fractal surfaces[J].
 Journal of Applied Physics, 1998, 84(7): 3617–3624
%}
clc;
clear;
%% 给定人工模拟粗糙表面的参数，生成3D粗糙表面
% 以um为单位定义数值的小数点位置
M=10;
G=1.36*10^-5;      % 1.36*10^-5um
D=2.4;
L=1;               % 1um
gamma=1.5;
[x,y,Z] = fractal_surf_3d(M,G,D,L,gamma);
%% 根据表面粗糙度和材料属性计算临界面积a_c_prime
% <<a_c_prime.jpg>>
% a_c_prime的面积单位与G的单位保持一致
G1=G;     %修正G的取值单位为um,使得a_c_prime的单位是um^2
E=72*10^9*10^-12;      %单位是N/um^2
mu=0.17;
E_star = (2*(1-mu^2)/E)^-1;
H=5.5*10^9*10^-12;     %单位是N/um^2
a_c_prime = a_c_function(D,E_star,H,G1,gamma);         %a_c_prime是临界接触面积，单位是um^2
%% 根据接触表面之间的距离统计得出接触表面的面积大小S_prime和a_L_prime，进而计算出总接触力F_t
% 定义表面间距的取值范围，以nm为单位
mssd=-50:0.0001:50; %单位是nm
S_prime=mssd;
a_L_prime=mssd;
F_e=mssd;
F_p=mssd;
S_e=mssd;
S_p=mssd;
S_base=L*L;     %单位是um^2
[width,height]=size(Z);
for ii=1:length(mssd)
    % 得到在从1-50的各个表面间距下，统计接触面面积，并存储在S_prime中
    S_prime(ii)=sum(sum(Z>mssd(ii)))/(width*height)*S_base; 
    a_L_prime(ii)=S_prime(ii)*(3-D)/(D-1);
    % 根据统计接触面面积S_prime和a_L_prime计算接触力和真实接触面积中的
    % 弹性成分F_e、S_e和塑性成分F_p、S_p
    if a_L_prime(ii)<=a_c_prime
        F_e(ii)=0;
        F_p(ii)=H*S_prime(ii);   %单位是N
        S_e(ii)=0;
        S_p(ii)=S_prime(ii);            %单位是um^2
    else
        F_e(ii)=2^((11-2*D)/2)/(3*pi^((4-D)/2))*(D-1)/(5-2*D)*(log(gamma))^(1/2)*G^(D-2)...
                *E_star*a_L_prime(ii)^((4-D)/2)*(1-(a_c_prime/a_L_prime(ii))^((5-2*D)/2));
        F_p(ii)=((D-1)/(3-D)) * H * a_L_prime(ii) * (a_c_prime/a_L_prime(ii))^((3-D)/2);
        S_e(ii)=(D-1)/(6-2*D)*(1-(a_c_prime/a_L_prime(ii))^((3-D)/2))*a_L_prime(ii);
        S_p(ii)=(D-1)/(3-D)*(a_c_prime/a_L_prime(ii))^((3-D)/2)*a_L_prime(ii);
    end
end
F_t=F_e+F_p;    %单位是N
S_t=S_e+S_p;    %单位是um^2
P_t=F_t./S_base;   %单位是N/um^2
P_t=P_t*10^6;   %单位是MPa(N/mm^2)
P_e=F_e./S_base*10^6;
P_p=F_p*10^6;
% P_t1=cumsum(P_t,'reverse'); %这部分内容大概是不需要累积的吧！之前的想法大概是错误的！210426
% P_e1=cumsum(P_e,'reverse');
% P_p1=cumsum(P_p,'reverse');
%% 绘制输出结果
% 绘制接触面的接触载荷P_t(MPa)相对于表面间距mssd(nm)的变化曲线
figure('name','接触载荷P_t(MPa)相对于表面间距mssd(nm)的变化曲线')
semilogy(mssd,P_t);
hold on;
semilogy(mssd,P_e);
hold on;
semilogy(mssd,P_p)
xlabel('Mean surface separation distance(nm)')
ylabel('Average contact pressure(MPa)')
title('平均接触载荷P_t(MPa)相对于表面间距mssd(nm)的变化曲线')
legend 'Total' 'Elastic component' 'Plastic component'
% xlim([20,50]);

% 输出S_prime(um^2)相对于表面间距mssd(nm)的变化曲线
figure('name','真实接触面积S_prime(um^2)相对于表面间距mssd(nm)的变化曲线')
plot(mssd,S_t/S_base);
hold on;
plot(mssd,S_e/S_base);
hold on;
plot(mssd,S_p/S_base);
xlabel('Mean surface separation distance(nm)')
ylabel('Real/Apparent contact areas')
title({'真实接触面积S与理论接触面积S_ prime的比值';'相对于表面间距mssd(nm)的变化曲线'})
legend 'Total' 'Elastic component' 'Plastic component'
% xlim([20,50]);
%% 















