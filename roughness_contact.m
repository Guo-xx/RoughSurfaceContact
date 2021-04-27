%%
%{
W-M�����˹�ģ��ֲڱ��棬�����нӴ�����غɷ������ο��������£�
[1]YAN W, KOMVOPOULOS K. Contact analysis of elastic-plastic fractal surfaces[J].
 Journal of Applied Physics, 1998, 84(7): 3617�C3624
%}
clc;
clear;
%% �����˹�ģ��ֲڱ���Ĳ���������3D�ֲڱ���
% ��umΪ��λ������ֵ��С����λ��
M=10;
G=1.36*10^-5;      % 1.36*10^-5um
D=2.4;
L=1;               % 1um
gamma=1.5;
[x,y,Z] = fractal_surf_3d(M,G,D,L,gamma);
%% ���ݱ���ֲڶȺͲ������Լ����ٽ����a_c_prime
% <<a_c_prime.jpg>>
% a_c_prime�������λ��G�ĵ�λ����һ��
G1=G;     %����G��ȡֵ��λΪum,ʹ��a_c_prime�ĵ�λ��um^2
E=72*10^9*10^-12;      %��λ��N/um^2
mu=0.17;
E_star = (2*(1-mu^2)/E)^-1;
H=5.5*10^9*10^-12;     %��λ��N/um^2
a_c_prime = a_c_function(D,E_star,H,G1,gamma);         %a_c_prime���ٽ�Ӵ��������λ��um^2
%% ���ݽӴ�����֮��ľ���ͳ�Ƶó��Ӵ�����������СS_prime��a_L_prime������������ܽӴ���F_t
% ����������ȡֵ��Χ����nmΪ��λ
mssd=-50:0.0001:50; %��λ��nm
S_prime=mssd;
a_L_prime=mssd;
F_e=mssd;
F_p=mssd;
S_e=mssd;
S_p=mssd;
S_base=L*L;     %��λ��um^2
[width,height]=size(Z);
for ii=1:length(mssd)
    % �õ��ڴ�1-50�ĸ����������£�ͳ�ƽӴ�����������洢��S_prime��
    S_prime(ii)=sum(sum(Z>mssd(ii)))/(width*height)*S_base; 
    a_L_prime(ii)=S_prime(ii)*(3-D)/(D-1);
    % ����ͳ�ƽӴ������S_prime��a_L_prime����Ӵ�������ʵ�Ӵ�����е�
    % ���Գɷ�F_e��S_e�����Գɷ�F_p��S_p
    if a_L_prime(ii)<=a_c_prime
        F_e(ii)=0;
        F_p(ii)=H*S_prime(ii);   %��λ��N
        S_e(ii)=0;
        S_p(ii)=S_prime(ii);            %��λ��um^2
    else
        F_e(ii)=2^((11-2*D)/2)/(3*pi^((4-D)/2))*(D-1)/(5-2*D)*(log(gamma))^(1/2)*G^(D-2)...
                *E_star*a_L_prime(ii)^((4-D)/2)*(1-(a_c_prime/a_L_prime(ii))^((5-2*D)/2));
        F_p(ii)=((D-1)/(3-D)) * H * a_L_prime(ii) * (a_c_prime/a_L_prime(ii))^((3-D)/2);
        S_e(ii)=(D-1)/(6-2*D)*(1-(a_c_prime/a_L_prime(ii))^((3-D)/2))*a_L_prime(ii);
        S_p(ii)=(D-1)/(3-D)*(a_c_prime/a_L_prime(ii))^((3-D)/2)*a_L_prime(ii);
    end
end
F_t=F_e+F_p;    %��λ��N
S_t=S_e+S_p;    %��λ��um^2
P_t=F_t./S_base;   %��λ��N/um^2
P_t=P_t*10^6;   %��λ��MPa(N/mm^2)
P_e=F_e./S_base*10^6;
P_p=F_p*10^6;
% P_t1=cumsum(P_t,'reverse'); %�ⲿ�����ݴ���ǲ���Ҫ�ۻ��İɣ�֮ǰ���뷨����Ǵ���ģ�210426
% P_e1=cumsum(P_e,'reverse');
% P_p1=cumsum(P_p,'reverse');
%% ����������
% ���ƽӴ���ĽӴ��غ�P_t(MPa)����ڱ�����mssd(nm)�ı仯����
figure('name','�Ӵ��غ�P_t(MPa)����ڱ�����mssd(nm)�ı仯����')
semilogy(mssd,P_t);
hold on;
semilogy(mssd,P_e);
hold on;
semilogy(mssd,P_p)
xlabel('Mean surface separation distance(nm)')
ylabel('Average contact pressure(MPa)')
title('ƽ���Ӵ��غ�P_t(MPa)����ڱ�����mssd(nm)�ı仯����')
legend 'Total' 'Elastic component' 'Plastic component'
% xlim([20,50]);

% ���S_prime(um^2)����ڱ�����mssd(nm)�ı仯����
figure('name','��ʵ�Ӵ����S_prime(um^2)����ڱ�����mssd(nm)�ı仯����')
plot(mssd,S_t/S_base);
hold on;
plot(mssd,S_e/S_base);
hold on;
plot(mssd,S_p/S_base);
xlabel('Mean surface separation distance(nm)')
ylabel('Real/Apparent contact areas')
title({'��ʵ�Ӵ����S�����۽Ӵ����S_ prime�ı�ֵ';'����ڱ�����mssd(nm)�ı仯����'})
legend 'Total' 'Elastic component' 'Plastic component'
% xlim([20,50]);
%% 















