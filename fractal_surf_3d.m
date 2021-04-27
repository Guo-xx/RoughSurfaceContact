function [x,y,Z] = fractal_surf_3d(M1,G1,D1,L1,gamma1)
    % ��ά�Ĵֲڱ���������ļ����ģ�⣬��΢��Ϊ��λ������ֵ��С����λ��
    G = G1;        %�߶����Ų�����Ҳ�з��δֲڶ�,
    D = D1;        %����ά��2<D<3
    gamma = gamma1;     %�����ռ�Ƶ�ʦ�>1
    L = L1;          %ȡ�����ȣ�=70umʱ����ȡ�������´�����Ҫ32minģ������һ���ֲڱ���
    Ls = 6*0.543*10^-3;     %��ֹ���ȣ�ȡΪ6�����Ͼ�����������Ϲ�ľ�����Ϊ=0.543nm
    M = M1;         %����������ص���
    phi_mn = 2*pi*rand;  %�����λ��[0,2*pi]
    n1 = 0;         %W-M���������Ƶ����������=0
    n_max = floor(log(L/Ls)/log(gamma)); %���Ƶ��������

    Z=zeros(10,10);
    index_x=1;
    index_y=1;
    step=0.0005;
    for x = 0.01:step:L
        for y = 0.01:step:L
            Z(index_x,index_y) = 0;
            for m = 1:1:M
                for n = n1:1:n_max
                    phi_mn = 2*pi*rand;
                    Z(index_x,index_y)= Z(index_x,index_y)+L*(G/L)^(D-2)*(log(gamma)/M)^0.5*gamma^((D-3)*n)*...
                        (cos(phi_mn)-cos((2*pi*gamma^n*(x^2+y^2)^0.5)*cos(atan(y/x)-pi*m/M)+phi_mn));
                end
            end
            index_y = index_y+1;
        end
        index_y=1;
        index_x = index_x+1;
    end
    x=0.01:step:L;
    [X, Y] = meshgrid(x);
    Z = Z-(max(max(Z))+min(min(Z)))/2;
    Z = Z*1000;
    T = delaunay(X,Y);
    trisurf(T,X,Y,Z)
    shading interp; % ʹ3άͼ���ӹ⻬�Ĳ�ֵ����
%     x=0.01:0.001:L;
%     y=0.01:0.001:L;
%     Z=max(Z,-0.02);
%     Z=Z*10^3;
%     mesh(x,y,Z,'.')
    xlabel('x(um)')
    ylabel('y(um)')
    zlabel('z(nm)')
    title('3D�ֲڱ�����η�ģ��')
end

