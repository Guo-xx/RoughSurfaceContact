function [x,y,Z] = fractal_surf_3d(M1,G1,D1,L1,gamma1)
    % 三维的粗糙表面分析法的计算机模拟，以微米为单位定义数值的小数点位置
    G = G1;        %高度缩放参数，也叫分形粗糙度,
    D = D1;        %分形维数2<D<3
    gamma = gamma1;     %轮廓空间频率γ>1
    L = L1;          %取样长度，=70um时，该取样长度下大致需要32min模拟生成一个粗糙表面
    Ls = 6*0.543*10^-3;     %截止长度，取为6个材料晶格常数，如材料硅的晶格常数为=0.543nm
    M = M1;         %曲面褶皱的重叠数
    phi_mn = 2*pi*rand;  %随机相位，[0,2*pi]
    n1 = 0;         %W-M函数的最低频率序列数，=0
    n_max = floor(log(L/Ls)/log(gamma)); %最高频率序列数

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
    shading interp; % 使3维图更加光滑的插值命令
%     x=0.01:0.001:L;
%     y=0.01:0.001:L;
%     Z=max(Z,-0.02);
%     Z=Z*10^3;
%     mesh(x,y,Z,'.')
    xlabel('x(um)')
    ylabel('y(um)')
    zlabel('z(nm)')
    title('3D粗糙表面分形法模拟')
end

