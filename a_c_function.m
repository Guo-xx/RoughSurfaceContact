function [ac_prime] = a_c_function(D,E,H,G,gama)
%UNTITLED ���ݴֲڶȺͲ������Լ���ac_prime
%   �˴���ʾ��ϸ˵��
ac_prime = (2^(11-2*D)/(9*pi^(4-D))*G^(2*D-4)*(E/H)^2*log(gama))^(1/(D-2));
end

