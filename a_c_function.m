function [ac_prime] = a_c_function(D,E,H,G,gama)
%UNTITLED 根据粗糙度和材料属性计算ac_prime
%   此处显示详细说明
ac_prime = (2^(11-2*D)/(9*pi^(4-D))*G^(2*D-4)*(E/H)^2*log(gama))^(1/(D-2));
end

