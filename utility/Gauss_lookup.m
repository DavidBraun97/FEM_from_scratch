function [GP,GW] = Gauss_lookup()
% function gets Gauss points and weights from table
GP_table = [ -1/sqrt(3.); +1/sqrt(3)];
GW_table = [ 1; 1];

GP = GP_table;
GW = GW_table;
end

