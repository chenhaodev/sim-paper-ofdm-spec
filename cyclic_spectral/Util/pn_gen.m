function [pn] = pn_gen(N1, N2)
% 1. function [pn] = pn_gen(N):	generate 1*N pseudo random sequence from {-1, +1}
% 2. function [pn] = pn_gen(N1, N2): generate N1*N2 pseudo random matrix from {-1, +1}
% example: pn_gen(10) => [-1 1 -1 1 1 ...]

if nargin == 1
	seq_valid = 1; 
	mtx_valid = 0;
elseif nargin == 2
	seq_valid = 0; 
	mtx_valid = 1;
else
	seq_valid = 0; 
	mtx_valid = 0;
	disp('parameter error');
end

if seq_valid == 1 && mtx_valid == 0 
	pn = zeros(1, N1); 
	temp = randi([0 1], 1, N1);
	index = find(temp == 0);
	temp(index) = -1;
	pn = temp;
elseif seq_valid == 0 && mtx_valid == 1 
	pn = zeros(N1, N2); 
	temp = randi([0 1], N1, N2);
	index = find(temp == 0);
	temp(index) = -1;
	pn = temp;
else
	disp('check your input parameter');
end
