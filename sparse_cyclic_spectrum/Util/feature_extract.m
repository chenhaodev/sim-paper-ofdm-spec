function [out] = feature_extract(in, f, f_range, a, a_range) 
% extract the feature with boundary
% e.g. [out] = feature_extract(Cyclic_Spectrum, f, [+0.2], a, [+0.2]) % 0 ~ 0.2f, 0 ~ 0.2a

	f_len = length(f);
	f_r = round(f_len*f_range(1));
	f_index_range = 1 : f_r;
	a_len = length(a);
	a_r = round(a_len*a_range(1));
	a_index_range = 1 : a_r;
	% feature extraction
	out = in(a_index_range, f_index_range);
end
