function [out] = feature_test(in, f, f_range, a, a_range) 
% extract the feature with boundary
% e.g. [out] = feature_test(Cyclic_Spectrum, f, [-0.25 +0.25], a, [-0.15 +0.15]) 

	f_mid = round(length(f)/2); 
	f_bios_l = round(f_mid*f_range(1));
	f_bios_r = round(f_mid*f_range(2));
	f_index_range = (f_mid + f_bios_l) : (f_mid + f_bios_r);
	a_mid = round(length(a)/2);
	a_bios_l = round(a_mid*a_range(1));
	a_bios_r = round(a_mid*a_range(2));
	a_index_range = (a_mid + a_bios_l) : (a_mid + a_bios_r);
	% feature extraction
	out = in(a_index_range, f_index_range);
end
