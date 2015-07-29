%generate a ricean distributed flag
%s: power in two dimensions
%k: factor linear
%normalize: 0 don't normalize 1=>normalize such that E{x*x'}=1;
%x resulting ricean vector.
%
%s can be a vector. in this case the x contains length(s) i.i.d RV
%if length(s)==length(k) then k is taken element per element
%length(k)==1 then k applies only to the first element
function x=gen_rice(s,k,normalize)
% %make row vectors and complete k to s by adding zeros
s=s(:)';
k=[k(:)' zeros(1,length(s)-length(k))];
n=length(s);
%variance (per one dimension) of scatter part
sigm2=s./(k+1)/2;



%power of ricean part
a2=s.*k./(k+1);


%generate vector
x=sqrt(a2).*exp(j*2*pi*rand(1,n))+ sqrt(sigm2).*([1 j]*randn(2,n));


%normalize if necessary
if (normalize),
%profile
prof=[a2+2*sigm2];
norm_fac= sqrt(prof*prof');
x=x/norm_fac;
end

