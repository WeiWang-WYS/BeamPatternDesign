function   [M_matrix w_vector] = HybridStructuralization(BeamformingVectorRef)

BeamformingVectorRef_normalized = 2*BeamformingVectorRef/max(abs(BeamformingVectorRef));
for ix = 1:1:length(BeamformingVectorRef_normalized)
    z = BeamformingVectorRef_normalized(ix);
    y_re_1 = 1/2*(real(z) + sqrt(-imag(z)^2 + 4*imag(z)^2/(imag(z)^2+real(z)^2))); 
	y_im_1 = 1/2*( imag(z)  - real(z)/imag(z)* sqrt(-imag(z)^2 + 4*imag(z)^2/(imag(z)^2+real(z)^2)) );
 
    y_re_2 = 1/2*(real(z) - sqrt(-imag(z)^2 + 4*imag(z)^2/(imag(z)^2+real(z)^2)));
	y_im_2 = 1/2*(imag(z) + real(z)/imag(z)* sqrt(-imag(z)^2 + 4*imag(z)^2/(imag(z)^2+real(z)^2)));
    M_matrix(ix,1) = y_re_1 + y_im_1*j; 
    M_matrix(ix,2) = y_re_2 + y_im_2*j;
end


w_vector = inv(M_matrix'*M_matrix)*M_matrix'*BeamformingVectorRef;
