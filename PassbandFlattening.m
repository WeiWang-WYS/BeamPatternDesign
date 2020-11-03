function   [x_vector SDR_PerIndicator] = PassbandFlattening(AoD_range, Fmatrix, Nt, AngleStep, TransitionBand)

SteeringVector=@(theta) 1/sqrt(Nt)*exp(theta*j*pi*(0:1:Nt-1)');

DicMatrix = [];
for angle =  AoD_range(1) + TransitionBand  : AngleStep : AoD_range(2) - TransitionBand 
    DicMatrix = [DicMatrix SteeringVector(angle)];
end

[RowNum BasisNum] = size(Fmatrix);
 
[DicRow DicCol] = size(DicMatrix);
 
cvx_begin quiet
    variable X(BasisNum,BasisNum) hermitian
    variable r   
    maximize r
    subject to
    for ii=1:1:DicCol
     	Gmatrix = (Fmatrix'*DicMatrix(:,ii)*DicMatrix(:,ii)'*Fmatrix);  
        real(trace(Gmatrix*X)) >= r;
    end
    trace(X) <= 1;
    X==hermitian_semidefinite(BasisNum);
cvx_end  
 

AT=chol(X);
 
val_min_max = 0;

for iRan = 1:1:10000 
    x_vector_temp = AT'*(1/sqrt(2)*randn(BasisNum,1)+1/sqrt(2)*randn(BasisNum,1)*j);
    x_vector_temp = x_vector_temp/norm(x_vector_temp,2);
   	for ii=1:1:DicCol
     	Gmatrix = (Fmatrix'*DicMatrix(:,ii)*DicMatrix(:,ii)'*Fmatrix);  
        val_temp(ii) =  real(x_vector_temp'*Gmatrix*x_vector_temp);
    end
    val = min(val_temp);
    if val>val_min_max
        val_min_max = val; 
        x_vector_min_max = x_vector_temp;
    end
end 
x_vector = x_vector_min_max;
val_max = val_min_max;

 
SDR_PerIndicator = [val_max r];