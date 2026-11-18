function y=getDampedPseudoInverseLeft(M,nu)
%Returns the damped left pseudoinverse
%Uses pinv so a parameter with no effect on the data (a zero column in M,
%e.g. the volume fraction of a single domain) gets a zero step instead of NaN.
   MtM=M'*M;   
   wMtM=MtM+nu*diag(diag(MtM));
   iwMtM=pinv(wMtM);
   Mt=M';
   y=iwMtM*Mt;
end
