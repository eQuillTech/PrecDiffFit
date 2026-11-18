function y=getCovarianceInverseM(M)
%Returns the inverse of the covariance matrix
%A parameter with no effect on the data (zero column in M) is reported with
%infinite variance instead of breaking the inversion.
   MtM=M'*M;   
   y=pinv(MtM);
   degenerate=(diag(MtM)==0);
   y(degenerate,degenerate)=diag(Inf(sum(degenerate),1));
end
