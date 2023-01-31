% Calculation of information measures of bivariate model

function Results=Bivariate_aux_fun(Complexity,nsim,nscales,index_proc)
for n=1:nsim
        for is=1:nscales
            Results.MSEall(is,n)=Complexity{is,n}.glob;
            Results.MSE_j_j(is,n)=Complexity{is,n}.ind(index_proc(1));
            Results.MSE_j_ji(is,n)=Complexity{is,n}.ind_ij(index_proc(1),index_proc(2));
        end
end

Results.Ti_j=Results.MSE_j_j-Results.MSE_j_ji;
end