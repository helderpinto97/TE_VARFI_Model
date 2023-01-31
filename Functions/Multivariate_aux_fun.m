% Calculation of information measures for a trivariate model
function Results=Multivariate_aux_fun(Complexity,nsim,nscales,index_proc)

for n=1:nsim
        for is=1:nscales
            Results.MSEall(is,n)=Complexity{is,n}.glob;
            Results.MSE_j_j(is,n)=Complexity{is,n}.ind(index_proc(1));
            Results.MSE_j_ji(is,n)=Complexity{is,n}.ind_ij(index_proc(1),index_proc(2));
            Results.MSE_j_jk(is,n)=Complexity{is,n}.ind_ij(index_proc(1),index_proc(3));
            Results.MSE_j_jki(is,n)=Complexity{is,n}.ind_all(index_proc(1));
            Results.Ti_j(is,n)=Results.MSE_j_j(is,n)-Results.MSE_j_ji(is,n);
            Results.Tl_j(is,n)=Results.MSE_j_j(is,n)-Results.MSE_j_jk(is,n);
            Results.Til_j(is,n)=Results.MSE_j_j(is,n)-Results.MSE_j_jki(is,n);
            Results.Iil_j(is,n)=Results.Til_j(is,n)-Results.Ti_j(is,n)-Results.Tl_j(is,n);
           
            % Calculation of PID Measures (Redundancy, Sinergy, Unique Information)
            Results.Ril_j(is,n)=min(Results.Ti_j(is,n),Results.Tl_j(is,n));
            Results.Ui_j(is,n)=Results.Ti_j(is,n)-Results.Ril_j(is,n); % Unique transfer
            Results.Ul_j(is,n)=Results.Tl_j(is,n)-Results.Ril_j(is,n); % Unique transfer
            Results.Sil_j(is,n)=Results.Til_j(is,n)-Results.Ui_j(is,n)-Results.Ul_j(is,n)-Results.Ril_j(is,n);
        end
end

end
