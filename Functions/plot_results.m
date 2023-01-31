function [f1,f2,f3,l]=plot_results(Results,tau_vector)

% Create legend for figures
l=cell(size(Results.d_matrix,1),1);
for index=1:size(Results.d_matrix,1)
    l{index,1}=strcat('d=[', regexprep(num2str(Results.d_matrix(index,:)),'\s+',','),']');
end

% Figure to Complexiy Measures
f1=figure('Name','Complexity Measures','Position', get(0, 'Screensize'));
subplot(3,2,[1 2]);
plot(tau_vector,Results.MSEall);
title('Global Complexity');
legend(l,'Location','best');
subplot(3,2,3);
plot(tau_vector,Results.MSE_j_j);
title('C(j|{j})');
legend(l,'Location','best');
subplot(3,2,4);
plot(tau_vector,Results.MSE_j_ji)
title('C(j|{ji})');
legend(l,'Location','best');
subplot(3,2,5);
plot(tau_vector,Results.MSE_j_jk)
title('C(j|{jk})');
legend(l,'Location','best');
subplot(3,2,6);
plot(tau_vector,Results.MSE_j_jki)
title('C(j|{ijk})');
legend(l,'Location','best');


% Figure to Transfer Entropy
f2=figure('Name','Transfer Entropy','Position', get(0, 'Screensize'));
subplot(2,2,1);
plot(tau_vector,Results.Ti_j);
title('T_{i \rightarrow j}');
legend(l,'Location','best');
subplot(2,2,2);
plot(tau_vector,Results.Tl_j);
title('T_{k \rightarrow j}');
legend(l,'Location','best');
subplot(2,2,[3 4]);
plot(tau_vector,Results.Til_j)
title('T_{ik \rightarrow j}');
legend(l,'Location','best');
sgtitle('Transfer Entropy','FontWeight','bold')


% Figure to PID
f3=figure('Name','PID Decomposition','Position', get(0, 'Screensize'));
subplot(2,2,1);
plot(tau_vector,Results.Ui_j);
title('U_{i \rightarrow j}');
legend(l,'Location','best');
subplot(2,2,2);
plot(tau_vector,Results.Ul_j);
title('U_{k \rightarrow j}');
legend(l,'Location','best');
subplot(2,2,3);
plot(tau_vector,Results.Ril_j)
title('R_{ik \rightarrow j}');
legend(l,'Location','best');
subplot(2,2,4);
plot(tau_vector,Results.Sil_j)
title('S_{ik \rightarrow j}');
legend(l,'Location','best');
sgtitle({'Partial Information Decomposition',['j=RR' ' ' 'i=RESP' ' ' 'k=ABP']},'FontWeight','bold')

end
