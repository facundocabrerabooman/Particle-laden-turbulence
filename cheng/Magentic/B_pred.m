function [B_pred_indv,B_pred0] = B_pred(Colis_Param,z) 
    for i =1:size(Colis_Param.I,2)
        B_pred_indv(:,i) = 2*pi*1e-7*Colis_Param.R(i)^2*Colis_Param.N(i)*Colis_Param.I(i)./sqrt(((z - Colis_Param.Z(i)).^2 + Colis_Param.R(i)^2).^3);
    end
    B_pred_indv = B_pred_indv*1e4; % Gauss
    B_pred0 = sum(B_pred_indv,2);
end