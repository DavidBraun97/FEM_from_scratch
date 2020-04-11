function [] = visualizeresults(Xbar,dbar_e,stresses)
% This function determines the approximated stresses and strains from the
% given discrete deflection D
% Inputs:
% Xbar          array holding the absolut position of each node
% dbar_e        nodal position after deflection D
% stresses      approximated (discrete) stresses

figure('Name','resulting deflection','NumberTitle','off');
[DOF_e,numele] = size(Xbar);

% calculate nodal position after deflection D
Xdef=dbar_e+Xbar;
% plot results (patch color indicating local stress intensity)
v_stress_e = zeros(numele,1);
for e=1:numele
    xbar_e = Xbar(:,e);
    stress_e_mises = 0;
    for i=1:4
        sigma11 = stresses(i+(i-1)*2,e);
        sigma22 = stresses(i+(i-1)*2+1,e);
        sigma12 = stresses(i+(i-1)*2+2,e);
        stress_e_mises = stress_e_mises +...
                         sqrt(sigma11^2-sigma11*sigma22+sigma22^2+3*sigma12^2);
    end
    v_stress_e(e) = stress_e_mises;
    patch(xbar_e([1:2:8]),xbar_e([2:2:8]),[0.9 0.9 0.9],'EdgeColor','k','Marker','o','MarkerFaceColor','k');

    hold on
    axis equal
end
v_stress_e = normalize(v_stress_e,'range',[0.6,1]);
for e=1:numele
    xdef_e = Xdef(:,e);
    r = v_stress_e(e);
    col = [r (1-r*0.7) 0];
    patch(xdef_e([1:2:8]),xdef_e([2:2:8]),col,'EdgeColor','r','Marker','o','MarkerFaceColor','k');
    hold on
    axis equal
end
end

