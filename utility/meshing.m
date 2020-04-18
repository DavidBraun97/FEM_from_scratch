function [Xbar,EDOF,GDOF,flag] = meshing(h,a,b,flag)
% This function takes a rectangular geometry as input and performs a simple
% meshing algorithm.
% Inputs:
% h     element size
% a,b   geometric parameter of rect. geometry
% flag  structured array containing properties
% Outputs:
% Xbar  array holding the absolut position of each node
% EDOF  array holding the element degrees of freedom
% GDOF  array holding the global degrees of freedom
% flag  structured array containing properties

%% MESHING + VISUALIZATION
% Currently this function exclusively supports 2D-bilinear elements.
% Via Else IF you could add different elements (e.g. triangular..)
if flag.type == "2D-bilinear"
    figure('Name','discretized geometry','NumberTitle','off');
    x0=0;
    y0=0;
    xbar = {};
    edof = {};
    e = 0;
    % How will our meshing algorithm work?
    % Using nested for loops we will be iterating over all finite elements. 
    
    % At each element we will be gathering:
    %  all elemental DOF 
    %  all nodal positions

    
    %     --------a-------
    %     | ..           |
    %     | ..           |
    %     | 12-11  16-15 b
    %     | 9 -10  13-14 |
    %     | 4 - 3  8 - 7 |
    %     | 1 - 2  5 - 6 |
    %     ----------------
    %     XXXXXXXXXXXXXXXX     <- beam clamped
    
    
    for j = 1:b/h        % loop in vertical direction
        for i = 1:a/h    % loop in horizontal direction
            e = e + 1;
            % nodal position of element are calculated and added to xbar
            xbar_e = [x0,y0,x0+h,y0,x0+h,y0+h,x0,y0+h];
            xbar{e} = xbar_e';
            edof{e} = [1+(e-1)*8:1:1+(e-1)*8+7]';
            
            % Visualize element
            patch(xbar_e([1:2:8]),xbar_e([2:2:8]),[0.9 0.9 0.9],'EdgeColor','k','Marker','o','MarkerFaceColor','k');
            hold on
            
            % The following loop serves the purpose of visualizing all
            % element DOF using the annotation feature of MATLAB
            for j=1:2
                xy_hor = [x0+(j-1)*h,y0,0.2,0];
                xy_ver = [x0+(j-1)*h,y0,0,0.2];
                ar = annotation('arrow');
                set(ar,'parent', gca, ...
                   'position', xy_hor, ...
                   'Color', [0.4 0.1 0.8], 'LineWidth',1 );
                ar = annotation('arrow');
                set(ar,'parent', gca, ...
                   'position', xy_ver, ...
                   'Color', [0.4 0.1 0.8], 'LineWidth', 1);
                
                hold on
                xy_hor = [x0+(j-1)*h,y0+h,0.2,0];
                xy_ver = [x0+(j-1)*h,y0+h,0,0.2];
                ar = annotation('arrow');
                set(ar,'parent', gca, ...
                    'position', xy_hor, ...
                    'Color', [0.4 0.1 0.8], 'LineWidth', 1);
                ar = annotation('arrow');
                set(ar,'parent', gca, ...
                    'position', xy_ver, ...
                    'Color', [0.4 0.1 0.8], 'LineWidth', 1);
            end
            axis equal 
            grid on
            x0 = x0+h;
            y0=y0;
        end
        x0 = 0;
        y0 = y0+h;
    end  
    % save number of elements to structured array
    flag.numele = e;
    % postprocess output
    Xbar = cell2mat(xbar);
    EDOF = cell2mat(edof);
%% Determine global DOF
    [DOF_e,numele] = size(Xbar);
    explored = reshape(Xbar,[DOF_e*numele,1]);
    explored = reshape(explored,[2,(DOF_e*numele)/2])';
    GDOF = zeros(DOF_e,numele);
    gdof=1;
    n=1;
    for e=1:numele
        for i=1:2:DOF_e
            tuple = explored(n,:);
            
            [id,idx_tuple] = ismember(tuple,explored(1:n-1,:),'rows');
            if ~id
                GDOF(i,e) = gdof;
                GDOF(i+1,e) = gdof+1;
                gdof = gdof+2;
            end
            if id
                % ceil() rounds up to the next adjacent integer
                e_tuple = ceil(idx_tuple/(DOF_e/2));  
                GDOF(i,e) = GDOF((2*idx_tuple-1)-(e_tuple-1)*DOF_e,e_tuple); 
                GDOF(i+1,e) = GDOF(2*idx_tuple-(e_tuple-1)*DOF_e,e_tuple);
            end
            n=n+1;
        end
    end
    formatSpec = 'Mashed the geometry into %i elements.';
    fprintf(formatSpec,numele)
end
end