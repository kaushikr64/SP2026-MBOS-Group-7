function [eVecs,eVals] = SortEigsPhi(Phi,pri)
% This function takes the monodromy matrix of a periodic orbit (the STM at
% exactly one period), solves for the eigenvalues and eigenvectors, and
% finaly orders them. If there is no prior information given by Phi_1k then
% the function will sort the Eigenvalues in order such that the two values
% of unity pair are put first, the 

[eVecs_init,eVals_init] = eig(Phi);

    
eVecs = zeros(6,6);
eVals = zeros(6,6);

% If there is no prior knowledge of the EigenValue order
if ~exist('pri','var')
    sort_mtrx = (eVals_init*ones(6,1))';
    sort_mtrx(2,:) = (abs(imag(sort_mtrx))>1e-4);
    num_cmplx = sum(sort_mtrx(2,:));
    
    unt_indx = 1;
    if num_cmplx == 0
        real_indx = 3;
        cmplx_indx = 0;
    elseif num_cmplx == 2
        real_indx = 3;
        cmplx_indx = 5;
    else
        real_indx = 0;
        cmplx_indx = 3;
    end
    % vpa(sort_mtrx)
    % real_indx
    % cmplx_indx
    for k = [1,3,5]
        % Finding Unity Pairs
        % Does the EigenValue lie on the real axis?
        if abs(imag(eVals_init(k,k))) < 1e-4
            % Is the EigenValue Unity?
            if abs(real(eVals_init(k,k))-1) < 1e-3
                eVals(unt_indx,unt_indx) = eVals_init(k,k);
                eVals(unt_indx+1,unt_indx+1) = eVals_init(k+1,k+1);
                eVecs(:,unt_indx) = eVecs_init(:,k);
                eVecs(:,unt_indx+1) = eVecs_init(:,k+1);
            % EigenValues is real but not unity
            else 
                eVals(real_indx,real_indx) = eVals_init(k,k);
                eVals(real_indx+1,real_indx+1) = eVals_init(k+1,k+1);
                eVecs(:,real_indx) = eVecs_init(:,k);
                eVecs(:,real_indx+1) = eVecs_init(:,k+1);
                real_indx = real_indx+2;
            end
        % Eigenvalue not unity and is complex
        else
            eVals(cmplx_indx,cmplx_indx) = eVals_init(k,k);
            eVals(cmplx_indx+1,cmplx_indx+1) = eVals_init(k+1,k+1);
            eVecs(:,cmplx_indx) = eVecs_init(:,k);
            eVecs(:,cmplx_indx+1) = eVecs_init(:,k+1);
            cmplx_indx = cmplx_indx+2;
        end
    end

% If we have some prior organization of EigenValues/Vectors
else
    % for i = 1:6
    %     vec_align = real(eVecs_init'*pri(:,i));
    %     [~,max_indx] = max(vec_align);
    %     eVecs(:,i) = eVecs_init(:,max_indx);
    %     eVals(i,i) = eVals_init(max_indx,max_indx);
    %     eVecs_init(:,max_indx) = zeros(6,1);
    % end

    for i = 1:6
        val_align = abs(eVals_init*ones(6,1)-pri(i,i));
        [~,min_indx] = min(val_align);
        eVecs(:,i) = eVecs_init(:,min_indx);
        eVals(i,i) = eVals_init(min_indx,min_indx);
        eVals_init(min_indx,min_indx) = 1e8+1e8i;
    end
    
    % Try again
    % algn_vals = pri.eVals*ones(6,1);
    % algn_vals = algn_vals([1,3,5]);
    % for i = [1,3,5]
    %     clc_vals = eVals_init*ones(6,1);
    %     clc_vals = clc_vals([1,3,5]);
    %     [~,min_indx] = min(abs(clc_vals-algn_vals(k)));
    %     min_indx = min_indx*2-1;
    % 
    %     eVecs(:,i) = eVecs_init(:,max_indx);
    %     eVecs(:,i+1) = eVecs_init(:,max_indx+1);
    %     eVals(i,i) = eVals_init(max_indx,max_indx);
    %     eVals(i+1,i+1) = eVals_init(max_indx+1,max_indx+1);

        
        

end
end