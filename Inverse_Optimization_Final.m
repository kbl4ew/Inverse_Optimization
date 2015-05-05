%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  FILE:        Final Research Project - Inverse Optimization
%
%  Data:        Tuna Sales Analytics Data collected via R's bayesm package
%               Data has been modified and compressed into myDataFull3.csv
%
%  PURPOSE:     
%               Main Code for the Final Research Project
%
%  AUTHOR:      Kevin Li (kbl4ew@berkeley.edu)
%               Department of IEOR
%               University of California, Berkeley 
%
%  Project:     Inverse Optimization - Applied on Tuna Data
%
%  Reference:   IEOR 290A - Lecture 27 & 28:
%               Estimating an individual utility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Load myDataFull3.csv via data import
% Step 2: Run the following code
u = myData(:, 2:7); % u is a vector containing the 6 input variables
[u_r, u_c] = size(u);
% Decision Variables (Observed as data)
d = myData(:, 8:21); % d is a matrix containing the 14 decision variables
[d_r, d_c] = size(d);
% Finding the indices that correspond to 0% and 100% Display Advertising
display_ad = myData(:,8:14)';
Indices = find(display_ad ~= 0 & display_ad ~= 1); % finding constraints 
                                                   % with slack
% Identifying indices of the lambda matrix to be set to 0
column_no = int64(mod(Indices, 7)+1); 
row_no = int64(floor(Indices ./7));
real_indices = [column_no, row_no];
% Matrices of the constraints
I_1 = eye(7);
I_2 = eye(6);
z_1 = zeros(7);
z_2 = zeros(7,6);

A = [-I_1 z_1; I_1 z_1; z_1 -I_1; zeros(6,7) zeros(6,7)];
B = [z_2; z_2; z_2; -I_2];
% Invoking cvx
cvx_begin sdp
    % Definition of Beta
    variable Q(14,14); 
    variable F(6,14);
    variable k(14, 1);
    variable lambda(27,u_r);
    minimize(0);
    subject to 
        %2*Q*d'- F'*u'- repmat(k, 1, u_r) + A'*lambda == zeros(14, u_r);
        2*Q*d'- F'*u'- repmat(k, 1, u_r) + A'*lambda == zeros(14, u_r);
        %lambda(real_indices)== 0;
        lambda(:) >= 0;
        for i = 1:length(column_no)
            lambda(column_no(i), row_no(i)) == 0;
        end
        Q >= eye(14);
cvx_end