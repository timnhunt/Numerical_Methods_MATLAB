function [L,U,P] = LU_Decomposition(A)
%LU_Decomposition(A)
%Created by Tim Hunt
%Created on 10/25/19
%This function decomposes a square matrix, A, into its L, U, and P
%components.
%INPUTS: A - The square matrix to be decomposed
%OUTPUTS: L - The lower triangular matrix
%         U - The upper triangular matrix
%         P - The position matrix

%Chack that A is a square matrix
[m,n] = size(A);
if m ~= n
    error('A must be a square matrix');
end

%Create L and P as identity matrices with the same size as A
P = eye(m);
L = P;

%Create U equal to A
U = A;
tic
%Iterate one less than the number of rows in A
for iter = 1:m-1
    
    %Create a column matrix with the values that will determine the
    %partial pivoting of this iteration
    pivot_col = U(iter:m,iter);
    %Find the index of the maximum value in the pivot column
    [~,max_row] = max(abs(pivot_col));
    %Pivot if the maximum is not already at the pivot position
    if max_row ~= 1
        %Swap the row with the max into the row with the pivot position for
        %U and P
        U = pivot(U,max_row + iter - 1,iter,1:m);
        P = pivot(P,max_row + iter - 1,iter,1:m);
        %Swap the necessary values in L
        if iter > 1
            L = pivot(L,max_row + iter - 1,iter,1:iter - 1);
        end
    end
    
    %Perform Gauss elimination
    for row = iter + 1:m
        elim_co=U(row,iter)/U(iter,iter);
        U(row,:)=U(row,:)-elim_co*U(iter,:);
        L(row,iter)=elim_co;
    end
end
   toc
end

%This function swaps rows in a square matrix
function A_new = pivot(A_in,row_from,row_to,col)
A_new = A_in;
A_new(row_to,col) = A_in(row_from,col);
A_new(row_from,col) = A_in(row_to,col);
end