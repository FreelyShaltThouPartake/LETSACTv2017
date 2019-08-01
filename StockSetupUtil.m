function [AstockUtil, BstockUtil] = StockSetupUtil(N_Util, Tcrit_Util, Y, initfleetUtil, P_Util)

% Old Util >= 0 constraint

% rowindex = constraint #
% colindex = variable #
% squish converts 3D or 2D index position to equivalent vector position

% AoldUtil = spalloc(Y*N_Util*Tcrit_Util, length(X0Util),nonzeros);
AstockUtil = sparse(Y*N_Util*Tcrit_Util, N_Util*Y*(1+Tcrit_Util));
BstockUtil = zeros(Y*N_Util*Tcrit_Util, 1);

for k=1:Y
    k
    for i=1:N_Util
        for j=1:Tcrit_Util
            % constaint # given by rowindex
            rowindex = i + N_Util*(j-1) + N_Util*Tcrit_Util*(k-1);

            if j>=k
                temp = 1;
                for m=j-k+1:j
                    temp = temp*P_Util(i,m);
                end
                BstockUtil(rowindex) = BstockUtil(rowindex) - initfleetUtil(i,j-k+1)*temp;

                startindex = N_Util*Y;
                for m=1:k-1
                    temp = 1;
                    for n=m+j-k+1:j
                        temp = temp*P_Util(i,n);
                    end
                    colindex = startindex + squish(i,m+j-k,m,N_Util,Tcrit_Util);
                    AstockUtil(rowindex, colindex) = AstockUtil(rowindex, colindex) - temp;
                end

                AstockUtil(rowindex, startindex+squish(i,j,k,N_Util,Tcrit_Util)) = AstockUtil(rowindex, startindex+squish(i,j,k,N_Util,Tcrit_Util)) - 1;

            elseif j<k
                
                startindex = 0;
                temp = 1;
                for m=1:j
                    temp = temp*P_Util(i,m);
                end
                AstockUtil(rowindex, startindex+squish(i,k-j,1,N_Util,Tcrit_Util)) = AstockUtil(rowindex, startindex+squish(i,k-j,1,N_Util,Tcrit_Util)) + temp;

                startindex = N_Util*Y;
                for m=1:j-1
                    temp = 1;
                    for n=m+1:j
                        temp = temp*P_Util(i,n);
                    end
                    colindex = startindex+squish(i,m,m+k-j,N_Util,Tcrit_Util);
                    AstockUtil(rowindex, colindex) = AstockUtil(rowindex, colindex) - temp;
                end

                AstockUtil(rowindex, startindex+squish(i,j,k,N_Util,Tcrit_Util)) = AstockUtil(rowindex, startindex+squish(i,j,k,N_Util,Tcrit_Util)) - 1;
            end
        end
    end
end