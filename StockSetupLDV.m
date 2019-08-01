function [AstockLDV, BstockLDV] = StockSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon_LP, initoldLDV, initnewLDV, P_LDV)
% Sets up the matrix and vector for calculating vehicle stock

% rowindex = constraint #
% colindex = variable #
% squish converts 3D or 2D index position to equivalent vector position

AstockLDV = sparse(analysis_horizon_LP*N_LDV*Tcrit_LDV, N_LDV*analysis_horizon_LP*(Tcrit_LDV+1));
% AstockLDV = spalloc(Y*N_LDV*Tcrit_LDV, length(X0LDV),49280);
BstockLDV = zeros(analysis_horizon_LP*N_LDV*Tcrit_LDV, 1);
newLDVend_index = N_LDV*analysis_horizon_LP;

for k=1:analysis_horizon_LP
    for i=1:N_LDV
        for j=1:Tcrit_LDV
            % constaint # given by rowindex
            rowindex = i + N_LDV*(j-1) + N_LDV*Tcrit_LDV*(k-1);

            if j>k
                temp = 1;
                for m=j-k+1:j
                    temp = temp*P_LDV(i,m);
                end
                BstockLDV(rowindex) = BstockLDV(rowindex) - initoldLDV(i,j-k)*temp;

                startindex = newLDVend_index;
                for m=j-k+1:j-1
                    temp = 1;
                    for n=m+1:j
                        temp = temp*P_LDV(i,n);
                    end
                    colindex = startindex+squish(i,m,m-j+k,N_LDV,Tcrit_LDV);
                    AstockLDV(rowindex, colindex) = AstockLDV(rowindex, colindex) - temp;
                end

                AstockLDV(rowindex, startindex+squish(i,j,k,N_LDV,Tcrit_LDV)) = AstockLDV(rowindex, startindex+squish(i,j,k,N_LDV,Tcrit_LDV)) - 1;

            elseif j<k
                
                startindex = 0;
                temp = 1;
                for m=1:j
                    temp = temp*P_LDV(i,m);
                end
                AstockLDV(rowindex, startindex+squish(i,k-j,1,N_LDV,Tcrit_LDV)) = AstockLDV(rowindex, startindex+squish(i,k-j,1,N_LDV,Tcrit_LDV)) + temp;

                startindex = newLDVend_index;
                for m=1:j-1
                    
                    temp = 1;
                    for n=m+1:j
                        temp = temp*P_LDV(i,n);
                    end
                    colindex = startindex+squish(i,m,k-j+m,N_LDV,Tcrit_LDV);
                    AstockLDV(rowindex, colindex) = AstockLDV(rowindex, colindex) - temp;
                end

                AstockLDV(rowindex, startindex+squish(i,j,k,N_LDV,Tcrit_LDV)) = AstockLDV(rowindex, startindex+squish(i,j,k,N_LDV,Tcrit_LDV)) - 1;

            elseif j==k
                
                temp = 1;
                for m=1:j
                    temp = temp*P_LDV(i,m);
                end
                BstockLDV(rowindex) = BstockLDV(rowindex) - initnewLDV(i)*temp;

                startindex = newLDVend_index;
                for m=1:j-1
                    temp=1;
                    for n=m+1:j
                        temp = temp*P_LDV(i,n);
                    end
                    colindex = startindex+squish(i,m,m,N_LDV,Tcrit_LDV);
                    AstockLDV(rowindex, colindex) = AstockLDV(rowindex, colindex) - temp;
                end

                AstockLDV(rowindex, startindex+squish(i,j,k,N_LDV,Tcrit_LDV)) = AstockLDV(rowindex, startindex+squish(i,j,k,N_LDV,Tcrit_LDV)) - 1;
            end
        end
    end
end

end

