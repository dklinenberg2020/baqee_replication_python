
function shock = input_shock(data,shock,shock_index,i,ngrid,intensity)

C = data.C;
N = data.N;
CN = data.CN;
CF= data.CF;

dlogt = zeros(C+CN,CN+CF);
dlogtau = zeros(C+CN,CN+CF);

if shock_index == 1 % universal change in iceberg trade cost
    
    dlogtau(:,1:CN) = (log(1+(i/ngrid)*intensity/100)-log(1+((i-1)/ngrid)*intensity/100))*ones(C+CN,CN);
    for c=1:C
        dlogtau(c,(c-1)*N+1:c*N) = 0;
        dlogtau(C+(c-1)*N+1:C+c*N,(c-1)*N+1:c*N) = 0;
    end
    
elseif shock_index == 2 % universal change in tariff
    
    dlogt(:,1:CN) = (log(1+(i/ngrid)*intensity/100)-log(1+((i-1)/ngrid)*intensity/100))*ones(C+CN,CN);
    for c=1:C
        dlogt(c,(c-1)*N+1:c*N) = 0;
        dlogt(C+(c-1)*N+1:C+c*N,(c-1)*N+1:c*N) = 0;
    end
end
    
    dX = sum(data.Omega_total_tilde(1:C+CN,C+1:end).*(dlogt+dlogtau),2); % Matrix of size (C+CN,1)
    
    shock.dlogt = dlogt;
    shock.dlogtau = dlogtau;
    shock.dX = dX;

end