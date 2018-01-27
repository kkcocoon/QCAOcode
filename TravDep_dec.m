function [blktree_dec,ci] = TravDep_dec(blktree_dec,j,blkcode,ci,Tk)

jStack = zeros(32,1,'int32');  % stack
si = 1;
jStack(si) = j;

tlen = length(blktree_dec);

while si>0
    j = jStack(si); si = si-1; 

    if abs(blktree_dec(j))>=2*Tk
        if (4*j+1<=tlen)  
            si=si+1; jStack(si)=4*j+1;
            si=si+1; jStack(si)=4*j;
            si=si+1; jStack(si)=4*j-1;
            si=si+1; jStack(si)=4*j-2;             
        end
        continue;
    end
    
    if j>1 & mod(j,4)==1 & abs(blktree_dec(j-1))<Tk & abs(blktree_dec(j-2))<Tk & abs(blktree_dec(j-3))<Tk    
        blktree_dec(j) = Tk;
        if (4*j+1<=tlen)  
            si=si+1; jStack(si)=4*j+1;
            si=si+1; jStack(si)=4*j;
            si=si+1; jStack(si)=4*j-1;
            si=si+1; jStack(si)=4*j-2;    
        else
            sn = (-1)^blkcode(ci);  ci = ci + 1;
            blktree_dec(j) = int32(sn)*(Tk + Tk/2);
        end

    elseif blkcode(ci)==1

        ci = ci + 1;
        blktree_dec(j) = Tk;

        if (4*j+1<=tlen) 
            si=si+1; jStack(si)=4*j+1;
            si=si+1; jStack(si)=4*j;
            si=si+1; jStack(si)=4*j-1;
            si=si+1; jStack(si)=4*j-2;                
        else
            sn = (-1)^blkcode(ci);  ci = ci + 1;
            blktree_dec(j) = int32(sn)*(Tk + Tk/2);
        end
    else
        ci = ci + 1;
    end
end