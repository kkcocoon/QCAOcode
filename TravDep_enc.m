function [blkcodej,blkimg_rec] = TravDep_enc(blktree,j,Tk,tlev,blkimg_rec)

blkcodej = int8([]);   

jStack = zeros(32,1,'int32'); 
si = 1;
jStack(si) = j;

blksize2 = length(blkimg_rec);
tlen = length(blktree);
leafb = tlen-blksize2+1-1;

while si>0
    j = jStack(si); si = si-1; 

    if abs(blktree(j))>=2*Tk
        if (4*j+1<=tlen)  
            si=si+1; jStack(si)=4*j+1;
            si=si+1; jStack(si)=4*j;
            si=si+1; jStack(si)=4*j-1;
            si=si+1; jStack(si)=4*j-2;             
        end
        continue;
    end

    if j>1 & mod(j,4)==1 & abs(blktree(j-1))<Tk & abs(blktree(j-2))<Tk & abs(blktree(j-3))<Tk
        if (4*j+1<=tlen) 
            si=si+1; jStack(si)=4*j+1;
            si=si+1; jStack(si)=4*j;
            si=si+1; jStack(si)=4*j-1;
            si=si+1; jStack(si)=4*j-2;        
        else
            blkcodej = [blkcodej, int8((sign(blktree(j))~=1))]; 
            blkimg_rec(j-leafb) = sign(blktree(j))*(Tk + Tk/2);      
        end
    elseif abs(blktree(j)) >= Tk
        
        blkcodej = [blkcodej, 1];
        
        if (4*j+1<=tlen)  
            si=si+1; jStack(si)=4*j+1;
            si=si+1; jStack(si)=4*j;
            si=si+1; jStack(si)=4*j-1;
            si=si+1; jStack(si)=4*j-2;              
        else
            blkcodej = [blkcodej, int8((sign(blktree(j))~=1))]; 
            blkimg_rec(j-leafb) = sign(blktree(j))*(Tk + Tk/2);      
        end
    else
        blkcodej = [blkcodej, 0];
    end
end