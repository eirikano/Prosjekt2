function Cnext=Cnum(v,r,x,Cstart)

Cnext=zeros(length(Cstart),1);
Cnext(1,1)=Cstart(1,1); %boundary cond.
for i=2:length(x)-1
    Cnext(i,1)=Cstart(i,1) + r*(Cstart(i+1,1)-2*Cstart(i,1)+Cstart(i-1,1));
    if i==length(x)-1
        Cnext(i+1,1)=v.C_0; %boundary cond.
    end
end


