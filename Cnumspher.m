function Cnext=Cnumspher(v,dt,D,x,Cstart)
dx=x(2)-x(1);

Cnext=zeros(length(Cstart),1);
Cnext(1,1)=Cstart(1,1); %boundary cond.
for i=2:length(x)-1
    Cnext(i,1)=Cstart(i,1) + (dt*D/dx)*((2/x(i))*(Cstart(i+1)-Cstart(i))+(1/dx)*(Cstart(i+1)-2*Cstart(i)+Cstart(i-1)));
    if i==length(x)-1
        Cnext(i+1,1)=v.C_0; %boundary cond.
    end
end


