function m=selection_roulette_wheel(p)
N=length(p);
p=p/sum(p);
p=cumsum(p);
r=rand;
id=find(p>r);
m=min(id);
if sum(m(:))==0
    m =100;
end