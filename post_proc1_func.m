function [G, I,box_out] = post_proc1_func(fn, v, M)
    flag_plot=0;
    best=5;
    id=find(M(:,1)==fn);
    id_sort=sort(id,'ascend');
    id_sort10=id_sort(1:best);
    box_out=M(id_sort10,3:6);
    v.CurrentTime=fn/30;
    I=readFrame(v);
    Nlin=size(I,1);
    Ncol=size(I,2);
    if flag_plot==1
        figure; imshow(I)
        hold on
        cc=gray;
    end
    kentro_col=round((M(id_sort10,4)+M(id_sort10,6))/2);
    kentro_lin=round((M(id_sort10,3)+M(id_sort10,5))/2);
    [jj,ii]=meshgrid(1:Ncol,1:Nlin);
    dcol=abs(M(id_sort10,4)-M(id_sort10,6));
    dlin=abs(M(id_sort10,5)-M(id_sort10,3));
    sigcol=round(dcol/4);
    siglin=round(dlin/4);
    iistart=kentro_lin-3*siglin;
    iistart=max([ones(size(iistart,1),1),iistart],[],2);

    iistop=kentro_lin+3*siglin;
    iistop=min([Nlin*ones(size(iistop,1),1),iistop],[],2);

    jjstart=kentro_col-3*sigcol;
    jjstart=max([ones(size(jjstart,1),1),jjstart],[],2);

    jjstop=kentro_col+3*sigcol;
    jjstop=min([Ncol*ones(size(jjstop,1),1),jjstop],[],2);
    G=zeros(Nlin,Ncol);
    for k=1:best
        [ii,jj]=meshgrid(iistart(k):iistop(k),jjstart(k):jjstop(k));
        factor=M(id_sort10(k),2);
        g=factor*exp(-((ii-kentro_lin(k)).^2+(jj-kentro_col(k)).^2)/(2*sigcol(k)*siglin(k)));
        g = g';
        G(iistart(k):iistop(k),jjstart(k):jjstop(k))= G(iistart(k):iistop(k),jjstart(k):jjstop(k))+factor*g;
        a=[M(id_sort10(k),4),M(id_sort10(k),6),M(id_sort10(k),6),M(id_sort10(k),4),M(id_sort10(k),4)];
        b=[M(id_sort10(k),3),M(id_sort10(k),3),M(id_sort10(k),5),M(id_sort10(k),5),M(id_sort10(k),3)];

        if flag_plot==1
            hp = plot(a,b,'o-');
            hp.MarkerFaceColor=cc(max(1,size(cc,1)-k+1),:);
            plot(kentro_col,kentro_lin,'*');
        end
    end
    kentro_lin=(M(id_sort10(k),3)+M(id_sort10(k),5))/2;
end