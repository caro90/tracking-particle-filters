clear all
close all
vid=0;
vid_in=VideoReader('Data/cards_courtyard_B_T.mp4');
net_out=csvread('Data/Coordinates_my_test3.csv');
kalm=0;
alpha=0.9;

sigma_x=2;
sigma_y=sigma_x;
sigma_vx=1;
sigma_vy=sigma_vx;
speed_max=5;

method=1;
vmin=0;
vmax=255;
nb=64;  % 255;
dd=(vmax-vmin)/nb;
cc=linspace(vmin+dd/2,vmax-dd/2,nb);

if -1>0
    pos=[436-16,347-16];
    psize1=18;
    psize2=8;
end

fig1=figure;
fig1.Position=[742 387 532 250];
fig2=figure;
fig2.Position=[742 70 532 250];

psize1=8;
psize2=8;

[G, I]=post_proc1_func(1520, vid_in, net_out);
Nlin=size(I,1);
Ncol=size(I,2);
%% PF init
Npf=500;
Ndim=5;  % [x,y,vx,vy,scale]
PF=zeros(Ndim,Npf);  % [x,y,vx,vy,
PF(1,:)=(Ncol-1)*rand(1,Npf)+1;
PF(2,:)=(Nlin-1)*rand(1,Npf)+1;
% PF(1,:)=psize2*randn(1,Npf)+pos(2);
% PF(2,:)=psize1*randn(1,Npf)+pos(1);

% PF(1,:)=pos(1)+5*sigma_x*randn(1,Npf);
% PF(2,:)=pos(2)+5*sigma_y*randn(1,Npf);

id=find(PF(1,:)>Ncol | PF(1,:)<1 | PF(2,:)>Nlin | PF(2,:)<1);
PF(1,id)=(Ncol-1)*rand(1,length(id))+1;
PF(2,id)=(Nlin-1)*rand(1,length(id))+1;


PF(3:4,:)=2*speed_max*rand(2,Npf)-speed_max;
PF(5,:)=rand(1,Npf)+0.5;
W=ones(1,Npf)/Npf;

figure(fig1); imshow(I,[]);
hold on;
%plot(pos(2),pos(1),'o');

resampled=1;
iter=0;
for fn=1520:1530  % 180
    [G, I]=post_proc1_func(fn, vid_in, net_out);
    Nlin=size(I,1);
    Ncol=size(I,2);
    figure(fig1); imshow(uint8(I));
    iter=iter+1;
    N=size(I,1);
    M=size(I,2);
  
    %% sample PFs
    cv2=sum((W-1/Npf).^2);
    ess=Npf/(1+cv2);
    if mod(iter,1)==0
        resampled=1;
        for i=1:Npf
            d(i)=selection_roulette_wheel(W);
            PF_temp(:,i)=PF(:,d(i));
        end
        PF=PF_temp;
    else
        resampled=0;
    end
   
    %% evolve PFs
    % position 1-2: x,y
    PFnew(1,:)=PF(1,:)+PF(3,:)+sigma_x*randn(1,Npf);
    PFnew(2,:)=PF(2,:)+PF(4,:)+sigma_y*randn(1,Npf);
    PFnew(3,:)=PF(3,:)+sigma_vx*randn(1,Npf);
    PFnew(4,:)=PF(4,:)+sigma_vy*randn(1,Npf);
    PFnew(5,:)=PF(5,:);
    
    %% evaluate PFs
    for i=1:Npf
        ipf=round(PF(2,i));
        jpf=round(PF(1,i));

        if ipf-psize1>=1 && ipf+psize1<=Nlin && jpf-psize2>=1 && jpf+psize2<=Ncol
            Bhat(i) = G(ipf, jpf);    
        else       
            Bhat(i)=0;
        end
     end
     if resampled==1
         W=Bhat;
     else
         W=W.*Bhat;
     end
     W=W/sum(W);
     
   

    pf_xavg=sum(PF(1,:).*W);
    pf_yavg=sum(PF(2,:).*W);
    plot(pf_xavg,pf_yavg,'o','MarkerFaceColor','g');
    quiver(pf_xavg,pf_yavg,sum(PF(3,:).*W),sum(PF(4,:).*W));
    if -1>0
    if ~exist('hplot_pfavg','var')
        hplot_pfavg=plot(pf_xavg,pf_yavg,'o','MarkerFaceColor','g');
    else
        hplot_pfavg.XData=pf_xavg;
        hplot_pfavg.YData=pf_yavg;
    end

    if ~exist('hquiv_pfavg','var')
        hquiv_pfavg=quiver(pf_xavg,pf_yavg,sum(PF(3,:).*W),sum(PF(4,:).*W));
    else
        hquiv_pfavg.XData=sum(PF(1,:).*W);
        hquiv_pfavg.YData=sum(PF(2,:).*W);
        hquiv_pfavg.UData=sum(PF(3,:).*W);
        hquiv_pfavg.VData=sum(PF(4,:).*W);
    end

    end
    drawnow
    if 1>0
       hq1=quiver(PF(1,:),PF(2,:),PF(3,:),PF(4,:));
    end

    drawnow;
    PF=PFnew;

    drawnow
    if  vid==1
            c1=fr;
        fr1 = getframe(fig1);
        fr2 = getframe(fig2);
        if c1<10
            fname=['.\jpegs\f0000',int2str(c1),'.jpg'];
        elseif c1<100
            fname=['.\jpegs\f000',int2str(c1),'.jpg'];
        elseif c1<1000
            fname=['.\jpegs\f00',int2str(c1),'.jpg'];
        elseif c1<10000
            fname=['.\jpegs\f0',int2str(c1),'.jpg'];
        else
            fname=[pth,'.\jpegs\f',int2str(c1),'.jpg'];
        end

        imwrite([fr1.cdata,fr2.cdata],fname);
    end
end


