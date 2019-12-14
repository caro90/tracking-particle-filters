clear all
close all
vid=0;
vid_in=VideoReader('cards_courtyard_B_T.mp4');
net_out=csvread('Coordinates_my_test3.csv');
kalm=0;
alpha=0.9;
% I=double(imread('D:\costas\data\video\dataset2014\PETS2006\input\in000088.jpg'));
% pth='D:\costas\data\video\dataset2014\PETS2006\input\';
% pth='D:\costas\data\video\dataset2014\pedestrians\input\';
% pth='D:\costas\data\video\MAMo\Malaya2\';
% pth='D:\costas\data\video\tracking\httpcvlab.hanyang.ac.krtracker_benchmarkdatasets.html\Walking\img\';

sigma_x=2;
sigma_y=sigma_x;
sigma_vx=1;
sigma_vy=sigma_vx;
speed_max=5;

% if -1>0
%     fr_start=448;
%     fr=fr_start;
%     if fr<10
%         fname=['in00000',num2str(fr),'.jpg'];
%     elseif fr<100
%         fname=['in0000',num2str(fr),'.jpg'];
%     elseif fr<1000
%         fname=['in000',num2str(fr),'.jpg'];
%     end
% end
% 
% if 1>0
%     fr_start=2;
%     fr=fr_start;
%     pth='D:\costas\data\video\tracking\BoBot\Vid_B_cup\';
%     if fr<10
%         fname=['img000',num2str(fr),'.jpg'];
%     elseif fr<100
%         fname=['img00',num2str(fr),'.jpg'];
%     elseif fr<1000
%         fname=['img0',num2str(fr),'.jpg'];
%     end
% end
% 
% if -1>0
%     fr_start=2;
%     fr=fr_start;
%     if fr<10
%         fname=['000',num2str(fr),'.jpg'];
%     elseif fr<100
%         fname=['00',num2str(fr),'.jpg'];
%     elseif fr<1000
%         fname=['0',num2str(fr),'.jpg'];
%     end
% end


%I=double(imread([pth,fname]));
%[xx,yy]=getpts;
% pos(2)=round((xx(1)+xx(end))/2);
% pos(1)=round((yy(1)+yy(end))/2);
% psize1=round((max(yy)-min(yy))/2);
% psize2=round((max(xx)-min(xx))/2);

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

% pat=double(I(pos(1)-psize1:pos(1)+psize1,pos(2)-psize2:pos(2)+psize2,:));
% patR=pat(:,:,1);
% patG=pat(:,:,2);
% patB=pat(:,:,3);
% H=zeros(nb,1);
% [Hr0]=hist(patR(:),cc);
% [Hg0]=hist(patG(:),cc);
% [Hb0]=hist(patB(:),cc);
% Hr0=Hr0/sum(Hr0(:));
% Hg0=Hg0/sum(Hg0(:));
% Hb0=Hb0/sum(Hb0(:));

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



%%

% nbins=50;
% vmin=0;
% vmax=255;
% dd=(vmax-vmin)/nbins;
% cc=linspace(vmin+dd/2,vmax-dd/2,nbins);


resampled=1;
iter=0;
for fn=1520:1530  % 180
    [G, I]=post_proc1_func(fn, vid_in, net_out);
    Nlin=size(I,1);
    Ncol=size(I,2);
    % I=round(255*rgb2hsv(I));
    figure(fig1); imshow(uint8(I));

    
    iter=iter+1;
    if -1>0
    if fr<10
        fname=['in00000',num2str(fr),'.jpg'];
    elseif fr<100
        fname=['in0000',num2str(fr),'.jpg'];
    elseif fr<1000
        fname=['in000',num2str(fr),'.jpg'];
    end
    end
    
    if -1>0
    if fr<10
        fname=['img000',num2str(fr),'.jpg'];
    elseif fr<100
        fname=['img00',num2str(fr),'.jpg'];
    elseif fr<1000
        fname=['img0',num2str(fr),'.jpg'];
    end
    end    

    if -1>0
        if fr<10
            fname=['000',num2str(fr),'.jpg'];
        elseif fr<100
            fname=['00',num2str(fr),'.jpg'];
        elseif fr<1000
            fname=['0',num2str(fr),'.jpg'];
        end
    end
%     I=imread([pth,fname]);
% %     I=round(255*rgb2hsv(I));
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
               
                if -1>0
                    Ir_pf=double(I(ipf-psize1:ipf+psize1,jpf-psize2:jpf+psize2,1));
                    Ig_pf=double(I(ipf-psize1:ipf+psize1,jpf-psize2:jpf+psize2,2));
                    Ib_pf=double(I(ipf-psize1:ipf+psize1,jpf-psize2:jpf+psize2,3));
                    [Hr_pf]=hist(Ir_pf(:),cc);
                    [Hg_pf]=hist(Ig_pf(:),cc);
                    [Hb_pf]=hist(Ib_pf(:),cc);
                    Hr_pf=Hr_pf/sum(Hr_pf(:));
                    Hg_pf=Hg_pf/sum(Hg_pf(:));
                    Hb_pf=Hb_pf/sum(Hb_pf(:));
                    Bhat1=sum(sqrt(Hr0.*Hr_pf));
                    Bhat2=sum(sqrt(Hg0.*Hg_pf));
                    Bhat3=sum(sqrt(Hb0.*Hb_pf));
                    Bhat(i)=Bhat1*Bhat2*Bhat3;
                end
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
     
   
% %     id=find(PFnew(1,:)>Ncol);
%     
% 
%     I1=zeros(N,M);
%     I2=zeros(N,M);
%     I3=zeros(N,M);
%     
%     % calc 
%     if method==1        
%         Ir=double(I(:,:,1));
%         Ig=double(I(:,:,2));
%         Ib=double(I(:,:,3));
% %         [Hr]=hist(patR(:),cc);
% %         [Hg]=hist(patG(:),cc);
% %         [Hb]=hist(patB(:),cc);
% %         Hr=Hr/sum(Hr(:));
% %         Hg=Hg/sum(Hg(:));
% %         Hb=Hb/sum(Hb(:));
%         for v=0:255
%             id1=find(Ir==v);
%             id2=find(Ig==v);
%             id3=find(Ib==v);
%             [temp,n]=min(abs(v-cc));
%             I1(id1)=Hr0(n);
%             I2(id2)=Hg0(n);
%             I3(id3)=Hb0(n);
%         end
%            
%         if -1>0
%         for i=1:N
%             for j=1:M
%                 v1=Ir(i,j);
%                 v2=Ig(i,j);
%                 v3=Ib(i,j);
%                 [temp,n]=min(abs(v1-cc));
%                 I1(i,j)=Hr(n);
%                 [temp,n]=min(abs(v2-cc));
%                 I2(i,j)=Hg(n);
%                 [temp,n]=min(abs(v3-cc));
%                 I3(i,j)=Hb(n);
%             end
%         end
%         end
%         
%         Ip=I1.*I2.*I3;
%         figure(fig1); imshow(I,[]);
%         hold on;
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
        %hold off
%         figure(fig2);
%         if iter==1
%             him=imshow(Ip,[]);
%         else
%             ch=findobj(allchild(fig2.CurrentAxes),'Type','Image');
%             ch.CData=Ip;
%         end
%         hold on;

        if 1>0
        %plot(pos(2),pos(1),'o');
%         if ~exist('hp1','var')
%             hp1=plot(PF(1,:),PF(2,:),'.y','MarkerFaceColor','g','Markersize',4);
%         else
%             hp1.XData=PF(1,:);
%             hp1.YData=PF(2,:);
%         end
%         if ~exist('hq1','var')
            hq1=quiver(PF(1,:),PF(2,:),PF(3,:),PF(4,:));
%         else
%             hq1.XData=PF(1,:);
%             hq1.YData=PF(2,:);
%             hq1.UData=PF(1,:);
%             hq1.VData=PF(2,:);
%         end
        end
        %hold off;
        drawnow;
        PF=PFnew;
        
        
        %%
        
        drawnow
        if  vid==1

%             c1=c1+1;
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
%         fr = getframe(gcf);
        imwrite([fr1.cdata,fr2.cdata],fname);
        %                         aviobj = addframe(aviobj,fr);
    end
end

%end


if -1>0
pat=double(I(436:436+31,347:347+31,:));
[H2d,cc]=hist2D(pat,0,255, 256);
% I1=project_hist_fun(I,H2d,cc1,cc2);
I1=project_hist_fun(double(I),H2d,cc,cc);

% [H,C]=hist(I(:),20);


if xroma==3
    R=I(:,:,1);
    G=I(:,:,2);
    B=I(:,:,3);
    maskR=I(436:436+31,347:347+31,1);
    maskR=I(436:436+31,347:347+31,2);
    maskR=I(436:436+31,347:347+31,3);
else
    I=double(rgb2gray(I));
    mask=I(436:436+31,347:347+31);
end


Ksize=10;
rmax=Ksize*sqrt(2);
i0=429;
j0=356;
hold on;
A=I(i0-Ksize:i0+Ksize,i0-Ksize:i0+Ksize);
for i=-Ksize:Ksize
    for j=-Ksize:Ksize
        plot(j,i,'.')
        v=I(i0+i,j0+j);
        r=sqrt((i-Ksize)^2+(j-Ksize)^2)/rmax;
        [temp,n]=min(abs(v-cc));
        H(n)=H(n)+abs(1-r)/Ksize;
    end
end




for i=1:N
    for j=1:M
        v1=I(i,j);
         [temp,n]=min(abs(v1-cc));
         I1(i,j)=H(n);
    end
end
    
end

if -1>0
   if fr<10
        fname=['in0000',num2str(fr),'.jpg'];
    elseif fr<100
        fname=['in000',num2str(fr),'.jpg'];
    elseif fr<1000
        fname=['in00',num2str(fr),'.jpg'];
    end 
end


        