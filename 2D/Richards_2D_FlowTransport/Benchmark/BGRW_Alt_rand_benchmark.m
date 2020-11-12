function [c]=BGRW_Alt_rand_benchmark(c0,cBC,dthtc,I,J,dx,dy,i0,j0,dtt,Lc,Vx,Vy,D1,D2)
%% BGRW function (random) - transport step in coupled flow & trnasport 2D

%% Initial condition
N=10^24; 
n0=round(c0*N); nn=zeros(J,I); nBC=round(cBC*N);
%% BGRW solution
v=Vy*dtt/Lc/dx; u=Vx*dtt/Lc/dy;
ru=2*D1*dtt/Lc/dx^2*ones(J,I);
rv=2*D2*dtt/Lc/dy^2*ones(J,I);
n=n0;
restr=0; restjump=0; restjumpx=0; restjumpy=0; restf=0;
    for y=1:J 
        for x=1:I 
            if n(y,x) > 0
                rx=ru(y,x); ry=rv(y,x); r=rx+ry;
                restr=n(y,x)*(1-r)+restr; nsta=floor(restr);
                restr=restr-nsta; njump=n(y,x)-nsta;
                nn(y,x)=nn(y,x)+nsta;
                restjump=njump*ry/r+restjump;
                njumpy=floor(restjump); restjump=restjump-njumpy;
                njumpx=njump-njumpy;
                if(njumpy)>0
                    restjumpy=njumpy*0.5*(1-v(y,x)/ry)+restjumpy;
                    nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                    nj(2)=njumpy-nj(1);
                    if y==1
                        nn(2,x)=nn(2,x)+nj(2); 
                    elseif y==J
                        nn(J-1,x)=nn(J-1,x)+nj(1); 
                    else
                        for i=1:2
                            yd=y+(2*i-3);
                            nn(yd,x)=nn(yd,x)+nj(i);
                        end
                    end
                end                
                if(njumpx)>0
                    restjumpx=njumpx*0.5*(1-u(y,x)/rx)+restjumpx;
                    nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                    nj(2)=njumpx-nj(1);
                    if x==1
                        nn(y,2)=nn(y,2)+nj(2); 
                    elseif x==I
                        nn(y,I-1)=nn(y,I-1)+nj(1); 
                    else
                        for i=1:2
                            xd=x+(2*i-3);
                            nn(y,xd)=nn(y,xd)+nj(i);
                        end
                    end
                end
            end
        end
    end
    %% Boundary conditions - concentration
    %%%% BCYLeft/Right
    nn(:,1)=nn(:,2); % no flux on \Gamma_{N}
    nn(1:j0,I)= nBC(1:j0,I);%zeros(j0,1); % 0.001*ones(j0,1); % nn(1:j0,I-1); % nBC(1:j0,I)/10;  % \Gamma_{D2}
    nn(j0+1:J,I)=nn(j0+1:J,I-1); % no flux on \Gamma_{N}
    %%%% BCXBottom/Upper
    nn(1,2:I-1)=nn(2,2:I-1); % no flux on \Gamma_{N}
    nn(J,1:i0)=nBC(J,1:i0);  %\Gamma_{D1}
    nn(J,i0+1:I-1)=nn(J-1,i0+1:I-1);
    %% Source term concentration
    restf = dthtc(2:J-1,2:I-1)*N +1e-3*nn(2:J-1,2:I-1)./(1+nn(2:J-1,2:I-1))*N*dtt/Lc+restf; 
    nf=floor(restf); 
    nn(2:J-1,2:I-1)=nn(2:J-1,2:I-1)+nf;
    n=nn; 
    c=n/N; 
