clearvars -except ut ur psyk vort M1 M2

ct=cputime;

sr=1;
br=20*sr;
vel=10;
Re=200;
T=10;                                      
w=0.5;

%Grid Generation and independent variable setup

rsize=(br-sr)/20;
asize=2*pi*sr/60;   
tsize=T/500;
    
theta=0:asize:(2*pi-asize);
rad=sr:rsize:br;
tim=0:tsize:T;

m=length(rad);
n=length(theta);
o=length(tim);

psi=20*ones(m*n,1);
psin=zeros(m,n+1);
coord=zeros(m,n+1);
omg=zeros(m*n,o);
omgn=zeros(m,n);

L=zeros(m*n,m*n);
U=zeros(m*n,m*n);
D=zeros(m*n,m*n);
B=zeros(m*n,1);
R=zeros(m*n,m*n);
D1=zeros(m*n,m*n);

uthe=zeros(length(rad),length(theta));
urad=zeros(length(rad),length(theta));

for i=1:length(theta)
    uthe(length(rad),i)=-vel*sin(theta(i));
    urad(length(rad),i)=vel*cos(theta(i));
end

for i=1:m
    for j=1:n
        coord(i,j)=rad(i)*(cos(theta(j))+1i*sin(theta(j)));
    end
end

for i=1:m
    coord(i,n+1)=rad(i)*(cos(2*pi)+1i*sin(2*pi));
end

%Boundary Condition Setup

for i=1:n
    psi(i)=20;
    psi(n*(m-1)+i)=vel*br*sin(theta(i))+20;
end

% L and U matrix generation

for i=n+1:(m-1)*n
    for j=n+1:n*(m-1)
        if(j>i)
            if(mod(i-1,n)==0)
                if(j==(i+1))
                    U(i,j)=1/(rad((j-mod(j,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    U(i,j)=(1/(2*rad((j-mod(j,n))/n)*rsize)+1/rsize^2);
                elseif(j==i-n)
                    U(i,j)=(-1/(2*rad((j-mod(j,n))/n+2)*rsize)+1/rsize^2);
                elseif((j-(i-1)-n)==0)
                    U(i,j)=1/(rad((j-mod(j,n))/n)*asize)^2;
                else
                    U(i,j)=0;
                end
            elseif(mod(i,n)==0)
                if(j==(i-1))
                    U(i,j)=1/(rad((j+1)/n)*asize)^2;
                elseif(j==(i-n+1))
                    U(i,j)=1/(rad((j-mod(j,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    U(i,j)=(1/(2*rad(j/n-1)*rsize)+1/rsize^2);
                elseif(j==(i-n))
                    U(i,j)=(-1/(2*rad(j/n+1)*rsize)+1/rsize^2);
                else
                    U(i,j)=0;
                end            
            else
                if(j==(i+1))
                    U(i,j)=1/(rad((j-1-mod(j-1,n))/n+1)*asize)^2;
                elseif(i==(j+1))
                    U(i,j)=1/(rad((j-1-mod(j-1,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    U(i,j)=(1/(2*rad((j-mod(j,n))/n)*rsize)+1/rsize^2);
                elseif(j==i-n)
                    U(i,j)=(-1/(2*rad((j-mod(j,n))/n+2)*rsize)+1/rsize^2);
                else
                    U(i,j)=0;
                end
            end
        else
            if(mod(i-1,n)==0)
                if(j==(i+1))
                    L(i,j)=1/(rad((j-mod(j,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    L(i,j)=(1/(2*rad((j-mod(j,n))/n)*rsize)+1/rsize^2);
                elseif(j==i-n)
                    L(i,j)=(-1/(2*rad((j-mod(j,n))/n+2)*rsize)+1/rsize^2);
                elseif((j-(i-1)-n)==0)
                    L(i,j)=1/(rad((j-mod(j,n))/n)*asize)^2;
                else
                    L(i,j)=0;
                end
            elseif(mod(i,n)==0)
                if(j==(i-1))
                    L(i,j)=1/(rad((j+1)/n)*asize)^2;
                elseif(j==(i-n+1))
                    L(i,j)=1/(rad((j-mod(j,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    L(i,j)=(1/(2*rad(j/n-1)*rsize)+1/rsize^2);
                elseif(j==(i-n))
                    L(i,j)=(-1/(2*rad(j/n+1)*rsize)+1/rsize^2);
                else
                    L(i,j)=0;
                end            
            else
                if(j==(i+1))
                    L(i,j)=1/(rad((j-1-mod(j-1,n))/n+1)*asize)^2;
                elseif(i==(j+1))
                    L(i,j)=1/(rad((j-1-mod(j-1,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    L(i,j)=(1/(2*rad((j-mod(j,n))/n)*rsize)+1/rsize^2);
                elseif(j==i-n)
                    L(i,j)=(-1/(2*rad((j-mod(j,n))/n+2)*rsize)+1/rsize^2);
                else
                    L(i,j)=0;
                end
            end
        end
    end
end

% D matrix Generation

for i=(n+1):(m-1)*n
    if(mod(i,n)==0)
        D(i,i)=-2*(1/rsize^2+1/(rad((i-mod(i,n))/n)*asize)^2);
    else
        D(i,i)=-2*(1/rsize^2+1/(rad((i-mod(i,n))/n+1)*asize)^2);
    end
end

for re=2:o

    % B Vector

    for i=n+1:(m-1)*n
        if(i>=n+1 && i<=2*n)
            B(i)=-(-1/(2*rad(2)*rsize)+1/rsize^2)*psi(i-n)-omg(i,re-1);
        elseif(i>=((m-2)*n+1) && i<=((m-1)*n))
            B(i)=-(1/(2*rad(m-1)*rsize)+1/rsize^2)*psi(i+n)-omg(i,re-1);
        else
            B(i)=-omg(i,re-1);
        end
    end

    % Solution of PSI

    res=B((n+1):(m-1)*n)-(L((n+1):(m-1)*n,(n+1):(m-1)*n)+D((n+1):(m-1)*n,(n+1):(m-1)*n)+U((n+1):(m-1)*n,(n+1):(m-1)*n))*psi((n+1):(m-1)*n);
    tem1=0;
    tem2=0;
    psitem=zeros(m*n,1);

    while(norm(res)>(0.0005))
        for i=n+1:(m-1)*n
            for j=i+1:(m-1)*n
                tem2=tem2+U(i,j)*psi(j);
            end
            for j=n+1:i-1
                tem1=tem1+L(i,j)*psitem(j);
            end
            psitem(i)=(B(i)-tem1-tem2)/D(i,i);
            tem1=0;
            tem2=0;
        end
        psi((n+1):(m-1)*n)=psitem((n+1):(m-1)*n);
        res=B((n+1):(m-1)*n)-(L((n+1):(m-1)*n,(n+1):(m-1)*n)+D((n+1):(m-1)*n,(n+1):(m-1)*n)+U((n+1):(m-1)*n,(n+1):(m-1)*n))*psi((n+1):(m-1)*n);
    end
    
    psidec(:,re)=psi;

    for i=1:m
        for j=1:n
            psin(i,j)=psi(n*(i-1)+j);
        end
    end

    for i=1:m
        psin(i,n+1)=psi(n*(i-1)+1);
    end

    % Velocity and Pressure Extraction

    for j=1:length(rad)-1
        for i=1:length(theta)
            uthe(j,i)=-(psin(j+1,i)-psin(j,i))/rsize;
            if(i==1)
                urad(j,i)=(psin(j,1)-psin(j,n))/(rad(j)*asize);
            else
                urad(j,i)=(psin(j,i)-psin(j,i-1))/(rad(j)*asize);
            end
        end
    end

    % R matrix generation

    for i=n+1:(m-1)*n
        for j=n+1:n*(m-1)
            if(mod(i-1,n)==0)
                if(j==(i+1))
                    R(i,j)=1/(Re*rad((j-mod(j,n))/n+1)^2*asize^2)-uthe((j-mod(j,n))/n+1,1)/(2*rad((j-mod(j,n))/n+1)*asize)+abs(uthe((j-mod(j,n))/n+1,1)/rad((j-mod(j,n))/n+1))/(2*asize); %1/(Re*rad((j-mod(j,n))/n+1)^2*asize^2);             %1/(rad((j-mod(j,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    R(i,j)=1/(Re*rsize^2)+1/(2*Re*rad((j-mod(j,n))/n)*rsize)-urad((j-mod(j,n))/n,1)/(2*rsize)+abs(urad((j-mod(j,n))/n,1))/(2*rsize);       %(1/(2*rad((j-mod(j,n))/n)*rsize)+1/rsize^2);
                elseif(j==i-n)
                    R(i,j)=1/(Re*rsize^2)-1/(2*Re*rad((j-mod(j,n))/n+2)*rsize)+urad((j-mod(j,n))/n+2,1)/(2*rsize)+abs(urad((j-mod(j,n))/n+2,1))/(2*rsize);  %1/(Re*rsize^2)+urad((j-mod(j,n))/n+2,1)/rsize-1/(rad((j-mod(j,n))/n+2)*rsize*Re);
                elseif((j-(i-1)-n)==0)
                    R(i,j)=1/(Re*rad((j-mod(j,n))/n)^2*asize^2)+uthe((j-mod(j,n))/n,1)/(2*rad((j-mod(j,n))/n)*asize)+abs(uthe((j-mod(j,n))/n,1)/rad((j-mod(j,n))/n))/(2*asize);    %uthe((j-mod(j,n))/n,1)/(rad((j-mod(j,n))/n)*asize)+1/(Re*rad((j-mod(j,n))/n)^2*asize^2);          %1/(rad((j-mod(j,n))/n)*asize)^2;
                else
                    R(i,j)=0;
                end
            elseif(mod(i,n)==0)
                if(j==(i-1))
                    R(i,j)=1/(Re*rad((j+1)/n)^2*asize^2)+uthe((j+1)/n,n)/(2*rad((j+1)/n)*asize)+abs(uthe((j+1)/n,n)/rad((j+1)/n))/(2*asize);    %uthe((j+1)/n,n)/(rad((j+1)/n)*asize)+1/(Re*rad((j+1)/n)^2*asize^2);                        %1/(rad((j+1)/n)*asize)^2;
                elseif(j==(i-n+1))
                    R(i,j)=1/(Re*rad((j-mod(j,n))/n+1)^2*asize^2)-uthe((j-mod(j,n))/n+1,n)/(2*rad((j-mod(j,n))/n+1)*asize)+abs(uthe((j-mod(j,n))/n+1,n)/rad((j-mod(j,n))/n+1))/(2*asize);          %1/(Re*rad((j-mod(j,n))/n+1)^2*asize^2);                        %1/(rad((j-mod(j,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    R(i,j)=1/(Re*rsize^2)+1/(2*Re*rad(j/n-1)*rsize)-urad(j/n-1,n)/(2*rsize)+abs(urad(j/n-1,n))/(2*rsize);                        %(1/(2*rad(j/n-1)*rsize)+1/rsize^2);
                elseif(j==(i-n))
                    R(i,j)=1/(Re*rsize^2)-1/(2*Re*rad(j/n+1)*rsize)+urad(j/n+1,n)/(2*rsize)+abs(urad(j/n+1,n))/(2*rsize);   %1/(Re*rsize^2)+urad(j/n+1,n)/rsize-1/(rad(j/n+1)*rsize*Re);                        %(-1/(2*rad(j/n+1)*rsize)+1/rsize^2);
                else
                    R(i,j)=0;
                end
            else
                if(j==(i+1))
                    R(i,j)=1/(Re*rad((j-1-mod(j-1,n))/n+1)^2*asize^2)-uthe((j-1-mod(j-1,n))/n+1,mod(i,n))/(2*rad((j-1-mod(j-1,n))/n+1)*asize)+abs(uthe((j-1-mod(j-1,n))/n+1,mod(i,n))/rad((j-1-mod(j-1,n))/n+1))/(2*asize);    %1/(Re*rad((j-1-mod(j-1,n))/n+1)^2*asize^2);                                %1/(rad((j-1-mod(j-1,n))/n+1)*asize)^2;
                elseif(i==(j+1))
                    R(i,j)=1/(Re*rad((j-1-mod(j-1,n))/n+1)^2*asize^2)+uthe((j-1-mod(j-1,n))/n+1,mod(i,n))/(2*rad((j-1-mod(j-1,n))/n+1)*asize)+abs(uthe((j-1-mod(j-1,n))/n+1,mod(i,n))/rad((j-1-mod(j-1,n))/n+1))/(2*asize);      %uthe((j-1-mod(j-1,n))/n+1,mod(i,n))/(rad((j-1-mod(j-1,n))/n+1)*asize)+1/(Re*rad((j-1-mod(j-1,n))/n+1)^2*asize^2);                                %1/(rad((j-1-mod(j-1,n))/n+1)*asize)^2;
                elseif(j==i+n)
                    R(i,j)=1/(Re*rsize^2)+1/(2*Re*rad((j-mod(j,n))/n)*rsize)-urad((j-mod(j,n))/n,mod(i,n))/(2*rsize)+abs(urad((j-mod(j,n))/n,mod(i,n)))/(2*rsize);                            %(1/(2*rad((j-mod(j,n))/n)*rsize)+1/rsize^2);
                elseif(j==i-n)
                    R(i,j)=1/(Re*rsize^2)-1/(2*Re*rad((j-mod(j,n))/n+2)*rsize)+urad((j-mod(j,n))/n+2,mod(i,n))/(2*rsize)+abs(urad((j-mod(j,n))/n+2,mod(i,n)))/(2*rsize);     %1/(Re*rsize^2)+urad((j-mod(j,n))/n+2,mod(i,n))/rsize-1/(rad((j-mod(j,n))/n+2)*rsize*Re);                                %(-1/(2*rad((j-mod(j,n))/n+2)*rsize)+1/rsize^2);
                else
                    R(i,j)=0;
                end
            end
        end
    end
    
    for i=1:n
        R(n+i,i)=1/(Re*rsize^2)-1/(2*Re*rad(2)*rsize)+urad(2,i)/(2*rsize)+abs(urad(2,i))/(2*rsize);                  %1/(Re*rsize^2)+urad(2,i)/rsize-1/(rad(2)*rsize*Re);
    end
    
    for i=1:n
        R((m-2)*n+i,(m-1)*n+i)=1/(Re*rsize^2)+1/(2*Re*rad(m-1)*rsize)-urad(m-1,i)/(2*rsize)+abs(urad(m-1,i))/(2*rsize);
    end

    for i=(n+1):(m-1)*n
        if(mod(i,n)==0)
            D1(i,i)=1/tsize-2/(Re*rsize^2)-2/(Re*rad((i-mod(i,n))/n)^2*asize^2)-abs(uthe((i-mod(i,n))/n,n)/rad((i-mod(i,n))/n))/asize-abs(urad((i-mod(i,n))/n,n))/rsize;                          %1/tsize-urad((i-mod(i,n))/n,n)/rsize-uthe((i-mod(i,n))/n,n)/(rad((i-mod(i,n))/n)*asize)+1/(rad((i-mod(i,n))/n)*rsize*Re)-2/(rsize^2*Re)-2/(Re*rad((i-mod(i,n))/n)^2*rsize);                                       %-2*(1/rsize^2+1/(rad((i-mod(i,n))/n)*asize)^2);
        else
            D1(i,i)=1/tsize-2/(Re*rsize^2)-2/(Re*rad((i-mod(i,n))/n+1)^2*asize^2)-abs(uthe((i-mod(i,n))/n+1,mod(i,n))/rad((i-mod(i,n))/n+1))/asize-abs(urad((i-mod(i,n))/n+1,mod(i,n)))/rsize;                  %1/tsize-urad((i-mod(i,n))/n+1,mod(i,n))/rsize-uthe((i-mod(i,n))/n+1,mod(i,n))/(rad((i-mod(i,n))/n+1)*asize)+1/(rad((i-mod(i,n))/n+1)*rsize*Re)-2/(rsize^2*Re)-2/(Re*rad((i-mod(i,n))/n+1)^2*rsize);                           %-2*(1/rsize^2+1/(rad((i-mod(i,n))/n+1)*asize)^2);
        end
    end

    R=tsize.*(R+D1);
    
    % Solution of OMEGA
    
    for i=1:n
        omg(n*(m-1)+i,re)=0;
        omg(i,re)=2*(psi(i)-psi(i+n))/rsize^2;
    end
    
    omg((n+1):(m-1)*n,re)=R((n+1):(m-1)*n,:)*omg(:,re-1);
    
    for i=1:m
        for j=1:n
            omgn(i,j)=omg(n*(i-1)+j,re);
        end
    end
    
    for i=1:m
        omgn(i,n+1)=omg(n*(i-1)+1,re);
    end
    
    figure(1)
    contour(real(coord),imag(coord),psin,100);
    colormap(jet)
    colorbar()
    M1(re)=getframe;
    
    
    figure(2)
    contour(real(coord),imag(coord),omgn,100);
    colormap(jet)
    colorbar()
    M2(re)=getframe;
    
end

disp(cputime-ct);