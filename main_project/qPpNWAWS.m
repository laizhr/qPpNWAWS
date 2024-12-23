function [y,meanlbd_iter,totaliter,y_total,cost_total,sinstar,yita] = qPpNWAWS(X,p,q,maxiter,tol,medIn)
%{
This file is part of codes for the following paper.
For any usage of this file, please cite the following paper as reference:

[1] Zhao-Rong Lai, Xiaotian Wu, Liangda Fang, Ziliang Chen, Cheng Li, 
"De-singularity Subgradient for the $q$-th-Powered $\ell_p$-Norm Weber Location Problem", 
the 39th AAAI Conference on Artificial Intelligence (AAAI, main track), 2025.

At the same time, it is encouraged to cite the following papers with previous related works:

[2] Zhao-Rong Lai, Xiaotian Wu, Liangda Fang, Ziliang Chen, "A De-singularity 
Subgradient Approach for the Extended Weber Location Problem", the 33rd 
International Joint Conference on Artificial Intelligence (IJCAI, main track), 2024.

[3] Yehuda Vardi and Cun-Hui Zhang. The multivariate L1-median and associated data depth.
Proceedings of the National Academy of Sciences of the United States of America, 
97(4):1423¨C1426, Feb. 2000.

[4] Khurrum Aftab, Richard Hartley, and Jochen Trumpf. Generalized weiszfeld 
algorithms for lq optimization. IEEE Transactions on Pattern Analysis and
Machine Intelligence, 37(4):728¨C745, Apr. 2015.

[5] Dingjiang Huang, Junlong Zhou, Bin Li, Steven C. H. Hoi, and Shuigeng Zhou. 
Robust median reversion strategy for online portfolio selection. IEEE Transactions 
on Knowledge and Data Engineering, 28(9):2480-2493, Sep. 2016.

[6] Bin Li, Doyen Sahoo, and Steven C.H. Hoi. OLPS: a toolbox for on-line portfolio 
selection. Journal of Machine Learning Research, 17(1):1242-1246, 2016.

[7] Zhao-Rong Lai and Haisheng Yang. A survey on gaps between mean-variance approach 
and exponential growth rate approach for portfolio optimization. ACM Computing Surveys, 
55(2):1¨C36, Mar. 2023. Article No. 25.
%}

    y_total=[];
    cost_total=[];
    lbd_iter1=0;
    meanlbd_iter=0;
    sinstar=0;
    yita = 0;

    [m,d] = size(X);
    shorten=0.1;

    if (nargin > 6)
        error ('Too many input arguments.') ;
    elseif (nargin < 6)
        medIn = X(1,:)';
            if(nargin < 5)
                tol = 1e-4 ;
                    if (nargin < 4)
                        maxiter = 100 ;
                            if (nargin < 3)
                                error ('The parameters p,q and the data matrix X are missing.') ;
                            end
                    end
            end      
    end
    
    xi=ones(m,1);
    iterdis = 1;
    tolf = 1e-14;
    iterdisf =1;
    iter = 1;
    y = medIn;
    costf=cost(y,X,p,q);
    
    %Begin the iteration
    while  (iter <= maxiter) && (iterdis >0) && (iterdisf>0)
        dist=zeros(m,d);
        lbd_iter=0;
        T=zeros(d,1);
        Tnum=zeros(d,1);
        Tden=zeros(d,1);
        R = zeros(d,1);
        xyd=zeros(m,1);
        a=zeros(d,1);
        xc=10;
        REa=0;
        
        for i=1:m
            for t=1:d
                dist(i,t) = X(i,t)-y(t);
            end
        end
        
        if  any(any(dist==0))
            disp('singular!');
            yita = yita+1;
          for t=1:d
              VT=find(dist(:,t)~=0);
              VTP=find(dist(:,t)==0);
            for i=VT'
                R(t,:)=R(t,:)+q*xi(i)*(norm(X(i,:)'-y,p))^(q-p)*abs(y(t)-X(i,t))^(p-2)*(y(t)-X(i,t));
            end
            if q==1 && p==1
                a(t)=sum(xi(VTP));
                if abs(R(t,:))>a(t)
                    REa=REa+1;
                end
            end
          end
          if norm(R)==0
              sinstar=sinstar+1;
              break
          end
          if q==1 && p==1 && REa==0
              sinstar=sinstar+1;
              break
          end
          if q==1 && 1<p
              for ii=1:m
                xm=X(ii,:)';
                xyd(ii)=norm(y-xm);
                xc=min(xyd);
              end
              if xc==0 
                  if norm(R,p/(p-1))<=xi(xyd==xc)
                    sinstar=sinstar+1;
                    break 
                  else 
                    R=sign(R).*abs(R).^(1/(p-1));
                  end
              end
          end         
          lambda=min(norm(R,p),1);
          cost_cur=cost(y,X,p,q);
          cost_total=[cost_total;cost_cur];
          while(1)
              lbd_iter=lbd_iter+1;
              ytmp=y-lambda*R/norm(R);
              cost_new=cost(ytmp,X,p,q);
                if cost_new<cost_cur || lbd_iter>20
                    meanlbd_iter=lbd_iter;
                    lbd_iter1=lbd_iter1+lbd_iter;  
                    Ty=ytmp;
                    break;                
                end
              lambda=lambda*shorten;
           end                
       else
          for t=1:d
            for i=1:m
                Tnum(t,:)=Tnum(t,:)+xi(i)*(norm(X(i,:)'-y,p))^(q-p)*abs(y(t)-X(i,t))^(p-2)*X(i,t);
                Tden(t,:)=Tden(t,:)+xi(i)*(norm(X(i,:)'-y,p))^(q-p)*abs(y(t)-X(i,t))^(p-2);
            end
            T(t,:)=Tnum(t,:)/Tden(t,:);
          end 
          Ty=T;           
       end
        iterdis = norm((Ty-y),1) - tol*norm(y,1);
        costf1=costf;
        costf=cost(Ty,X,p,q);
        iterdisf = norm((costf-costf1),1) - tolf*norm(costf1,1);
        iter = iter + 1;
        y=Ty;
        y_total=[y_total y];
    end    
    totaliter=iter+meanlbd_iter;
    y=y';
end
    
      
    
        
