function [check,ar,e,arcoeffs,e_check,ar_theory]=ARTimeseries3(n,sigma,alpha)
    % cho hệ số tự tương quan riêng , xây dựng chuỗi {x_t}
    % dự đoán hệ số tự hồi quy
    % n là số lượng điểm muốn sinh ra
    % alpha là hệ số tự tương quan riêng
    p=numel(alpha);
    b=zeros(p);% ma trận 2 chiều lưu các giá trị a(k,p)
    for i=1:p
        b(i,i)=alpha(i);
    end
    for k=1:p-1
        for i=1:k
            b(k+1,i)=b(k,i)-b(k+1,k+1)*b(k,k-i+1);
        end
    end
    Alpha = zeros(p,1);
    Alpha(p) = b(p,p);
     mm= 1;
     for u = 1:p
           mm = mm*(1-Alpha(u)^2);
     end
     T=zeros(n,1); %khai báo mảng chứa dữ liệu sinh ngẫu nhiên
     T(1)=sqrt(sigma/mm)*randn(1);
     for u=1:p
         m=1;
         sum=0;
         for k=1:u
             sum =sum + b(u,k)*T(u+1-k);
        m=m*(1-Alpha(k)^2);
        end
        T(u+1)= sum + sqrt(sigma*m/mm)*randn(1);
        end
        for u=p+1:n-1
        sum=0;
        for k=1:p
        sum =sum + b(p,k)*T(u+1-k);
        end
        T(u+1)= sum + sqrt(sigma)*randn(1);
        end
        [r_,lg] = xcorr(T,'biased');
        r_(lg<0) = [];
        [ar,e] = levinson(r_,p);
        [~,~,k] = levinson(r_,100);
        stem(k,'filled') 
        conf = sqrt(2)*erfinv(0.95)/sqrt(15000);
        hold on
        [X,Y] = ndgrid(xlim,conf*[-1 1]);
        %[X,Y] = ndgrid(xlim,T);
        plot(X,Y,'--r')
        hold off
        ar(1)=[];
        check=zeros(10,10);
        j=1;
        for i=1:10 
        check(i,:)=T(j:j+9,1);
        j=j+10; 
        end
        [arcoeffs,e_check] = arburg(T,p);
        arcoeffs(1)=[];
        ar_theory=zeros(1,p);
        ar_theory(1,:)=b(p,:);
        
end
