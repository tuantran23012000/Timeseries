function [ar,T,e]=ARTimeseries3(n,sigma,alpha)
    % cho hệ số tự tương quan riêng , xây dựng chuỗi {x_t}
    % dự đoán hệ số tự hồi quy
    % n là số lượng điểm muốn sinh ra
    % alpha là hệ số tự tương quan riêng
    p=numel(alpha);
    r=zeros(p,1); % ma trận 1 chiều lưu các tự tương quan mẫu
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
    ar(1)=[];
    end