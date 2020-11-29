function [ar,T] = ARTimeseries1(n,a,sigma)
    % tìm hệ số tự hồi quy ước lượng được từ bộ hệ số tự hồi quy gốc cho trước
    % xây dựng dãy {x_t}
     p = numel(a);
     Alpha = zeros(p,1);
     Alpha(p) = -a(p);
     b = zeros(p,p);
     for u = 1:p
        b(p,u) = a(u);
     end
     for u= p:-1:2
        for k =1:u-1
            b(u-1,k) = (b(u,k) + Alpha(p)*b(u,u-k))/(1-Alpha(p)^2);
            Alpha(u-1) = -b(u-1,u-1);
        end
     end
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
        T(u+1)= -sum + sqrt(sigma*m/mm)*randn(1);
     end
     for u=p+1:n-1
         sum=0;
         for k=1:p
             sum =sum + b(p,k)*T(u+1-k);
         end
        T(u+1)= -sum + sqrt(sigma)*randn(1);
     end
     [r,lg] = xcorr(T,'biased');
     r(lg<0) = [];
    [ar,e] = levinson(r,p);
    ar(1)=[];
end