function [ar,T]=ARTimeseries2(ro,num_n,sigma,dg)
    % tạo dãy quan sát {x_t} và dùng nó để dự đoán hệ số tự hồi quy 
    % p là tự tương quan cho trước
    % bậc của đa thức tự hồi quy :dg
    %num_n : số lượng điểm muốn sinh ra
    function x = block_levinson(y, L)
    %BLOCK_LEVINSON Block Levinson recursion for efficiently solving 
    %symmetric block Toeplitz matrix equations.
    %   BLOCK_LEVINSON(Y, L) solves the matrix equation T * x = y, where T 
    %   is a symmetric matrix with block Toeplitz structure, and returns the 
    %   solution vector x. The matrix T is never stored in full (because it 
    %   is large and mostly redundant), so the input parameter L is actually 
    %   the leftmost "block column" of T (the leftmost d columns where d is 
    %   the block dimension).
    
    %   Author: Keenan Pepper
    %   Last modified: 2007-12-23
    
    %   References:
    %     [1] Akaike, Hirotugu (1973). "Block Toeplitz Matrix Inversion".
    %     SIAM J. Appl. Math. 24 (2): 234-241
    
    s = size(L);
    d = s(2)    ;             % Block dimension
    N = s(1) / d  ;           % Number of blocks
    
    B = reshape(L, [d,N,d]);  % This is just to get the bottom block row B
    B = permute(B, [1,3,2]);  % from the left block column L
    B = flipdim(B, 3);
    B = reshape(B, [d,N*d]);
    
    f = L(1:d,:)^-1;          % "Forward" block vector
    b = f;                    % "Backward" block vector
    x = f * y(1:d);           % Solution vector
    
    for n = 2:N
        ef = B(:,(N-n)*d+1:N*d) * [f;zeros(d)];
        eb = L(1:n*d,:)' * [zeros(d);b];
        ex = B(:,(N-n)*d+1:N*d) * [x;zeros(d,1)];
        A = [eye(d),eb;ef,eye(d)]^-1;
        fn = [[f;zeros(d)],[zeros(d);b]] * A(:,1:d);
        bn = [[f;zeros(d)],[zeros(d);b]] * A(:,d+1:end);
        f = fn;
        b = bn;
        x = [x;zeros(d,1)] + b * (y((n-1)*d+1:n*d) - ex);
    end
    end
    function [ar,T] = ARTimeseries(n,a,sigma)
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
     n
     for u=p+1:n-1
         sum=0;
         for k=1:p
             sum =sum + b(p,k)*T(u+1-k);
         end
        T(u+1)= -sum + sqrt(sigma)*randn(1);
     end
     [r,lg] = xcorr(T,'biased');
     r(lg<0) = [];
    [ar,~] = levinson(r,p);
    ar(1)=[];
    end
    if dg==1
        L=toeplitz([1]);
    else
        L=toeplitz([1 ro(1:(dg-1))]);
    end
     % ma trận toeplitz biểu diễn các tự tương quan 
    y=ro';
    a = block_levinson(y, L)% tìm được hệ số tự hồi quy
    [ar,T] = ARTimeseries(num_n,a,sigma);
    
end